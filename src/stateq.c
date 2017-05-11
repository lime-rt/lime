/*
 *  stateq.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2017 The LIME development team
 *
 */

#include "lime.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_errno.h>


/*....................................................................*/
void
getFixedMatrix(molData *md, int ispec, struct grid *gp, int id, gsl_matrix *colli, configInfo *par){
  int ipart,k,l,ti;

  /* Initialize matrix with zeros */
  if(md[ispec].nlev<=0){
    if(!silent) bail_out("Matrix initialization error in stateq");
    exit(1);
  }
  gsl_matrix_set_zero(colli);

  /* Populate matrix with collisional transitions */
  for(ipart=0;ipart<md[ispec].npart;ipart++){
    struct cpData part = md[ispec].part[ipart];
    double *downrates = part.down;
    int di = md[ispec].part[ipart].densityIndex;
    if (di<0) continue;

    for(ti=0;ti<part.ntrans;ti++){
      int coeff_index = ti*part.ntemp + gp[id].mol[ispec].partner[ipart].t_binlow;
      double down = downrates[coeff_index]\
                  + gp[id].mol[ispec].partner[ipart].interp_coeff*(downrates[coeff_index+1]\
                  - downrates[coeff_index]);
      double up = down*md[ispec].gstat[part.lcu[ti]]/md[ispec].gstat[part.lcl[ti]]\
                *exp(-HCKB*(md[ispec].eterm[part.lcu[ti]]-md[ispec].eterm[part.lcl[ti]])/gp[id].t[0]);

      gsl_matrix_set(colli, part.lcu[ti], part.lcl[ti], gsl_matrix_get(colli, part.lcu[ti], part.lcl[ti]) - down*gp[id].dens[di]);
      gsl_matrix_set(colli, part.lcl[ti], part.lcu[ti], gsl_matrix_get(colli, part.lcl[ti], part.lcu[ti]) - up*gp[id].dens[di]);
    }

  }

  /* Does this work with >1 coll. part? */
  double *ctot  = malloc(sizeof(double)*md[ispec].nlev);
  for(k=0;k<md[ispec].nlev;k++){     
    ctot[k]=0.0;
    for(l=0;l<md[ispec].nlev;l++)
      ctot[k] += gsl_matrix_get(colli,k,l);
    gsl_matrix_set(colli,k,k,gsl_matrix_get(colli,k,k) - ctot[k]);
  }
  free(ctot);

  double *girtot = malloc(sizeof(double)*md[ispec].nlev);
  if(par->girdatfile!=NULL){
    for(k=0;k<md[ispec].nlev;k++){
      girtot[k] = 0;
      for(l=0;l<md[ispec].nlev;l++)
        girtot[k] += md[ispec].gir[k*md[ispec].nlev+l];
    }
    for(k=0;k<md[ispec].nlev;k++){
      gsl_matrix_set(colli,k,k,gsl_matrix_get(colli,k,k)+girtot[k]);
      for(l=0;l<md[ispec].nlev;l++){
        if(k!=l){
          if(par->girdatfile!=NULL)
            gsl_matrix_set(colli,k,l,gsl_matrix_get(colli,k,l)-md[ispec].gir[l*md[ispec].nlev+k]);
        }
      }
    }
  }
  free(girtot);

  /* Someone who is not this lazy could fix the loops instead of using transpose. */
  gsl_matrix_transpose(colli);
}

/*....................................................................*/
void
getMatrix(gsl_matrix *matrix, molData *md, int ispec, gridPointData *mp, gsl_matrix *colli){
  int k,l,li;

  /* Initialize matrix by copying the fixed part */
  gsl_matrix_memcpy(matrix, colli);

  /* Populate matrix with radiative transitions */
  for(li=0;li<md[ispec].nline;li++){
    k=md[ispec].lau[li];
    l=md[ispec].lal[li];
    gsl_matrix_set(matrix, k, k, gsl_matrix_get(matrix, k, k)+md[ispec].beinstu[li]*mp[ispec].jbar[li]+md[ispec].aeinst[li]);
    gsl_matrix_set(matrix, l, l, gsl_matrix_get(matrix, l, l)+md[ispec].beinstl[li]*mp[ispec].jbar[li]);
    gsl_matrix_set(matrix, k, l, gsl_matrix_get(matrix, k, l)-md[ispec].beinstl[li]*mp[ispec].jbar[li]);
    gsl_matrix_set(matrix, l, k, gsl_matrix_get(matrix, l, k)-md[ispec].beinstu[li]*mp[ispec].jbar[li]-md[ispec].aeinst[li]);
  }
}

/*....................................................................*/
void
stateq(int id, struct grid *gp, molData *md, const int ispec, configInfo *par\
  , struct blendInfo blends, int nextMolWithBlend, gridPointData *mp\
  , double *halfFirstDs, _Bool *luWarningGiven){

  int t,s,iter,status;
  double *opop,*oopop,*tempNewPop=NULL;
  double diff;
  char errStr[80];

  gsl_matrix *colli  = gsl_matrix_alloc(md[ispec].nlev, md[ispec].nlev);
  gsl_matrix *matrix = gsl_matrix_alloc(md[ispec].nlev, md[ispec].nlev);
  gsl_vector *newpop = gsl_vector_alloc(md[ispec].nlev);
  gsl_vector *rhVec  = gsl_vector_alloc(md[ispec].nlev);
  gsl_permutation *p = gsl_permutation_alloc (md[ispec].nlev);

  opop       = malloc(sizeof(*opop)      *md[ispec].nlev);
  oopop      = malloc(sizeof(*oopop)     *md[ispec].nlev);
  tempNewPop = malloc(sizeof(*tempNewPop)*md[ispec].nlev);

  for(t=0;t<md[ispec].nlev;t++){
    opop[t]=0.;
    oopop[t]=0.;
    gsl_vector_set(rhVec,t,0.);
  }
  gsl_vector_set(rhVec,md[ispec].nlev-1,1.);
  diff=1;
  iter=0;

  getFixedMatrix(md,ispec,gp,id,colli,par);

  while((diff>TOL && iter<MAXITER) || iter<5){
    getjbar(id,md,gp,ispec,par,blends,nextMolWithBlend,mp,halfFirstDs);

    getMatrix(matrix,md,ispec,mp,colli);

    /* this could also be done in getFixedMatrix */ 
    for(s=0;s<md[ispec].nlev;s++){
      gsl_matrix_set(matrix,md[ispec].nlev-1,s,1.);
    }

    status = gsl_linalg_LU_decomp(matrix,p,&s);
    if(status){
      if(!silent){
        sprintf(errStr, "LU decomposition failed for point %d, iteration %d (GSL error %d).", id, iter, status);
        bail_out(errStr);
      }
      exit(1);
    }

    status = gsl_linalg_LU_solve(matrix,p,rhVec,newpop);
    if(status){
      if(!silent && !(*luWarningGiven)){
        *luWarningGiven = 1;
        sprintf(errStr, "LU solver failed for point %d, iteration %d (GSL error %d).", id, iter, status);
        warning(errStr);
        warning("Doing LSE for this point. NOTE that no further warnings will be issued.");
      }
      lteOnePoint(md, ispec, gp[id].t[0], tempNewPop);
      for(s=0;s<md[ispec].nlev;s++)
        gsl_vector_set(newpop,s,tempNewPop[s]);
    }

    diff=0.;
    for(t=0;t<md[ispec].nlev;t++){
      gsl_vector_set(newpop,t,gsl_max(gsl_vector_get(newpop,t),1e-30));
      oopop[t]=opop[t];
      opop[t]=gp[id].mol[ispec].pops[t];

#pragma omp critical
      {
        gp[id].mol[ispec].pops[t]=gsl_vector_get(newpop,t);
      }

      if(gsl_min(gp[id].mol[ispec].pops[t],gsl_min(opop[t],oopop[t]))>minpop){
        diff=gsl_max(fabs(gp[id].mol[ispec].pops[t]-opop[t])/gp[id].mol[ispec].pops[t]\
            ,gsl_max(fabs(gp[id].mol[ispec].pops[t]-oopop[t])/gp[id].mol[ispec].pops[t],diff));
      }
    }
    iter++;
  }

  gsl_matrix_free(colli);
  gsl_matrix_free(matrix);
  gsl_vector_free(rhVec);
  gsl_vector_free(newpop);
  gsl_permutation_free(p);
  free(tempNewPop);
  free(opop);
  free(oopop);
}


