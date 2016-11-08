/*
 *  stateq.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2016 The LIME development team
 *
 */

#include "lime.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_errno.h>


void
stateq(int id, struct grid *gp, molData *md, const int ispec, configInfo *par\
  , struct blendInfo blends, int nextMolWithBlend, gridPointData *mp\
  , double *halfFirstDs, _Bool *luWarningGiven){

  int t,s,iter,status;
  double *opop,*oopop,*tempNewPop=NULL;
  double diff;
  char errStr[80];

  gsl_matrix *matrix = gsl_matrix_alloc(md[ispec].nlev+1, md[ispec].nlev+1);
  gsl_matrix *reduc  = gsl_matrix_alloc(md[ispec].nlev, md[ispec].nlev);
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

  while((diff>TOL && iter<MAXITER) || iter<5){
    getjbar(id,md,gp,ispec,par,blends,nextMolWithBlend,mp,halfFirstDs);

    getmatrix(id,matrix,md,gp,ispec,mp);
    for(s=0;s<md[ispec].nlev;s++){
      for(t=0;t<md[ispec].nlev-1;t++){
        gsl_matrix_set(reduc,t,s,gsl_matrix_get(matrix,t,s));
      }
      gsl_matrix_set(reduc,md[ispec].nlev-1,s,1.);
    }

    status = gsl_linalg_LU_decomp(reduc,p,&s);
    if(status){
      if(!silent){
        sprintf(errStr, "LU decomposition failed for point %d, iteration %d (GSL error %d).", id, iter, status);
        bail_out(errStr);
      }
      exit(1);
    }

    status = gsl_linalg_LU_solve(reduc,p,rhVec,newpop);
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

  gsl_matrix_free(matrix);
  gsl_matrix_free(reduc);
  gsl_vector_free(rhVec);
  gsl_vector_free(newpop);
  gsl_permutation_free(p);
  free(tempNewPop);
  free(opop);
  free(oopop);
}

void
getmatrix(int id, gsl_matrix *matrix, molData *md, struct grid *gp, int ispec, gridPointData *mp){
  int ti,k,l,li,ipart,di;
  struct getmatrix {
    double *ctot;
    gsl_matrix * colli;
  } *partner;

  partner = malloc(sizeof(struct getmatrix)*md[ispec].npart);

  /* Initialize matrix with zeros */
  for(ipart=0;ipart<md[ispec].npart;ipart++){
    partner[ipart].colli = gsl_matrix_alloc(md[ispec].nlev+1,md[ispec].nlev+1);
    if(md[ispec].nlev>0) partner[ipart].ctot  = malloc(sizeof(double)*md[ispec].nlev);
    else {
      if(!silent)bail_out("Matrix initialization error in stateq");
      exit(0);
    }
    for(k=0;k<md[ispec].nlev+1;k++){
      for(l=0;l<md[ispec].nlev+1;l++){
        gsl_matrix_set(matrix, k, l, 0.);
        gsl_matrix_set(partner[ipart].colli, k, l, 0.);
      }
    }
  }

  /* Populate matrix with radiative transitions */
  for(li=0;li<md[ispec].nline;li++){
    k=md[ispec].lau[li];
    l=md[ispec].lal[li];
    gsl_matrix_set(matrix, k, k, gsl_matrix_get(matrix, k, k)+md[ispec].beinstu[li]*mp[ispec].jbar[li]+md[ispec].aeinst[li]);
    gsl_matrix_set(matrix, l, l, gsl_matrix_get(matrix, l, l)+md[ispec].beinstl[li]*mp[ispec].jbar[li]);
    gsl_matrix_set(matrix, k, l, gsl_matrix_get(matrix, k, l)-md[ispec].beinstl[li]*mp[ispec].jbar[li]);
    gsl_matrix_set(matrix, l, k, gsl_matrix_get(matrix, l, k)-md[ispec].beinstu[li]*mp[ispec].jbar[li]-md[ispec].aeinst[li]);
  }

  /* Populate matrix with collisional transitions */
  for(ipart=0;ipart<md[ispec].npart;ipart++){
    struct cpData part = md[ispec].part[ipart];
    double *downrates = part.down;
    for(ti=0;ti<part.ntrans;ti++){
      int coeff_index = ti*part.ntemp + gp[id].mol[ispec].partner[ipart].t_binlow;
      double down = downrates[coeff_index]\
                  + gp[id].mol[ispec].partner[ipart].interp_coeff*(downrates[coeff_index+1]\
                  - downrates[coeff_index]);
      double up = down*md[ispec].gstat[part.lcu[ti]]/md[ispec].gstat[part.lcl[ti]]\
                *exp(-HCKB*(md[ispec].eterm[part.lcu[ti]]-md[ispec].eterm[part.lcl[ti]])/gp[id].t[0]);
      gsl_matrix_set(partner[ipart].colli, part.lcu[ti], part.lcl[ti], down);
      gsl_matrix_set(partner[ipart].colli, part.lcl[ti], part.lcu[ti], up);
    }

    for(k=0;k<md[ispec].nlev;k++){
      partner[ipart].ctot[k]=0.;
      for(l=0;l<md[ispec].nlev;l++)
        partner[ipart].ctot[k] += gsl_matrix_get(partner[ipart].colli,k,l);
    }
  }

  for(k=0;k<md[ispec].nlev;k++){
    for(ipart=0;ipart<md[ispec].npart;ipart++){
      di = md[ispec].part[ipart].densityIndex;
      if(di>=0)
        gsl_matrix_set(matrix,k,k,gsl_matrix_get(matrix,k,k)+gp[id].dens[di]*partner[ipart].ctot[k]);
    }
    for(l=0;l<md[ispec].nlev;l++){
      if(k!=l){
        for(ipart=0;ipart<md[ispec].npart;ipart++){
          di = md[ispec].part[ipart].densityIndex;
          if(di>=0)
            gsl_matrix_set(matrix,k,l,gsl_matrix_get(matrix,k,l)-gp[id].dens[di]*gsl_matrix_get(partner[ipart].colli,l,k));
        }
      }
    }
    gsl_matrix_set(matrix, md[ispec].nlev, k, 1.);
    gsl_matrix_set(matrix, k, md[ispec].nlev, 0.);
  }

  for(ipart=0;ipart<md[ispec].npart;ipart++){
    gsl_matrix_free(partner[ipart].colli);
    free(partner[ipart].ctot);
  }
  free(partner);
}


