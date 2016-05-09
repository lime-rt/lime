/*
 *  stateq.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015 The LIME development team
 *
 */

#include "lime.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_errno.h>


void
stateq(int id, struct grid *g, molData *m, int ispec, inputPars *par\
  , gridPointData *mp, double *halfFirstDs, _Bool *luWarningGiven){

  int t,s,iter,status;
  double *opop,*oopop,*tempNewPop=NULL;
  double diff;
  gsl_error_handler_t *defaultErrorHandler=NULL;
  char errStr[80];

  gsl_matrix *matrix = gsl_matrix_alloc(m[ispec].nlev+1, m[ispec].nlev+1);
  gsl_matrix *reduc  = gsl_matrix_alloc(m[ispec].nlev, m[ispec].nlev);
  gsl_vector *newpop = gsl_vector_alloc(m[ispec].nlev);
  gsl_vector *rhVec  = gsl_vector_alloc(m[ispec].nlev);
  gsl_matrix *svv    = gsl_matrix_alloc(m[ispec].nlev, m[ispec].nlev);
  gsl_vector *svs    = gsl_vector_alloc(m[ispec].nlev);
  gsl_vector *work   = gsl_vector_alloc(m[ispec].nlev);
  gsl_permutation *p = gsl_permutation_alloc (m[ispec].nlev);

  opop       = malloc(sizeof(*opop)      *m[ispec].nlev);
  oopop      = malloc(sizeof(*oopop)     *m[ispec].nlev);
  tempNewPop = malloc(sizeof(*tempNewPop)*m[ispec].nlev);

  for(t=0;t<m[ispec].nlev;t++){
    opop[t]=0.;
    oopop[t]=0.;
    gsl_vector_set(rhVec,t,0.);
  }
  gsl_vector_set(rhVec,m[ispec].nlev-1,1.);
  diff=1;
  iter=0;

  defaultErrorHandler = gsl_set_error_handler_off();
  /* While this is off, the gsl_matrix_* etc calls will not exit if they encounter a problem. However the usual problem they would have is out-of-range indices; the respective code is simple enough though that the likelihood of this sort of problem seems low. */

  while((diff>TOL && iter<MAXITER) || iter<5){
    getjbar(id,m,g,par,mp,halfFirstDs);
    getmatrix(id,matrix,m,g,ispec,mp);
    for(s=0;s<m[ispec].nlev;s++){
      for(t=0;t<m[ispec].nlev-1;t++){
        gsl_matrix_set(reduc,t,s,gsl_matrix_get(matrix,t,s));
      }
      gsl_matrix_set(reduc,m[ispec].nlev-1,s,1.);
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
      lteOnePoint(par, m, ispec, g[id].t[0], tempNewPop);
      for(s=0;s<m[ispec].nlev;s++)
        gsl_vector_set(newpop,s,tempNewPop[s]);
    }

    diff=0.;
    for(t=0;t<m[ispec].nlev;t++){
      gsl_vector_set(newpop,t,gsl_max(gsl_vector_get(newpop,t),1e-30));
      oopop[t]=opop[t];
      opop[t]=g[id].mol[ispec].pops[t];

#pragma omp critical
      {
        g[id].mol[ispec].pops[t]=gsl_vector_get(newpop,t);
      }

      if(gsl_min(g[id].mol[ispec].pops[t],gsl_min(opop[t],oopop[t]))>minpop){
        diff=gsl_max(fabs(g[id].mol[ispec].pops[t]-opop[t])/g[id].mol[ispec].pops[t],gsl_max(fabs(g[id].mol[ispec].pops[t]-oopop[t])/g[id].mol[ispec].pops[t],diff));
      }
    }
    iter++;
  }

  gsl_set_error_handler(defaultErrorHandler);

  gsl_matrix_free(matrix);
  gsl_matrix_free(reduc);
  gsl_matrix_free(svv);
  gsl_vector_free(rhVec);
  gsl_vector_free(newpop);
  gsl_vector_free(svs);
  gsl_vector_free(work);
  gsl_permutation_free(p);
  free(tempNewPop);
  free(opop);
  free(oopop);
}


void
getmatrix(int id, gsl_matrix *matrix, molData *m, struct grid *g, int ispec, gridPointData *mp){
  int p,t,k,l,ipart;
  struct getmatrix {
    double *ctot;
    gsl_matrix * colli;
  } *partner;

  partner= malloc(sizeof(struct getmatrix)*m[ispec].npart);

  /* Initialize matrix with zeros */
  for(ipart=0;ipart<m[ispec].npart;ipart++){
    partner[ipart].colli = gsl_matrix_alloc(m[ispec].nlev+1,m[ispec].nlev+1);
    if(m[ispec].nlev>0) partner[ipart].ctot  = malloc(sizeof(double)*m[ispec].nlev);
    else {
      if(!silent)bail_out("Matrix initialization error in stateq");
      exit(0);
    }
    for(t=0;t<m[ispec].nlev+1;t++){
      for(p=0;p<m[ispec].nlev+1;p++){
        gsl_matrix_set(matrix, t, p, 0.);
        gsl_matrix_set(partner[ipart].colli, t, p, 0.);
      }
    }
  }

  /* Populate matrix with radiative transitions */
  for(t=0;t<m[ispec].nline;t++){
    k=m[ispec].lau[t];
    l=m[ispec].lal[t];
    gsl_matrix_set(matrix, k, k, gsl_matrix_get(matrix, k, k)+m[ispec].beinstu[t]*mp[ispec].jbar[t]+m[ispec].aeinst[t]);
    gsl_matrix_set(matrix, l, l, gsl_matrix_get(matrix, l, l)+m[ispec].beinstl[t]*mp[ispec].jbar[t]);
    gsl_matrix_set(matrix, k, l, gsl_matrix_get(matrix, k, l)-m[ispec].beinstl[t]*mp[ispec].jbar[t]);
    gsl_matrix_set(matrix, l, k, gsl_matrix_get(matrix, l, k)-m[ispec].beinstu[t]*mp[ispec].jbar[t]-m[ispec].aeinst[t]);
  }

  /* Populate matrix with collisional transitions */
  for(ipart=0;ipart<m[ispec].npart;ipart++){
    for(t=0;t<m[ispec].ntrans[ipart];t++){
      gsl_matrix_set(partner[ipart].colli, m[ispec].lcu[t], m[ispec].lcl[t], g[id].mol[ispec].partner[ipart].down[t]);
      gsl_matrix_set(partner[ipart].colli, m[ispec].lcl[t], m[ispec].lcu[t], g[id].mol[ispec].partner[ipart].up[t]);
    }

    for(p=0;p<m[ispec].nlev;p++){
      partner[ipart].ctot[p]=0.;
      for(t=0;t<m[ispec].nlev;t++) partner[ipart].ctot[p]+=gsl_matrix_get(partner[ipart].colli,p,t);
    }
  }

  for(p=0;p<m[ispec].nlev;p++){
    for(ipart=0;ipart<m[ispec].npart;ipart++){
      gsl_matrix_set(matrix,p,p,gsl_matrix_get(matrix,p,p)+g[id].dens[ipart]*partner[ipart].ctot[p]);
    }
    for(t=0;t<m[ispec].nlev;t++){
      if(p!=t){
        for(ipart=0;ipart<m[ispec].npart;ipart++){
          gsl_matrix_set(matrix,p,t,gsl_matrix_get(matrix,p,t)-g[id].dens[ipart]*gsl_matrix_get(partner[ipart].colli,t,p));
        }
      }
    }
    gsl_matrix_set(matrix, m[ispec].nlev, p, 1.);
    gsl_matrix_set(matrix, p, m[ispec].nlev, 0.);
  }

  for(ipart=0;ipart<m[ispec].npart;ipart++){
    gsl_matrix_free(partner[ipart].colli);
    free(partner[ipart].ctot);
  }
  free(partner);
}


