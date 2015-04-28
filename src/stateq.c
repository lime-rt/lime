/*
 *  stateq.c
 *  LIME, The versatile 3D line modeling environment
 *
 *  Created by Christian Brinch on 15/11/06.
 *  Copyright 2006-2014, Christian Brinch,
 *  <brinch@nbi.dk>
 *  Niels Bohr institutet
 *  University of Copenhagen
 *	All rights reserved.
 *
 */

#include "lime.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>


void
stateq(int id, struct grid *g, molData *m, double *pstate, int ispec, inputPars *par){
  int t,s,iter;
  double *opop, *oopop;
  double diff;

  gsl_matrix *matrix = gsl_matrix_alloc(m[ispec].nlev+1, m[ispec].nlev+1);
  gsl_matrix *reduc  = gsl_matrix_alloc(m[ispec].nlev, m[ispec].nlev);
  gsl_vector *newpop = gsl_vector_alloc(m[ispec].nlev);
  gsl_vector *oldpop = gsl_vector_alloc(m[ispec].nlev);
  gsl_matrix *svv    = gsl_matrix_alloc(m[ispec].nlev, m[ispec].nlev);
  gsl_vector *svs    = gsl_vector_alloc(m[ispec].nlev);
  gsl_vector *work   = gsl_vector_alloc(m[ispec].nlev);
  gsl_permutation *p = gsl_permutation_alloc (m[ispec].nlev);

  opop	 = malloc(sizeof(double)*m[ispec].nlev);
  oopop	 = malloc(sizeof(double)*m[ispec].nlev);

  for(t=0;t<m[ispec].nlev;t++){
    opop[t]=0.;
    oopop[t]=0.;
    gsl_vector_set(oldpop,t,0.);
  }
  gsl_vector_set(oldpop,m[ispec].nlev-1,1.);
  diff=1;
  iter=0;

  while((diff>TOL && iter<MAXITER) || iter<5){
    getjbar(id,m,g,par);
    getmatrix(id,matrix,m,g,ispec);
    for(s=0;s<m[ispec].nlev;s++){
      for(t=0;t<m[ispec].nlev-1;t++){
        gsl_matrix_set(reduc,t,s,gsl_matrix_get(matrix,t,s));
      }
      gsl_matrix_set(reduc,m[ispec].nlev-1,s,1.);
    }

    gsl_linalg_LU_decomp(reduc,p,&s);
    if(gsl_linalg_LU_det(reduc,s) == 0){
      gsl_linalg_SV_decomp(reduc,svv, svs, work);
      gsl_linalg_SV_solve(reduc, svv, svs, oldpop, newpop);
      if(!silent) warning("Matrix is singular. Switching to SVD.");
    } else gsl_linalg_LU_solve(reduc,p,oldpop,newpop);

    diff=0.;
    for(t=0;t<m[ispec].nlev;t++){
      gsl_vector_set(newpop,t,gsl_max(gsl_vector_get(newpop,t),1e-30));
      oopop[t]=opop[t];
      opop[t]=g[id].mol[ispec].pops[t];
      g[id].mol[ispec].pops[t]=gsl_vector_get(newpop,t);
      if(gsl_min(g[id].mol[ispec].pops[t],gsl_min(opop[t],oopop[t]))>minpop){
        diff=gsl_max(fabs(g[id].mol[ispec].pops[t]-opop[t])/g[id].mol[ispec].pops[t],gsl_max(fabs(g[id].mol[ispec].pops[t]-oopop[t])/g[id].mol[ispec].pops[t],diff));
      }
    }
    iter++;
  }
  if(diff>TOL) *pstate=diff;
  gsl_matrix_free(matrix);
  gsl_matrix_free(reduc);
  gsl_matrix_free(svv);
  gsl_vector_free(oldpop);
  gsl_vector_free(newpop);
  gsl_vector_free(svs);
  gsl_vector_free(work);
  gsl_permutation_free(p);
  free(opop);
  free(oopop);
}


void
getmatrix(int id, gsl_matrix *matrix, molData *m, struct grid *g, int ispec){
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
    gsl_matrix_set(matrix, k, k, gsl_matrix_get(matrix, k, k)+m[ispec].beinstu[t]*m[ispec].jbar[t]+m[ispec].aeinst[t]);
    gsl_matrix_set(matrix, l, l, gsl_matrix_get(matrix, l, l)+m[ispec].beinstl[t]*m[ispec].jbar[t]);
    gsl_matrix_set(matrix, k, l, gsl_matrix_get(matrix, k, l)-m[ispec].beinstl[t]*m[ispec].jbar[t]);
    gsl_matrix_set(matrix, l, k, gsl_matrix_get(matrix, l, k)-m[ispec].beinstu[t]*m[ispec].jbar[t]-m[ispec].aeinst[t]);
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


