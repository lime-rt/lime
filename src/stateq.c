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


void
stateq(int id, struct grid *g, molData *m, int ispec, inputPars *par, gridPointData *mp, double *halfFirstDs){
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

  opop	 = malloc(sizeof(*opop)*m[ispec].nlev);
  oopop	 = malloc(sizeof(*oopop)*m[ispec].nlev);

  for(t=0;t<m[ispec].nlev;t++){
    opop[t]=0.;
    oopop[t]=0.;
    gsl_vector_set(oldpop,t,0.);
  }
  gsl_vector_set(oldpop,m[ispec].nlev-1,1.);
  diff=1;
  iter=0;

  while((diff>TOL && iter<MAXITER) || iter<5){
    getjbar(id,m,g,par,mp,halfFirstDs);
    getmatrix(id,matrix,m,g,ispec,mp);
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
getmatrix(int id, gsl_matrix *matrix, molData *m, struct grid *g, int ispec, gridPointData *mp){
  int ti,k,l,li,ipart,di;
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
    for(k=0;k<m[ispec].nlev+1;k++){
      for(l=0;l<m[ispec].nlev+1;l++){
        gsl_matrix_set(matrix, k, l, 0.);
        gsl_matrix_set(partner[ipart].colli, k, l, 0.);
      }
    }
  }

  /* Populate matrix with radiative transitions */
  for(li=0;li<m[ispec].nline;li++){
    k=m[ispec].lau[li];
    l=m[ispec].lal[li];
    gsl_matrix_set(matrix, k, k, gsl_matrix_get(matrix, k, k)+m[ispec].beinstu[li]*mp[ispec].jbar[li]+m[ispec].aeinst[li]);
    gsl_matrix_set(matrix, l, l, gsl_matrix_get(matrix, l, l)+m[ispec].beinstl[li]*mp[ispec].jbar[li]);
    gsl_matrix_set(matrix, k, l, gsl_matrix_get(matrix, k, l)-m[ispec].beinstl[li]*mp[ispec].jbar[li]);
    gsl_matrix_set(matrix, l, k, gsl_matrix_get(matrix, l, k)-m[ispec].beinstu[li]*mp[ispec].jbar[li]-m[ispec].aeinst[li]);
  }

  /* Populate matrix with collisional transitions */
  for(ipart=0;ipart<m[ispec].npart;ipart++){
    struct cpData part = m[ispec].part[ipart];
    double *downrates = part.down;
    for(ti=0;ti<part.ntrans;ti++){
      int coeff_index = ti*part.ntemp + g[id].mol[ispec].partner[ipart].t_binlow;
      double down = downrates[coeff_index] + g[id].mol[ispec].partner[ipart].interp_coeff*(downrates[coeff_index+1] - downrates[coeff_index]);
      double up = down*m[ispec].gstat[part.lcu[ti]]/m[ispec].gstat[part.lcl[ti]]*exp(-HCKB*(m[ispec].eterm[part.lcu[ti]]-m[ispec].eterm[part.lcl[ti]])/g[id].t[0]);
      gsl_matrix_set(partner[ipart].colli, part.lcu[ti], part.lcl[ti], down);
      gsl_matrix_set(partner[ipart].colli, part.lcl[ti], part.lcu[ti], up);
    }

    for(k=0;k<m[ispec].nlev;k++){
      partner[ipart].ctot[k]=0.;
      for(l=0;l<m[ispec].nlev;l++)
        partner[ipart].ctot[k] += gsl_matrix_get(partner[ipart].colli,k,l);
    }
  }

  for(k=0;k<m[ispec].nlev;k++){
    for(ipart=0;ipart<m[ispec].npart;ipart++){
      di = m[ispec].part[ipart].densityIndex;
      if(di>=0)
        gsl_matrix_set(matrix,k,k,gsl_matrix_get(matrix,k,k)+g[id].dens[di]*partner[ipart].ctot[k]);
    }
    for(l=0;l<m[ispec].nlev;l++){
      if(k!=l){
        for(ipart=0;ipart<m[ispec].npart;ipart++){
          di = m[ispec].part[ipart].densityIndex;
          if(di>=0)
            gsl_matrix_set(matrix,k,l,gsl_matrix_get(matrix,k,l)-g[id].dens[di]*gsl_matrix_get(partner[ipart].colli,l,k));
        }
      }
    }
    gsl_matrix_set(matrix, m[ispec].nlev, k, 1.);
    gsl_matrix_set(matrix, k, m[ispec].nlev, 0.);
  }

  for(ipart=0;ipart<m[ispec].npart;ipart++){
    gsl_matrix_free(partner[ipart].colli);
    free(partner[ipart].ctot);
  }
  free(partner);
}


