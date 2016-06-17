/*
 *  velospline.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015 The LIME development team
 *
 */

#include "lime.h"

void
getVelosplines(inputPars *par, struct grid *g){
  int i,k,j,l,s;
  double v[5], vel[3], x[3], d;
  gsl_matrix *matrix = gsl_matrix_alloc(5,5);
  gsl_vector *a = gsl_vector_alloc(5);
  gsl_vector *y = gsl_vector_alloc(5);
  gsl_permutation *p = gsl_permutation_alloc(5);
  
  for(i=0;i<par->pIntensity;i++){
    g[i].a0=malloc(g[i].numNeigh*sizeof(double));
    g[i].a1=malloc(g[i].numNeigh*sizeof(double));
    g[i].a2=malloc(g[i].numNeigh*sizeof(double));
    g[i].a3=malloc(g[i].numNeigh*sizeof(double));
    g[i].a4=malloc(g[i].numNeigh*sizeof(double));	
    
    for(k=0;k<g[i].numNeigh;k++){
      for(j=0;j<3;j++) x[j]=g[i].x[j];		
      for(l=0;l<5;l++){
        velocity(x[0],x[1],x[2],vel);			
        v[l]=veloproject(g[i].dir[k].xn,vel);
        for(j=0;j<3;j++) x[j]=x[j]+(g[i].dir[k].xn[j]*g[i].ds[k])/4.;
      }
      for(j=0;j<5;j++){
        d=g[i].ds[k]/4*j;
        gsl_matrix_set(matrix,j,0,d*d*d*d);		
        gsl_matrix_set(matrix,j,1,d*d*d);		
        gsl_matrix_set(matrix,j,2,d*d);		
        gsl_matrix_set(matrix,j,3,d);		
        gsl_matrix_set(matrix,j,4,1.);		
      }
      for(j=0;j<5;j++) gsl_vector_set(y,j,v[j]);	
      gsl_linalg_LU_decomp (matrix, p, &s);
      gsl_linalg_LU_solve (matrix, p, y, a);
      g[i].a0[k]=gsl_vector_get(a,4);
      g[i].a1[k]=gsl_vector_get(a,3);
      g[i].a2[k]=gsl_vector_get(a,2);
      g[i].a3[k]=gsl_vector_get(a,1);
      g[i].a4[k]=gsl_vector_get(a,0);
    }		
  }
  
  for(i=par->pIntensity;i<par->ncell;i++){
    g[i].a0=malloc(g[i].numNeigh*sizeof(double));
    g[i].a1=malloc(g[i].numNeigh*sizeof(double));
    g[i].a2=malloc(g[i].numNeigh*sizeof(double));
    g[i].a3=malloc(g[i].numNeigh*sizeof(double));
    g[i].a4=malloc(g[i].numNeigh*sizeof(double));
    for(j=0;j<g[i].numNeigh;j++){
      g[i].a0[j]=0.;
      g[i].a1[j]=0.;
      g[i].a2[j]=0.;
      g[i].a3[j]=0.;
      g[i].a4[j]=0.;
    }
  }
  
  gsl_permutation_free (p);
  gsl_vector_free (a);
  gsl_vector_free (y);
  gsl_matrix_free (matrix);
}


void
getVelosplines_lin(inputPars *par, struct grid *g){
  int i,k,j;
  double v[2];
  
  for(i=0;i<par->pIntensity;i++){
    g[i].a0=malloc(g[i].numNeigh*sizeof(double));
    g[i].a1=malloc(g[i].numNeigh*sizeof(double));		
    for(k=0;k<g[i].numNeigh;k++){
      v[0]=veloproject(g[i].dir[k].xn,g[i].vel);
      v[1]=veloproject(g[i].dir[k].xn,g[i].neigh[k]->vel);
      g[i].a1[k]=(v[0]-v[1])/(0-g[i].ds[k]);
      g[i].a0[k]=v[0];
    }		
  }
  
  for(i=par->pIntensity;i<par->ncell;i++){
    g[i].a0=malloc(g[i].numNeigh*sizeof(double));
    g[i].a1=malloc(g[i].numNeigh*sizeof(double));
    for(j=0;j<g[i].numNeigh;j++){
      g[i].a0[j]=0.;
      g[i].a1[j]=0.;
    }
  }	
}
