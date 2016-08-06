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
getVelocities(configInfo *par, struct grid *g){
  int i,k,j,l;
  double vel[3], x[3];
  
  for(i=0;i<par->pIntensity;i++){

    g[i].v1=malloc(3*g[i].numNeigh*sizeof(double));
    g[i].v2=malloc(3*g[i].numNeigh*sizeof(double));
    g[i].v3=malloc(3*g[i].numNeigh*sizeof(double));

    velocity(g[i].x[0],g[i].x[1],g[i].x[2],g[i].vel);
    
    for(k=0;k<g[i].numNeigh;k++){
      for(j=0;j<3;j++) x[j]=g[i].x[j];		
      for(l=0;l<5;l++){
        velocity(x[0],x[1],x[2],vel);	

        if (l==1) {
	  g[i].v1[3*k]=vel[0]; g[i].v1[3*k+1]=vel[1]; g[i].v1[3*k+2]=vel[2];
        }
        if (l==2) {
          g[i].v2[3*k]=vel[0]; g[i].v2[3*k+1]=vel[1]; g[i].v2[3*k+2]=vel[2];
        }
        if (l==3) {
          g[i].v3[3*k]=vel[0]; g[i].v3[3*k+1]=vel[1]; g[i].v3[3*k+2]=vel[2];
        }
		
        for(j=0;j<3;j++) x[j]=x[j]+(g[i].dir[k].xn[j]*g[i].ds[k])/4.;
      }
    }		
  }
  
  for(i=par->pIntensity;i<par->ncell;i++){
    /* Set velocity values also for sink points (otherwise Delaunay ray-tracing has problems) */
    velocity(g[i].x[0],g[i].x[1],g[i].x[2],g[i].vel);

    g[i].v1=malloc(3*g[i].numNeigh*sizeof(double));
    g[i].v2=malloc(3*g[i].numNeigh*sizeof(double));
    g[i].v3=malloc(3*g[i].numNeigh*sizeof(double));

    for(j=0;j<g[i].numNeigh;j++){
      g[i].v1[3*j]=0.; g[i].v1[3*j+1]=0.; g[i].v1[3*j+2]=0.;
      g[i].v2[3*j]=0.; g[i].v2[3*j+1]=0.; g[i].v2[3*j+2]=0.;
      g[i].v3[3*j]=0.; g[i].v3[3*j+1]=0.; g[i].v3[3*j+2]=0.;
    }
  }
}

void
getVelocities_pregrid(configInfo *par, struct grid *g){
  int i,j;
  /* Is this needed? */ 
  for(i=par->pIntensity;i<par->ncell;i++){
    for(j=0;j<3;j++) g[i].vel[j]=0.;    
  }
}
