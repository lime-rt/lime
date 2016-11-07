/*
 *  velospline.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2016 The LIME development team
 *
 */

#include "lime.h"

void
getVelocities(configInfo *par, struct grid *gp){
  int i,k,j,l;
  double vel[3], x[3];
  
  for(i=0;i<par->pIntensity;i++){
    gp[i].v1=malloc(3*gp[i].numNeigh*sizeof(double));
    gp[i].v2=malloc(3*gp[i].numNeigh*sizeof(double));
    gp[i].v3=malloc(3*gp[i].numNeigh*sizeof(double));

    velocity(gp[i].x[0],gp[i].x[1],gp[i].x[2],gp[i].vel);
    
    for(k=0;k<gp[i].numNeigh;k++){
      for(j=0;j<3;j++) x[j]=gp[i].x[j];		
      for(l=0;l<5;l++){
        velocity(x[0],x[1],x[2],vel);	

        if (l==1) {
	  gp[i].v1[3*k]=vel[0]; gp[i].v1[3*k+1]=vel[1]; gp[i].v1[3*k+2]=vel[2];
        }
        if (l==2) {
          gp[i].v2[3*k]=vel[0]; gp[i].v2[3*k+1]=vel[1]; gp[i].v2[3*k+2]=vel[2];
        }
        if (l==3) {
          gp[i].v3[3*k]=vel[0]; gp[i].v3[3*k+1]=vel[1]; gp[i].v3[3*k+2]=vel[2];
        }
		
        for(j=0;j<3;j++) x[j]=x[j]+(gp[i].dir[k].xn[j]*gp[i].ds[k])/4.;
      }
    }		
  }

  for(i=par->pIntensity;i<par->ncell;i++){
    /* Set velocity values also for sink points (otherwise Delaunay ray-tracing has problems) */
    velocity(gp[i].x[0],gp[i].x[1],gp[i].x[2],gp[i].vel);

    gp[i].v1=malloc(3*gp[i].numNeigh*sizeof(double));
    gp[i].v2=malloc(3*gp[i].numNeigh*sizeof(double));
    gp[i].v3=malloc(3*gp[i].numNeigh*sizeof(double));

    for(j=0;j<gp[i].numNeigh;j++){
      gp[i].v1[3*j]=0.; gp[i].v1[3*j+1]=0.; gp[i].v1[3*j+2]=0.;
      gp[i].v2[3*j]=0.; gp[i].v2[3*j+1]=0.; gp[i].v2[3*j+2]=0.;
      gp[i].v3[3*j]=0.; gp[i].v3[3*j+1]=0.; gp[i].v3[3*j+2]=0.;
    }
  }
}

void
getVelocities_pregrid(configInfo *par, struct grid *gp){
  int i,j;
  /* Is this needed? */ 
  for(i=par->pIntensity;i<par->ncell;i++){
    for(j=0;j<3;j++) gp[i].vel[j]=0.;    
  }
}

