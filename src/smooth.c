/*
 *  smooth.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2016 The LIME development team
 *
 */

#include "lime.h"


/* Based on Lloyds Algorithm (Lloyd, S. IEEE, 1982) */	
void
smooth(configInfo *par, struct grid *gp){
  double mindist;	/* Distance to closest neighbor				*/
  int k=0,j,i;		/* counters									*/
  int sg;		/* counter for smoothing the grid			*/
  int cn;
  double move[3];	/* Auxillary array for smoothing the grid	*/
  double dist;		/* Distance to a neighbor					*/
  struct cell *dc=NULL; /* Not used at present. */
  unsigned long numCells;

  for(sg=0;sg<N_SMOOTH_ITERS;sg++){
    delaunay(DIM, gp, (unsigned long)par->ncell, 0, 0, &dc, &numCells);
    distCalc(par, gp);

    for(i=0;i<par->pIntensity;i++){
      mindist=1e30;
      cn=-1;
      for(k=0;k<gp[i].numNeigh;k++){
        if(gp[i].neigh[k]->sink==0){
          if(gp[i].ds[k]<mindist){
            mindist=gp[i].ds[k];
            cn=k;
          }
        }
      }

      if(par->radius-sqrt(gp[i].x[0]*gp[i].x[0] + gp[i].x[1]*gp[i].x[1] + gp[i].x[2]*gp[i].x[2])<mindist) cn=-1;

      if(cn>-1) {
        for(k=0;k<DIM;k++){
          move[k] = gp[i].x[k] - gp[i].dir[cn].x[k]*0.20;
        }			  
        if((move[0]*move[0]+move[1]*move[1]+move[2]*move[2])<par->radiusSqu &&
           (move[0]*move[0]+move[1]*move[1]+move[2]*move[2])>par->minScaleSqu){
          for(k=0;k<DIM;k++) gp[i].x[k]=move[k];
        }
      }
    }
		
    for(i=par->pIntensity;i<par->ncell;i++){
      mindist=1e30;
      cn=-1;
      for(k=0;k<gp[i].numNeigh;k++){
        if(gp[i].neigh[k]->sink==1){
          if(gp[i].ds[k]<mindist){
            mindist=gp[i].ds[k];
            cn=k;
          }
        }
      }
      j=gp[i].neigh[cn]->id;
      for(k=0;k<DIM;k++){
        gp[i].x[k] = gp[i].x[k] - (gp[j].x[k]-gp[i].x[k]) * 0.15;
      }			
      dist=par->radius/sqrt(gp[i].x[0]*gp[i].x[0] + gp[i].x[1]*gp[i].x[1] + gp[i].x[2]*gp[i].x[2]);	
      for(k=0;k<DIM;k++){
        gp[i].x[k] *= dist;
      }	
    }
		
    if(!silent) progressbar((double)(sg+1)/(double)N_SMOOTH_ITERS, 5);	
    free(dc);
  }	
}


