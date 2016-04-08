/*
 *  smooth.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015 The LIME development team
 *
 */

#include "lime.h"


/* Based on Lloyds Algorithm (Lloyd, S. IEEE, 1982) */	
void
smooth(inputPars *par, struct grid *g){
  double mindist;	/* Distance to closest neighbor				*/
  int k=0,j,i;		/* counters									*/
  int sg;		/* counter for smoothing the grid			*/
  int cn;
  int smooth=20;	/* Amount of grid smoothing					*/
  double move[3];	/* Auxillary array for smoothing the grid	*/
  double dist;		/* Distance to a neighbor					*/
  struct cell *dc=NULL; /* Not used at present. */
  unsigned long numCells;

  for(sg=0;sg<smooth;sg++){
    for(i=0;i<par->ncell && !g[i].sink;i++){
      mindist=1e30;
      cn=-1;
      for(k=0;k<g[i].numNeigh;k++){
        if(g[i].neigh[k]->sink==0){
          if(g[i].ds[k]<mindist){
            mindist=g[i].ds[k];
            cn=k;
          }
        }
      }

      if(par->radius-sqrt(g[i].x[0]*g[i].x[0] + g[i].x[1]*g[i].x[1] + g[i].x[2]*g[i].x[2])<mindist) cn=-1;
      			
      if(cn>-1) {
        for(k=0;k<DIM;k++){
          move[k] = g[i].x[k] - g[i].dir[cn].x[k]*0.20;
        }			  
        if((move[0]*move[0]+move[1]*move[1]+move[2]*move[2])<par->radiusSqu &&
           (move[0]*move[0]+move[1]*move[1]+move[2]*move[2])>par->minScaleSqu){
          for(k=0;k<DIM;k++) g[i].x[k]=move[k];
        }
      }
    }
		
    for(i=par->pIntensity;i<par->ncell;i++){
      mindist=1e30;
      cn=-1;
      for(k=0;k<g[i].numNeigh;k++){
        if(g[i].neigh[k]->sink==1){
          if(g[i].ds[k]<mindist){
            mindist=g[i].ds[k];
            cn=k;
          }
        }
      }
      j=g[i].neigh[cn]->id;
      for(k=0;k<DIM;k++){
        g[i].x[k] = g[i].x[k] - (g[j].x[k]-g[i].x[k]) * 0.15;
      }			
      dist=par->radius/sqrt(g[i].x[0]*g[i].x[0] + g[i].x[1]*g[i].x[1] + g[i].x[2]*g[i].x[2]);	
      for(k=0;k<DIM;k++){
        g[i].x[k] *= dist;
      }	
    }
		
    delaunay(DIM, g, (unsigned long)par->ncell, 0, &dc, &numCells);
    distCalc(par, g);	    
    if(!silent) progressbar((double)(sg+1)/(double)smooth, 5);	
    if(dc!=NULL) free(dc);
  }	
}


