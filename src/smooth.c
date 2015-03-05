/*
 *  smooth.c
 *  LIME, The versatile 3D line modeling environment 
 *
 *  Created by Christian Brinch on 08/07/06.
 *  Copyright 2006-2014, Christian Brinch, 
 *  <brinch@nbi.dk>
 *  Niels Bohr institutet
 *  University of Copenhagen
 *	All rights reserved.
 *
 */

#include <qhull_a.h>
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

      // if(par->radius-sqrt(pow(g[i].x[0],2)+pow(g[i].x[1],2)+pow(g[i].x[2],2))<mindist) cn=-1;
      if(par->radius-sqrt(g[i].x[0]*g[i].x[0] + g[i].x[1]*g[i].x[1] + g[i].x[2]*g[i].x[2])<mindist) cn=-1;
      			
      if(cn>-1) {
        for(k=0;k<DIM;k++){
          move[k] = g[i].x[k] - g[i].dir[cn].x[k]*0.20;
        }			  
        // if(sqrt(move[0]*move[0]+move[1]*move[1]+move[2]*move[2])<par->radius &&
        //    sqrt(move[0]*move[0]+move[1]*move[1]+move[2]*move[2])>par->minScale){
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
      // dist=par->radius/sqrt(pow(g[i].x[0],2)+pow(g[i].x[1],2)+pow(g[i].x[2],2));	
      dist=par->radius/sqrt(g[i].x[0]*g[i].x[0] + g[i].x[1]*g[i].x[1] + g[i].x[2]*g[i].x[2]);	
      for(k=0;k<DIM;k++){
        g[i].x[k] *= dist;
      }	
    }
		
    qhull(par, g);	
    distCalc(par, g);	    
    if(!silent) progressbar((double)(sg+1)/(double)smooth, 5);	
  }	
}


