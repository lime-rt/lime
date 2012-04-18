/*
 *  fit.c
 *  LIME, The versatile 3D line modeling environment 
 *
 *  Created by Marco Padovani on 06/04/10.
 *  Copyright 2006-2012, Marco Padovani 
 *
 *  Institut de Ciències de l'Espai (IEEC-CSIC) - 
 *  UAB - Catalunya, Spain. 
 *  All rights reserved.
 *
 */

#include "lime.h"

void
getclosest(double x, double y, double z, long *best, long *size, double *rx, double *ry, double *rz){

	long i, AA=0, BB=0;
	double dist=1e30;
	
	*best = -1;
	
	if (x>=0 && y>=0 && z>=0)	{ AA=0; BB=size[0]; }
	if (x>=0 && y>=0 && z< 0)	{ AA=size[0]; BB=size[1]; }
	if (x>=0 && y< 0 && z>=0)	{ AA=size[1]; BB=size[2]; }
	if (x>=0 && y< 0 && z< 0)	{ AA=size[2]; BB=size[3]; }
	if (x< 0 && y>=0 && z>=0)	{ AA=size[3]; BB=size[4]; }
	if (x< 0 && y>=0 && z< 0)	{ AA=size[4]; BB=size[5]; }
	if (x< 0 && y< 0 && z>=0)	{ AA=size[5]; BB=size[6]; }
	if (x< 0 && y< 0 && z< 0)	{ AA=size[6]; BB=size[7]; }
	
	for(i=AA;i<BB;i++){
		if(sqrt(pow(rx[i]-x,2)+pow(ry[i]-y,2)+pow(rz[i]-z,2)) <= dist){
			dist=sqrt(pow(rx[i]-x,2)+pow(ry[i]-y,2)+pow(rz[i]-z,2));
			*best=i;
		}
	}

	if (*best == -1) {
	  if(!silent) bail_out("Error in getclosest routine");
	  exit(1);
	}

}	
