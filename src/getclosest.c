/*
 *  getclosest.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2017 The LIME development team
 *
 */

#include "lime.h"

void
getclosest(double x, double y, double z, long *best, long *size, double *rx, double *ry, double *rz){

	long i, AA=0, BB=0;
	double distSquared, trialD2;
	
	if (x>=0 && y>=0 && z>=0)	{ AA=0; BB=size[0]; }
	if (x>=0 && y>=0 && z< 0)	{ AA=size[0]; BB=size[1]; }
	if (x>=0 && y< 0 && z>=0)	{ AA=size[1]; BB=size[2]; }
	if (x>=0 && y< 0 && z< 0)	{ AA=size[2]; BB=size[3]; }
	if (x< 0 && y>=0 && z>=0)	{ AA=size[3]; BB=size[4]; }
	if (x< 0 && y>=0 && z< 0)	{ AA=size[4]; BB=size[5]; }
	if (x< 0 && y< 0 && z>=0)	{ AA=size[5]; BB=size[6]; }
	if (x< 0 && y< 0 && z< 0)	{ AA=size[6]; BB=size[7]; }
	
	i=AA;
        distSquared = (rx[i]-x)*(rx[i]-x)+(ry[i]-y)*(ry[i]-y)+(rz[i]-z)*(rz[i]-z);
	*best = i;
	for(i=AA+1;i<BB;i++){
                trialD2 = (rx[i]-x)*(rx[i]-x)+(ry[i]-y)*(ry[i]-y)+(rz[i]-z)*(rz[i]-z);
		if(trialD2 <= distSquared){
			distSquared=trialD2;
			*best=i;
		}
	}
}	
