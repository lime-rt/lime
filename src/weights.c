/*
 *  weights.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015 The LIME development team
 *
 */

#include "lime.h"

int
pointEvaluation(inputPars *par,double ran, double x, double y, double z){
  double weight1, weight2, val[99],normalizer=0.0,totalDensity=0.0;
  int i;

  density(par->minScale,par->minScale,par->minScale,val);
  for (i=0;i<par->numDensities;i++) normalizer += val[i];
  if (normalizer<=0.){
    if(!silent) bail_out("Sum of reference densities equals 0");
    exit(1);
  }
  //abundance(par->minScale,par->minScale,par->minScale,val2);
  density(x,y,z,val);
  for (i=0;i<par->numDensities;i++) totalDensity += val[i];
  //abundance(x,y,z,val2);
  weight1=pow(totalDensity/normalizer,0.2);

  weight2=0.;

  if(ran < weight1 || ran < weight2) return 1;
  else return 0;
}
