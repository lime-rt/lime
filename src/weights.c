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

#ifdef CAVITY_WALLS
  density(par->minScale,par->minScale,par->minScale,-1.0,val);
#else
  density(par->minScale,par->minScale,par->minScale,val);
#endif
  for (i=0;i<par->collPart;i++) normalizer += val[i];
  if (normalizer<=0.){
    if(!silent) bail_out("Error: Sum of reference densities equals 0");
    exit(1);
  }
  //abundance(par->minScale,par->minScale,par->minScale,val2);
#ifdef CAVITY_WALLS
  density(x,y,z,-1.0,val);
#else
  density(x,y,z,val);
#endif
  for (i=0;i<par->collPart;i++) totalDensity += val[i];
  //abundance(x,y,z,val2);
  weight1=pow(totalDensity/normalizer,0.2);

  weight2=0.;

  if(ran < weight1 || ran < weight2) return 1;
  else return 0;
}

int
pointEvaluationW(inputPars *par,double ran, double x, double y, double z){
  int flag;
  double weight1,val[9];

  weight1 = angletocavity(x,y,z);
  weight1 = sqrt(weight1);

  density(x,y,z,-1.0,val);
  if(val[0]+val[1]<1e0) weight1=0.;

  flag=0;
  if(ran<weight1) flag=1;

  return flag;
}
