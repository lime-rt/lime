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
  double weight1, weight2, val[99],n0;

  density(par->minScale,par->minScale,par->minScale,val);
  //abundance(par->minScale,par->minScale,par->minScale,val2);
  n0=val[0];
  density(x,y,z,val);
  //abundance(x,y,z,val2);
  weight1=pow(val[0]/n0,0.2);


  weight2=0.;

  if(ran < weight1 || ran < weight2) return 1;
  else return 0;
}
