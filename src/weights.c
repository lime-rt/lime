/*
 *  weights.c
 *  LIME, The versatile 3D line modeling environment 
 *
 *  Created by Christian Brinch on 11/16/06.
 *  Copyright 2006-2012, Christian Brinch, 
 *  <brinch@nbi.dk>
 *  Niels Bohr institutet
 *  University of Copenhagen
 *	All rights reserved.
 *
 */

#include "lime.h"

int
pointEvaluation(inputPars *par,double ran, double x, double y, double z){
  double weight1, weight2, val[99],n0,val2[99];

  density(par->minScale,par->minScale,par->minScale,val);
  //abundance(par->minScale,par->minScale,par->minScale,val2);
  n0=val[0];
  density(x,y,z,val);
  //abundance(x,y,z,val2);
  weight1=pow(val[0]/n0,0.5);


  weight2=0.;

  if(ran < weight1 || ran < weight2) return 1;
  else return 0;
}
