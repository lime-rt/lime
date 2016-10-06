/*
 *  weights.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2016 The LIME development team
 *
 */

#include "lime.h"

int pointEvaluation(configInfo *par, const double uniformRandom, double *r, double *weighting_n0, double weighting_w){
  double fracDensity;

  fracDensity = gridDensity(par, r, weighting_n0, weighting_w);

  if(uniformRandom < fracDensity) return 1;
  else return 0;
}

