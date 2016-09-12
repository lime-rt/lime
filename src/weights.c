/*
 *  weights.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015 The LIME development team
 *
 */

#include "lime.h"

int pointEvaluation(configInfo *par, const double uniformRandom, double *r, double *nref){
  double fracDensity;

  fracDensity = gridDensity(par, r, nref);

  if(uniformRandom < fracDensity) return 1;
  else return 0;
}

