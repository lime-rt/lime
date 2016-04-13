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
pointEvaluation(inputPars *par, double ran, double x, double y, double z){
  double fracDensity;

  gridDensity(par, x, y, z, &fracDensity);

  if(ran < fracDensity) return 1;
  else return 0;
}

