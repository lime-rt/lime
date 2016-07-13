/*
 *  defaults.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015 The LIME development team
 *
 */

#include "lime.h"
#include "assert.h"


void __attribute__((weak))
    density(double x, double y, double z, double *density){
      if(!silent) bail_out("Density is not defined in model.c but is needed by LIME!");
      exit(1);
    }

void __attribute__((weak))
    temperature(double x, double y, double z, double *temperature){
      if(!silent) bail_out("Temperature is not defined in model.c but is needed by LIME!");
      exit(1);
    }

void __attribute__((weak))
    abundance(double x, double y, double z, double *abundance){
      if(!silent) bail_out("Abundance is not defined in model.c but is needed by LIME!");
      exit(1);
    }

void __attribute__((weak))
    doppler(double x, double y, double z, double *doppler){
      if(!silent) bail_out("Doppler velocity is not defined in model.c but is needed by LIME!");
      exit(1);
    }

void __attribute__((weak))
    velocity(double x, double y, double z, double *vel){
      if(!silent) bail_out("Velocity field is not defined in model.c but is needed by LIME!");
      exit(1);
    }

void __attribute__((weak))
    magfield(double x, double y, double z, double *B){
      B[0]=0.0;
      B[1]=0.0;
      B[2]=0.0;
    }

void __attribute__((weak))
    gasIIdust(double x, double y, double z, double *gas2dust){
      *gas2dust=100.;
    }
