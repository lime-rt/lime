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
    }

void __attribute__((weak))
    gasIIdust(double x, double y, double z, double *gas2dust){
      *gas2dust=100.;
    }

void __attribute__((weak))
    gridDensity(inputPars par, double x, double y, double z, double *fracDensity){
      /* The grid points within the model are chosen randomly via the rejection method with a probability distribution which the present function is intended to provide. The user may supply their own version of this within model.c; the default here implements the grid-point selection function used in LIME<=1.5.
      */
      double val[99],totalDensity=0.0;
      int i;
      static double maxTotalDensity=0.0;

      if(maxTotalDensity==0.0){
        density(par.minScale,par.minScale,par.minScale,val);
        for (i=0;i<par.collPart;i++)
          maxTotalDensity += val[i];

        if (maxTotalDensity<=0.){
          if(!silent) bail_out("Error: Sum of reference densities equals 0");
          exit(1);
        }
      }

      density(x,y,z,val);
      for (i=0;i<par.collPart;i++) totalDensity += val[i];
      *fracDensity = pow(totalDensity/maxTotalDensity,0.2);
    }

void __attribute__((weak))
    gridDensity_old(inputPars *par, double x, double y, double z, double *fracDensity){
      /* The grid points within the model are chosen randomly via the rejection method with a probability distribution which the present function is intended to provide. The user may supply their own version of this within model.c; the default here implements the grid-point selection function used in LIME<=1.5.
      */
      double val[99],totalDensity=0.0;
      int i;
      static double maxTotalDensity=0.0;

      if(maxTotalDensity==0.0){
        density(par->minScale,par->minScale,par->minScale,val);
        for (i=0;i<par->collPart;i++)
          maxTotalDensity += val[i];

        if (maxTotalDensity<=0.){
          if(!silent) bail_out("Error: Sum of reference densities equals 0");
          exit(1);
        }
      }

      density(x,y,z,val);
      for (i=0;i<par->collPart;i++) totalDensity += val[i];
      *fracDensity = pow(totalDensity/maxTotalDensity,0.2);
    }

