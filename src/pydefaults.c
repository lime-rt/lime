/*
 *  defaults.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2016 The LIME development team
 *
 */

#include "lime.h"


void
    density_default(double x, double y, double z, double *density){
      if(!silent) bail_out("Density is not defined in model.c but is needed by LIME!");
      exit(1);
    }

void
    temperature_default(double x, double y, double z, double *temperature){
      if(!silent) bail_out("Temperature is not defined in model.c but is needed by LIME!");
      exit(1);
    }

void
    abundance_default(double x, double y, double z, double *abundance){
      if(!silent) bail_out("Abundance is not defined in model.c but is needed by LIME!");
      exit(1);
    }

void
    doppler_default(double x, double y, double z, double *doppler){
      if(!silent) bail_out("Doppler velocity is not defined in model.c but is needed by LIME!");
      exit(1);
    }

void
    velocity_default(double x, double y, double z, double *vel){
      if(!silent) bail_out("Velocity field is not defined in model.c but is needed by LIME!");
      exit(1);
    }

void
    magfield_default(double x, double y, double z, double *B){
      B[0]=0.0;
      B[1]=0.0;
      B[2]=0.0;
    }

void
    gasIIdust_default(double x, double y, double z, double *gas2dust){
      *gas2dust=100.;
    }

double
    gridDensity_default(configInfo *par, double *r){
      /* The grid points within the model are chosen randomly via the rejection method with a probability distribution which the present function is intended to provide. The user may supply their own version of this within model.c; the default here implements the grid-point selection function used in LIME<=1.5.
      */
      double val[99],totalDensity=0.0,rSquared=0.0,fracDensity=0.0;
      int i;
      static double maxTotalDensity=0.0;

      rSquared = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
      if(rSquared>=par->radiusSqu)
        return 0.0;

      if(maxTotalDensity==0.0){
        density(par->minScale,par->minScale,par->minScale,val);
        for (i=0;i<par->numDensities;i++)
          maxTotalDensity += val[i];

        if (maxTotalDensity<=0.){
          if(!silent) bail_out("Error: Sum of reference densities equals 0");
          exit(1);
        }
      }

      density(r[0],r[1],r[2],val);
      for (i=0;i<par->numDensities;i++) totalDensity += val[i];
      fracDensity = pow(totalDensity/maxTotalDensity,0.2);

      return fracDensity;
    }

