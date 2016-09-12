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

double __attribute__((weak))
    gridDensity(configInfo *par, double *r, double *nref){
      /* The grid points within the model are chosen randomly via the rejection method with a probability distribution which the present function is intended to provide. The user may supply their own version of this within model.c; the default here implements the grid-point selection function used in LIME<=1.5.
      */
      double val[99],totalDensity=0.0,rSquared=0.0,fracDensity=0.0;
      int i;
      static double referenceDensity=0.0;

      rSquared = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
      if(rSquared>=par->radiusSqu)
        return 0.0;

      if(referenceDensity==0.0){
          density(par->minScale,par->minScale,par->minScale,val);
          for (i=0;i<par->numDensities;i++)
              if(nref[i] != -1){
                  referenceDensity += nref[i];
              }
              else{
                  referenceDensity += val[i];
              }

          if (referenceDensity<=0.){
              if(!silent) bail_out("Error: Sum of reference densities equals 0");
              exit(1);
          }
      }

      density(r[0],r[1],r[2],val);
      for (i=0;i<par->numDensities;i++) totalDensity += val[i];
      fracDensity = pow(totalDensity/referenceDensity,0.2);

      return fracDensity;
    }

