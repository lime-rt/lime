/*
 *  defaults.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#include <math.h>
#include "defaults.h" /* includes lime_config.h which defines configInfo */
#include "ufunc_types.h" /* for the USERFUNC_* macros */

void
default_density(double x, double y, double z, double *density){
  density[0] = 0.0;
  defaultFuncFlags |= (1 << USERFUNC_density);
}

void
default_temperature(double x, double y, double z, double *temperature){
  temperature[0] = 0.0;
  temperature[1] = 0.0;
  defaultFuncFlags |= (1 << USERFUNC_temperature);
}

void
default_abundance(double x, double y, double z, double *abundance){
  abundance[0] = -1.0;
  defaultFuncFlags |= (1 << USERFUNC_abundance);
}

void
default_molNumDensity(double x, double y, double z, double *dummy){
  dummy[0] = -1.0;
  defaultFuncFlags |= (1 << USERFUNC_molNumDensity);
}

void
default_doppler(double x, double y, double z, double *doppler){
  *doppler = 0.0;
  defaultFuncFlags |= (1 << USERFUNC_doppler);
}

void
default_velocity(double x, double y, double z, double *vel){
  vel[0] = 0.0;
  vel[1] = 0.0;
  vel[2] = 0.0;
//*** probably not good to hard-wire DIM to 3 in this way.
  defaultFuncFlags |= (1 << USERFUNC_velocity);
}

void
default_magfield(double x, double y, double z, double *B){
  B[0]=0.0;
  B[1]=0.0;
  B[2]=0.0;
  defaultFuncFlags |= (1 << USERFUNC_magfield);
}

void
default_gasIIdust(double x, double y, double z, double *gas2dust){
  *gas2dust=100.;
  defaultFuncFlags |= (1 << USERFUNC_gasIIdust);
}

double
default_gridDensity(configInfo *par, double *r, void (*density)(double x, double y, double z, double *val)){
  /*
The grid points within the model are chosen randomly via the rejection method with a probability distribution which the present function is intended to provide.

Notes:
  - The present function is interpreted by LIME as giving *relative* probabilities, the ultimate normalization being set by the desired number of grid points conveyed to the task via par->pIntensity.
  - If par->samplingAlgorithm is chosen to be zero (the current default value), further manipulations to the probability distribution are performed according to the set value of par->sampling.
  - The user may supply their own version of the present function within model.c; the default here implements the grid-point selection function used in LIME<=1.5.
  */
  double val[99],totalDensity=0.0,rSquared=0.0,fracDensity=0.0;
  int i;

  rSquared = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
  if(rSquared>=par->radiusSqu)
    return 0.0;

  density(r[0],r[1],r[2],val);
  for (i=0;i<par->numDensities;i++) totalDensity += val[i];
  fracDensity = pow(totalDensity,defaultDensyPower)/par->gridDensGlobalMax;

  defaultFuncFlags |= (1 << USERFUNC_gridDensity);

  return fracDensity;
}

