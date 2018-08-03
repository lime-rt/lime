/*
 *  lime_defaults.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#include "lime.h" /* includes defaults.h */

void __attribute__((weak))
density(double x, double y, double z, double *density){
  default_density(x,y,z, density);
}

void __attribute__((weak))
temperature(double x, double y, double z, double *temperature){
  default_temperature(x,y,z, temperature);
}

/* One of the following two must be defined by the user: */
void __attribute__((weak))
abundance(double x, double y, double z, double *abun){
  default_abundance(x,y,z, abun);
}

void __attribute__((weak))
molNumDensity(double x, double y, double z, double *nDensity){
  default_molNumDensity(x,y,z, nDensity);
}

void __attribute__((weak))
doppler(double x, double y, double z, double *doppler){
  default_doppler(x,y,z, doppler);
}

void __attribute__((weak))
velocity(double x, double y, double z, double *vel){
  default_velocity(x,y,z, vel);
}

void __attribute__((weak))
magfield(double x, double y, double z, double *B){
  default_magfield(x,y,z, B);
}

void __attribute__((weak))
gasIIdust(double x, double y, double z, double *gas2dust){
  default_gasIIdust(x,y,z, gas2dust);
}

double __attribute__((weak))
gridDensity(configInfo *par, double *r){
  return default_gridDensity(par, r, density);
}

