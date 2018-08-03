/*
 *  defaults.h
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
Files which include this must
  - define a struct named configInfo which has members as follows:
      double radiusSqu;
      double gridDensGlobalMax;
      int numDensities;

  - define macros USERFUNC_* as used in defaults.c.
 */

#ifndef DEFAULTS_H
#define DEFAULTS_H

#include "lime_config.h" /* for configInfo. */

#define DENSITY_POWER           0.2

extern int defaultFuncFlags;
extern double defaultDensyPower;

void	default_density(double x, double y, double z, double *density);
void	default_temperature(double x, double y, double z, double *temperature);
void	default_abundance(double x, double y, double z, double *abundance);
void	default_molNumDensity(double x, double y, double z, double *dummy);
void	default_doppler(double x, double y, double z, double *doppler);
void	default_velocity(double x, double y, double z, double *vel);
void	default_magfield(double x, double y, double z, double *B);
void	default_gasIIdust(double x, double y, double z, double *gas2dust);
double	default_gridDensity(configInfo *par, double *r, void (*density)(double x, double y, double z, double *val));

#endif /* DEFAULTS_H */

