/*
 *  ufunc_types.h
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#ifndef UFUNC_TYPES_H
#define UFUNC_TYPES_H

#include "lime_config.h" /* for configInfo */

#define USERFUNC_density       0
#define USERFUNC_temperature   1
#define USERFUNC_abundance     2
#define USERFUNC_molNumDensity 3
#define USERFUNC_doppler       4
#define USERFUNC_velocity      5
#define USERFUNC_magfield      6
#define USERFUNC_gasIIdust     7
#define USERFUNC_gridDensity   8

void	density(double,double,double,double *);
void	temperature(double,double,double,double *);
void	abundance(double,double,double,double *);
void	molNumDensity(double,double,double,double *);
void	doppler(double,double,double, double *);
void	velocity(double,double,double,double *);
void	magfield(double,double,double,double *);
void	gasIIdust(double,double,double,double *);
double	gridDensity(configInfo*, double*);

#endif /* UFUNC_TYPES_H */

