/*
 *  model1a.c -- Benchmark model for LIME
 *
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2016 The LIME development team
 *
 */

#include "lime.h"

void
input (inputPars * par, image * img)
{
  par->radius = 7.8e16;
  par->minScale = 1e13;
  par->pIntensity = 4000;
  par->sinkPoints = 3000;
  par->dust = NULL;
  par->moldatfile[0] = "molec.dat";
  par->antialias = 1;
  par->sampling = 2;
  par->outputfile = "model1a.pop";

  img[0].nchan = 32;
  img[0].velres = 100.;
  img[0].trans = 0;
  img[0].pxls = 128;
  img[0].imgres = 0.1;
  img[0].theta = 0.0;
  img[0].distance = 140 * PC;
  img[0].source_vel = 0;
  img[0].unit = 0;
  img[0].filename = "model1a.fits";
}

void
density (double x, double y, double z, double *density)
{
  double r;

  r = sqrt (x * x + y * y + z * z);
  if ((r >= 1e13) && (r <= 7.8e16))
    density[0] = 2e7 * pow (r / 1e13, -2) * 1e6;
  else
    density[0] = 0.;
}

void
temperature (double x, double y, double z, double *temperature)
{
  temperature[0] = 20.;
}

void
abundance (double x, double y, double z, double *abundance)
{
  abundance[0] = 1e-8;
}

void
doppler (double x, double y, double z, double *doppler)
{
  /*
     The width a=0.15 km/s given by Zadelhoff et al. is a total
     thermal+turbulent linewidth, so we need to remove the thermal
     component.
  */
  doppler[0] = sqrt(150 * 150 - 2 * KBOLTZ / AMU);
}

void
velocity (double x, double y, double z, double *vel)
{
  vel[0] = 0;
  vel[1] = 0;
  vel[2] = 0;
}
