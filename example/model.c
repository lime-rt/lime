/*
 *  model.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015 The LIME development team
 *
 */

#include "lime.h"

/******************************************************************************/

void
input(inputPars *par, image *img){
  int i;

  /*
   * Basic parameters. See cheat sheet for details.
   */
  par->radius                   = 2000*AU;
  par->minScale                 = 0.5*AU;
  par->pIntensity               = 4000;
  par->sinkPoints               = 3000;
  par->dust                     = "jena_thin_e6.tab";
  par->moldatfile[0]            = "hco+@xpol.dat";
  par->antialias                = 4;
  par->sampling                 = 2; // log distr. for radius, directions distr. uniformly on a sphere.

  par->outputfile               = "populations.pop";
  par->binoutputfile            = "restart.pop";
  par->gridfile                 = "grid.vtk";

  /*
   * Definitions for image #0. Add blocks with successive values of i for additional images.
   */
  i=0;
  img[i].nchan                  = 61;             // Number of channels
  img[i].velres                 = 500.;           // Channel resolution in m/s
  img[i].trans                  = 3;              // zero-indexed J quantum number
  img[i].pxls                   = 100;            // Pixels per dimension
  img[i].imgres                 = 0.1;            // Resolution in arc seconds
  img[i].theta                  = 0.0;            // 0: face-on, pi/2: edge-on
  img[i].distance               = 140*PC;         // source distance in m
  img[i].source_vel             = 0;              // source velocity in m/s
  img[i].unit                   = 0;              // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau
  img[i].filename               = "image0.fits";  // Output filename
}

/******************************************************************************/

void
density(double x, double y, double z, double *density){
  /*
   * Define variable for radial coordinate
   */
  double r, rToUse;
  const double rMin = 0.1*AU; /* This cutoff should be chosen smaller than par->minScale but greater than zero (to avoid a singularity at the origin). */

  /*
   * Calculate radial distance from origin
   */
  r=sqrt(x*x+y*y+z*z);
  /*
   * Calculate a spherical power-law density profile
   * (Multiply with 1e6 to go to SI-units)
   */
  if(r>rMin)
    rToUse = r;
  else
    rToUse = rMin;

  density[0] = 1.5e6*pow(rToUse/(300*AU),-1.5)*1e6;
}

/******************************************************************************/

void
temperature(double x, double y, double z, double *temperature){
  int i,x0=0;
  double r;
  /*
   * Array containing temperatures as a function of radial
   * distance from origin (this is an example of a tabulated model)
   */
  double temp[2][10] = {
      {2.0e13, 5.0e13, 8.0e13, 1.1e14, 1.4e14, 1.7e14, 2.0e14, 2.3e14, 2.6e14, 2.9e14},
      {44.777, 31.037, 25.718, 22.642, 20.560, 19.023, 17.826, 16.857, 16.050, 15.364}
  };
  /*
   * Calculate radial distance from origin
   */
  r=sqrt(x*x+y*y+z*z);
  /*
   * Linear interpolation in temperature input
   */
  if(r > temp[0][0] && r<temp[0][9]){
    for(i=0;i<9;i++){
      if(r>temp[0][i] && r<temp[0][i+1]) x0=i;
    }
  }
  if(r<temp[0][0])
    temperature[0]=temp[1][0];
  else if (r>temp[0][9])
    temperature[0]=temp[1][9];
  else
    temperature[0]=temp[1][x0]+(r-temp[0][x0])*(temp[1][x0+1]-temp[1][x0])/(temp[0][x0+1]-temp[0][x0]);
}

/******************************************************************************/

void
abundance(double x, double y, double z, double *abundance){
  /*
   * Here we use a constant abundance. Could be a
   * function of (x,y,z).
   */
  abundance[0] = 1.e-9;
}

/******************************************************************************/

void
doppler(double x, double y, double z, double *doppler){
  /*
   * 200 m/s as the doppler b-parameter. This
   * can be a function of (x,y,z) as well.
   * Note that *doppler is a pointer, not an array.
   * Remember the * in front of doppler.
   */
  *doppler = 200.;
}

/******************************************************************************/

void
velocity(double x, double y, double z, double *vel){
  double r, rToUse, ffSpeed;
  const double rMin = 0.1*AU; /* This cutoff should be chosen smaller than par->minScale but greater than zero (to avoid a singularity at the origin). */

  /*
   * Calculate radial distance from origin
   */
  r = sqrt(x*x+y*y+z*z);
  if(r>rMin)
    rToUse = r;
  else
    rToUse = rMin;

  /*
   * Free-fall velocity in the radial direction onto a central 
   * mass of 1.0 solar mass
   */  
  ffSpeed = sqrt(2*GRAV*1.989e30/rToUse);

  vel[0] = -x*ffSpeed/rToUse;
  vel[1] = -y*ffSpeed/rToUse;
  vel[2] = -z*ffSpeed/rToUse;
}

/******************************************************************************/


