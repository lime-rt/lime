/*
 *  model.c
 *  LIME, The versatile 3D line modeling tool
 *
 *  Created by Christian Brinch on 11/05/07.
 *  Copyright 2006-2013, Christian Brinch,
 *  <brinch@nbi.dk>
 *  Niels Bohr institutet
 *  University of Copenhagen
 *	All rights reserved.
 *
 */

#include "lime.h"

/******************************************************************************/

void
input(inputPars *par, image *img){
  /*
   * Basic parameters. See cheat sheet for details.
   */
  par->radius			= 2000*AU;
  par->minScale	   		= 0.5*AU;
  par->pIntensity    	= 4000;
  par->sinkPoints    	= 3000;
  par->dust				= "jena_thin_e6.tab";
  par->moldatfile[0] 	= "hco+@xpol.dat";
  par->antialias		= 4;
  par->sampling			= 2; // log distr. for radius, directions distr. uniformly on a sphere.

  par->outputfile 		= "populations.pop";
  par->binoutputfile 	= "restart.pop";
  par->gridfile			= "grid.vtk";

  /*
   * Definitions for image #0. Add blocks for additional images.
   */
  img[0].nchan			= 61;		  // Number of channels
  img[0].velres			= 500.;       // Channel resolution in m/s
  img[0].trans			= 3;          // zero-indexed J quantum number
  img[0].pxls			= 100;	      // Pixels per dimension
  img[0].imgres			= 0.1;		  // Resolution in arc seconds
  img[0].theta			= 0.0;		  // 0: face-on, pi/2: edge-on
  img[0].distance		= 140*PC;	  // source distance in m
  img[0].source_vel		= 0;          // source velocity in m/s
  img[0].unit			= 0;		  // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau
  img[0].filename		= "image0.fits";	// Output filename
}

/******************************************************************************/

void
density(double x, double y, double z, double *density){
  /*
   * Define variable for radial coordinate
   */
  double r;
  /*
   * Calculate radial distance from origin
   */
  r=sqrt(x*x+y*y+z*z);
  /*
   * Calculate a spherical power-law density profile
   * (Multiply with 1e6 to go to SI-units)
   */
  density[0] = 1.5e6*pow(r/(300*AU),-1.5)*1e6;
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
   * Calculate coordinate distance from origin
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
  double R, r;

  R=sqrt(x*x+y*y+z*z);
  if (R>0.){
/*
 * Free-fall velocity in the radial direction onto a central 
 * mass of 1.0 solar mass
 */  
    r=-sqrt(2*6.67e-11*1.989e30/R);

    vel[0]=x*r/R;
    vel[1]=y*r/R;
    vel[2]=z*r/R;
  } else {
    vel[0]=0.0;
    vel[1]=0.0;
    vel[2]=0.0;
  }
}

/******************************************************************************/


