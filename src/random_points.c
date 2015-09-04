/*
 *  random_points.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2015 The LIME development team
 *
 */

#include "lime.h"

void randomsViaRejection(inputPars *par, unsigned int desiredNumPoints, gsl_rng *randGen\
  , double (*numberDensyFunc)(locusType), locusType *outRandLocations, double *outRandDensities){

  extern double densityNormalizer;
  double lograd; // The logarithm of the model radius.
  double logmin; // Logarithm of par->minScale.
  double r,theta,phi,sinPhi,z,semiradius;
  double uniformRandom, density, acceptanceChance;
  int k,i,di;
  int pointIsAccepted;
  locusType location;

  lograd=log10(par->radius);
  logmin=log10(par->minScale);

  /* Sample pIntensity number of points */
  for(k=0;k<desiredNumPoints;k++){
    uniformRandom=gsl_rng_uniform(randGen);
    pointIsAccepted=0;
    /* Pick a point and check if we like it or not */
    do{
      if(par->sampling==0){
        r=pow(10,logmin+gsl_rng_uniform(randGen)*(lograd-logmin));
        theta=2.*PI*gsl_rng_uniform(randGen);
        phi=PI*gsl_rng_uniform(randGen);
        sinPhi=sin(phi);
        location.x[0]=r*cos(theta)*sinPhi;
        location.x[1]=r*sin(theta)*sinPhi;
        if(DIM==3) location.x[2]=r*cos(phi);
      } else if(par->sampling==1){
        location.x[0]=(2*gsl_rng_uniform(randGen)-1)*par->radius;
        location.x[1]=(2*gsl_rng_uniform(randGen)-1)*par->radius;
        if(DIM==3) location.x[2]=(2*gsl_rng_uniform(randGen)-1)*par->radius;
      } else if(par->sampling==2){
        r=pow(10,logmin+gsl_rng_uniform(randGen)*(lograd-logmin));
        theta=2.*PI*gsl_rng_uniform(randGen);
        if(DIM==3) {
          z=2*gsl_rng_uniform(randGen)-1.;
          semiradius=r*sqrt(1.-z*z);
          z*=r;
          location.x[2]=z;
        } else {
          semiradius=r;
        }
        location.x[0]=semiradius*cos(theta);
        location.x[1]=semiradius*sin(theta);
      } else {
        if(!silent) bail_out("Don't know how to sample model");
        exit(1);
      }
      density = numberDensyFunc(location);
      if(density>0.0){
        acceptanceChance = pow(density/densityNormalizer, DENSITY_POWER);
        if (uniformRandom<acceptanceChance) pointIsAccepted=1;
      }
    } while(!pointIsAccepted);

    for(di=0;di<DIM;di++)
      outRandLocations[k].x[di]=location.x[di];
    outRandDensities[k] = density;

    if(!silent) progressbar((double) k/((double)desiredNumPoints-1), 4);
  }
}

double densityFunc3D(locusType location){
  extern double modelRadiusSquared;
  double rSquared=0.0;
  int i;
  double vals[99],totalDensity=0.0; //**** define MAX_N_COLL_PARTNERS in lime.h and dimension vals to that rather than 99.
  extern int numCollisionPartners;

  for(i=0;i<DIM;i++)
    rSquared += location.x[i]*location.x[i];

  if(rSquared>=modelRadiusSquared)
    return 0.0;

  density(location.x[0],location.x[1],location.x[2],vals);//****** Unsafe to assume DIM==3 like this.
  for (i=0;i<numCollisionPartners;i++) totalDensity += vals[i];
  return totalDensity;
}


