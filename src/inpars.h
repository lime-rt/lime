/*
 *  inpars.h
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#ifndef INPARS_H
#define INPARS_H

#include "dims.h"

/* input parameters */
typedef struct {
  double radius,minScale,tcmb,*nMolWeights,*dustWeights;
  double (*gridDensMaxLoc)[DIM],*gridDensMaxValues,*collPartMolWeights;
  int sinkPoints,pIntensity,blend,*collPartIds,traceRayAlgorithm,samplingAlgorithm;
  int sampling,lte_only,init_lte,antialias,polarization,nThreads,nSolveIters;
  char **girdatfile,**moldatfile,**collPartNames;
  char *outputfile,*binoutputfile,*gridfile,*pregrid,*restart,*dust;
  char *gridInFile,**gridOutFiles;
  _Bool resetRNG,doSolveRTE;
} inputPars;

/* Image information */
typedef struct {
  int nchan,trans,molI;
  double velres;
  double imgres;
  int pxls;
  int unit;
  char *units;
  double freq,bandwidth;
  char *filename;
  double source_vel;
  double theta,phi,incl,posang,azimuth;
  double distance;
  _Bool doInterpolateVels;
} image;

#endif /* INPARS_H */
