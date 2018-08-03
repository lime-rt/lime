/*
 *  lime_config.h
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#ifndef LIME_CONFIG_H
#define LIME_CONFIG_H

#include "dims.h"

#define MAX_NSPECIES            100
#define MAX_NIMAGES             100
#define NUM_GRID_STAGES         5
#define MAX_N_COLL_PART         20
#define DENSITY_POWER           0.2
#define MAX_N_HIGH              10

typedef struct {
  /* Elements also present in struct inpars: */
  double radius,minScale,tcmb,*nMolWeights;
  double (*gridDensMaxLoc)[DIM],*gridDensMaxValues,*collPartMolWeights;
  int sinkPoints,pIntensity,blend,*collPartIds,traceRayAlgorithm,samplingAlgorithm;
  int sampling,lte_only,init_lte,antialias,polarization,nThreads,nSolveIters;
  int collPartUserSetFlags;
  char **girdatfile,**moldatfile,**collPartNames;
  char *outputfile,*binoutputfile,*gridfile,*pregrid,*restart,*dust;
  char *gridInFile,**gridOutFiles;
  _Bool resetRNG,doSolveRTE;

  /* New elements: */
  double radiusSqu,minScaleSqu,taylorCutoff,gridDensGlobalMax;
  int ncell,nImages,nSpecies,numDensities,doPregrid,numGridDensMaxima,numDims;
  int nLineImages,nContImages,dataFlags,nSolveItersDone;
  _Bool doInterpolateVels,useAbun,doMolCalcs;
  _Bool writeGridAtStage[NUM_GRID_STAGES],useVelFuncInRaytrace,edgeVelsAvailable;
  _Bool needToInitPops,needToInitSND,SNDhasBeenInit,popsHasBeenInit;
} configInfo;

struct spec {
  double *intense;
  double *tau;
  double stokes[3];
  int numRays;
};

/* Image information */
typedef struct {
  /* Elements also present in struct image: */
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

  /* New elements: */
  struct spec *pixel;
  int *imgunits;
  int numunits;
  double rotMat[3][3];
  int doline;
} imageInfo;

#endif /* LIME_CONFIG_H */

