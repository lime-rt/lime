/*
 *  lime_config.h
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2017 The LIME development team
 *
 */

#ifndef LIME_CONFIG_H
#define LIME_CONFIG_H

#include "dims.h"

#define SQRT_PI                 (sqrt(M_PI))           /* sqrt(pi)	*/
#define ARCSEC_TO_RAD           (M_PI/180.0/3600.0)
#define EPS                     1.0e-30                /* general use small number */

#define MAX_NSPECIES            100
#define MAX_NIMAGES             100
#define NUM_GRID_STAGES         5
#define MAX_N_COLL_PART         20
#define TYPICAL_ISM_DENS        1000.0
#define STR_LEN_0               80
#define DENSITY_POWER           0.2
#define MAX_N_HIGH              10


extern int silent;

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
} configInfo;

#endif /* LIME_CONFIG_H */

