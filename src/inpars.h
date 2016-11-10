/*
 *  inpars.h
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2016 The LIME development team
 *
 */

#ifndef INPARS_H
#define INPARS_H

/* input parameters */
typedef struct {
  double radius,minScale,tcmb,*nMolWeights,*dustWeights;
  int sinkPoints,pIntensity,blend,*collPartIds,traceRayAlgorithm,samplingAlgorithm;
  char *outputfile,*binoutputfile;
//  char *inputfile; unused at present.
  char *gridfile;
  char *pregrid;
  char *restart;
  char *dust;
  int sampling,lte_only,init_lte,antialias,polarization,nThreads;
  char **moldatfile;
  char *gridInFile,**gridOutFiles;
  int nSolveIters;
  double (*gridDensMaxLoc)[DIM], *gridDensMaxValues;
  _Bool resetRNG;
} inputPars;

#endif /* INPARS_H */
