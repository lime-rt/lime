/*
 *  tree_random.h
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#ifndef TREERAND_H
#define TREERAND_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_qrng.h>

#include "dims.h" /* Should define the macro N_DIMS which specifies the number of spatial dimensions. Note that this is a _maximum_ value, the actual value 'numDims' which is passed via the function interfaces may be <= N_DIMS. */


#define TREE_N_RANDOMS          10000	/* Arbitrary - experiment around a bit to find a good value. */
#define TREE_MAX_RECURSION      20	/* >20 is not safe with single-precision arithmetic. */
#define TREE_MAX_N_TRIALS       1000	/* Arbitrary - experiment to find a good value. */
#define TREE_DITHER             0.1	/* must be >=0, <1. */
#define TREE_leafBufLenI	10000
#define TREE_inRandBufLenI	10000
#define TREE_abstandFrac	0.1	/* Buffer to leave at the cell borders for quasi-random points, expressed as a fraction of the mean distance between points. Only used if treeRandConstType.doQuasiRandom==1. */

#define TREE_MIN_ARA_RANGE	0.25
#define TREE_MAX_N_EXTRAS	10000
#define TREE_STRLEN		80

#define ERR_NO_GOOD_POINTS	1
#define ERR_NO_POINTS_DESIRED	2
#define ERR_TOO_FEW_POINTS	3

#define TREE_MSG_MESSAGE	0
#define TREE_MSG_WARN		1
#define TREE_MSG_ERROR		2

typedef struct {
  configInfo par;
  gsl_rng_type *randGenType;
  unsigned long int randSeed;
  gsl_qrng_type *quasiRandGenType;
  int numDims,numInRandoms,verbosity,totalNumHighPoints;
  int maxRecursion,maxNumTrials,leafBufLenI,inRandBufLenI;
  double abstandFrac,dither,wholeFieldOrigin[N_DIMS],wholeFieldWidth[N_DIMS];
  double (*allHighPointLoc)[N_DIMS];
  double *allHighPointDensy;
  unsigned int desiredNumPoints;
  _Bool doShuffle,doQuasiRandom;
  void (*monitorFunc)(const int numDims, const int cellI, double fieldOrigin[N_DIMS]\
    , double fieldWidth[N_DIMS], unsigned int desiredNumPoints, double (*outRandLocations)[N_DIMS]\
    , unsigned int firstPointI, unsigned int actualNumPoints);
} treeRandConstType; /* Fields of this struct should be set before generateRandoms() is called. */ 

typedef struct {
  int numSubFields;
  double maxNumTrialsDbl,(*inRandLocations)[N_DIMS];
  gsl_rng *randGen; /* Random number generator - should be the value returned by gsl_rng_alloc(). */
  gsl_qrng *quasiRandGen; /* Quasi-random number generator - should be the value returned by gsl_qrng_alloc(). */
} treeRandInternalType; /* Fields of this struct are constant but are set internally. */

typedef struct {
  int numHighPoints;
  int axisIndices[N_DIMS];
  double fieldOrigin[N_DIMS],fieldWidth[N_DIMS],axisSigns[N_DIMS],absRanAcceptRange[N_DIMS];
  double expectedDesNumPoints,sumDensity,maxDensity,densityIntegral;
  double (*highPointLocations)[N_DIMS];
  double *highPointDensities;
} subCellType;

typedef struct{
  subCellType *leaves;
  int lastLeafI,maxLeafI;
} treeType;

void setConstDefaults(treeRandConstType*);
void treeGenerateRandoms(treeRandConstType*, double (*numberDensyFunc)(configInfo*, double*), double (*)[], double*);
void freeRinc(treeRandConstType);

#endif /* TREERAND_H */


