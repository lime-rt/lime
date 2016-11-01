#ifndef TREERAND_H
#define TREERAND_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_qrng.h>

#define TREE_N_RANDOMS          1000	/* Arbitrary - experiment around a bit to find a good value. */
#define TREE_MAX_RECURSION      20	/* >20 is not safe with single-precision arithmetic. */
#define TREE_MAX_N_TRIALS       1000	/* Arbitrary - experiment to find a good value. */
#define TREE_DITHER             0.1	/* must be >=0, <1. */

/*
Note that these structs and functions require the following definitions:
  configInfo
  DIM (the number of spatial dimensions)
*/

typedef struct {
  configInfo par;
  gsl_rng *randGen;
  int numSubFields,numInRandoms,verbosity;
  double (*inRandLocations)[DIM];
  int maxRecursion,maxNumTrials;
  double dither,maxNumTrialsDbl,overSamplingRatio,abstandFrac;
  int leafBufLenI,inRandBufLenI;
  unsigned int desiredNumPoints;
  _Bool doShuffle,doQuasiRandom;
} treeRandConstType;

typedef struct{
  double maxDensity,densityIntegral,fieldOrigin[DIM],fieldWidth[DIM];
} leafType;

typedef struct{
  leafType *leaves;
  int lastLeafI,maxLeafI;
} treeType;

typedef struct {
  double fieldOrigin[DIM],fieldWidth[DIM],fieldVolume;
  double expectedDesNumPoints;
  int numHighPoints;
  double sumDensity,maxDensity,densityIntegral;
  double (*highPointLocations)[DIM];
  double *highPointDensities;
} treeRandVarType;

void	calcSubCellValues(treeRandConstType*, treeRandVarType*, double (*numberDensyFunc)(configInfo*, double*));
void	initializeTree(treeRandConstType*, treeRandVarType*, double (*numberDensyFunc)(configInfo*, double*), treeType*);
void	constructTheTree(treeRandConstType*, treeRandVarType*, int, double (*numberDensyFunc)(configInfo*, double*), treeType*);
void	fillTheTree(treeRandConstType*, treeType*, double (*numberDensyFunc)(configInfo*, double*), double(*)[DIM], double*);
void	freeRinc(treeRandConstType);
void	freeRinv(treeRandVarType);

#endif /* TREERAND_H */


