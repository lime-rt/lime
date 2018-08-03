/*
 *  tree_random.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#include "lime.h"
#include "tree_random.h"

/*
The Treegrid algorithm attempts to generate random points within an orthotopic working space which follow a distribution given by some user-supplied number-density function. (An orthotope is a generalization to N dimensions of the concept of a rectangle.) The algorithm is designed to avoid the lengthy run times which can plague the (otherwise robust) rejection method when there are significant variations in the value of the point number density function within the working space. The algorithm makes use of a 2^D tree approach (where D is the number of spatial dimensions) to split the working space into smaller and smaller segments until a space is reached within which the variation of the number density function is within acceptable bounds.

The algorithm is performed in two sections: the first (implemented in function _constructTheTree) decides on the structure of the tree; the second (implemented in function _fillTheTree) decides on the number of points required for each sub-cell, then generates these points via the standard rejection method.

The rejection method takes a starting set of points which have a constant (i.e. flat) probability distribution, then tests each one of these, rejecting points based on the ratio between the value of the desired number density at that point and its maximum value in the relevant subcell. The user has the choice to use either completely random points for the starting set, or points which have a flat probability distribution but which are chosen from a quasi-random sequence, the default being the Halton sequence:

	Halton, J. (1964), Algorithm 247: Radical-inverse quasi-random point sequence, ACM, p. 701

One reason for choosing a quasi-random sequence is to avoid pairs of points which are close together. This is useful if the points are used for the vertices of a finite-element mesh, because it is inefficient to solve these sorts of problems if the FE cells are very divergent in size and shape.

Points chosen from a quasi-random sequence require no seed and thus are always the same set. This has obvious implications for aliasing-type problems if the locations of the starting points at least have the same fractional positions in each of the sub-cells. Two modifications have been made to tree_random to ameliorate such problems, as given below:

	- The division between subcells is not at exactly 1/2 the distance along each dimension axis of the parent cell: some dither is introduced.

	- The axes of each new sub-cell are randomly permuted, and the signs of displacement along each axis are chosen randomly.

One final adjustment which was needed for good results with quasi-random sequences was to introduce a small buffer zone near sub-cell boundaries within which points are not allowed. Without this, points in abutting sub-cells could approach each other arbitrarily closely, even though points within the same sub-cell cannot.
*/

/*....................................................................*/
void __attribute__((weak))
treePrintMessage(const int status, const char message[TREE_STRLEN]){

  if(     status==TREE_MSG_MESSAGE)
    printf("%s\n", message);
  else if(status==TREE_MSG_WARN)
    printf("Warning: %s\n", message);
  else if(status==TREE_MSG_ERROR)
    printf("Error: %s\n", message);
  else{
    printf("Error! Message status %d not understood.\n", status);
exit(1);
  }

}

/*....................................................................*/
void setConstDefaults(treeRandConstType *rinc){
  rinc->randGenType         = (gsl_rng_type *)gsl_rng_ranlxs2;
  rinc->randSeed            = 342972;
  rinc->quasiRandGenType    = (gsl_qrng_type *)gsl_qrng_halton;
  rinc->numDims             = N_DIMS;
  rinc->numInRandoms        = TREE_N_RANDOMS;
  rinc->verbosity           = 0;
  rinc->totalNumHighPoints  = 0;
  rinc->maxRecursion        = TREE_MAX_RECURSION;
  rinc->maxNumTrials        = TREE_MAX_N_TRIALS;
  rinc->leafBufLenI         = TREE_leafBufLenI;
  rinc->inRandBufLenI       = TREE_inRandBufLenI;
  rinc->abstandFrac         = TREE_abstandFrac;
  rinc->dither              = TREE_DITHER;
  rinc->wholeFieldOrigin[0] = 0.0;
  rinc->wholeFieldOrigin[1] = 0.0;
  rinc->wholeFieldWidth[0]  = 1.0;
  rinc->wholeFieldWidth[1]  = 1.0;
  rinc->allHighPointLoc     = NULL;
  rinc->allHighPointDensy   = NULL;
  rinc->desiredNumPoints    = 0;
  rinc->doShuffle           = 1;
  rinc->doQuasiRandom       = 1;
  rinc->monitorFunc         = NULL;
}

/*....................................................................*/
void _freeSubCell(subCellType *subCell){
  free(subCell->highPointLocations);
  free(subCell->highPointDensities);
}

/*....................................................................*/
void _copySubCell(subCellType *inObject, subCellType *outObject){
  int i,di;

  outObject->numHighPoints        = inObject->numHighPoints;
  outObject->expectedDesNumPoints = inObject->expectedDesNumPoints;
  outObject->sumDensity           = inObject->sumDensity;
  outObject->maxDensity           = inObject->maxDensity;
  outObject->densityIntegral      = inObject->densityIntegral;

  for(di=0;di<N_DIMS;di++){
    outObject->axisIndices[di]       = inObject->axisIndices[di];
    outObject->fieldOrigin[di]       = inObject->fieldOrigin[di];
    outObject->fieldWidth[di]        = inObject->fieldWidth[di];
    outObject->axisSigns[di]         = inObject->axisSigns[di];
    outObject->absRanAcceptRange[di] = inObject->absRanAcceptRange[di];
  }

  if(outObject->numHighPoints>0){
    outObject->highPointLocations = malloc(sizeof(*(outObject->highPointLocations))*outObject->numHighPoints);
    outObject->highPointDensities = malloc(sizeof(*(outObject->highPointDensities))*outObject->numHighPoints);
    for(i=0;i<outObject->numHighPoints;i++){
      for(di=0;di<N_DIMS;di++){
        outObject->highPointLocations[i][di] = inObject->highPointLocations[i][di];
      }
      outObject->highPointDensities[i] = inObject->highPointDensities[i];
    }
  }else{
    outObject->highPointLocations = NULL;
    outObject->highPointDensities = NULL;
  }
}

/*....................................................................*/
void _calcPointLocations(double axesMid[N_DIMS], subCellType *cell\
  , double inRandLocation[N_DIMS], const int numDims, double r[N_DIMS]){
  /*
inRandLocation values should be in the range [-0.5,0.5).
  */
  int di;

  for(di=0;di<numDims;di++)
    r[di] = axesMid[di] + cell->axisSigns[di]*cell->fieldWidth[di]*inRandLocation[cell->axisIndices[di]];
}

/*....................................................................*/
double _calcDensityWithMask(treeRandConstType *rinc, treeRandInternalType *rini\
  , subCellType *cell, double inRandLocation[N_DIMS]\
  , double (*numberDensyFunc)(configInfo*, double*), double axesMid[N_DIMS], _Bool *pointIsInRange, double r[N_DIMS]){
  /*
inRandLocation values should be in the range [-0.5,0.5).
  */

  double density=0.0;
  int di;

  *pointIsInRange = 1; /* the default. */

  if(rinc->doQuasiRandom){ /* only this, because all points are accepted if they are pure random. */
    for(di=0;di<rinc->numDims;di++){
      if(abs(inRandLocation[cell->axisIndices[di]]) > cell->absRanAcceptRange[di]){
        *pointIsInRange = 0;
        break;
      }
    }
  }

  _calcPointLocations(axesMid, cell, inRandLocation, rinc->numDims, r);

  if(*pointIsInRange)
    density = numberDensyFunc(&(rinc->par), r);

  return density;
}

/*....................................................................*/
void _shuffleSubCell(treeRandInternalType *rini, const int numDims, subCellType *cell){
  /*
This chooses a random sign and a random permutation of the axes to apply to the 'universal' set of random points.
  */

  int di;

  for(di=0;di<numDims;di++)
    if(gsl_rng_uniform(rini->randGen)>0.5) cell->axisSigns[di] = -1.0;

  gsl_ran_shuffle(rini->randGen, cell->axisIndices, numDims, sizeof(int));
}

/*....................................................................*/
int _calcSubCellDensities(treeRandConstType *rinc, treeRandInternalType *rini\
  , subCellType *cell, double (*numberDensyFunc)(configInfo*, double*)){
  /*
See general remarks about the algorithm at the head of the present module.

Notes:
  - The function expects rini->inRandLocations to be malloc'd by the calling routine.
  */

  int di,ppi,numGoodPoints=0,status=0;
  double density,axesMid[N_DIMS],subCellVolume,dummyR[N_DIMS];
  _Bool pointIsInRange;

  subCellVolume = 1.0;
  for(di=0;di<rinc->numDims;di++){
    axesMid[di] = cell->fieldOrigin[di] + cell->fieldWidth[di]*0.5;
    subCellVolume *= cell->fieldWidth[di];
  }

  ppi = 0;
  cell->sumDensity = _calcDensityWithMask(rinc, rini, cell, rini->inRandLocations[ppi]\
    , numberDensyFunc, axesMid, &pointIsInRange, dummyR);

  cell->maxDensity = cell->sumDensity;
  if(pointIsInRange) numGoodPoints++;

  for(ppi=1;ppi<rinc->numInRandoms;ppi++){
    density = _calcDensityWithMask(rinc, rini, cell, rini->inRandLocations[ppi]\
      , numberDensyFunc, axesMid, &pointIsInRange, dummyR);

    cell->sumDensity += density;
    if(pointIsInRange) numGoodPoints++;
    if (density>cell->maxDensity)
      cell->maxDensity = density;
  }


  if(numGoodPoints<1)
    return ERR_NO_GOOD_POINTS;

  if(cell->highPointDensities!=NULL){
    for(ppi=0;ppi<cell->numHighPoints;ppi++){
      if (cell->highPointDensities[ppi]>(cell->maxDensity)) /**** check that the location falls within the mask?? */
        cell->maxDensity = cell->highPointDensities[ppi];
    }
  }

  cell->densityIntegral = (cell->sumDensity/(double)numGoodPoints)*subCellVolume;

  return status;
}

/*....................................................................*/
int _initializeTree(treeRandConstType *rinc, treeRandInternalType *rini\
  , subCellType *cell, double (*numberDensyFunc)(configInfo*, double*), treeType *tree){
  /*
See general remarks about the algorithm at the head of the present module.

Notes:
  - The function expects rini->randGen to be correctly allocated by the calling routine via the standard GSL call.
  - The calling routine should free rini->inRandLocations after use.
  */

  int di,i,ppi,status=0;
  double point[rinc->numDims];

  rini->randGen = gsl_rng_alloc(rinc->randGenType);
  gsl_rng_set(rini->randGen,rinc->randSeed);

  if(rinc->doQuasiRandom)
    rini->quasiRandGen = gsl_qrng_alloc(rinc->quasiRandGenType, (unsigned int)rinc->numDims);
  else
    rini->quasiRandGen = NULL;

  cell->numHighPoints = rinc->totalNumHighPoints;
  for(di=0;di<rinc->numDims;di++){
    cell->axisIndices[di] = di;
    cell->axisSigns[di]   = 1.0;
    cell->fieldOrigin[di] = rinc->wholeFieldOrigin[di];
    cell->fieldWidth[di]  = rinc->wholeFieldWidth[di];
    cell->absRanAcceptRange[di] = 0.5;
  }
  cell->expectedDesNumPoints = (double)rinc->desiredNumPoints;

  if(cell->numHighPoints>0){
    cell->highPointLocations = malloc(sizeof(*(cell->highPointLocations))*cell->numHighPoints);
    cell->highPointDensities = malloc(sizeof(*(cell->highPointDensities))*cell->numHighPoints);
    for(i=0;i<cell->numHighPoints;i++){
      for(di=0;di<rinc->numDims;di++){
        cell->highPointLocations[i][di] = rinc->allHighPointLoc[i][di];
      }
      cell->highPointDensities[i] = rinc->allHighPointDensy[i];
    }
  }else{
    cell->highPointLocations = NULL;
    cell->highPointDensities = NULL;
  }

  rini->maxNumTrialsDbl = (double)rinc->maxNumTrials;

  rini->numSubFields = 1;
  for(di=0;di<rinc->numDims;di++)
    rini->numSubFields *= 2;

  rini->inRandLocations = malloc(sizeof(*(rini->inRandLocations))*rinc->numInRandoms);
  if(rinc->doQuasiRandom){
    for(ppi=0;ppi<rinc->numInRandoms;ppi++){
      gsl_qrng_get(rini->quasiRandGen, point);
      for(di=0;di<rinc->numDims;di++)
        rini->inRandLocations[ppi][di] = point[di] - 0.5;
    }
  }else{
    for(ppi=0;ppi<rinc->numInRandoms;ppi++){
      for(di=0;di<rinc->numDims;di++)
        rini->inRandLocations[ppi][di] = gsl_rng_uniform(rini->randGen) - 0.5;
    }
  }

  if(rinc->doShuffle) _shuffleSubCell(rini, rinc->numDims, cell);
  status = _calcSubCellDensities(rinc, rini, cell, numberDensyFunc);

/*** check dither is in range? desired num pts > 0? */

  tree->lastLeafI = 0;
  tree->maxLeafI = rinc->leafBufLenI;
  tree->leaves = malloc(sizeof(*tree->leaves)*tree->maxLeafI);

  return status;
}

/*....................................................................*/
int _constructTheTree(treeRandConstType *rinc, treeRandInternalType *rini\
  , subCellType *cell, int levelI, double (*numberDensyFunc)(configInfo*, double*)\
  , treeType *tree){
  /*
See general remarks about the algorithm at the head of the present module.

The present function is that one that constructs the tree, the relevant information being returned in the list (i.e. pointer) in argument 'tree'. The space is subdivided until the variation of the supplied point distribution function within each sub-cell is judged small enough to allow the rejection method to follow it efficiently.

Notes:
  - The function expects rini->randGen to be correctly allocated by the calling routine via the standard GSL call.
  - The function expects tree->leaves to be allocated by the calling routine to size tree->maxLeafI.
  */
  int di,sci,ppi,status=0;
  _Bool doSubdivide,pointAllocated,subFieldExcluded;
  int *subFieldThisPoint=NULL,*indices=NULL;
  double meanDensity,expectedNumTrialsDbl,dither,xLo,xHi, totalDensityIntegral;
  double sideMaxPoints[rinc->numDims],sideMidPoints[rinc->numDims],r[rinc->numDims];
  subCellType subCell[rini->numSubFields];
  char spaces[2*levelI+1],message[TREE_STRLEN];

  if (rinc->verbosity>1){
    memset(spaces,' ',2*levelI);
    spaces[2*levelI]='\0';
    sprintf(message, "\n%sEntering _constructTheTree, level %d", spaces, levelI);
    treePrintMessage(TREE_MSG_MESSAGE, message);
  }

  if (cell->expectedDesNumPoints<=0.0 && cell->sumDensity>0) { /* should test this before entering the function, but just to be sure... */
    if (rinc->verbosity>1){
      sprintf(message, "%sdesiredNumPoints==0!\n%sLeaving randomsViaTree, level %d", spaces, spaces, levelI);
      treePrintMessage(TREE_MSG_ERROR, message);
    }
    return ERR_NO_POINTS_DESIRED;
  }

  /*. . . . . . . . . . . . . .  . . . . . . . . . . . . . . . . . . . .*/
  /* Decide whether to sub-divide.
  */
  if(rinc->verbosity>3){
    sprintf(message, "%sDecide whether to subdivide.", spaces);
    treePrintMessage(TREE_MSG_MESSAGE, message);
  }

  doSubdivide = 0; /* the default */
  if(levelI<rinc->maxRecursion && cell->sumDensity>0.0){
    meanDensity = cell->sumDensity/(double)rinc->numInRandoms;
    expectedNumTrialsDbl = cell->expectedDesNumPoints*cell->maxDensity/meanDensity;
    if (expectedNumTrialsDbl>rini->maxNumTrialsDbl)
      doSubdivide = 1;
  }else if(cell->sumDensity<=0.0 && cell->maxDensity>0.0){
    /* This means there is a concentration of density (as indicated by the non-zero maximum) which is still too compact for the algorithm to find (hence the zero-valued sum). We just have to subdivide until we zoom in enough to see it, or run out of levels.
    */
    doSubdivide = 1;
  }

  /*. . . . . . . . . . . . . .  . . . . . . . . . . . . . . . . . . . .*/
  if(doSubdivide){
    if(rinc->verbosity>2){
      sprintf(message, "%sSubdividing.", spaces);
      treePrintMessage(TREE_MSG_MESSAGE, message);
    }
    if(rinc->verbosity>3){
      sprintf(message, "%sCalculate origin and dimensions.", spaces);
      treePrintMessage(TREE_MSG_MESSAGE, message);
    }

    /* Calculate the origin and dimensions for the sub-cubes:
    */
    for(di=0;di<rinc->numDims;di++){
      sideMaxPoints[di] = cell->fieldOrigin[di] + cell->fieldWidth[di];
      if (rinc->dither>0.0){
        dither = 0.5*rinc->dither*(gsl_rng_uniform(rini->randGen)-0.5)*cell->fieldWidth[di];
        sideMidPoints[di] = cell->fieldOrigin[di] + cell->fieldWidth[di]*0.5 + dither;
/**** check it is within range? Nah just check dither is in [0,1). Also maybe warn if dither too close to 0 or 1? */
      }else
        sideMidPoints[di] = cell->fieldOrigin[di] + cell->fieldWidth[di]*0.5;
    }

    for(sci=0;sci<rini->numSubFields;sci++){
      for(di=0;di<rinc->numDims;di++){
        if (sci&(1<<di)){
          /*
Relies on the fact that each dimension is divided in 2. Take an example in which rinc->numDims==3. Here there are 3**2=8 subFields. The 6th subField has number 5 in the range 0 to 7; 5 in binary is 011; thus we take this to be the subField which is in the upper half of the X axis, the upper half of the Y axis, but the lower half of the Z axis (1, 1 and 0, starting from the LSB).
          */
          xLo = sideMidPoints[di];
          xHi = sideMaxPoints[di];
        }else{
          xLo = cell->fieldOrigin[di];
          xHi = sideMidPoints[di];
        }

        subCell[sci].axisIndices[di] = di;
        subCell[sci].axisSigns[di]   = 1.0;
        subCell[sci].fieldOrigin[di] = xLo;
        subCell[sci].fieldWidth[di] = xHi - xLo;
        subCell[sci].absRanAcceptRange[di] = 0.5;
      }
      subCell[sci].expectedDesNumPoints = 0;
    }

    /* Allocate the high points between subfields:
    */
    if (rinc->verbosity>3){
      sprintf(message, "%sAllocate density maxima among the subfields.", spaces);
      treePrintMessage(TREE_MSG_MESSAGE, message);
    }

    for(sci=0;sci<rini->numSubFields;sci++)
      subCell[sci].numHighPoints = 0;

    if (cell->numHighPoints==0 || cell->highPointLocations==NULL || cell->highPointDensities==NULL){
      for(sci=0;sci<rini->numSubFields;sci++){
        subCell[sci].highPointDensities = NULL;
        subCell[sci].highPointLocations = NULL;
      }

    } else {
      subFieldThisPoint = malloc(sizeof(int)*rinc->numInRandoms);
      for(ppi=0;ppi<cell->numHighPoints;ppi++){
        if (rinc->verbosity>4){
          sprintf(message, "%s  Density maximum %d", spaces, ppi);
          treePrintMessage(TREE_MSG_MESSAGE, message);
        }

        for(di=0;di<rinc->numDims;di++)
          r[di] = cell->highPointLocations[ppi][di];

        pointAllocated = 0;

        /* Loop over all the subfields and stop when the one which contains the point is found.
        */
        sci = 0;
        while (!pointAllocated && sci<rini->numSubFields){
          subFieldExcluded = 0;

          /* Loop through the dimensions and, for each one, check to see if the respective coordinate of the point is outside the range for that subfield. If so, exit the dimensions loop, because the point can't be in the subfield: there is no point checking any more dimensions.
          */ 
          di = 0;
          while (!subFieldExcluded && di<rinc->numDims){
            if ((r[di]  < subCell[sci].fieldOrigin[di])\
            ||  (r[di] >= subCell[sci].fieldOrigin[di] + subCell[sci].fieldWidth[di]))
              subFieldExcluded = 1;

            di++;
          }
          if (!subFieldExcluded){
            pointAllocated = 1;
            subCell[sci].numHighPoints++;
            subFieldThisPoint[ppi] = sci;
          }
          sci++;
        }
/**** complain if !pointAllocated - bug */
      }

      for(sci=0;sci<rini->numSubFields;sci++){
        if (subCell[sci].numHighPoints>0){
          subCell[sci].highPointDensities = malloc(sizeof(double                            )*subCell[sci].numHighPoints);
          subCell[sci].highPointLocations = malloc(sizeof(*(subCell[sci].highPointLocations))*subCell[sci].numHighPoints);
        } else {
          subCell[sci].highPointDensities = NULL; /*** but this is set already..? */
          subCell[sci].highPointLocations = NULL; /*** but this is set already..? */
        }
      }

      indices = malloc(sizeof(int)*rini->numSubFields);
      memset(indices, 0, sizeof(int)*rini->numSubFields);
      for(ppi=0;ppi<cell->numHighPoints;ppi++){
        sci = subFieldThisPoint[ppi];
        subCell[sci].highPointDensities[indices[sci]] = cell->highPointDensities[ppi];
        for(di=0;di<rinc->numDims;di++)
          subCell[sci].highPointLocations[indices[sci]][di] = cell->highPointLocations[ppi][di];
        indices[sci]++;
      }
      free(indices);
      free(subFieldThisPoint);
    }

    /* Find the sum and max of the density values in each sub-cube, using the offset-and-scaled rini->inRandLocations:
    */
    if (rinc->verbosity>3){
      sprintf(message, "%sCalculate max and sum densitities.", spaces);
      treePrintMessage(TREE_MSG_MESSAGE, message);
    }

    for(sci=0;sci<rini->numSubFields;sci++){
      if(rinc->doShuffle) _shuffleSubCell(rini, rinc->numDims, &subCell[sci]);
      status = _calcSubCellDensities(rinc, rini, &(subCell[sci]), numberDensyFunc);
      if(status!=0) return status;
    }

    /* The following gives us a better estimate of the total density integral because it used 2^N times as many sample points.
    */
    totalDensityIntegral = 0.0;
    for(sci=0;sci<rini->numSubFields;sci++)
      totalDensityIntegral += subCell[sci].densityIntegral;

    /* Calculate an expectation value for the desired number of points for each subcube.
    */
    if (rinc->verbosity>3){
      sprintf(message, "%sCalculate approx. desired number of points.", spaces);
      treePrintMessage(TREE_MSG_MESSAGE, message);
    }

    for(sci=0;sci<rini->numSubFields;sci++)
      subCell[sci].expectedDesNumPoints = cell->expectedDesNumPoints\
                                        * subCell[sci].densityIntegral/totalDensityIntegral;

    /* Now call the function for each valid subcube.
    */
    for(sci=0;sci<rini->numSubFields;sci++){
      if (rinc->verbosity>4){
        sprintf(message, "%sapprox %e points required for sub-field %d."\
          , spaces, subCell[sci].expectedDesNumPoints, sci);
        treePrintMessage(TREE_MSG_MESSAGE, message);
      }
      if (rinc->verbosity>3){
        sprintf(message, "%sCall the function for the %dth sub-field.", spaces, sci);
        treePrintMessage(TREE_MSG_MESSAGE, message);
      }

      status = _constructTheTree(rinc, rini, &subCell[sci], levelI+1, numberDensyFunc, tree);
      if(status!=0) return status;
    }

    /* Free up everything else malloc'd in this block:
    */
    if (rinc->verbosity>3){
      sprintf(message, "%sFree everything.", spaces);
      treePrintMessage(TREE_MSG_MESSAGE, message);
    }

    for(sci=0;sci<rini->numSubFields;sci++)
      _freeSubCell(&(subCell[sci]));

  } else { /* don't divide into sub-cubes, just record the specs in the list of leaves. */
    if (rinc->verbosity>2){
      sprintf(message, "%sNot subdividing.", spaces);
      treePrintMessage(TREE_MSG_MESSAGE, message);
    }

    _copySubCell(cell, &(tree->leaves[tree->lastLeafI]));

    /* Increment the leaf index.
    */
    tree->lastLeafI++;
    if(tree->lastLeafI>=tree->maxLeafI){
      tree->maxLeafI += rinc->leafBufLenI;
      tree->leaves = realloc(tree->leaves, sizeof(*tree->leaves)*tree->maxLeafI);
    }
  }
  if (rinc->verbosity>1){
    sprintf(message, "%sLeaving _constructTheTree, level %d\n", spaces, levelI);
    treePrintMessage(TREE_MSG_MESSAGE, message);
  }

  return status;
}

/*....................................................................*/
_Bool _pointIsAccepted(treeRandConstType *rinc, treeRandInternalType *rini\
  , subCellType *cell, double inRandLocation[N_DIMS]\
  , double (*numberDensyFunc)(configInfo*, double*), double axesMid[rinc->numDims]\
  , double r[N_DIMS], double *density){

  double acceptanceChance;
  _Bool pointIsInRange;

  *density = _calcDensityWithMask(rinc, rini, cell, inRandLocation\
    , numberDensyFunc, axesMid, &pointIsInRange, r);

  if(!pointIsInRange)
    return 0; /* We'd get density and thus acceptanceChance ==0 in any case but no harm in getting out early. */

  acceptanceChance = *density/cell->maxDensity;
  /* Note that we need to screen out cell->maxDensity==0 beforehand. */

  if (acceptanceChance>1.0) acceptanceChance = 1.0; /* Can happen if numberDensyFunc has narrow maxima which our scatter-gun sampling has missed. */

  /* which should suffice to exclude acceptance of density<=0, because this will give (1-p)>=1; but gsl_rng_uniform is always <1. */
  if (gsl_rng_uniform(rini->randGen)>(1.0-acceptanceChance))
    return 1;
  else
    return 0;
}

/*....................................................................*/
int _fillTheTree(treeRandConstType *rinc, treeRandInternalType *rini, treeType *tree\
  , double (*numberDensyFunc)(configInfo*, double*)\
  , double (*outRandLocations)[rinc->numDims], double *outRandDensities){
  /*
See general remarks about the algorithm at the head of the present module.

The present function takes the tree constructed in function _constructTheTree() and fills it with random points. What matters in the tree is not the branches, but the leaves. Each 'leaf' is a sub-cell which the _constructTheTree() algorithm decided it was unnecessary to further subdivide. The entire list of leaves tiles the whole of the initial working space.

The present function first goes through the list of leaves and calculates the desired number of points required for each. These values should follow a binomial distribution, since the total desired number of points is a fixed input, but the expectation value of the number desired for any sub-cell is clearly the fraction of (estimated) total density contained in that cell times the total desired number of points.

The second loop through the cells (or leaves) generates the points. This is done on a cell-by-cell basis by supplying a list of evenly-distributed random locations and testing them against the density function in that cell via the standard rejection method. The algorithm uses the same set of points for each cell (with locations scaled and offset as appropriate to fit neatly into the cell), which is faster than generating a new set for each cell.

It is difficult a priori to estimate the number of input evenly-distributed randoms to supply. If we denote this number by N_in and the desired number by N_out, then N_out/N_in is expected to be proportional to I/(d_max*V), where I is the integral of the density function within the cell, d_max is the maximum in the cell of this function and V is the volume of the cell. Since point selection is a random process however, there will be scatter about this value. Facility is therefore provided for adding additional evenly-distributed randoms to the input list.

The random locations in the input list are optionally quasi-random, which avoids the occasional close neighbours which can occur in a set of entirely random points, and gives rise to more even-sized Delaunay cells without the necessity to undertake a time-consuming smoothing of the entire set of points after generation.

Notes:
  - The function expects rini->randGen to be correctly allocated by the calling routine via the standard GSL call.
  */
  double totalIntegral,probabilities[tree->lastLeafI],axesMid[rinc->numDims],point[rinc->numDims];
  int i,di,ppi,status=0;
  double r[rinc->numDims],density;
  double inRandLocation[N_DIMS];
  unsigned int desiredNumsPoints[tree->lastLeafI];
  unsigned int totalNumPoints,numPointsThisCell,currentPointI;
  gsl_qrng *qrSeqGen = NULL;
  double standoffs[tree->lastLeafI];
  char message[TREE_STRLEN];

  if (rinc->verbosity>1){
    sprintf(message, "  Entering _fillTheTree\n");
    treePrintMessage(TREE_MSG_MESSAGE, message);
  }

  totalIntegral = 0.0;
  for(i=0;i<tree->lastLeafI;i++){
    totalIntegral += tree->leaves[i].densityIntegral;
    standoffs[i] = 0.0;
  }

  /*
Now we calculate the number of points we want for each leaf. This is a two-step algorithm: the first step calculates the expectation value, the second step generates a random number from the multinomial distribution which has that expectation value.
  */
  for(i=0;i<tree->lastLeafI;i++)
    probabilities[i] = tree->leaves[i].densityIntegral/totalIntegral;

  gsl_ran_multinomial(rini->randGen, (size_t)tree->lastLeafI, rinc->desiredNumPoints, probabilities, desiredNumsPoints);

  /*
For quasi-random points we will need a second pass at this, because the final list of locations will be drawn from a set no one of which lies nearer than a standoff value from the subcell boundary. However, since the standoff value depends on the number of points desired, we don't know it a priori; thus the first pass calculates the standoff values for each leaf, the second uses these in a refinement of the calculation of the expectation values and thus the desired number of points.
  */
  if(rinc->doQuasiRandom){
    for(i=0;i<tree->lastLeafI;i++){
      if(desiredNumsPoints[i]<=0)
    continue; /* because increasing the standoff can only decrease the desired number of points. */

      /*
We want to leave a little bit of padding at the edge of the cell. The quasi-random points avoid each other inside the cell ok, but we want to also make sure that the likelihood of points in adjacent cells being close is low. The buffer is scaled to the approx. mean separation between nearest neighbours <NN>, which we here approximate (for a cell of unit side length) by

	                1
	<NN> ~ ---------------------
	        N^{1/rinc->numDims}

where N is the number of points in the cell.
      */
      standoffs[i] = rinc->abstandFrac*pow((double)desiredNumsPoints[i], -1.0/(double)rinc->numDims);

      for(di=0;di<rinc->numDims;di++){
        tree->leaves[i].absRanAcceptRange[di] = 0.5*((tree->leaves[i].fieldWidth[di] - 2.0*standoffs[i])/tree->leaves[i].fieldWidth[di]);
        if(tree->leaves[i].absRanAcceptRange[di] < TREE_MIN_ARA_RANGE)
          tree->leaves[i].absRanAcceptRange[di] = TREE_MIN_ARA_RANGE;
/**** should also raise a status */
      }

      status = _calcSubCellDensities(rinc, rini, &(tree->leaves[i]), numberDensyFunc);
      if(status!=0) return status;
    }

    totalIntegral = 0.0;
    for(i=0;i<tree->lastLeafI;i++){
      totalIntegral += tree->leaves[i].densityIntegral;
    }

    for(i=0;i<tree->lastLeafI;i++)
      probabilities[i] = tree->leaves[i].densityIntegral/totalIntegral;

    gsl_ran_multinomial(rini->randGen, (size_t)tree->lastLeafI, rinc->desiredNumPoints, probabilities, desiredNumsPoints);
  }

  /* Now fill in the points.
  */
  totalNumPoints = 0;
  for(i=0;i<tree->lastLeafI;i++){
    if (rinc->verbosity>2){
      sprintf(message, "    Leaf %04d/%04d", i+1, tree->lastLeafI);
      treePrintMessage(TREE_MSG_MESSAGE, message);
    }

    if(desiredNumsPoints[i]<=0)
  continue;

    for(di=0;di<rinc->numDims;di++)
      axesMid[di] = tree->leaves[i].fieldOrigin[di] + tree->leaves[i].fieldWidth[di]*0.5;

    numPointsThisCell = 0;
    for(ppi=0;ppi<rinc->numInRandoms;ppi++){
      /* Note that we cannot arrive here if desiredNumsPoints[i]<=0, which is the only circumstance in which tree->leaves[i].maxDensity could ==0. */
      if(_pointIsAccepted(rinc, rini, &(tree->leaves[i]), rini->inRandLocations[ppi], numberDensyFunc, axesMid, r, &density)){
        currentPointI = totalNumPoints + numPointsThisCell;
        for(di=0;di<rinc->numDims;di++)
          outRandLocations[currentPointI][di] = r[di];
        outRandDensities[currentPointI] = density;

        numPointsThisCell++;

        if(numPointsThisCell>=desiredNumsPoints[i])
    break;
      }
    }

    if(numPointsThisCell<desiredNumsPoints[i]){ /* the standard set of positions was not enough - generate some extra points. */
      if(rinc->doQuasiRandom)
        qrSeqGen = gsl_qrng_clone(rini->quasiRandGen); /* because we want new points to avoid the old, which means they have to start in the QRG sequence where the old ones left off. */

      for(ppi=0;ppi<TREE_MAX_N_EXTRAS;ppi++){
        /* Generate a new random location. */
        if(rinc->doQuasiRandom){
          gsl_qrng_get(qrSeqGen, point);
          for(di=0;di<rinc->numDims;di++)
            inRandLocation[di] = point[di] - 0.5;

        }else{
          for(di=0;di<rinc->numDims;di++)
            inRandLocation[di] = gsl_rng_uniform(rini->randGen) - 0.5;
        }

        /* Note that we cannot arrive here if desiredNumsPoints[i]<=0, which is the only circumstance in which tree->leaves[i].maxDensity could ==0. */
        if(_pointIsAccepted(rinc, rini, &(tree->leaves[i]), inRandLocation, numberDensyFunc, axesMid, r, &density)){
          currentPointI = totalNumPoints + numPointsThisCell;
          for(di=0;di<rinc->numDims;di++)
            outRandLocations[currentPointI][di] = r[di];
          outRandDensities[currentPointI] = density;

          numPointsThisCell++;

          if(numPointsThisCell>=desiredNumsPoints[i])
      break;
        }
      }

      if(numPointsThisCell<desiredNumsPoints[i]){
        if (rinc->verbosity>2){
          sprintf(message, "Couldn't make enough points.");
          treePrintMessage(TREE_MSG_ERROR, message);
        }
        return ERR_TOO_FEW_POINTS;
      }
    }

    if(rinc->monitorFunc!=NULL)
      rinc->monitorFunc(rinc->numDims, i, tree->leaves[i].fieldOrigin\
        , tree->leaves[i].fieldWidth, desiredNumsPoints[i]\
        , outRandLocations, totalNumPoints, numPointsThisCell);

    totalNumPoints += numPointsThisCell;
  } /* end loop over leaves.*/

  if(rinc->doQuasiRandom)
    gsl_qrng_free(qrSeqGen);

  if (rinc->verbosity>1){
    sprintf(message, "  Leaving _fillTheTree\n");
    treePrintMessage(TREE_MSG_MESSAGE, message);
  }

  return status;
}

/*....................................................................*/
void treeGenerateRandoms(treeRandConstType *rinc, double (*numberDensyFunc)(configInfo*, double*)\
  , double (*outRandLocations)[rinc->numDims], double *outRandDensities){
  /*
Before calling this function you should do:

  outRandDensities = malloc(sizeof(*outRandDensities)*<desired number of points>);
  outRandLocations = malloc(sizeof(*outRandLocations)*<desired number of points>);

The appropriate values of rinc should also be set. Call the function then as

  generateRandoms(&rinc, numberDensyFunc, outRandLocations, outRandDensities);

After you've finished with the results, do:

  free(outRandLocations);
  free(outRandDensities);

  */

  subCellType cell;
  treeRandInternalType rini;
  treeType tree;
  int levelI=0,status=0;

  if(rinc->verbosity>0) treePrintMessage(TREE_MSG_MESSAGE, "+++++ Call _initializeTree()");
  status = _initializeTree(rinc, &rini, &cell, numberDensyFunc, &tree);
  if(status!=0) exit(1);

  if(rinc->verbosity>0) treePrintMessage(TREE_MSG_MESSAGE, "+++++ Call _constructTheTree()");
  status = _constructTheTree(rinc, &rini, &cell, levelI, numberDensyFunc, &tree);
  if(status!=0) exit(1);

  if(rinc->verbosity>0) treePrintMessage(TREE_MSG_MESSAGE, "+++++ Call _fillTheTree()");
  status = _fillTheTree(rinc, &rini, &tree, numberDensyFunc, outRandLocations, outRandDensities);
  if(status!=0) exit(1);

  if(rinc->verbosity>0) treePrintMessage(TREE_MSG_MESSAGE, "+++++ Clean up");
  free(tree.leaves);
  free(rini.inRandLocations);
  gsl_rng_free(rini.randGen);
  gsl_qrng_free(rini.quasiRandGen);
  _freeSubCell(&cell);
}


