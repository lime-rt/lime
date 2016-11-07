/*
 *  tree_random.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2016 The LIME development team
 *
 */

#include "lime.h"
#include "tree_random.h"


/*
The Treegrid algorithm attempts to generate random points within an orthotopic working space which follow a distribution given by some user-supplied number-density function. (An orthotope is a generalization to N dimensions of the concept of a rectangle.) The algorithm is designed to avoid the lengthy run times which can plague the (otherwise robust) rejection method when there are significant variations in the value of the point number density function within the working space. The algorithm makes use of a 2^D tree approach (where D is the number of spatial dimensions) to split the working space into smaller and smaller segments until a space is reached within which the variation of the number density function is within acceptable bounds.
*/

/*....................................................................*/
void calcSubCellValues(treeRandConstType *rinc, treeRandVarType *rinv\
  , double (*numberDensyFunc)(configInfo*, double*)){
  /*
See general remarks about the algorithm at the head of the present module.

Notes:
  - The function expects rinc->inRandLocations to be malloc'd by the calling routine.
  */

  int di,ppi;
  double r[DIM],density;

  rinv->fieldVolume = 1.0;
  for(di=0;di<DIM;di++)
    rinv->fieldVolume *= rinv->fieldWidth[di];

  rinv->sumDensity = 0.0;
  ppi = 0;
  for(di=0;di<DIM;di++)
    r[di] = rinc->inRandLocations[ppi][di]*rinv->fieldWidth[di] + rinv->fieldOrigin[di];
  rinv->sumDensity = numberDensyFunc(&(rinc->par), r);
  rinv->maxDensity = rinv->sumDensity;

  for(ppi=1;ppi<rinc->numInRandoms;ppi++){
    for(di=0;di<DIM;di++)
      r[di] = rinc->inRandLocations[ppi][di]*rinv->fieldWidth[di] + rinv->fieldOrigin[di];

    density = numberDensyFunc(&(rinc->par), r);
    rinv->sumDensity += density;
    if (density>rinv->maxDensity)
      rinv->maxDensity = density;
  }

  if(rinv->highPointDensities!=NULL){
    for(ppi=0;ppi<rinv->numHighPoints;ppi++){
      if (rinv->highPointDensities[ppi]>(rinv->maxDensity))
        rinv->maxDensity = rinv->highPointDensities[ppi];
    }
  }

  rinv->densityIntegral = rinv->sumDensity*rinv->fieldVolume/(double)rinc->numInRandoms;
}


/*....................................................................*/
void initializeTree(treeRandConstType *rinc, treeRandVarType *rinv\
  , double (*numberDensyFunc)(configInfo*, double*), treeType *tree){
  /*
See general remarks about the algorithm at the head of the present module.

Notes:
  - The function expects rinc->randGen to be correctly allocated by the calling routine via the standard GSL call.
  - The calling routine should free rinc->inRandLocations after use.
  */

  int di, ppi;

  rinc->numSubFields = 1;
  for(di=0;di<DIM;di++)
    rinc->numSubFields *= 2;

  rinc->inRandLocations = malloc(sizeof(*(rinc->inRandLocations))*rinc->numInRandoms);
  for(ppi=0;ppi<rinc->numInRandoms;ppi++){
    for(di=0;di<DIM;di++)
      rinc->inRandLocations[ppi][di] = gsl_rng_uniform(rinc->randGen);
  }

  calcSubCellValues(rinc, rinv, numberDensyFunc);

  rinc->leafBufLenI = 10000;
  rinc->inRandBufLenI = 10000;
  rinc->abstandFrac = 0.1; /* Buffer to leave at the cell borders for quasi-random points, expressed as a fraction of the mean distance between points. */

//*** check dither is in range? desired num pts > 0?

  tree->lastLeafI = 0;
  tree->maxLeafI = rinc->leafBufLenI;
  tree->leaves = malloc(sizeof(*tree->leaves)*tree->maxLeafI);
}

/*....................................................................*/
void constructTheTree(treeRandConstType *rinc, treeRandVarType *rinv, int levelI\
  , double (*numberDensyFunc)(configInfo*, double*), treeType *tree){
  /*
See general remarks about the algorithm at the head of the present module.

The present function is that one that constructs the tree, the relevant information being returned in the list (i.e. pointer) in argument 'tree'. The space is subdivided until the variation of the supplied point distribution function within each sub-cell is judged small enough to allow the rejection method to follow it efficiently.

Notes:
  - The function expects rinc->randGen to be correctly allocated by the calling routine via the standard GSL call.
  - The function expects tree->leaves to be allocated by the calling routine to size tree->maxLeafI.
  */
  int di,sci,ppi;
  _Bool doSubdivide,pointAllocated,subFieldExcluded;
  int *subFieldThisPoint=NULL,*indices=NULL;
  double meanDensity,expectedNumTrialsDbl,dither,xLo,xHi, totalDensityIntegral;
  double sideMaxPoints[DIM],sideMidPoints[DIM],r[DIM];
  treeRandVarType subRinv[rinc->numSubFields];
  char spaces[2*levelI+1];

  if (rinc->verbosity>0){
    memset(spaces,' ',2*levelI);
    spaces[2*levelI]='\0';
    printf("\n%sEntering constructTheTree, level %d\n", spaces, levelI);
  }

  if (rinv->expectedDesNumPoints<=0.0 && rinv->sumDensity>0) { /* should test this before entering the function, but just to be sure... */
    if (rinc->verbosity>0) printf("%sdesiredNumPoints==0!\n", spaces);
    if (rinc->verbosity>0) printf("%sLeaving randomsViaTree, level %d\n", spaces, levelI);
    exit(1);
  }

  /*. . . . . . . . . . . . . .  . . . . . . . . . . . . . . . . . . . .*/
  /* Decide whether to sub-divide.
  */
  if(rinc->verbosity>2) printf("%sDecide whether to subdivide.\n", spaces);
  doSubdivide = 0; /* the default */
  if(levelI<rinc->maxRecursion && rinv->sumDensity>0.0){
    meanDensity = rinv->sumDensity/(double)rinc->numInRandoms;
    expectedNumTrialsDbl = rinv->expectedDesNumPoints*rinv->maxDensity/meanDensity;
    if (expectedNumTrialsDbl>rinc->maxNumTrialsDbl)
      doSubdivide = 1;
  }else if(rinv->sumDensity<=0.0 && rinv->maxDensity>0.0){
    /* This means there is a concentration of density (as indicated by the non-zero maximum) which is still too compact for the algorithm to find (hence the zero-valued sum). We just have to subdivide until we zoom in enough to see it, or run out of levels.
    */
    doSubdivide = 1;
  }

  /*. . . . . . . . . . . . . .  . . . . . . . . . . . . . . . . . . . .*/
  if(doSubdivide){
    if(rinc->verbosity>1) printf("%sSubdividing.\n", spaces);
    if(rinc->verbosity>2) printf("%sCalculate origin and dimensions.\n", spaces);

    /* Calculate the origin and dimensions for the sub-cubes:
    */
    for(di=0;di<DIM;di++){
      sideMaxPoints[di] = rinv->fieldOrigin[di] + rinv->fieldWidth[di];
      if (rinc->dither>0.0){
        dither = 0.5*rinc->dither*(gsl_rng_uniform(rinc->randGen)-0.5)*rinv->fieldWidth[di];
        sideMidPoints[di] = rinv->fieldOrigin[di] + rinv->fieldWidth[di]/2.0 + dither;
//**** check it is within range? Nah just check dither is in [0,1). Also maybe warn if dither too close to 0 or 1?
      }else
        sideMidPoints[di] = rinv->fieldOrigin[di] + rinv->fieldWidth[di]/2.0;
    }

    for(sci=0;sci<rinc->numSubFields;sci++){
      for(di=0;di<DIM;di++){
        if (sci&(1<<di)){
          /*
Relies on the fact that each dimension is divided in 2. Take an example in which DIM==3. Here there are 3**2=8 subFields. The 6th subField has number 5 in the range 0 to 7; 5 in binary is 011; thus we take this to be the subField which is in the upper half of the X axis, the upper half of the Y axis, but the lower half of the Z axis (1, 1 and 0, starting from the LSB).
          */
          xLo = sideMidPoints[di];
          xHi = sideMaxPoints[di];
        }else{
          xLo = rinv->fieldOrigin[di];
          xHi = sideMidPoints[di];
        }

        subRinv[sci].fieldOrigin[di] = xLo;
        subRinv[sci].fieldWidth[di] = xHi - xLo;
      }
    }

    /* Allocate the high points between subfields:
    */
    if (rinc->verbosity>2) printf("%sAllocate density maxima among the subfields.\n", spaces);

    for(sci=0;sci<rinc->numSubFields;sci++)
      subRinv[sci].numHighPoints = 0;

    if (rinv->numHighPoints==0 || rinv->highPointLocations==NULL || rinv->highPointDensities==NULL){
      for(sci=0;sci<rinc->numSubFields;sci++){
        subRinv[sci].highPointDensities = NULL;
        subRinv[sci].highPointLocations = NULL;
      }

    } else {
      subFieldThisPoint = malloc(sizeof(int)*rinc->numInRandoms);
      for(ppi=0;ppi<rinv->numHighPoints;ppi++){
        if (rinc->verbosity>3) printf("%s  Density maximum %d\n", spaces, ppi);
        for(di=0;di<DIM;di++)
          r[di] = rinv->highPointLocations[ppi][di];

        pointAllocated = 0;

        /* Loop over all the subfields and stop when the one which contains the point is found.
        */
        sci = 0;
        while (!pointAllocated && sci<rinc->numSubFields){
          subFieldExcluded = 0;

          /* Loop through the dimensions and, for each one, check to see if the respective coordinate of the point is outside the range for that subfield. If so, exit the dimensions loop, because the point can't be in the subfield: there is no point checking any more dimensions.
          */ 
          di = 0;
          while (!subFieldExcluded && di<DIM){
            if ((r[di]  < subRinv[sci].fieldOrigin[di])\
            ||  (r[di] >= subRinv[sci].fieldOrigin[di] + subRinv[sci].fieldWidth[di]))
              subFieldExcluded = 1;

            di++;
          }
          if (!subFieldExcluded){
            pointAllocated = 1;
            subRinv[sci].numHighPoints++;
            subFieldThisPoint[ppi] = sci;
          }
          sci++;
        }
//**** complain if !pointAllocated - bug
      }

      for(sci=0;sci<rinc->numSubFields;sci++){
        if (subRinv[sci].numHighPoints>0){
          subRinv[sci].highPointDensities = malloc(sizeof(double                            )*subRinv[sci].numHighPoints);
          subRinv[sci].highPointLocations = malloc(sizeof(*(subRinv[sci].highPointLocations))*subRinv[sci].numHighPoints);
        } else {
          subRinv[sci].highPointDensities = NULL; //*** but this is set already..?
          subRinv[sci].highPointLocations = NULL; //*** but this is set already..?
        }
      }

      indices = malloc(sizeof(int)*rinc->numSubFields);
      memset(indices, 0, sizeof(int)*rinc->numSubFields);
      for(ppi=0;ppi<rinv->numHighPoints;ppi++){
        sci = subFieldThisPoint[ppi];
        subRinv[sci].highPointDensities[indices[sci]] = rinv->highPointDensities[ppi];
        for(di=0;di<DIM;di++)
          subRinv[sci].highPointLocations[indices[sci]][di] = rinv->highPointLocations[ppi][di];
        indices[sci]++;
      }
      free(indices);
      free(subFieldThisPoint);
    }

    /* Find the sum and max of the density values in each sub-cube, using the offset-and-scaled rinc->inRandLocations:
    */
    if (rinc->verbosity>2) printf("%sCalculate max and sum densitities.\n", spaces);

    for(sci=0;sci<rinc->numSubFields;sci++)
      calcSubCellValues(rinc, &subRinv[sci], numberDensyFunc);

    totalDensityIntegral = 0.0;
    for(sci=0;sci<rinc->numSubFields;sci++)
      totalDensityIntegral += subRinv[sci].densityIntegral;

    /* Calculate an expectation value for the desired number of points for each subcube.
    */
    if (rinc->verbosity>2) printf("%sCalculate approx. desired number of points.\n", spaces);
    for(sci=0;sci<rinc->numSubFields;sci++)
      subRinv[sci].expectedDesNumPoints = rinv->expectedDesNumPoints\
      *subRinv[sci].densityIntegral/totalDensityIntegral;

    /* Now call the function for each valid subcube.
    */
    for(sci=0;sci<rinc->numSubFields;sci++){
      if (rinc->verbosity>3)\
        printf("%sapprox %e points required for sub-field %d.\n"\
        , spaces, subRinv[sci].expectedDesNumPoints, sci);
      if (rinc->verbosity>2) printf("%sCall the function for the %dth sub-field.\n", spaces, sci);

      constructTheTree(rinc, &subRinv[sci], levelI+1, numberDensyFunc, tree);
    }

    /* Free up everything else malloc'd in this block:
    */
    if (rinc->verbosity>2) printf("%sFree everything.\n", spaces);
    for(sci=0;sci<rinc->numSubFields;sci++)
      freeRinv(subRinv[sci]);

  } else { /* don't divide into sub-cubes, just record the specs in the list of leaves. */
    if (rinc->verbosity>1) printf("%sNot subdividing.\n", spaces);

    tree->leaves[tree->lastLeafI].densityIntegral = rinv->densityIntegral;
    tree->leaves[tree->lastLeafI].maxDensity = rinv->maxDensity;
    for(di=0;di<DIM;di++){
      tree->leaves[tree->lastLeafI].fieldOrigin[di] = rinv->fieldOrigin[di];
      tree->leaves[tree->lastLeafI].fieldWidth[di]  = rinv->fieldWidth[di];
    }

    /* Increment the leaf index.
    */
    tree->lastLeafI++;
    if(tree->lastLeafI>=tree->maxLeafI){
      tree->maxLeafI += rinc->leafBufLenI;
      tree->leaves = realloc(tree->leaves, sizeof(*tree->leaves)*tree->maxLeafI);
    }
  }
  if (rinc->verbosity>0) printf("%sLeaving constructTheTree, level %d\n\n", spaces, levelI);
}

/*....................................................................*/
void fillTheTree(treeRandConstType *rinc, treeType *tree\
  , double (*numberDensyFunc)(configInfo*, double*)\
  , double (*outRandLocations)[DIM], double *outRandDensities){
  /*
See general remarks about the algorithm at the head of the present module.

The present function takes the tree constructed in function constructTheTree() and fills it with random points. What matters in the tree is not the branches, but the leaves. Each 'leaf' is a sub-cell which the constructTheTree() algorithm decided it was unnecessary to further subdivide. The entire list of leaves tiles the whole of the initial working space.

The present function first goes through the list of leaves and calculates the desired number of points required for each. These values should follow a binomial distribution, since the total desired number of points is a fixed input, but the expectation value of the number desired for any sub-cell is clearly the fraction of (estimated) total density contained in that cell times the total desired number of points.

The second loop through the cells (or leaves) generates the points. This is done on a cell-by-cell basis by supplying a list of evenly-distributed random locations and testing them against the density function in that cell via the standard rejection method. The algorithm uses the same set of points for each cell (with locations scaled and offset as appropriate to fit neatly into the cell), which is faster than generating a new set for each cell.

It is difficult a priori to estimate the number of input evenly-distributed randoms to supply. If we denote this number by N_in and the desired number by N_out, then N_out/N_in is expected to be proportional to I/(d_max*V), where I is the integral of the density function within the cell, d_max is the maximum in the cell of this function and V is the volume of the cell. Since point selection is a random process however, there will be scatter about this value. Facility is therefore provided for adding additional evenly-distributed randoms to the input list.

The random locations in the input list are optionally quasi-random, which avoids the occasional close neighbours which can occur in a set of entirely random points, and gives rise to more even-sized Delaunay cells without the necessity to undertake a time-consuming smoothing of the entire set of points after generation.

Notes:
  - The function expects rinc->randGen to be correctly allocated by the calling routine via the standard GSL call.
  */
  double totalIntegral,probabilities[tree->lastLeafI],axesMid[DIM],scales[DIM],sign,point[DIM];
  int i,di,maxInRandI,lastInRandI,piIn,indices[DIM];
  double (*inRandLocations)[DIM]=NULL,r[DIM],density,acceptanceChance,shrinkFrac;
  unsigned int desiredNumsPoints[tree->lastLeafI];
  unsigned int piOutStart,piOut,lastPiOut;
  gsl_qrng *qrSeqGen = NULL;

  if(rinc->doQuasiRandom)
    qrSeqGen = gsl_qrng_alloc(gsl_qrng_halton, DIM);

  totalIntegral = 0.0;
  for(i=0;i<tree->lastLeafI;i++){
    totalIntegral += tree->leaves[i].densityIntegral;
  }

  for(i=0;i<tree->lastLeafI;i++)
    probabilities[i] = tree->leaves[i].densityIntegral/totalIntegral;

  gsl_ran_multinomial(rinc->randGen, (size_t)tree->lastLeafI, rinc->desiredNumPoints, probabilities, desiredNumsPoints);

  /* Do an initial malloc of the inRandLocations.
  */
  maxInRandI = rinc->inRandBufLenI;
  inRandLocations = malloc(sizeof(*inRandLocations)*maxInRandI);
  lastInRandI = 0;

  for(di=0;di<DIM;di++)
    indices[di] = di;

  /* Now fill in the points.
  */
  piOutStart = 0;
  for(i=0;i<tree->lastLeafI;i++){
    if(rinc->doShuffle){
      sign = -1.0;
      if(gsl_rng_uniform(rinc->randGen)>0.5) sign = 1.0;

      gsl_ran_shuffle (rinc->randGen, indices, DIM, sizeof(int));

    }else
      sign = 1.0;

    for(di=0;di<DIM;di++){
      axesMid[di] = tree->leaves[i].fieldOrigin[di] + tree->leaves[i].fieldWidth[di]*0.5;
      scales[di] = sign*tree->leaves[i].fieldWidth[di];
    }

    if(rinc->doQuasiRandom){
      /*
We want to leave a little bit of padding at the edge of the cell. The quasi-random points avoid each other inside the cell ok, but we want to also make sure that the likelihood of points in adjacent cells being close is low. The buffer is scaled to the approx. mean separation between nearest neighbours <NN>, which we here approximate (for a cell of unit side length) by

	            1
	<NN> ~ -----------
	        N^{1/DIM}

where N is the number of points in the cell.
      */
      shrinkFrac = 1.0 - 2.0*rinc->abstandFrac*pow((double)desiredNumsPoints[i], -1.0/(double)DIM);
      for(di=0;di<DIM;di++)
        scales[di] *= shrinkFrac;
    }

    piOut = piOutStart;
    lastPiOut = piOutStart + desiredNumsPoints[i];
    piIn = 0;
    while(piOut<lastPiOut){
      if(piIn>=lastInRandI){ /* Expect never to see '>' but include this case for rigour. */
        /* Generate a new random location. */
        if(rinc->doQuasiRandom){
          gsl_qrng_get(qrSeqGen, point);
          for(di=0;di<DIM;di++)
            inRandLocations[lastInRandI][di] = point[di] - 0.5;

        }else{
          for(di=0;di<DIM;di++)
            inRandLocations[lastInRandI][di] = gsl_rng_uniform(rinc->randGen) - 0.5;
        }

        lastInRandI++;
        if(lastInRandI>=maxInRandI){ /* need more room. */
          maxInRandI += rinc->inRandBufLenI;
          inRandLocations = realloc(inRandLocations, sizeof(*inRandLocations)*maxInRandI);
        }
      }

      /* Scale and offset the input location.
      */
      for(di=0;di<DIM;di++)
        r[di] = axesMid[di] + scales[di]*inRandLocations[piIn][indices[di]];

      density = numberDensyFunc(&(rinc->par), r);
      acceptanceChance = density/tree->leaves[i].maxDensity;
      /* Note that the 'while' loop will not be entered if desiredNumsPoints[i]<=0, which is the only circumstance in which tree->leaves[i].maxDensity could ==0. */

      if (acceptanceChance>1.0) acceptanceChance = 1.0; /* Can happen if numberDensyFunc has narrow maxima which our scatter-gun sampling has missed. */
      if (gsl_rng_uniform(rinc->randGen)>(1.0-acceptanceChance)){ /* which should suffice to exclude acceptance of density<=0, because this will give (1-p)>=1; but gsl_rng_uniform is always <1. */
        for(di=0;di<DIM;di++)
          outRandLocations[piOut][di] = r[di];
        outRandDensities[piOut] = density;

        piOut++;
      }
      piIn++;
    } /* End of point generation for cell i. */

    piOutStart += desiredNumsPoints[i];
  }

  free(inRandLocations);
}

/*....................................................................*/
void freeRinc(treeRandConstType rinc){
  free(rinc.inRandLocations);
  gsl_rng_free(rinc.randGen);
}

/*....................................................................*/
void freeRinv(treeRandVarType rinv){
  free(rinv.highPointLocations);
  free(rinv.highPointDensities);
}

