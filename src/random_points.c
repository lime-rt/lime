/*
 *  random_points.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2015 The LIME development team
 *
 */

#include "lime.h"

void initializeTree(double *fieldOrigin, double *fieldDimensions\
  , unsigned int desiredNumPoints, gsl_rng *randGen, double (*numberDensyFunc)(locusType)\
  , int *numSubFields, double *fieldVolume, double *sumDensity, double *maxDensity\
  , int numHighPoints, locusType *highPointLocations, double *highPointDensities\
  , int numInRandoms, locusType **inRandLocations, double **inRandDensities){

  int di, i;
  double x;
  locusType location;

  *fieldVolume = 1.0;
  *numSubFields = 1;
  for(di=0;di<DIM;di++){
    *fieldVolume *= fieldDimensions[di];
    *numSubFields *= 2;
  }

  /* Generate the evenly-distributed random points:
  */
  (*inRandDensities) = malloc(sizeof(double   )*numInRandoms);
  (*inRandLocations) = malloc(sizeof(locusType)*numInRandoms);

  *sumDensity = 0.0;
  for(i=0;i<numInRandoms;i++){
    for(di=0;di<DIM;di++){
      x = gsl_rng_uniform(randGen);
      (*inRandLocations)[i].x[di] = x;
      location.x[di] = fieldOrigin[di] + fieldDimensions[di]*x;
    }

    (*inRandDensities)[i] = numberDensyFunc(location);
    *sumDensity += (*inRandDensities)[i];
  }

  /* Find the maximum density:
  */
  i = 0;
  *maxDensity = (*inRandDensities)[i];
  for(i=1;i<numInRandoms;i++){
    if ((*inRandDensities)[i]>(*maxDensity))
      *maxDensity = (*inRandDensities)[i];
  }

  if(highPointDensities!=NULL){
    for(i=0;i<numHighPoints;i++)
      if (highPointDensities[i]>(*maxDensity))
        *maxDensity = highPointDensities[i];
  }
}

void randomsViaTree(int levelI, int numSubFields, double *fieldOrigin, double *fieldDimensions\
  , double fieldVolume, unsigned int desiredNumPoints, unsigned int startI\
  , int numHighPoints, locusType *highPointLocations, double *highPointDensities\
  , int numInRandoms, locusType *inRandLocations, double *inRandDensities\
  , double sumDensity, double maxDensity, gsl_rng *randGen, double (*numberDensyFunc)(locusType)\
  , void (*monitorFunc)(locusType*, double*, unsigned int, unsigned int, double*, double*)\
  , locusType *outRandLocations, double *outRandDensities, int verbosity){

  /*
This function is designed to be called recursively to choose random points having a distribution which follows some user-supplied number-density function, within an orthotopic space. (An orthotope is a generalization to N dimensions of the concept of a rectangle.)

Each function call is supplied with the dimensions of its orthotope, the required number of points, and a function which returns a number-density value for a location argument. Other arguments are also supplied which are described later. The function must decide whether to further divide the orthotope, cycling then through recursive calls for all these sub-orthotopes, or to go ahead and generate its required number of points via the rejection method.

The points are loaded into, and thus returned in, the appropriate series of entries of the pointers outRandLocations and outRandDensities.

The only criterion in deciding whether to subdivide or not is whether the expected amount of work to do in trying points exceeds a preset cutoff value. Essentially this is proportional to the number of points desired from the subField divided by the 'filling fraction' of the density function. A density function that is highly peaked in the present field will likely have a small filling fraction and thus likely require subdivision.

Explanation of the input arguments:

	- levelI: the current recursion depth.

	- numSubFields: should be set to 2**DIM by the calling program. It would be safer to calculate it anew with each function call, but a bit slower; that's why I am trusting the caller.

	- fieldOrigin: a vector of length DIM giving the origin of the present orthotope (aka field).

	- fieldDimensions: a vector of length DIM giving the side lengths of the present orthotope.

	- fieldVolume: the volume of the present orthotope. It would be safer (but slower) to calculate it anew with each function call.

	- desiredNumPoints: the number of output random points desired in the present field.

	- startI: this is the index of the first unfilled values in the vectors of points outRandLocations and outRandDensities (see below).

	- numHighPoints: the number of entries (can be 0) in the following:

	- highPointLocations, highPointDensities: these lists allow the user to specify maxima of the number-density function. This is useful in the case that these peaks in the function are so narrow that they are unlikely to be properly sampled by the algorithm, at least in the early stages of subField division. Setting the number-density maximum thus at its true peak value(s) ensures that fields into which the highPointLocations fall will correctly assess the need to subdivide.

	- numInRandoms: the length of the list of inRandLocations and inRandDensities.

	- inRandLocations: a list of random positions distributed evenly through a dimension-N orthotope of unit side lengths. The same set of points, appropriately scaled and offset, is used for all the sub-fields. This is probably quicker than calculating new random positions for each sub-field, and there doesn't seem to be any reason not to reuse the one set of positions.

	- inRandDensities: values of the number-density function evaluated at the inRandLocations. These are of course re-evaluated for each sub-field.

	- sumDensity: the sum of the values of inRandDensities. (Note that highPointDensities are not included.)

	- maxDensity: the maximum of inRandDensities (and here highPointDensities are also included, if non-NULL and numHighPoints>0).

	- randGen: random number generator.

	- numberDensyFunc: this is a user-supplied function which must take a locusType argument and return a double. Note that it is assumed that numberDensyFunc 'knows' the value of DIM a priori.

	- monitorFunc: this is optional: if not set to NULL, it should be a user-supplied function which takes the arguments shown. The purpose of this is essentially to aid in bug finding and/or logging.

	- outRandLocations, outRandDensities: the pointers which are filled by successive calls to the function randomsViaTree().


Note that the following are expected to be either macros defined in a header file or globally accessible variables:
	* MAX_RECURSION: limit on the depth of recursion.

	* DIM: dimensionality of the space being populated.

	* MAX_N_TRIALS_TREE: essentially a limit on the number of rejection-method discards permitted. Cells for which the estimated number of discards would be higher than this are subdivided.

	* TREE_DITHER: when subdividing, the field dimensions are divided in two, but with this dither applied. It must be less than unity.
  */

  int doSubdivide, expectedNumTrials, di, sci, ppi, pointAllocated;
  int subFieldExcluded, localStartI, piIn, piOut, lastPiOut, timeOut;
  int *subFieldNumHighPoints=NULL, *subFieldThisPoint=NULL, *indices=NULL;
  unsigned int remainingDesNum;
  unsigned int *subFieldDesNumPoints=NULL;
  double meanDensity, dither, xLo, xHi, totalDensityIntegral, remainingDensy;
  double probability, acceptanceChance, x, density;
  double *sideMidPoints=NULL, *sideMaxPoints=NULL, *subFieldVolumes=NULL;
  double *subFieldMaxDensities=NULL, *subFieldSumDensities=NULL;
  double *subFieldDensityIntegrals=NULL;
  double **subFieldOrigins=NULL, **subFieldLengths=NULL, **subFieldDensities=NULL;
  double **subFieldHighDensities=NULL;
  locusType location, **subFieldHighLocations=NULL;

  char spaces[2*levelI+1];

  if (verbosity>0){
    memset(spaces,' ',2*levelI);
    spaces[2*levelI]='\0';
    printf("\n%sEntering randomsViaTree, level %d, startI=%d\n", spaces, levelI, startI);
  }

  if (desiredNumPoints<1) { // should test this before entering the function, but just to be sure...
    if (verbosity>0) printf("%sdesiredNumPoints==0!\n", spaces);
    if (verbosity>0) printf("%sLeaving randomsViaTree, level %d\n", spaces, levelI);
    exit(1);
  }

  /*. . . . . . . . . . . . . .  . . . . . . . . . . . . . . . . . . . .*/
  /* Decide whether to sub-divide.
  */
  if (verbosity>2) printf("%sDecide whether to subdivide.\n", spaces);
  doSubdivide = 0; // default
  if (levelI<MAX_RECURSION){
//****    if sumDensity<=0.0 raise a 'bug' exception. This should not happen if desiredNumPoints>0.
    meanDensity = sumDensity/(double)numInRandoms;
    expectedNumTrials = (double)(desiredNumPoints)*maxDensity/meanDensity;
    if (expectedNumTrials>MAX_N_TRIALS_TREE)
      doSubdivide = 1;
  }

  /*. . . . . . . . . . . . . .  . . . . . . . . . . . . . . . . . . . .*/
  if (doSubdivide){
    if (verbosity>1) printf("%sSubdividing.\n", spaces);
    if (verbosity>2) printf("%sCalculate origin and dimensions.\n", spaces);

    /* Calculate the origin and dimensions for the sub-cubes:
    */
    sideMidPoints   = malloc(sizeof(double)*DIM);
    sideMaxPoints   = malloc(sizeof(double)*DIM);
    subFieldVolumes = malloc(sizeof(double)  *numSubFields);
    subFieldOrigins = malloc(sizeof(double *)*numSubFields);
    subFieldLengths = malloc(sizeof(double *)*numSubFields);
    for(sci=0;sci<numSubFields;sci++){
      subFieldOrigins[sci] = malloc(sizeof(double)*DIM);
      subFieldLengths[sci] = malloc(sizeof(double)*DIM);
    }

    for(di=0;di<DIM;di++){
      sideMaxPoints[di] = fieldOrigin[di] + fieldDimensions[di];
      if (TREE_DITHER>0.0){
        dither = 0.5*TREE_DITHER*(gsl_rng_uniform(randGen)-0.5)*fieldDimensions[di];
        sideMidPoints[di] = fieldOrigin[di] + fieldDimensions[di]/2.0 + dither;
//**** check it is within range? Nah just check dither is in [0,1). Also maybe warn if dither too close to 0 or 1?
      }else
        sideMidPoints[di] = fieldOrigin[di] + fieldDimensions[di]/2.0;
    }

    for(sci=0;sci<numSubFields;sci++){
      subFieldVolumes[sci] = 1.0;
      for(di=0;di<DIM;di++){
        if (sci&(1<<di)){
          /* Relies on the fact that each dimension is divided in 2. Take an example in which DIM==3. Here there are 3**2=8 subFields. The 6th subField has number 5 in the range 0 to 7; 5 in binary is 011; thus we take this to be the subField which is in the upper half of the X axis, the upper half of the Y axis, but the lower half of the Z axis (1, 1 and 0, starting from the LSB). */
          xLo = sideMidPoints[di];
          xHi = sideMaxPoints[di];
        }else{
          xLo = fieldOrigin[di];
          xHi = sideMidPoints[di];
        }

        subFieldOrigins[sci][di] = xLo;
        subFieldLengths[sci][di] = xHi - xLo;
        subFieldVolumes[sci] *= subFieldLengths[sci][di];
      }
    }

    /* Allocate the high points between subFields:
    */
    if (verbosity>2) printf("%sAllocate density maxima among the subfields.\n", spaces);
    subFieldNumHighPoints = malloc(sizeof(int)        *numSubFields);
    subFieldHighDensities = malloc(sizeof(double    *)*numSubFields);
    subFieldHighLocations = malloc(sizeof(locusType *)*numSubFields);
    memset(subFieldNumHighPoints, 0, sizeof(int)*numSubFields);

    if (numHighPoints==0 || highPointLocations==NULL || highPointDensities==NULL){
      for(sci=0;sci<numSubFields;sci++){
        subFieldHighDensities[sci] = NULL;
        subFieldHighLocations[sci] = NULL;
      }

    } else {
      subFieldThisPoint = malloc(sizeof(int)*numInRandoms);
      for(ppi=0;ppi<numHighPoints;ppi++){
        if (verbosity>3) printf("%s  Density maximum %d\n", spaces, ppi);
        for(di=0;di<DIM;di++)
          location.x[di] = highPointLocations[ppi].x[di];

        pointAllocated = 0;

        /* Loop over all the subfields and stop when the one which contains the point is found.
        */
        sci = 0;
        while (!pointAllocated && sci<numSubFields){
          subFieldExcluded = 0;

          /* Loop through the dimensions and, for each one, check to see if the respective coordinate of the point is outside the range for that subfield. If so, exit the dimensions loop, because the point can't be in the subfield: there is no point checking any more dimensions.
          */ 
          di = 0;
          while (!subFieldExcluded && di<DIM){
            if ((location.x[di]  < subFieldOrigins[sci][di])\
            ||  (location.x[di] >= subFieldOrigins[sci][di] + subFieldLengths[sci][di]))
              subFieldExcluded = 1;

            di++;
          }
          if (!subFieldExcluded){
            pointAllocated = 1;
            subFieldNumHighPoints[sci]++;
            subFieldThisPoint[ppi] = sci;
          }
          sci++;
        }
//**** complain if !pointAllocated - bug
      }

      for(sci=0;sci<numSubFields;sci++){
        if (subFieldNumHighPoints[sci]>0){
          subFieldHighDensities[sci] = malloc(sizeof(double   )*subFieldNumHighPoints[sci]);
          subFieldHighLocations[sci] = malloc(sizeof(locusType)*subFieldNumHighPoints[sci]);
        } else {
          subFieldHighDensities[sci] = NULL;
          subFieldHighLocations[sci] = NULL;
        }
      }

      indices = malloc(sizeof(int)*numSubFields);
      memset(indices, 0, sizeof(int)*numSubFields);
      for(ppi=0;ppi<numHighPoints;ppi++){
        sci = subFieldThisPoint[ppi];
        subFieldHighDensities[sci][indices[sci]] = highPointDensities[ppi];
        for(di=0;di<DIM;di++)
          subFieldHighLocations[sci][indices[sci]].x[di] = highPointLocations[ppi].x[di];
        indices[sci]++;
      }
      free(indices);
      free(subFieldThisPoint);
    }

    /* Find the sum and max of the density values in each sub-cube, using the offset-and-scaled inRandLocations:
    */
    if (verbosity>2) printf("%sCalculate max and sum densitities.\n", spaces);
    subFieldDensities    = malloc(sizeof(double *)*numSubFields);
    subFieldMaxDensities = malloc(sizeof(double)  *numSubFields);
    subFieldSumDensities = malloc(sizeof(double)  *numSubFields);
    for(sci=0;sci<numSubFields;sci++){
      subFieldDensities[sci] = malloc(sizeof(double)*numInRandoms);
      subFieldSumDensities[sci] = 0.0;
      for(ppi=0;ppi<numInRandoms;ppi++){
        for(di=0;di<DIM;di++)
          location.x[di] = inRandLocations[ppi].x[di]*subFieldLengths[sci][di] + subFieldOrigins[sci][di];
        subFieldDensities[sci][ppi] = numberDensyFunc(location);
        subFieldSumDensities[sci] += subFieldDensities[sci][ppi];
      }

      ppi = 0;
      subFieldMaxDensities[sci] = subFieldDensities[sci][ppi];
      for(ppi=1;ppi<numInRandoms;ppi++)
        if (subFieldDensities[sci][ppi]>subFieldMaxDensities[sci])
          subFieldMaxDensities[sci] = subFieldDensities[sci][ppi];

      if(highPointDensities!=NULL){
        for(ppi=0;ppi<numHighPoints;ppi++)
          if (highPointDensities[ppi]>subFieldMaxDensities[sci])
            subFieldMaxDensities[sci] = highPointDensities[ppi];
      }
    }

    /* Estimate the integral of the density in each sub-cube:
    */
    if (verbosity>2) printf("%sEstimate density integrals.\n", spaces);
    subFieldDensityIntegrals = malloc(sizeof(double)*numSubFields);
    totalDensityIntegral = 0.0;
    for(sci=0;sci<numSubFields;sci++){
      subFieldDensityIntegrals[sci] = subFieldSumDensities[sci]\
        * subFieldVolumes[sci] / (double)numInRandoms;
      totalDensityIntegral += subFieldDensityIntegrals[sci];

      if (verbosity>3) printf("%sSummed density for sub-field %d is %e.\n", spaces, sci, subFieldSumDensities[sci]);
      if (verbosity>3) printf("%sIntegrated density for sub-field %d is %e.\n", spaces, sci, subFieldDensityIntegrals[sci]);
    }

//*** raise exception if totalDensityIntegral<=0

    /* Calculate a desired number of points for each subcube. We'll only bother with cubes with non-zero densities.
    */
    if (verbosity>2) printf("%sCalculate desired number of points.\n", spaces);
    subFieldDesNumPoints = malloc(sizeof(unsigned int)*numSubFields);
    remainingDesNum = desiredNumPoints;
    remainingDensy = totalDensityIntegral;
    for(sci=0;sci<numSubFields;sci++){
      /* Generate a random number from a binomial distribution for the desired number of points in the sub-cube:
      */
      probability = subFieldDensityIntegrals[sci]/remainingDensy;
      if (probability>1.0) probability = 1.0;
//*********** have a think what happens (or is supposed to happen) if this ==1.

      if (probability>0){
        subFieldDesNumPoints[sci] = gsl_ran_binomial(randGen, probability, remainingDesNum);
        remainingDesNum -= subFieldDesNumPoints[sci];
        remainingDensy -= subFieldDensityIntegrals[sci];
      }else
        subFieldDesNumPoints[sci] = 0;
    }

//*** raise exception if remainingDesNum>0

    /* Now call the function for each valid subcube. Ignore subcubes with desiredNumPoints<=0.
    */
    localStartI = startI;
    for(sci=0;sci<numSubFields;sci++){
      if (verbosity>3) printf("%s%d points required for sub-field %d.\n", spaces, subFieldDesNumPoints[sci], sci);
      if (subFieldDesNumPoints[sci]>0){
        if (verbosity>2) printf("%sCall the function for the %dth sub-field.\n", spaces, sci);

        randomsViaTree(levelI+1, numSubFields, subFieldOrigins[sci]\
          , subFieldLengths[sci], subFieldVolumes[sci], subFieldDesNumPoints[sci]\
          , localStartI, subFieldNumHighPoints[sci]\
          , subFieldHighLocations[sci], subFieldHighDensities[sci], numInRandoms, inRandLocations, subFieldDensities[sci]\
          , subFieldSumDensities[sci], subFieldMaxDensities[sci]\
          , randGen, numberDensyFunc, monitorFunc, outRandLocations, outRandDensities, verbosity);

        localStartI += subFieldDesNumPoints[sci];
      }
    }

    /* Free up everything else mallocd in this block:
    */
    if (verbosity>2) printf("%sFree everything.\n", spaces);
    free(subFieldDesNumPoints);
    free(subFieldDensityIntegrals);
    free(subFieldSumDensities);
    free(subFieldMaxDensities);
    for(sci=0;sci<numSubFields;sci++){
      if (subFieldDensities[sci]!=NULL)
        free(subFieldDensities[sci]);
      if (subFieldHighDensities[sci]!=NULL)
        free(subFieldHighDensities[sci]);
      if (subFieldHighLocations[sci]!=NULL){
        free(subFieldHighLocations[sci]);
      }
      if (subFieldOrigins[sci]!=NULL)
        free(subFieldOrigins[sci]);
      if (subFieldLengths[sci]!=NULL)
        free(subFieldLengths[sci]);
    }
    free(subFieldDensities);
    free(subFieldHighDensities);
    free(subFieldHighLocations);
    free(subFieldNumHighPoints);
    free(subFieldOrigins);
    free(subFieldLengths);
    free(subFieldVolumes);
    free(sideMaxPoints);
    free(sideMidPoints);

  } else { /* don't divide into sub-cubes, just generate the required points for this cube. */
    if (verbosity>1) printf("%sNot subdividing.\n", spaces);
    /* We already have a bunch of evenly-distributed randoms with associated density values. So first let's loop through these, accepting any that pass the density criterion; we'll exit the loop (and the function) if we get enough points before we are finished.
    */
    piOut = startI;
    lastPiOut = startI + desiredNumPoints;
    piIn = 0;
    while(piIn<numInRandoms && piOut<lastPiOut){
      acceptanceChance = inRandDensities[piIn]/maxDensity;
      if (acceptanceChance>1.0) acceptanceChance = 1.0;
      if (gsl_rng_uniform(randGen)>(1.0-acceptanceChance)){ // which should suffice to exclude acceptance of density<=0, because this will give (1-p)>=1; but gsl_rng_uniform is always <1.
        for(di=0;di<DIM;di++){
          outRandLocations[piOut].x[di] = fieldOrigin[di] + fieldDimensions[di]*inRandLocations[piIn].x[di];
        }
        outRandDensities[piOut] = inRandDensities[piIn];

        piOut++;
      }
      piIn++;
    }

    timeOut = 0;
    while (piOut<lastPiOut || timeOut){
      /* If we entered this loop, it means our input random points were not sufficient to make up the desired number. Generate more until we have enough.
      */
      for(di=0;di<DIM;di++)
        location.x[di] = fieldOrigin[di] + fieldDimensions[di]*gsl_rng_uniform(randGen);

      density = numberDensyFunc(location);
      acceptanceChance = density/maxDensity;
      if (acceptanceChance>1.0) acceptanceChance = 1.0;
      if (gsl_rng_uniform(randGen)>(1.0-acceptanceChance)){
        for(di=0;di<DIM;di++)
          outRandLocations[piOut].x[di] = location.x[di];
        outRandDensities[piOut] = density;

        piOut++;
      }
      //**** set timeOut.
    }

    if (verbosity>2) printf("%spiOut=%d, %d expected\n", spaces, piOut, lastPiOut);

    if (monitorFunc!=NULL)
      monitorFunc(outRandLocations, outRandDensities, desiredNumPoints, startI, fieldOrigin, fieldDimensions);

  }
  if (verbosity>0) printf("%sLeaving randomsViaTree, level %d\n\n", spaces, levelI);
}

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
  extern double modelRadiusSquared, minDensity;
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

  if(totalDensity<minDensity)
    return minDensity;
  else
    return totalDensity;
}


