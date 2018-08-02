/*
 *  init.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#include "lime.h"

/*....................................................................*/
int
copyInpars(const inputPars inpars, image *inimg, const int nImages\
  , configInfo *par, imageInfo **img){
  /*
The parameters visible to the user have now been strictly confined to members of the structs 'inputPars' and 'image', both of which are defined in inpars.h. There are however further internally-set values which is is convenient to bundle together with the user-set ones in configuration structs. At present we have a fairly clunky arrangement in which the user-set values are copied member-by-member from the user-dedicated structs to the generic internal structs 'configInfo' and 'imageInfo'. This is done in the present function. This seems like duplicated effort but it is a good principle to separate the two sorts of structs, for several reasons, as follows. (i) We will usually want more config parameters than user-settable ones. The separation leaves it clearer which things the user needs to (or can) set. (ii) The separation allows checking and screening out of impossible combinations of parameters. (iii) We can adopt new names (for clarity) for config parameters without bothering the user with a changed interface.
  */
  int i,j,numGirDatFiles,nSpecies=0;

  par->radius            = inpars.radius;
  par->minScale          = inpars.minScale;
  par->pIntensity        = inpars.pIntensity;
  par->sinkPoints        = inpars.sinkPoints;
  par->samplingAlgorithm = inpars.samplingAlgorithm;
  par->sampling          = inpars.sampling;
  par->tcmb              = inpars.tcmb;
  par->lte_only          = inpars.lte_only;
  par->init_lte          = inpars.init_lte;
  par->blend             = inpars.blend;
  par->antialias         = inpars.antialias;
  par->polarization      = inpars.polarization;
  par->nThreads          = inpars.nThreads;
  par->nSolveIters       = inpars.nSolveIters;
  par->traceRayAlgorithm = inpars.traceRayAlgorithm;
  par->resetRNG          = inpars.resetRNG;
  par->doSolveRTE        = inpars.doSolveRTE;

  /* Somewhat more carefully copy over the strings:
  */
  copyInparStr(inpars.dust,          &(par->dust));
  copyInparStr(inpars.outputfile,    &(par->outputfile));
  copyInparStr(inpars.binoutputfile, &(par->binoutputfile));
  copyInparStr(inpars.restart,       &(par->restart));
  copyInparStr(inpars.gridfile,      &(par->gridfile));
  copyInparStr(inpars.pregrid,       &(par->pregrid));
  copyInparStr(inpars.gridInFile,    &(par->gridInFile));

  par->gridOutFiles = malloc(sizeof(char *)*NUM_GRID_STAGES);
  for(i=0;i<NUM_GRID_STAGES;i++)
    copyInparStr(inpars.gridOutFiles[i], &(par->gridOutFiles[i]));

  for(i=0;i<NUM_GRID_STAGES;i++)
    par->writeGridAtStage[i] = 0;

  /* If the user has provided a list of moldatfile names, the corresponding elements of par->moldatfile will be non-NULL. Thus we can deduce the number of files (species) from the number of non-NULL elements.
  */
  nSpecies=0;
  while(nSpecies<MAX_NSPECIES && !charPtrIsNullOrEmpty(inpars.moldatfile[nSpecies]))
  nSpecies++;

  numGirDatFiles=0;
  while(numGirDatFiles<MAX_NSPECIES && !charPtrIsNullOrEmpty(inpars.girdatfile[numGirDatFiles]))
    numGirDatFiles++;

  if(numGirDatFiles<=0)
    par->girdatfile = NULL;
  else if(numGirDatFiles!=nSpecies){
    if(!silent) bail_out("Number of girdatfiles different from number of species.");
exit(1);
  }else{
    par->girdatfile=malloc(sizeof(char *)*nSpecies);
    for(i=0;i<nSpecies;i++){
      copyInparStr(inpars.girdatfile[i], &(par->girdatfile[i]));
    }
  }

  /* Copy over the moldatfiles.
  */
  if(nSpecies <= 0){
    par->moldatfile = NULL; /* This will be tested for all line images, so we can never get par->nLineImages>0 if no moldata files have been supplied. */

  } else {
    par->moldatfile=malloc(sizeof(char *)*nSpecies);
    for(i=0;i<nSpecies;i++){
      copyInparStr(inpars.moldatfile[i], &(par->moldatfile[i]));
    }
  }

  /* Copy over the collision-partner pointers:
  */
  par->collPartIds  = malloc(sizeof(int   )*MAX_N_COLL_PART);
  for(i=0;i<MAX_N_COLL_PART;i++) par->collPartIds[i] = inpars.collPartIds[i];
  par->nMolWeights  = malloc(sizeof(double)*MAX_N_COLL_PART);
  for(i=0;i<MAX_N_COLL_PART;i++) par->nMolWeights[i] = inpars.nMolWeights[i];
  par->collPartNames = malloc(sizeof(char*)*MAX_N_COLL_PART);
  for(i=0;i<MAX_N_COLL_PART;i++) copyInparStr(inpars.collPartNames[i], &(par->collPartNames[i]));
  par->collPartMolWeights = malloc(sizeof(double)*MAX_N_COLL_PART);
  for(i=0;i<MAX_N_COLL_PART;i++) par->collPartMolWeights[i] = inpars.collPartMolWeights[i];

  /* Copy over the grid-density-maximum data and find out how many were set:
  */
  par->gridDensMaxValues = malloc(sizeof(*(par->gridDensMaxValues))*MAX_N_HIGH);
  par->gridDensMaxLoc    = malloc(sizeof(*(par->gridDensMaxLoc))*MAX_N_HIGH);
  for(i=0;i<MAX_N_HIGH;i++){
    par->gridDensMaxValues[i] = inpars.gridDensMaxValues[i];
    for(j=0;j<DIM;j++) par->gridDensMaxLoc[i][j] = inpars.gridDensMaxLoc[i][j];
  }

  /* Copy over user-set image values to the generic struct:
  */
  if(nImages>0){
    *img = malloc(sizeof(**img)*nImages);
    for(i=0;i<nImages;i++){
      (*img)[i].nchan      = inimg[i].nchan;
      (*img)[i].trans      = inimg[i].trans;
      (*img)[i].molI       = inimg[i].molI;
      (*img)[i].velres     = inimg[i].velres;
      (*img)[i].imgres     = inimg[i].imgres;
      (*img)[i].pxls       = inimg[i].pxls;
      (*img)[i].unit       = inimg[i].unit;
      copyInparStr(inimg[i].units, &((*img)[i].units));
      (*img)[i].freq       = inimg[i].freq;
      (*img)[i].bandwidth  = inimg[i].bandwidth;
      copyInparStr(inimg[i].filename, &((*img)[i].filename));
      (*img)[i].source_vel = inimg[i].source_vel;
      (*img)[i].theta      = inimg[i].theta;
      (*img)[i].phi        = inimg[i].phi;
      (*img)[i].incl       = inimg[i].incl;
      (*img)[i].posang     = inimg[i].posang;
      (*img)[i].azimuth    = inimg[i].azimuth;
      (*img)[i].distance   = inimg[i].distance;
      (*img)[i].doInterpolateVels = inimg[i].doInterpolateVels; // This is only accessed if par->traceRayAlgorithm==1.
    }
  }

  return nSpecies;
}

/*....................................................................*/
void
setOtherEasyConfigValues(const int nImages, configInfo *par, imageInfo **img){
  /*
Already set by this time should be all the user-settable parameters, plus:
	nSpecies

Harder things which we will calculate and set later:
	gridDensGlobalMax
	nLineImages
	nContImages
	doInterpolateVels
	doMolCalcs
	useVelFuncInRaytrace
	edgeVelsAvailable
	dataFlags
  */

  int i,j;
  _Bool doThetaPhi;
  double cos_pa,sin_pa,cosPhi,sinPhi,cos_incl,sin_incl,cosTheta,sinTheta,cos_az,sin_az;
  double tempRotMat[3][3],auxRotMat[3][3];
  int row,col;

  i = 0;
  while(i<MAX_N_HIGH && par->gridDensMaxValues[i]>=0) i++;
  par->numGridDensMaxima = i;

  /* Now set the additional values in par. (Note that some of these can be redefined if we read grid points from a file.) */
  par->ncell = par->pIntensity + par->sinkPoints;
  par->radiusSqu = par->radius*par->radius;
  par->minScaleSqu=par->minScale*par->minScale;
  par->doPregrid = (par->pregrid==NULL)?0:1;
  par->nSolveItersDone = 0; /* This can be set to some non-zero value if the user reads in a grid file at dataStageI==5. */
  par->useAbun = 1; /* Can be unset within readOrBuildGrid(). */
  par->dataFlags = 0; /* default */
  par->numDensities = 0; /* default */

  for(i=0;i<NUM_GRID_STAGES;i++){
    if(par->gridOutFiles[i] != NULL)
      par->writeGridAtStage[i] = 1;
  };

  par->nImages = nImages;
  par->numDims = DIM;

  /*
Now we need to calculate the cutoff value used in calcSourceFn(). The issue is to decide between

  y = (1 - exp[-x])/x

or the approximation using the Taylor expansion of exp(-x), which to 3rd order is

  y ~ 1 - x/2 + x^2/6.

The cutoff will be the value of abs(x) for which the error in the exact expression equals the next unused Taylor term, which is x^3/24. This error can be shown to be given for small |x| by epsilon/|x|, where epsilon is the floating-point precision of the computer. Hence the cutoff evaluates to

  |x|_cutoff = (24*epsilon)^{1/4}.

  */
  par->taylorCutoff = pow(24.*DBL_EPSILON, 0.25);

  for(i=0;i<par->nImages;i++){
    (*img)[i].imgres=(*img)[i].imgres*ARCSEC_TO_RAD;
    (*img)[i].pixel = malloc(sizeof(*((*img)[i].pixel))*(*img)[i].pxls*(*img)[i].pxls);
    /*
The image rotation matrix is used within traceray() to transform the coordinates of a vector (actually two vectors - the ray direction and its starting point) as initially specified in the observer-anchored frame into the coordinate frame of the model. In linear algebra terms, the model-frame vector v_mod is related to the vector v_obs as expressed in observer- (or image-) frame coordinates via the image rotation matrix R by

    v_mod = R * v_obs,				1

the multiplication being the usual matrix-vector style. Note that the ith row of R is the ith axis of the model frame with coordinate values expressed in terms of the observer frame.

The matrix R can be broken into a sequence of several (3 at least are needed for full degrees of freedom) simpler rotations. Since these constituent rotations are usually easier to conceive in terms of rotations of the model in the observer framework, it is convenient to invert equation (1) to give

    v_obs = R^T * v_mod,			2

where ^T here denotes transpose. Supposing now we rotate the model in a sequence R_3^T followed by R_2^T followed by R_1^T, equation (2) can be expanded to give

    v_obs = R_1^T * R_2^T * R_3^T * v_mod.	3

Inverting everything to return to the format of equation (1), which is what we need, we find

    v_mod = R_3 * R_2 * R_1 * v_obs.		4

LIME provides two different schemes of {R_1, R_2, R_3}: {PA, phi, theta} and {PA, inclination, azimuth}. As an example, consider phi, which is a rotation of the model from the observer Z axis towards the X. The matching obs->mod rotation matrix is therefore

            ( cos(ph)  0  -sin(ph) )
            (                      )
    R_phi = (    0     0     1     ).
            (                      )
            ( sin(ph)  0   cos(ph) )

    */

    doThetaPhi = (((*img)[i].incl<-900.)||((*img)[i].azimuth<-900.)||((*img)[i].posang<-900.))?1:0;

    if(doThetaPhi){
      /* For the present PA is not implemented for the theta/phi scheme. Thus we just load the identity matrix at present.
      */
      (*img)[i].rotMat[0][0] = 1.0;
      (*img)[i].rotMat[0][1] = 0.0;
      (*img)[i].rotMat[0][2] = 0.0;
      (*img)[i].rotMat[1][0] = 0.0;
      (*img)[i].rotMat[1][1] = 1.0;
      (*img)[i].rotMat[1][2] = 0.0;
      (*img)[i].rotMat[2][0] = 0.0;
      (*img)[i].rotMat[2][1] = 0.0;
      (*img)[i].rotMat[2][2] = 1.0;
    }else{
      /* Load PA rotation matrix R_PA:
      */
      cos_pa   = cos((*img)[i].posang);
      sin_pa   = sin((*img)[i].posang);
      (*img)[i].rotMat[0][0] =  cos_pa;
      (*img)[i].rotMat[0][1] = -sin_pa;
      (*img)[i].rotMat[0][2] =  0.0;
      (*img)[i].rotMat[1][0] =  sin_pa;
      (*img)[i].rotMat[1][1] =  cos_pa;
      (*img)[i].rotMat[1][2] =  0.0;
      (*img)[i].rotMat[2][0] =  0.0;
      (*img)[i].rotMat[2][1] =  0.0;
      (*img)[i].rotMat[2][2] =  1.0;
    }

    if(doThetaPhi){
      /* Load phi rotation matrix R_phi:
      */
      cosPhi   = cos((*img)[i].phi);
      sinPhi   = sin((*img)[i].phi);
      auxRotMat[0][0] =  cosPhi;
      auxRotMat[0][1] =  0.0;
      auxRotMat[0][2] = -sinPhi;
      auxRotMat[1][0] =  0.0;
      auxRotMat[1][1] =  1.0;
      auxRotMat[1][2] =  0.0;
      auxRotMat[2][0] =  sinPhi;
      auxRotMat[2][1] =  0.0;
      auxRotMat[2][2] =  cosPhi;
    }else{
      /* Load inclination rotation matrix R_inc:
      */
      cos_incl = cos((*img)[i].incl + M_PI);
      sin_incl = sin((*img)[i].incl + M_PI);
      auxRotMat[0][0] =  cos_incl;
      auxRotMat[0][1] =  0.0;
      auxRotMat[0][2] = -sin_incl;
      auxRotMat[1][0] =  0.0;
      auxRotMat[1][1] =  1.0;
      auxRotMat[1][2] =  0.0;
      auxRotMat[2][0] =  sin_incl;
      auxRotMat[2][1] =  0.0;
      auxRotMat[2][2] =  cos_incl;
    }

    for(row=0;row<3;row++){
      for(col=0;col<3;col++){
        tempRotMat[row][col] = 0.0;
        for(j=0;j<3;j++)
          tempRotMat[row][col] += auxRotMat[row][j]*(*img)[i].rotMat[j][col];
      }
    }

    if(doThetaPhi){
      /* Load theta rotation matrix R_theta:
      */
      cosTheta = cos((*img)[i].theta);
      sinTheta = sin((*img)[i].theta);
      auxRotMat[0][0] =  1.0;
      auxRotMat[0][1] =  0.0;
      auxRotMat[0][2] =  0.0;
      auxRotMat[1][0] =  0.0;
      auxRotMat[1][1] =  cosTheta;
      auxRotMat[1][2] =  sinTheta;
      auxRotMat[2][0] =  0.0;
      auxRotMat[2][1] = -sinTheta;
      auxRotMat[2][2] =  cosTheta;
    }else{
      /* Load azimuth rotation matrix R_az:
      */
      cos_az   = cos((*img)[i].azimuth + M_PI/2.0);
      sin_az   = sin((*img)[i].azimuth + M_PI/2.0);
      auxRotMat[0][0] =  cos_az;
      auxRotMat[0][1] = -sin_az;
      auxRotMat[0][2] =  0.0;
      auxRotMat[1][0] =  sin_az;
      auxRotMat[1][1] =  cos_az;
      auxRotMat[1][2] =  0.0;
      auxRotMat[2][0] =  0.0;
      auxRotMat[2][1] =  0.0;
      auxRotMat[2][2] =  1.0;
    }

    for(row=0;row<3;row++){
      for(col=0;col<3;col++){
        (*img)[i].rotMat[row][col] = 0.0;
        for(j=0;j<3;j++)
          (*img)[i].rotMat[row][col] += auxRotMat[row][j]*tempRotMat[j][col];
      }
    }
  }
}

/*....................................................................*/
void
parChecks(configInfo *par){
  int id;
  FILE *fp;
  char message[STR_LEN_1];

  /* Check if moldat files exist. */
  for(id=0;id<par->nSpecies;id++){
    if((fp=fopen(par->moldatfile[id], "r"))==NULL){
      if(!silent){
        snprintf(message, STR_LEN_1, "Moldat file %s not found locally - fetching it from LAMDA", par->moldatfile[id]);
        printMessage(message);
      }
      openSocket(par->moldatfile[id]);
    } else {
      checkFirstLineMolDat(fp, par->moldatfile[id]);
      fclose(fp);
    }
  }

  if(par->gridInFile==NULL){
    /* Check that the mandatory parameters now have 'sensible' settings (i.e., that they have been set at all). Raise an exception if not. */
    if (par->radius<=0){
      if(!silent) bail_out("You must define the radius parameter.");
exit(1);
    }
    if (par->minScale<=0){
      if(!silent) bail_out("You must define the minScale parameter.");
exit(1);
    }
    if (par->pIntensity<=0){
      if(!silent) bail_out("You must define the pIntensity parameter.");
exit(1);
    }
    if (par->sinkPoints<=0){
      if(!silent) bail_out("You must define the sinkPoints parameter.");
exit(1);
    }
  }

}

/*....................................................................*/
int
checkUserFunctions(configInfo *par, _Bool checkForSingularities){
  /*
Run through all the user functions and set flags in the global defaultFuncFlags for those which have defaulted. NOTE however that we will not check which of these functions the user has provided until readOrBuildGrid(), because this will depend on the appropriate data being present or not in any grid file we read in. There are two exceptions to this:

	- The velocity() function, because this is not only called within readOrBuildGrid() to calculate velocities at the grid node positions and at sample locations along the edges between grids, but also potentially within raytrace() to calculate velocities along ray paths through cells. Therefore we perform extra tests for the presence of a user-supplied velocity function near the end of the present function.

	- The density() function, because we need to know the number of densities early on, in case we need to call the default gridDensity() function. Thus we test for this below.
  */

  double x,y,z,tempGridDensGlobalMax;
  double *dummyDens = malloc(sizeof(double)*MAX_N_COLL_PART);
  double *dummyAbun = malloc(sizeof(double)*1);
  double *dummyNmol = malloc(sizeof(double)*1);
  double dummyT[2],dummyTurbDop,dummyVel[DIM],dummyB[3],dummyG2d,dummyR[3],dummyNdens;
  int numDensities=0,i,tempNumDensities;

//**** should give them all values because you don't know what some screwy user routine might return.

  /* just to stop compiler warnings because this return value is currently unused. */
//  (void)dummyDens;
  (void)dummyAbun;
  (void)dummyNmol;
  (void)dummyT;
  (void)dummyTurbDop;
  (void)dummyVel;
  (void)dummyB;
  (void)dummyG2d;
  (void)dummyR;
  (void)dummyNdens;

  tempGridDensGlobalMax = par->gridDensGlobalMax;
  tempNumDensities = par->numDensities;
  par->gridDensGlobalMax = 1.0; /* Dummy value, needed in case the default gridDensity() is called. */
  par->numDensities = 0; /* Ditto. */

  if(checkForSingularities){
    x = 0.0;
    y = 0.0;
    z = 0.0;
  }else{
    x = par->minScale;
    y = par->minScale;
    z = par->minScale;
  }

  dummyR[0] = x;
  dummyR[1] = y;
  dummyR[2] = z;

  for(i=0;i<MAX_N_COLL_PART;i++) dummyDens[i] = -1.0; /* We expect that no real function will return such values, so this may be used as an indicator for the number of values returned. */

  density(      x, y, z, dummyDens);
  temperature(  x, y, z, dummyT);
  abundance(    x, y, z, dummyAbun);
  molNumDensity(x, y, z, dummyNmol);
  doppler(      x, y, z, &dummyTurbDop);
  velocity(     x, y, z, dummyVel);
  magfield(     x, y, z, dummyB);
  gasIIdust(    x, y, z, &dummyG2d);

  if(checkForSingularities){ /* The default functions have no singularities so we needn't bother with them. */
    if(!bitIsSet(defaultFuncFlags, USERFUNC_density))       reportInfsAtOrigin(1, dummyDens, "density");
    if(!bitIsSet(defaultFuncFlags, USERFUNC_temperature))   reportInfsAtOrigin(1, dummyT, "temperature");
    if(!bitIsSet(defaultFuncFlags, USERFUNC_abundance))     reportInfsAtOrigin(1, dummyAbun, "abundance");
    if(!bitIsSet(defaultFuncFlags, USERFUNC_molNumDensity)) reportInfsAtOrigin(1, dummyNmol, "molNumDensity");
    if(!bitIsSet(defaultFuncFlags, USERFUNC_doppler))       reportInfsAtOrigin(0, &dummyTurbDop, "doppler");
    if(!bitIsSet(defaultFuncFlags, USERFUNC_velocity))      reportInfsAtOrigin(3, dummyVel, "velocity");
    if(!bitIsSet(defaultFuncFlags, USERFUNC_magfield))      reportInfsAtOrigin(3, dummyB, "magfield");
    if(!bitIsSet(defaultFuncFlags, USERFUNC_gasIIdust))     reportInfsAtOrigin(0, &dummyG2d, "gasIIdust");
  }

  /* Calculate the number of density values returned:
  */
  for(i=0;i<MAX_N_COLL_PART;i++){
    if(dummyDens[i]<0){ /* Note that NaNs test negative for everything, thus even singular-valued density functions will return the right number of density elements set. */
      break;
    }
  }
  numDensities = i;

  dummyNdens = gridDensity(par, dummyR);

  free(dummyNmol);
  free(dummyAbun);
  free(dummyDens);

  par->gridDensGlobalMax = tempGridDensGlobalMax;
  par->numDensities = tempNumDensities;

  return numDensities;
}

/*....................................................................*/
void
calcGridDensGlobalMax(configInfo *par){
  /* In the case that we need to calculate grid point locations, we will need to call the function gridDensity(). The default one is ok, but this (i) needs the user to supply a density function, and (ii) requires par->gridDensGlobalMax etc to be calculated.
  */

  int i,j;
  double tempPointDensity,r[3];
  const gsl_rng_type *ranNumGenType = gsl_rng_ranlxs2;
  gsl_rng *randGen;
  _Bool foundGoodValue;

  if(bitIsSet(defaultFuncFlags, USERFUNC_gridDensity)){
    if(bitIsSet(defaultFuncFlags, USERFUNC_density)){
      if(!silent) bail_out("You need to supply either a density() function or a gridDensity() function which doesn't call density().");
exit(1);
    }

    par->gridDensGlobalMax = 1.0; /* We need some sort of >0 value for par->gridDensGlobalMax before we call the default gridDensity(). */

    /* First try gridDensity() at the origin of coordinates, where the density is often highest:
    */
    for(i=0;i<DIM;i++) r[i] = 0.0;
    tempPointDensity = gridDensity(par, r);

    if(isinf(tempPointDensity) || isnan(tempPointDensity)){
      if(!silent) warning("There is a singularity at the origin of the gridDensity() function.");
    }else if(tempPointDensity<=0.0){
      if(!silent) warning("The gridDensity() function returns zero at the origin.");
    }else if (tempPointDensity>par->gridDensGlobalMax)
      par->gridDensGlobalMax = tempPointDensity;

    if(isinf(tempPointDensity) || isnan(tempPointDensity) || tempPointDensity<=0.0){
      /* Try gridDensity() a little offset from the origin.
      */
      for(i=0;i<DIM;i++) r[i] = par->minScale;
      tempPointDensity = gridDensity(par, r);

      if(!isinf(tempPointDensity) && !isnan(tempPointDensity) && tempPointDensity>0.0)
        par->gridDensGlobalMax = tempPointDensity;

      else{
        /* Hmm ok, let's try a spread of random locations!
        */
        randGen = gsl_rng_alloc(ranNumGenType);	/* Random number generator */
        if(fixRandomSeeds)
          gsl_rng_set(randGen,140978);
        else
          gsl_rng_set(randGen,time(0));

        foundGoodValue = FALSE; /* the default. */
        for(j=0;j<NUM_RAN_DENS;j++){
          for(i=0;i<DIM;i++) r[i] = par->radius*(2.0*gsl_rng_uniform(randGen) - 1.0);
          tempPointDensity = gridDensity(par, r);
          if(!isinf(tempPointDensity) && !isnan(tempPointDensity) && tempPointDensity>0.0){
            foundGoodValue = TRUE;
        break;
          }
        }

        if(foundGoodValue){
          if(tempPointDensity>par->gridDensGlobalMax)
            par->gridDensGlobalMax = tempPointDensity;

        }else if(par->numGridDensMaxima>0){
          /* Test any maxima the user has provided:
          */
          par->gridDensGlobalMax = par->gridDensMaxValues[0];
          for(i=1;i<par->numGridDensMaxima;i++)
            if(par->gridDensMaxValues[i]>par->gridDensGlobalMax) par->gridDensGlobalMax = par->gridDensMaxValues[i];
        }else{
#ifdef KLUDGE_FOR_BAD_DENS
          /* This has been added under protest to cope with modellib's crappy density functions. */
          defaultDensyPower = 1.0;
          par->gridDensGlobalMax = 1.0;
#else
          if(!silent) bail_out("Can't find non-pathological values of the gridDensity() function.");
exit(1);
#endif
        }
      }
    }
  }
}

/*....................................................................*/
void
parseInputWhenIncompleteGridFile(configInfo *par, _Bool checkForSingularities){
  /*
Does a bunch of checks and initializations relating to the density() and gridDensity() functions.
  */
  int numFuncDensities;

  /* Checks which functions the user has supplied and which fall through to default (setting thereby the appropriate flag in defaultFuncFlags). It also checks for singularities and counts how many density values are returned by the density() function.
  */
  numFuncDensities = checkUserFunctions(par, checkForSingularities);

  /* Calculate par->numDensities if we have not read it from file.
  */
  if(!allBitsSet(par->dataFlags, DS_mask_density)){ /* this mask is set only if a density column has been read from file. */
    if(bitIsSet(defaultFuncFlags, USERFUNC_density)){
      if(!silent) bail_out("You need to provide a density() function.");
exit(1);
    }

    par->numDensities = numFuncDensities;
  }

  if(par->numDensities<=0){
    if(!silent) bail_out("No density values returned.");
exit(1);
  }

  if(!anyBitSet(par->dataFlags, DS_mask_x)){ /* This should only happen if we did not read a file. */
    calcGridDensGlobalMax(par);
  }

#ifndef KLUDGE_FOR_BAD_DENS
  if(par->samplingAlgorithm==0)
    defaultDensyPower = DENSITY_POWER;
  else
    defaultDensyPower = TREE_POWER;
#endif
}

/*....................................................................*/
void
parseImagePars(configInfo *par, imageInfo **img){
  _Bool changedInterp;
  int i,j,id;
  char message[STR_LEN_1+1];
  char *pch_sep = " ,:_", *pch, *pch_end, *units_str;

  /* Allocate pixel space and parse image information.
  */
  for(i=0;i<par->nImages;i++){
    (*img)[i].imgunits = NULL;

    /* If user has not supplied a units string then use unit value (default 0) to maintain backwards compatibility */
    if((*img)[i].units == NULL){
      (*img)[i].numunits = 1;
      (*img)[i].imgunits = malloc(sizeof(*(*img)[i].imgunits));
      if((*img)[i].imgunits == NULL){
        if(!silent) bail_out("Error allocating memory for single unit");
      }
      (*img)[i].imgunits[0] = (*img)[i].unit;
    }else{
      /* Otherwise parse image units, populate imgunits array with appropriate image identifiers and track number
       * of units requested */
      copyInparStr((*img)[i].units, &(units_str));
      pch = strtok(units_str, pch_sep);
      j = 0;
      while(pch){
        j++;
        (*img)[i].imgunits = realloc((*img)[i].imgunits, sizeof(*(*img)[i].imgunits)*j);
        if((*img)[i].imgunits == NULL){
          if(!silent) bail_out("Error allocating memory for multiple units");
        }
        (*img)[i].imgunits[j-1] = (int)strtol(pch, &pch_end, 0);
        if(*pch_end){
          snprintf(message, STR_LEN_1, "Image %d: units string contains '%s' which could not be converted to an integer", i, pch_end);
          if(!silent) bail_out(message);
exit(1);
        }
        pch = strtok(NULL, pch_sep);
      }
      (*img)[i].numunits = j;
    } /* End parse of units. */
    
    if((*img)[i].nchan == 0 && (*img)[i].velres<0 ){ /* => user has set neither nchan nor velres. One of the two is required for a line image. */
      /* Assume continuum image */

      if(par->polarization){
        (*img)[i].nchan=3;

        if(!silent && bitIsSet(defaultFuncFlags, USERFUNC_magfield) && !allBitsSet(par->dataFlags, DS_mask_magfield)){
          warning("You need to supply a magfield function for a polarized image.");
        }
      }else
        (*img)[i].nchan=1;

      if((*img)[i].freq<0){
        if(!silent){
          snprintf(message, STR_LEN_1, "Image %d: you must set freq for a continuum image.", i);
          bail_out(message);
        }
exit(1);
      }

      if(!silent && ((*img)[i].trans>-1 || (*img)[i].bandwidth>-1.0)){
        snprintf(message, STR_LEN_1, "Image %d: bandwidth and trans are ignored for a continuum image.", i);
        warning(message);
      }
      (*img)[i].doline=0;

    }else{ /* => user has set one of either nchan or velres, or possibly both. */
      /* Assume line image */

      /*
For a valid line image, the user must set one of the following pairs:
  bandwidth, velres (if they also set nchan, this is overwritten)
  bandwidth, nchan (if they also set velres, this is overwritten)
  velres, nchan (if they also set bandwidth, this is overwritten)

The presence of one of these combinations at least is checked here, although the actual calculation is done in raytrace(), because it depends on moldata info which we have not yet got.
      */
      if((*img)[i].bandwidth > 0 && (*img)[i].velres > 0){
        if(!silent && (*img)[i].nchan > 0){
          snprintf(message, STR_LEN_1, "Image %d: your nchan value will be overwritten.", i);
          warning(message);
        }

      }else if((*img)[i].nchan <= 0 || ((*img)[i].bandwidth <= 0 && (*img)[i].velres <= 0)) {
        if(!silent){
          snprintf(message, STR_LEN_1, "Image %d: insufficient info to calculate nchan, velres and bandwidth.", i);
          bail_out(message);
        }
exit(1);
      }

      /* Check that we have keywords which allow us to calculate the image frequency (if necessary) after reading in the moldata file:
      */
      if((*img)[i].trans>-1){ /* => user has set trans, possibly also freq. */
        if(!silent && (*img)[i].freq > 0){
          snprintf(message, STR_LEN_1, "Image %d: you set trans, so I'm ignoring freq.", i);
          warning(message);
        }

        if((*img)[i].molI < 0){
          if(!silent && par->nSpecies>1){
            snprintf(message, STR_LEN_1, "Image %d: you did not set molI, so I'm assuming the 1st molecule.", i);
            warning(message);
          }
          (*img)[i].molI = 0;
        }
      }else if((*img)[i].freq<0){ /* => user has set neither trans nor freq. */
        if(!silent){
          snprintf(message, STR_LEN_1, "Image %d: you must set either freq or trans (plus optionally molI).", i);
          bail_out(message);
        }
exit(1);
      }/* else => the user has set freq. */

      (*img)[i].doline=1;
    } /* End check of line or continuum. */

    if((*img)[i].imgres<0.0){
      if(!silent){
        snprintf(message, STR_LEN_1, "Image %d: you must set imgres.", i);
        bail_out(message);
      }
exit(1);
    }

    if((*img)[i].pxls<0){
      if(!silent){
        snprintf(message, STR_LEN_1, "Image %d: you must set pxls.", i);
        bail_out(message);
      }
exit(1);
    }

    if((*img)[i].distance<0.0){
      if(!silent){
        snprintf(message, STR_LEN_1, "Image %d: you must set distance.", i);
        bail_out(message);
      }
exit(1);
    }

    for(id=0;id<((*img)[i].pxls*(*img)[i].pxls);id++){
      (*img)[i].pixel[id].intense = malloc(sizeof(double)*(*img)[i].nchan);
      (*img)[i].pixel[id].tau = malloc(sizeof(double)*(*img)[i].nchan);
    }
  }

  par->nLineImages = 0;
  par->nContImages = 0;
  par->doInterpolateVels = 0;
  changedInterp=FALSE;
  for(i=0;i<par->nImages;i++){
    if((*img)[i].doline){

#ifdef IS_PYTHON
      /* This is done because we want to avoid calling the velocity() function within raytrace(). */
      if(par->traceRayAlgorithm==1 && !(*img)[i].doInterpolateVels\
      && !bitIsSet(defaultFuncFlags, USERFUNC_velocity)\
      && par->nThreads>1){
        changedInterp = TRUE;
        (*img)[i].doInterpolateVels = TRUE;
      }
#endif

      if(par->traceRayAlgorithm==1 && !(*img)[i].doInterpolateVels\
      && bitIsSet(defaultFuncFlags, USERFUNC_velocity)){
        if(!silent) bail_out("par->traceRayAlgorithm==1 && !img[i].doInterpolateVels requires you to supply a velocity function.");
exit(1);
      }
      par->nLineImages++;
    }else
      par->nContImages++;

    if((*img)[i].doInterpolateVels)
      par->doInterpolateVels = 1;
  }

  if(!silent && changedInterp)
    warning("You cannot call a python velocity function when multi-threaded. Vels will be interpolated from grid values.");

}

/*....................................................................*/
void
furtherParChecks(configInfo *par){
  FILE *fp;
  char message[STR_LEN_1+1];

  if(par->nContImages>0){
    if(par->dust==NULL){
      if(!silent) bail_out("You must point par->dust to a dust opacity file for a continuum image.");
exit(1);
    }else{
      if((fp=fopen(par->dust, "r"))==NULL){
        if(!silent){
          snprintf(message, STR_LEN_1, "Couldn't open dust opacity data file %s", par->dust);
          bail_out(message);
        }
exit(1);
      }
      fclose(fp);
    }
  }

  par->useVelFuncInRaytrace = FALSE; /* the default */
  if(par->nLineImages>0 && par->traceRayAlgorithm==0){
    if(bitIsSet(defaultFuncFlags, USERFUNC_velocity)){
      if(!silent) warning("No velocity function supplied - raytracing will have lower precision.");
    }else
      par->useVelFuncInRaytrace = TRUE;
  }

#ifdef IS_PYTHON
  if(par->nThreads>1 && par->useVelFuncInRaytrace){
    par->useVelFuncInRaytrace = FALSE;
    if(!silent)
      warning("You cannot call a python velocity function when multi-threaded.");
  }
#endif

  par->edgeVelsAvailable=0; /* default value, this is set within getEdgeVelocities(). */

  if(!silent){
    if(par->lte_only && par->nSolveIters>0)
      warning("Requesting par->nSolveIters>0 will have no effect if LTE calculation is also requested.");

    if(allBitsSet(par->dataFlags, DS_mask_populations) && par->lte_only)
      warning("LTE calculation will overwrite the population values read from file.");
  }

  if(par->nSolveIters>par->nSolveItersDone || par->lte_only) /* To save the user having to set par->doSolveRTE as well as par->nSolveIters>0 or par->lte_only. */
    par->doSolveRTE = TRUE;

  if(allBitsSet(par->dataFlags, DS_mask_populations))
    par->popsHasBeenInit = TRUE;
  else
    par->popsHasBeenInit = FALSE;

  if(!par->popsHasBeenInit && par->doSolveRTE)
    par->needToInitPops = TRUE;
  else
    par->needToInitPops = FALSE;

  if((!par->lte_only && par->doSolveRTE) || par->nLineImages>0)
    par->needToInitSND = TRUE;
  else
    par->needToInitSND = FALSE;

  par->SNDhasBeenInit = FALSE; /* default */

  par->doMolCalcs = par->doSolveRTE || par->nLineImages>0;
  if(par->doMolCalcs && par->moldatfile==NULL){
    if(!silent) bail_out("You must point par->moldatfile to a data file.");
exit(1);
  }

  if(!silent && !par->doMolCalcs && par->init_lte)
    warning("Your choice of par->init_lte will have no effect.");

  if(par->nSpecies>0 && !par->doMolCalcs){
    if(!silent) bail_out("If you want only continuum calculations you must supply zero moldatfiles.");
exit(1);
  }
}

/*....................................................................*/
void gridPopsInit(configInfo *par, molData *md, struct grid *gp){
  int i,id,ilev;

  for(i=0;i<par->nSpecies;i++){
    /* Allocate space for populations etc */
    for(id=0;id<par->ncell; id++){
      free(gp[id].mol[i].pops);
      gp[id].mol[i].pops = malloc(sizeof(double)*md[i].nlev);
      gp[id].mol[i].pops[0] = 1.0;
      for(ilev=1;ilev<md[i].nlev;ilev++)
        gp[id].mol[i].pops[ilev] = 0.0;
    }
  }

  par->popsHasBeenInit = TRUE;
}

/*....................................................................*/
void specNumDensInit(configInfo *par, molData *md, struct grid *gp){
  int i,id,ilev;

  for(i=0;i<par->nSpecies;i++){
    /* Allocate space for populations etc */
    for(id=0;id<par->ncell; id++){
      free(gp[id].mol[i].specNumDens);
      gp[id].mol[i].specNumDens = malloc(sizeof(double)*md[i].nlev);
      gp[id].mol[i].specNumDens[0] = 0.0;
      for(ilev=1;ilev<md[i].nlev;ilev++)
        gp[id].mol[i].specNumDens[ilev] = 0.0;
    }
  }

  par->SNDhasBeenInit = TRUE;
}

