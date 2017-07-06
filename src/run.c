/*
 *  run.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2017 The LIME development team
 *
TODO:
 */

#include <locale.h>

#include "lime.h"
#include "gridio.h" /* For countDensityCols() */

int defaultFuncFlags = 0;
double defaultDensyPower = DENSITY_POWER;

/*....................................................................*/
void
reportInfAtOrigin(const double value, const char *funcName){
  char message[STR_LEN_0];

  if(isinf(value) && !silent){
    sprintf(message, "You have a singularity at the origin of your %s() function.", funcName);
    warning(message);
  }
}

/*....................................................................*/
void
reportInfsAtOrigin(const int numElements, const double *values, const char *funcName){
  int i;
  char message[STR_LEN_0];

  for(i=0;i<numElements;i++){
    if(isinf(values[i]) && !silent){
      sprintf(message, "You have a singularity at the origin in return %d of your %s() function.", i, funcName);
      warning(message);
    }
  }
}

/*....................................................................*/
void
parseInput(const inputPars inpar, image *inimg, const int nImages, configInfo *par, imageInfo **img, molData **md){
  /*
The parameters visible to the user have now been strictly confined to members of the structs 'inputPars' and 'image', both of which are defined in inpars.h. There are however further internally-set values which is is convenient to bundle together with the user-set ones. At present we have a fairly clunky arrangement in which the user-set values are copied member-by-member from the user-dedicated structs to the generic internal structs 'configInfo' and 'imageInfo'. This is done in the present function, along with some checks and other initialization.
  */

  int i,j,id,status=0,numGirDatFiles;
  double BB[3],normBSquared,dens[MAX_N_COLL_PART],r[DIM];
  FILE *fp;
  char message[80];
  _Bool doThetaPhi;
  double cos_pa,sin_pa,cosPhi,sinPhi,cos_incl,sin_incl,cosTheta,sinTheta,cos_az,sin_az;
  double tempRotMat[3][3],auxRotMat[3][3];
  int row,col;
  char *pch_sep = " ,:_", *pch, *pch_end, *units_str;

  double *dummyDens = malloc(sizeof(double)*1);
  double *dummyAbun = malloc(sizeof(double)*1);
  double *dummyNmol = malloc(sizeof(double)*1);
  double dummyT[2],dummyTurbDop,dummyVel[DIM],dummyB[3],dummyG2d,dummyR[]={0.0,0.0,0.0},dummyNdens;

  /* just to stop compiler warnings because this return value is currently unused. */
  (void)dummyDens;
  (void)dummyAbun;
  (void)dummyNmol;
  (void)dummyT;
  (void)dummyTurbDop;
  (void)dummyVel;
  (void)dummyB;
  (void)dummyG2d;
  (void)dummyR;
  (void)dummyNdens;

  /* Check that the mandatory parameters now have 'sensible' settings (i.e., that they have been set at all). Raise an exception if not. */
  if (inpar.radius<=0){
    if(!silent) bail_out("You must define the radius parameter.");
exit(1);
  }
  if (inpar.minScale<=0){
    if(!silent) bail_out("You must define the minScale parameter.");
exit(1);
  }
  if (inpar.pIntensity<=0){
    if(!silent) bail_out("You must define the pIntensity parameter.");
exit(1);
  }
  if (inpar.sinkPoints<=0){
    if(!silent) bail_out("You must define the sinkPoints parameter.");
exit(1);
  }

  /* Copy over user-set parameters to the configInfo versions. (This seems like duplicated effort but it is a good principle to separate the two structs, for several reasons, as follows. (i) We will usually want more config parameters than user-settable ones. The separation leaves it clearer which things the user needs to (or can) set. (ii) The separation allows checking and screening out of impossible combinations of parameters. (iii) We can adopt new names (for clarity) for config parameters without bothering the user with a changed interface.)
  */
  par->radius       = inpar.radius;
  par->minScale     = inpar.minScale;
  par->pIntensity   = inpar.pIntensity;
  par->sinkPoints   = inpar.sinkPoints;
  par->samplingAlgorithm=inpar.samplingAlgorithm;
  par->sampling     = inpar.sampling;
  par->tcmb         = inpar.tcmb;
  par->lte_only     = inpar.lte_only;
  par->init_lte     = inpar.init_lte;
  par->blend        = inpar.blend;
  par->antialias    = inpar.antialias;
  par->polarization = inpar.polarization;
  par->nThreads     = inpar.nThreads;
  par->nSolveIters  = inpar.nSolveIters;
  par->traceRayAlgorithm = inpar.traceRayAlgorithm;
  par->resetRNG     = inpar.resetRNG;
  par->doSolveRTE   = inpar.doSolveRTE;

  /* Somewhat more carefully copy over the strings:
  */
  copyInparStr(inpar.dust,          &(par->dust));
  copyInparStr(inpar.outputfile,    &(par->outputfile));
  copyInparStr(inpar.binoutputfile, &(par->binoutputfile));
  copyInparStr(inpar.restart,       &(par->restart));
  copyInparStr(inpar.gridfile,      &(par->gridfile));
  copyInparStr(inpar.pregrid,       &(par->pregrid));
  copyInparStr(inpar.gridInFile,    &(par->gridInFile));

  par->gridOutFiles = malloc(sizeof(char *)*NUM_GRID_STAGES);
  for(i=0;i<NUM_GRID_STAGES;i++)
    copyInparStr(inpar.gridOutFiles[i], &(par->gridOutFiles[i]));

  /* Now set the additional values in par. */
  par->ncell = inpar.pIntensity + inpar.sinkPoints;
  par->radiusSqu = inpar.radius*inpar.radius;
  par->minScaleSqu=inpar.minScale*inpar.minScale;
  par->doPregrid = (inpar.pregrid==NULL)?0:1;
  par->nSolveItersDone = 0; /* This can be set to some non-zero value if the user reads in a grid file at dataStageI==5. */
  par->useAbun = 1; /* Can be unset within readOrBuildGrid(). */

  for(i=0;i<NUM_GRID_STAGES;i++)
    par->writeGridAtStage[i] = 0;
  par->dataFlags = 0;

  if(!par->doPregrid && par->restart){
    par->nSpecies=0; /* This will get set during popsin(). */
    par->girdatfile = NULL;
  }else{
    /* If the user has provided a list of moldatfile names, the corresponding elements of par->moldatfile will be non-NULL. Thus we can deduce the number of files (species) from the number of non-NULL elements.
    */
    par->nSpecies=0;
    while(par->nSpecies<MAX_NSPECIES && inpar.moldatfile[par->nSpecies]!=NULL && strlen(inpar.moldatfile[par->nSpecies])>0)
    par->nSpecies++;

    numGirDatFiles=0;
    while(numGirDatFiles<MAX_NSPECIES && inpar.girdatfile[numGirDatFiles]!=NULL && strlen(inpar.girdatfile[numGirDatFiles])>0)
      numGirDatFiles++;

    if(numGirDatFiles<=0)
      par->girdatfile = NULL;
    else if(numGirDatFiles!=par->nSpecies){
      if(!silent) bail_out("Number of girdatfiles different from number of species.");
exit(1);
    }else{
      par->girdatfile=malloc(sizeof(char *)*par->nSpecies);
      for(id=0;id<par->nSpecies;id++){
        copyInparStr(inpar.girdatfile[id], &(par->girdatfile[id]));
      }
    }
  }

  /* Copy over the moldatfiles.
  */
  if(par->nSpecies <= 0){
    par->moldatfile = NULL; /* This will be tested for all line images, so we can never get par->nLineImages>0 if no moldata files have been supplied. */

  } else {
    par->moldatfile=malloc(sizeof(char *)*par->nSpecies);
    for(id=0;id<par->nSpecies;id++){
      copyInparStr(inpar.moldatfile[id], &(par->moldatfile[id]));
    }

    /* Check if files exist. */
    for(id=0;id<par->nSpecies;id++){
      if((fp=fopen(par->moldatfile[id], "r"))==NULL){
        if(!silent){
          sprintf(message, "Moldat file %s not found locally - fetching it from LAMDA", par->moldatfile[id]);
          printMessage(message);
        }
        openSocket(par->moldatfile[id]);
      } else {
        checkFirstLineMolDat(fp, par->moldatfile[id]);
        fclose(fp);
      }
    }
  }

  /* Copy over the collision-partner pointers:
  */
  par->collPartIds  = malloc(sizeof(int   )*MAX_N_COLL_PART);
  for(i=0;i<MAX_N_COLL_PART;i++) par->collPartIds[i] = inpar.collPartIds[i];
  par->nMolWeights  = malloc(sizeof(double)*MAX_N_COLL_PART);
  for(i=0;i<MAX_N_COLL_PART;i++) par->nMolWeights[i] = inpar.nMolWeights[i];
  par->collPartNames = malloc(sizeof(char*)*MAX_N_COLL_PART);
  for(i=0;i<MAX_N_COLL_PART;i++) copyInparStr(inpar.collPartNames[i], &(par->collPartNames[i]));
  par->collPartMolWeights = malloc(sizeof(double)*MAX_N_COLL_PART);
  for(i=0;i<MAX_N_COLL_PART;i++) par->collPartMolWeights[i] = inpar.collPartMolWeights[i];

  /* Copy over the grid-density-maximum data and find out how many were set:
  */
  par->gridDensMaxValues = malloc(sizeof(*(par->gridDensMaxValues))*MAX_N_HIGH);
  par->gridDensMaxLoc    = malloc(sizeof(*(par->gridDensMaxLoc))*MAX_N_HIGH);
  for(i=0;i<MAX_N_HIGH;i++){
    par->gridDensMaxValues[i] = inpar.gridDensMaxValues[i];
    for(j=0;j<DIM;j++) par->gridDensMaxLoc[i][j] = inpar.gridDensMaxLoc[i][j];
  }
  i = 0;
  while(i<MAX_N_HIGH && par->gridDensMaxValues[i]>=0) i++;
  par->numGridDensMaxima = i;

  /*
Run through all the user functions and set flags in the global defaultFuncFlags for those which have defaulted. NOTE however that we will not check which of these functions the user has provided until readOrBuildGrid(), because this will depend on the appropriate data being present or not in any grid file we read in. There are two exceptions to this:

	- The velocity() function, because this is not only called within readOrBuildGrid() to calculate velocities at the grid node positions and at sample locations along the edges between grids, but also potentially within raytrace() to calculate velocities along ray paths through cells. Therefore we perform extra tests for the presence of a user-supplied velocity function near the end of the present function.

	- The density() function, because we need to know the number of densities early on, in case we need to call the default gridDensity() function. Thus we test for this below.
  */
  density(      0.0,0.0,0.0, dummyDens);
  if(!bitIsSet(defaultFuncFlags, FUNC_BIT_density)) reportInfsAtOrigin(1, dummyDens, "density");
  temperature(  0.0,0.0,0.0, dummyT);
  abundance(    0.0,0.0,0.0, dummyAbun);
  molNumDensity(0.0,0.0,0.0, dummyNmol);
  if(!bitIsSet(defaultFuncFlags, FUNC_BIT_molNumDensity)) reportInfsAtOrigin(1, dummyNmol, "molNumDensity");
  doppler(      0.0,0.0,0.0, &dummyTurbDop);
  velocity(     0.0,0.0,0.0, dummyVel);
  if(!bitIsSet(defaultFuncFlags, FUNC_BIT_velocity)) reportInfsAtOrigin(3, dummyVel, "velocity");
  magfield(     0.0,0.0,0.0, dummyB);
  gasIIdust(    0.0,0.0,0.0, &dummyG2d);

  free(dummyNmol);
  free(dummyAbun);
  free(dummyDens);

  par->gridDensGlobalMax = 1.0; /* Dummy value, needed in case the default gridDensity() is called. */
  par->numDensities = 0; /* Ditto. */
  dummyNdens = gridDensity(par, dummyR);

  /* Calculate par->numDensities.
  */
  if(!(par->doPregrid || par->restart)){ /* These switches cause par->numDensities to be set in routines they call. */
    par->numDensities = 0; /* default. */
    if(par->gridInFile!=NULL){
      status = countDensityCols(par->gridInFile, &(par->numDensities));
      if (status){
        if(!silent){
          printf(message, "countDensityCols() status return %d", status);
          bail_out(message);
        }
exit(1);
      }
    }

    if(par->numDensities<=0){
      if(bitIsSet(defaultFuncFlags, FUNC_BIT_density)){
        if(!silent) bail_out("You need to provide a density() function.");
exit(1);
      }

      /* Find out how many density functions we have (which sets par->numDensities).
      */
      for(i=0;i<MAX_N_COLL_PART;i++) dens[i] = -1.0;
      density(0.0,0.0,0.0,dens);
      /* Testing for a singularity in the density function at the origin is now performed in readOrBuildGrid() */
      i = 0;
      while(i<MAX_N_COLL_PART && dens[i]>=0) i++;
      par->numDensities = i;

      if(par->numDensities<=0){
        if(!silent) bail_out("No density values returned.");
exit(1);
      }
    }
  }

  if(!(par->doPregrid || par->restart || par->gridInFile!=NULL)){
    /* In this case we will need to calculate grid point locations, thus we will need to call the function gridDensity(). The default one is ok, but this (i) needs the user to supply a density function, and (ii) requires par->gridDensGlobalMax etc to be calculated.
    */
    if(bitIsSet(defaultFuncFlags, FUNC_BIT_gridDensity)){
      if(bitIsSet(defaultFuncFlags, FUNC_BIT_density)){
        if(!silent) bail_out("You need to supply either a density() function or a gridDensity() function which doesn't call density().");
exit(1);
      }

      /* See if we can deduce a global maximum for the grid point number density function. Set the starting value from the unnormalized number density at the origin of coordinates:
      */
      par->gridDensGlobalMax = 1.0;
      for(i=0;i<DIM;i++) r[i] = 0.0;
      par->gridDensGlobalMax = gridDensity(par, r);
      if(isinf(par->gridDensGlobalMax)){
        if(!silent) bail_out("There is a singularity at the origin of the gridDensity() function.");
exit(1);
      }else if(par->gridDensGlobalMax<=0.0){
        if(!silent) bail_out("Zero values of the gridDensity() function are not permitted.");
exit(1);
      }

      /* Test now any maxima the user has provided:
      */
      for(i=0;i<par->numGridDensMaxima;i++)
        if(par->gridDensMaxValues[i]>par->gridDensGlobalMax) par->gridDensGlobalMax = par->gridDensMaxValues[i];
    }
  }

  for(i=0;i<NUM_GRID_STAGES;i++){
    if(par->gridOutFiles[i] != NULL)
      par->writeGridAtStage[i] = 1;
  };

  /*
Now we need to calculate the cutoff value used in calcSourceFn(). The issue is to decide between

  y = (1 - exp[-x])/x

or the approximation using the Taylor expansion of exp(-x), which to 3rd order is

  y ~ 1 - x/2 + x^2/6.

The cutoff will be the value of abs(x) for which the error in the exact expression equals the next unused Taylor term, which is x^3/24. This error can be shown to be given for small |x| by epsilon/|x|, where epsilon is the floating-point precision of the computer. Hence the cutoff evaluates to

  |x|_cutoff = (24*epsilon)^{1/4}.

  */
  par->taylorCutoff = pow(24.*DBL_EPSILON, 0.25);
  par->nImages = nImages;
  par->numDims = DIM;

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
      (*img)[i].doInterpolateVels = inimg[i].doInterpolateVels;
    }
  }

  /* Allocate pixel space and parse image information.
  */
  for(i=0;i<nImages;i++){
    (*img)[i].imgunits = NULL;

    /* If user has not supplied a units string then use unit value (default 0) to maintain backwards compatibility */
    if((*img)[i].units == NULL){
      (*img)[i].numunits = 1;
      (*img)[i].imgunits = malloc(sizeof(*(*img)[i].imgunits));
      if((*img)[i].imgunits == NULL){
        if(!silent) bail_out("Error allocating memory for single unit");
      }
      (*img)[i].imgunits[0] = inimg[i].unit;
    }
    else{
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
          sprintf(message, "Units string contains '%s' which could not be converted to an integer", pch_end);
          if(!silent) bail_out(message);
exit(1);
        }
        pch = strtok(NULL, pch_sep);
      }
      (*img)[i].numunits = j;
    }
    
    if((*img)[i].nchan == 0 && (*img)[i].velres<0 ){ /* => user has set neither nchan nor velres. One of the two is required for a line image. */
      /* Assume continuum image */

      if(par->polarization){
        (*img)[i].nchan=3;

        if(!silent){
          /* Do a sketchy check which might indicate if the user has forgotten to supply a magfield function, and warn if this comes up positive. Note: there is no really robust way at present to distinguish the default magfield function (which, if called, indicates that the user forgot to supply their own) from one the user has supplied but which happens to set the B field to 0 at the origin.
          */
          magfield(par->minScale,par->minScale,par->minScale,BB);
          normBSquared = BB[0]*BB[0] + BB[1]*BB[1] + BB[2]*BB[2];
          if(normBSquared <= 0.) warning("Zero B field - did you remember to supply a magfield function?");
        }
      }else
        (*img)[i].nchan=1;

      if((*img)[i].freq<0){
        if(!silent) bail_out("You must set image freq for a continuum image.");
exit(1);
      }

      if(par->dust==NULL){
        if(!silent) bail_out("You must point par.dust to a dust opacity file for a continuum image.");
exit(1);
      }else{
        if((fp=fopen(par->dust, "r"))==NULL){
          if(!silent){
            sprintf(message, "Couldn't open dust opacity data file %s", par->dust);
            bail_out(message);
          }
exit(1);
        }
        fclose(fp);
      }

      if((*img)[i].trans>-1 || (*img)[i].bandwidth>-1.)
        if(!silent) warning("Image bandwidth and trans are ignored for a continuum image.");

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
        if(!silent && (*img)[i].nchan > 0)
          warning("Your nchan value will be overwritten.");

      }else if((*img)[i].nchan <= 0 || ((*img)[i].bandwidth <= 0 && (*img)[i].velres <= 0)) {
        if(!silent) bail_out("Insufficient info to calculate nchan, velres and bandwidth.");
exit(1);
      }

      /* Check that we have keywords which allow us to calculate the image frequency (if necessary) after reading in the moldata file:
      */
      if((*img)[i].trans>-1){ /* => user has set trans, possibly also freq. */
        if(!silent && (*img)[i].freq > 0)
          warning("You set image trans, so I'm ignoring freq.");

        if((*img)[i].molI < 0){
          if(par->nSpecies>1 && !silent) warning("You did not set image molI, so I'm assuming the 1st molecule.");
          (*img)[i].molI = 0;
        }
      }else if((*img)[i].freq<0){ /* => user has set neither trans nor freq. */
        if(!silent) bail_out("You must set either freq or trans (plus optionally molI).");
exit(1);
      }/* else => the user has set freq. */

      (*img)[i].doline=1;
    }
    (*img)[i].imgres=(*img)[i].imgres/206264.806;
    (*img)[i].pixel = malloc(sizeof(*((*img)[i].pixel))*(*img)[i].pxls*(*img)[i].pxls);
    for(id=0;id<((*img)[i].pxls*(*img)[i].pxls);id++){
      (*img)[i].pixel[id].intense = malloc(sizeof(double)*(*img)[i].nchan);
      (*img)[i].pixel[id].tau = malloc(sizeof(double)*(*img)[i].nchan);
    }

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

  par->nLineImages = 0;
  par->nContImages = 0;
  par->doInterpolateVels = 0;
  for(i=0;i<nImages;i++){
    if((*img)[i].doline){
      if(par->traceRayAlgorithm==1 && !(*img)[i].doInterpolateVels && bitIsSet(defaultFuncFlags, FUNC_BIT_velocity)){
        if(!silent) bail_out("par->traceRayAlgorithm==1 && !img[i].doInterpolateVels requires you to supply a velocity function.");
exit(1);
      }
      par->nLineImages++;
    }else
      par->nContImages++;

    if((*img)[i].doInterpolateVels)
      par->doInterpolateVels = 1;
  }

  if(par->nLineImages>0 && par->traceRayAlgorithm==0 && !par->doPregrid){
    if(bitIsSet(defaultFuncFlags, FUNC_BIT_velocity)){
      par->useVelFuncInRaytrace = 0;
      if(!silent) warning("No velocity function supplied - raytracing will have lower precision.");
    }else
      par->useVelFuncInRaytrace = 1;
  }else
    par->useVelFuncInRaytrace = 0;

  par->edgeVelsAvailable=0; /* default value, this is set within getEdgeVelocities(). */

  par->doMolCalcs = par->doSolveRTE || par->nLineImages>0;
  if(par->doMolCalcs && par->moldatfile==NULL && (par->doPregrid || !par->restart)){
    if(!silent) bail_out("You must point par->moldatfile to a data file if you want the RTE solved.");
exit(1);
  }

  if(par->nSpecies>0 && !par->doMolCalcs){
    if(!silent) bail_out("If you want only continuum calculations you must supply zero moldatfiles.");
exit(1);
  }

  if(!silent && !par->doMolCalcs && par->init_lte)
    warning("Your choice of par.init_lte will have no effect.");
  if(!silent && !par->doMolCalcs && par->lte_only)
    warning("Your choice of par.lte_only will have no effect.");
  if(!silent && par->nSolveIters>0 && par->lte_only)
    warning("Requesting par.nSolveIters>0 will have no effect if LTE calculation is also requested.");

  /* Allocate moldata array.
  */
  if(par->nSpecies>0){
    (*md)=malloc(sizeof(molData)*par->nSpecies);
    for( i=0; i<par->nSpecies; i++ ){
      (*md)[i].nlev  = -1;
      (*md)[i].nline = -1;
      (*md)[i].npart = -1;
      (*md)[i].amass = -1.0;
      (*md)[i].part = NULL;
      (*md)[i].lal = NULL;
      (*md)[i].lau = NULL;
      (*md)[i].aeinst = NULL;
      (*md)[i].gir = NULL;
      (*md)[i].freq = NULL;
      (*md)[i].beinstu = NULL;
      (*md)[i].beinstl = NULL;
      (*md)[i].eterm = NULL;
      (*md)[i].gstat = NULL;
      (*md)[i].cmb = NULL;
    }
  } /* otherwise leave it at NULL - we will not be using it. */

  if(par->samplingAlgorithm==0)
    defaultDensyPower = DENSITY_POWER;
  else
    defaultDensyPower = TREE_POWER;
}

/*....................................................................*/
int
run(inputPars inpars, image *inimg, const int nImages){
  /* Run LIME with inpars and the output fits files specified.

     This routine may be used as an interface to LIME from external
     programs. In this case, inpars and img must be specified by the
     external program.
  */
  int i,gi,si,status=0;
  int initime=time(0);
  int popsdone=0;
  molData *md=NULL;
  configInfo par;
  imageInfo *img=NULL;
  struct grid *gp=NULL;
  char message[80];
  int nEntries=0;
  double *lamtab=NULL,*kaptab=NULL;

  /*Set locale to avoid trouble when reading files*/
  setlocale(LC_ALL, "C");

  if(!silent) greetings();
  if(!silent) screenInfo();

#ifdef FASTEXP
  calcExpTableEntries(FAST_EXP_MAX_TAYLOR, FAST_EXP_NUM_BITS);
#endif
  fillErfTable();

  parseInput(inpars, inimg, nImages, &par, &img, &md); /* Sets par.numDensities for !(par.doPregrid || par.restart) */

  if(!silent && par.nThreads>1){
    sprintf(message, "Number of threads used: %d", par.nThreads);
    printMessage(message);
  }

  if(par.doPregrid){
    mallocAndSetDefaultGrid(&gp, (size_t)par.ncell, (size_t)par.nSpecies);
    predefinedGrid(&par,gp); /* Sets par.numDensities */
    checkUserDensWeights(&par); /* Needs par.numDensities */
  }else if(par.restart){
    popsin(&par,&gp,&md,&popsdone);
  }else{
    checkUserDensWeights(&par); /* Needs par.numDensities */
    readOrBuildGrid(&par,&gp);
  }

  if(par.dust != NULL)
    readDustFile(par.dust, &lamtab, &kaptab, &nEntries);

  /* Make all the continuum images:
  */
  if(par.nContImages>0){
    for(i=0;i<par.nImages;i++){
      if(!img[i].doline){
        raytrace(i, &par, gp, md, img, lamtab, kaptab, nEntries);
        writeFitsAllUnits(i, &par, img);
      }
    }
  }

  if(par.doMolCalcs){
    if(!popsdone){
      molInit(&par, md);
      calcGridMolDoppler(&par, md, gp);
    }
    if(par.useAbun)
      calcGridMolDensities(&par, &gp);

    for(gi=0;gi<par.ncell;gi++){
      for(si=0;si<par.nSpecies;si++)
        gp[gi].mol[si].specNumDens = malloc(sizeof(double)*md[si].nlev);
    }

    if(!popsdone && ((par.lte_only && !allBitsSet(par.dataFlags, DS_mask_populations))\
                     || par.nSolveIters>par.nSolveItersDone))
      levelPops(md, &par, gp, &popsdone, lamtab, kaptab, nEntries);

    calcGridMolSpecNumDens(&par,md,gp);

    par.nSolveItersDone = par.nSolveIters;
  }

  if(onlyBitsSet(par.dataFlags & DS_mask_all_but_mag, DS_mask_3))
    writeGridIfRequired(&par, gp, NULL, 3);
  else if(onlyBitsSet(par.dataFlags & DS_mask_all_but_mag, DS_mask_5)){
    writeGridIfRequired(&par, gp, md, 5);
  }else if(!silent){
    sprintf(message, "Data flags %x match neither mask 3 %x (cont.) or 5 %x (line).", par.dataFlags, DS_mask_3, DS_mask_5);
    warning(message);
  }

  freeSomeGridFields((unsigned int)par.ncell, (unsigned short)par.nSpecies, gp);

  /* Now make the line images.   */

  if(par.nLineImages>0){
    for(i=0;i<par.nImages;i++){
      if(img[i].doline){
        raytrace(i, &par, gp, md, img, lamtab, kaptab, nEntries);
        writeFitsAllUnits(i, &par, img);
      }
    }
  }

  if(!silent){
    if(par.nImages>0) reportOutput(img[0].filename);
    goodnight(initime);
  }

  freeGrid((unsigned int)par.ncell, (unsigned short)par.nSpecies, gp);
  freeMolData(par.nSpecies, md);
  freeImgInfo(par.nImages, img);
  freeConfigInfo(&par);

  if(par.dust != NULL){
    free(kaptab);
    free(lamtab);
  }

  return status;
}


