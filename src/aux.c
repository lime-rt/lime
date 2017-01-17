/*
 *  aux.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2016 The LIME development team
 *
TODO:
  - The test to run photon() etc in levelPops just tests dens[0]. This is a bit sloppy.
 */

#include "lime.h"
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <float.h>


/*....................................................................*/
double
dotProduct3D(const double *vA, const double *vB){
  return vA[0]*vB[0] + vA[1]*vB[1] + vA[2]*vB[2];
}

/*....................................................................*/
void
mallocInputPars(inputPars *par){
  int i;

  par->moldatfile = malloc(sizeof(char *)*MAX_NSPECIES);

  par->collPartIds  = malloc(sizeof(int)*MAX_N_COLL_PART);
  for(i=0;i<MAX_N_COLL_PART;i++) par->collPartIds[i] = 0;
  par->nMolWeights  = malloc(sizeof(double)*MAX_N_COLL_PART);
  for(i=0;i<MAX_N_COLL_PART;i++) par->nMolWeights[i] = -1.0;
  par->dustWeights  = malloc(sizeof(double)*MAX_N_COLL_PART);
  for(i=0;i<MAX_N_COLL_PART;i++) par->dustWeights[i] = -1.0;
}

/*....................................................................*/
void
copyInparStr(const char *inStr, char **outStr){
  if(inStr==NULL || strlen(inStr)<=0 || strlen(inStr)>STR_LEN_0){
    *outStr = NULL;
  }else{
    *outStr = malloc(sizeof(**outStr)*(STR_LEN_0+1));
    strcpy(*outStr, inStr);
  }
}

/*....................................................................*/
void
parseInput(const inputPars inpar, image *inimg, const int nImages, configInfo *par, imageInfo **img, molData **md){
  /*
The parameters visible to the user have now been strictly confined to members of the structs 'inputPars' and 'image', both of which are defined in inpars.h. There are however further internally-set values which is is convenient to bundle together with the user-set ones. At present we have a fairly clunky arrangement in which the user-set values are copied member-by-member from the user-dedicated structs to the generic internal structs 'configInfo' and 'imageInfo'. This is done in the present function, along with some checks and other initialization.
  */

  int i,j,id;
  double BB[3],normBSquared,dens[MAX_N_COLL_PART],r[DIM];
  double dummyVel[DIM];
  FILE *fp;
  char message[80];
  _Bool doThetaPhi;
  double cos_pa,sin_pa,cosPhi,sinPhi,cos_incl,sin_incl,cosTheta,sinTheta,cos_az,sin_az;
  double tempRotMat[3][3],auxRotMat[3][3];
  int row,col;

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

  for(i=0;i<NUM_GRID_STAGES;i++)
    par->writeGridAtStage[i] = 0;
  par->dataFlags = 0;

  /* If the user has provided a list of moldatfile names, the corresponding elements of par->moldatfile will be non-NULL. Thus we can deduce the number of files (species) from the number of non-NULL elements.
  */
  par->nSpecies=0;
  while(par->nSpecies<MAX_NSPECIES && inpar.moldatfile[par->nSpecies]!=NULL && strlen(inpar.moldatfile[par->nSpecies])>0)
    par->nSpecies++;

  /* Copy over the moldatfiles.
  */
  if(par->nSpecies == 0){
    par->nSpecies = 1;/* ******************* should not confuse continuum and line in this way. */
    par->moldatfile = NULL;

  } else {
    par->moldatfile=malloc(sizeof(char *)*par->nSpecies);
    for(id=0;id<par->nSpecies;id++)
      copyInparStr(inpar.moldatfile[id], &(par->moldatfile[id]));

    /* Check if files exist. */
    for(id=0;id<par->nSpecies;id++){
      if((fp=fopen(par->moldatfile[id], "r"))==NULL) {
        openSocket(par->moldatfile[id]);
      } else {
        fclose(fp);
      }
    }
  }

  /* Copy over the collision-partner pointers:
  */
  par->collPartIds  = malloc(sizeof(int)*MAX_N_COLL_PART);
  for(i=0;i<MAX_N_COLL_PART;i++) par->collPartIds[i] = inpar.collPartIds[i];
  par->nMolWeights  = malloc(sizeof(double)*MAX_N_COLL_PART);
  for(i=0;i<MAX_N_COLL_PART;i++) par->nMolWeights[i] = inpar.nMolWeights[i];
  par->dustWeights  = malloc(sizeof(double)*MAX_N_COLL_PART);
  for(i=0;i<MAX_N_COLL_PART;i++) par->dustWeights[i] = inpar.dustWeights[i];

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

  /* Calculate par->numDensities.
  */
  if(!(par->doPregrid || par->restart)){ /* These switches cause par->numDensities to be set in routines they call. */
    /* Find out how many density functions we have (which sets par->numDensities).
    */
    for(i=0;i<MAX_N_COLL_PART;i++) dens[i] = -1.0;
    density(0.0,0.0,0.0,dens); /* Note that the example density function in LIME-1.5 generated a singularity at r==0! Such uglinesses should not be encouraged. I've fixed it now, thus I can use 0s here in (relative) safety. */
    i = 0;
    while(i<MAX_N_COLL_PART && dens[i]>=0) i++;
    par->numDensities = i;

    if(par->numDensities<=0){
      if(!silent) bail_out("No density values returned.");
      exit(1);
    }
  }

  /* See if we can deduce a global maximum for the grid point number density function. Set the starting value from the unnormalized number density at the origin of coordinates:
  */
  par->gridDensGlobalMax = 1.0;
  for(i=0;i<DIM;i++) r[i] = par->minScale;//****0.0;
  par->gridDensGlobalMax = gridDensity(par, r);

  /* Test now any maxima the user has provided:
  */
  for(i=0;i<par->numGridDensMaxima;i++)
    if(par->gridDensMaxValues[i]>par->gridDensGlobalMax) par->gridDensGlobalMax = par->gridDensMaxValues[i];

  if (par->gridDensGlobalMax<=0.0){
    if(!silent) bail_out("Cannot normalize the grid-point number density function.");
    exit(1);
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
  *img = malloc(sizeof(**img)*nImages);
  for(i=0;i<nImages;i++){
    (*img)[i].nchan      = inimg[i].nchan;
    (*img)[i].trans      = inimg[i].trans;
    (*img)[i].molI       = inimg[i].molI;
    (*img)[i].velres     = inimg[i].velres;
    (*img)[i].imgres     = inimg[i].imgres;
    (*img)[i].pxls       = inimg[i].pxls;
    (*img)[i].unit       = inimg[i].unit;
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
  }

  /* Allocate pixel space and parse image information.
  */
  for(i=0;i<nImages;i++){
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

      if(par->moldatfile==NULL){
        if(!silent) bail_out("You must point par->moldatfile to a data file for a line image.");
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
    (*img)[i].pixel = malloc(sizeof(spec)*(*img)[i].pxls*(*img)[i].pxls);
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
      cos_incl = cos((*img)[i].incl + PI);
      sin_incl = sin((*img)[i].incl + PI);
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
      cos_az   = cos((*img)[i].azimuth + PI/2.0);
      sin_az   = sin((*img)[i].azimuth + PI/2.0);
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
  for(i=0;i<nImages;i++){
    if((*img)[i].doline)
      par->nLineImages++;
    else
      par->nContImages++;
  }

  /* Check that the user has supplied the velocity function (needed in raytracing unless par->doPregrid). Note that the other previously mandatory functions (density, abundance, doppler and temperature) may not be necessary if the user reads in the appropriate values from a file. This is tested at the appropriate place in readOrBuildGrid().
  */
  if(par->nLineImages>0 && (!par->doPregrid || par->traceRayAlgorithm==1))
    velocity(0.0,0.0,0.0, dummyVel);

  /* Allocate moldata array.
  */
  (*md)=malloc(sizeof(molData)*par->nSpecies);
  for( i=0; i<par->nSpecies; i++ ){
    (*md)[i].part = NULL;
    (*md)[i].lal = NULL;
    (*md)[i].lau = NULL;
    (*md)[i].aeinst = NULL;
    (*md)[i].freq = NULL;
    (*md)[i].beinstu = NULL;
    (*md)[i].beinstl = NULL;
    (*md)[i].eterm = NULL;
    (*md)[i].gstat = NULL;
    (*md)[i].cmb = NULL;
  }
}

void checkUserDensWeights(configInfo *par){
  /*
This deals with three user-settable vectors: par->collPartIds, par->nMolWeights and par->dustWeights. We have to see if these (optional) parameters were set, do some basic checks on them, and make sure they have the same numbers of elements as the number of density values, which by this time should be stored in par->numDensities.
  */
  int i,j,numUserSetCPIds,numUserSetNMWs,numUserSetDWs;
  int *uniqueCPIds=NULL;
  double sum;

  /* Get the number of par->collPartIds set by the user:
  */
  i = 0;
  while(i<MAX_N_COLL_PART && par->collPartIds[i]>0) i++;
  numUserSetCPIds = i;

  if(numUserSetCPIds>0){
    /* Check that they are unique.
    */
    uniqueCPIds = malloc(sizeof(int)*numUserSetCPIds);
    for(i=0;i<numUserSetCPIds;i++){
      for(j=0;j<i;j++){
        if(par->collPartIds[i]==uniqueCPIds[j]){
          if(!silent) bail_out("Your list of par.collPartIds is not unique.");
          exit(1);
        }
      }
      uniqueCPIds[i] = par->collPartIds[i];
    }
    free(uniqueCPIds);
  }

  /* Get the number of par->nMolWeights set by the user:
  */
  i = 0;
  while(i<MAX_N_COLL_PART && par->nMolWeights[i]>=0.0) i++;
  numUserSetNMWs = i;

  if(numUserSetNMWs>0){
    /* Check that they do not sum to zero.
    */
    sum = 0.0;
    for(i=0;i<numUserSetNMWs;i++){
      sum += par->nMolWeights[i];
    }
    if(sum<=0.0){
      if(!silent) bail_out("At least some of your par.nMolWeights must be non-zero!");
      exit(1);
    }
  }

  /* Get the number of par->dustWeights set by the user:
  */
  i = 0;
  while(i<MAX_N_COLL_PART && par->dustWeights[i]>=0.0) i++;
  numUserSetDWs = i;

  if(numUserSetDWs>0){
    /* Check that they do not sum to zero.
    */
    sum = 0.0;
    for(i=0;i<numUserSetDWs;i++){
      sum += par->dustWeights[i];
    }
    if(sum<=0.0){
      if(!silent) bail_out("At least some of your par.dustWeights must be non-zero!");
      exit(1);
    }
  }

  /* Check if we have either 0 par->collPartIds or the same number as the number of density values.
  */
  if(numUserSetCPIds != par->numDensities){
    free(par->collPartIds);
    par->collPartIds = NULL;
    /* Note that in the present case we will (for a line-emission image) look for the collision partners listed in the moldatfiles and set par->collPartIds from them. For that to happen, we require the number of collision partners found in the files to equal par->numDensities. */

    /* numUserSetCPIds==0 is ok, this just means the user has not set the parameter at all, but for other values we should issue some warnings, because if the user sets any at all, they should set the same number as there are returns from density():
    */
    if(numUserSetCPIds > 0)
      if(!silent) warning("par.collPartIds will be ignored - there should be 1 for each density.");
  }else{
    par->collPartIds = realloc(par->collPartIds, sizeof(*(par->collPartIds))*par->numDensities);
  }

  /* Check if we have either 0 par->nMolWeights or the same number as the number of density values.
  */
  if(numUserSetNMWs != par->numDensities){
    free(par->nMolWeights);
    par->nMolWeights = NULL;
    /* Note that in the present case we will (for a line-emission image) look for the collision partners listed in the moldatfiles and set par->nMolWeights from them. */

    /* numUserSetNMWs==0 is ok, this just means the user has not set the parameter at all, but for other values we should issue some warnings, because if the user sets any at all, they should set the same number as there are returns from density():
    */
    if(numUserSetNMWs > 0)
      if(!silent) warning("par->nMolWeights will be ignored - there should be 1 for each density() return.");
  }else{
    par->nMolWeights = realloc(par->nMolWeights, sizeof(*(par->nMolWeights))*par->numDensities);
  }

  /* Check if we have either 0 par->dustWeights or the same number as the number of density values. Note that the treatment of the dust weights is stricter, since we need knu for the continuum case, in which we may not have access to collision partner information from moldat files.
  */
  if(numUserSetDWs != par->numDensities){
    if(numUserSetDWs == 0){
      /* This is ok, this just means the user has not set the parameter at all. Revert to the previous algorithm, but with a warning, because the previous algorithm is dangerous.
      */
      par->dustWeights = realloc(par->dustWeights, sizeof(*(par->dustWeights))*par->numDensities);
      par->dustWeights[0] = 1.0;
      for(i=1;i<par->numDensities;i++)
        par->dustWeights[i] = 0.0;

      if(!silent) warning("User didn't set par.dustWeights. Using the first density to calculate k_nu.");

    }else{
      if(!silent) bail_out("There must be 1 value of par.dustWeights for each density() return.");
      exit(1);
    }
  }else{
    par->dustWeights = realloc(par->dustWeights, sizeof(*(par->dustWeights))*par->numDensities);
  }
}

float
invSqrt(float x){
  /* The magic Quake(TM) fast inverse square root algorithm   */
  /* Can _only_ be used on 32-bit machine architectures       */
  float xhalf = 0.5f*x;
  int i = *(int*)&x;
  i = 0x5f3759df - (i>>1);
  x = *(float*)&i;
  x = x*(1.5f - xhalf*x*x);
  return x;
}

void checkGridDensities(configInfo *par, struct grid *gp){
  int i;
  static _Bool warningAlreadyIssued=0;
  char errStr[80];

  if(!silent){ /* Warn if any densities too low. */
    i = 0;
    while(i<par->pIntensity && !warningAlreadyIssued){
      if(gp[i].dens[0]<TYPICAL_ISM_DENS){
        warningAlreadyIssued = 1;
        sprintf(errStr, "gp[%d].dens[0] at %.1e is below typical values for the ISM (~%.1e).", i, gp[i].dens[0], TYPICAL_ISM_DENS);
        warning(errStr);
        warning("This could give you convergence problems. NOTE: no further warnings will be issued.");
      }
      i++;
    }
  }
}

void lineBlend(molData *m, configInfo *par, struct blendInfo *blends){
  /*
This obtains information on all the lines of all the radiating species which have other lines within some cutoff velocity separation.

A variable of type 'struct blendInfo' has a nested structure which can be illustrated diagrammaticaly as follows.

  Structs:	blendInfo		molWithBlends		lineWithBlends		blend

  Variables:	blends
		  .numMolsWithBlends     ____________________
		  .*mols--------------->|.molI               |
		                        |.numLinesWithBlends |   ___________
		                        |.*lines--------------->|.lineI     |
		                        |____________________|  |.numBlends |           ________
		                        |        etc         |  |.*blends------------->|.molJ   |
		                                                |___________|          |.lineJ  |
		                                                |    etc    |          |.deltaV |
		                                                                       |________|
		                                                                       |   etc  |

Pointers are indicated by a * before the attribute name and an arrow to the memory location pointed to.
  */
  int molI, lineI, molJ, lineJ;
  int nmwb, nlwb, numBlendsFound, li, bi;
  double deltaV;
  struct blend *tempBlends=NULL;
  struct lineWithBlends *tempLines=NULL;

  /* Dimension blends.mols first to the total number of species, then realloc later if need be.
  */
  (*blends).mols = malloc(sizeof(struct molWithBlends)*par->nSpecies);
  (*blends).numMolsWithBlends = 0;

  nmwb = 0;
  for(molI=0;molI<par->nSpecies;molI++){
    tempBlends = malloc(sizeof(struct blend)*m[molI].nline);
    tempLines  = malloc(sizeof(struct lineWithBlends)*m[molI].nline);

    nlwb = 0;
    for(lineI=0;lineI<m[molI].nline;lineI++){
      numBlendsFound = 0;
      for(molJ=0;molJ<par->nSpecies;molJ++){
        for(lineJ=0;lineJ<m[molJ].nline;lineJ++){
          if(!(molI==molJ && lineI==lineJ)){
            deltaV = (m[molJ].freq[lineJ] - m[molI].freq[lineI])*CLIGHT/m[molI].freq[lineI];
            if(fabs(deltaV)<maxBlendDeltaV){
              tempBlends[numBlendsFound].molJ   = molJ;
              tempBlends[numBlendsFound].lineJ  = lineJ;
              tempBlends[numBlendsFound].deltaV = deltaV;
              numBlendsFound++;
            }
          }
        }
      }

      if(numBlendsFound>0){
        tempLines[nlwb].lineI = lineI;
        tempLines[nlwb].numBlends = numBlendsFound;
        tempLines[nlwb].blends = malloc(sizeof(struct blend)*numBlendsFound);
        for(bi=0;bi<numBlendsFound;bi++)
          tempLines[nlwb].blends[bi] = tempBlends[bi];

        nlwb++;
      }
    }

    if(nlwb>0){
      (*blends).mols[nmwb].molI = molI;
      (*blends).mols[nmwb].numLinesWithBlends = nlwb;
      (*blends).mols[nmwb].lines = malloc(sizeof(struct lineWithBlends)*nlwb);
      for(li=0;li<nlwb;li++){
        (*blends).mols[nmwb].lines[li].lineI     = tempLines[li].lineI;
        (*blends).mols[nmwb].lines[li].numBlends = tempLines[li].numBlends;
        (*blends).mols[nmwb].lines[li].blends = malloc(sizeof(struct blend)*tempLines[li].numBlends);
        for(bi=0;bi<tempLines[li].numBlends;bi++)
          (*blends).mols[nmwb].lines[li].blends[bi] = tempLines[li].blends[bi];
      }

      nmwb++;
    }

    free(tempLines);
    free(tempBlends);
  }

  (*blends).numMolsWithBlends = nmwb;
  if(nmwb>0){
    if(!par->blend)
      if(!silent) warning("There are blended lines, but line blending is switched off.");

    (*blends).mols = realloc((*blends).mols, sizeof(struct molWithBlends)*nmwb);
  }else{
    if(par->blend)
      if(!silent) warning("Line blending is switched on, but no blended lines were found.");

    free((*blends).mols);
    (*blends).mols = NULL;
  }
}

/*....................................................................*/
void
levelPops(molData *md, configInfo *par, struct grid *gp, int *popsdone, double *lamtab, double *kaptab, const int nEntries){
  int id,iter,ilev,ispec,c=0,n,i,threadI,nVerticesDone,nItersDone,nlinetot;
  double percent=0.,*median,result1=0,result2=0,snr,delta_pop;
  int nextMolWithBlend;
  struct statistics { double *pop, *ave, *sigma; } *stat;
  const gsl_rng_type *ranNumGenType = gsl_rng_ranlxs2;
  struct blendInfo blends;
  _Bool luWarningGiven=0;
  gsl_error_handler_t *defaultErrorHandler=NULL;
  int RNG_seeds[par->nThreads];

  nlinetot = 0;
  for(ispec=0;ispec<par->nSpecies;ispec++)
    nlinetot += md[ispec].nline;

  gridPopsInit(par,md,gp);

  if(par->lte_only){
    LTE(par,gp,md);
    if(par->outputfile) popsout(par,gp,md);

  }else{ /* Non-LTE */
    stat=malloc(sizeof(struct statistics)*par->pIntensity);

    /* Random number generator */
    gsl_rng *ran = gsl_rng_alloc(ranNumGenType);
#ifdef TEST
    gsl_rng_set(ran, 1237106) ;
#else 
    gsl_rng_set(ran,time(0));
#endif

    gsl_rng **threadRans;
    threadRans = malloc(sizeof(gsl_rng *)*par->nThreads);

    for (i=0;i<par->nThreads;i++){
      threadRans[i] = gsl_rng_alloc(ranNumGenType);
      if (par->resetRNG==1) RNG_seeds[i] = (int)(gsl_rng_uniform(ran)*1e6);
      else gsl_rng_set(threadRans[i],(int)(gsl_rng_uniform(ran)*1e6));
    }

    calcGridCollRates(par,md,gp);
    calcGridLinesDustOpacity(par, md, lamtab, kaptab, nEntries, gp);

    /* Check for blended lines */
    lineBlend(md, par, &blends);

    if(par->init_lte) LTE(par,gp,md);

    for(id=0;id<par->pIntensity;id++){
      stat[id].pop=malloc(sizeof(double)*md[0].nlev*5);
      stat[id].ave=malloc(sizeof(double)*md[0].nlev);
      stat[id].sigma=malloc(sizeof(double)*md[0].nlev);
      for(ilev=0;ilev<md[0].nlev;ilev++) {
        for(iter=0;iter<5;iter++) stat[id].pop[ilev+md[0].nlev*iter]=gp[id].mol[0].pops[ilev];
      }
    }

    if(par->outputfile) popsout(par,gp,md);

    /* Initialize convergence flag */
    for(id=0;id<par->ncell;id++){
      gp[id].conv=0;
    }

    defaultErrorHandler = gsl_set_error_handler_off();
    /*
This is done to allow proper handling of errors which may arise in the LU solver within stateq(). It is done here because the GSL documentation does not recommend leaving the error handler at the default within multi-threaded code.

While this is off however, other gsl_* etc calls will not exit if they encounter a problem. We may need to pay some attention to trapping their errors.
    */

    nItersDone=0;
    while(nItersDone < par->nSolveIters){ /* Not a 'for' loop because we will probably later want to add a convergence criterion. */
      if(!silent) progressbar2(par, 0, nItersDone, 0, result1, result2);

      for(id=0;id<par->pIntensity;id++){
        for(ilev=0;ilev<md[0].nlev;ilev++) {
          for(iter=0;iter<4;iter++) stat[id].pop[ilev+md[0].nlev*iter]=stat[id].pop[ilev+md[0].nlev*(iter+1)];
          stat[id].pop[ilev+md[0].nlev*4]=gp[id].mol[0].pops[ilev];
        }
      }
      calcGridMolSpecNumDens(par,md,gp);

      nVerticesDone=0;
      omp_set_dynamic(0);
#pragma omp parallel private(id,ispec,threadI,nextMolWithBlend) num_threads(par->nThreads)
      {
        threadI = omp_get_thread_num();

        if (par->resetRNG==1) gsl_rng_set(threadRans[threadI],RNG_seeds[threadI]);
        /* Declare and allocate thread-private variables */
        gridPointData *mp;	// Could have declared them earlier
        double *halfFirstDs;	// and included them in private() I guess.
        mp=malloc(sizeof(gridPointData)*par->nSpecies);

#pragma omp for schedule(dynamic)
        for(id=0;id<par->pIntensity;id++){
#pragma omp atomic
          ++nVerticesDone;

          for (ispec=0;ispec<par->nSpecies;ispec++){
            mp[ispec].jbar = malloc(sizeof(double)*md[ispec].nline);
            mp[ispec].phot = malloc(sizeof(double)*md[ispec].nline*gp[id].nphot);
            mp[ispec].vfac = malloc(sizeof(double)*                gp[id].nphot);
            mp[ispec].vfac_loc = malloc(sizeof(double)*            gp[id].nphot);
          }
          halfFirstDs = malloc(sizeof(*halfFirstDs)*gp[id].nphot);

          if (threadI == 0){ /* i.e., is master thread. */
            if(!silent) progressbar(nVerticesDone/(double)par->pIntensity,10);
          }
          if(gp[id].dens[0] > 0 && gp[id].t[0] > 0){
            photon(id,gp,md,threadRans[threadI],par,nlinetot,blends,mp,halfFirstDs);
            nextMolWithBlend = 0;
            for(ispec=0;ispec<par->nSpecies;ispec++){
              stateq(id,gp,md,ispec,par,blends,nextMolWithBlend,mp,halfFirstDs,&luWarningGiven);
              if(par->blend && blends.mols!=NULL && ispec==blends.mols[nextMolWithBlend].molI)
                nextMolWithBlend++;
            }
          }
          if (threadI == 0){ /* i.e., is master thread */
            if(!silent) warning("");
          }
          freeGridPointData(par->nSpecies, mp);
          free(halfFirstDs);
        }
        free(mp);
      } /* end parallel block. */

      for(id=0;id<par->pIntensity;id++){
        snr=0;
        n=0;
        for(ilev=0;ilev<md[0].nlev;ilev++) {
          stat[id].ave[ilev]=0;
          for(iter=0;iter<5;iter++) stat[id].ave[ilev]+=stat[id].pop[ilev+md[0].nlev*iter];
          stat[id].ave[ilev]=stat[id].ave[ilev]/5.;
          stat[id].sigma[ilev]=0;
          for(iter=0;iter<5;iter++) {
            delta_pop = stat[id].pop[ilev+md[0].nlev*iter]-stat[id].ave[ilev];
            stat[id].sigma[ilev]+=delta_pop*delta_pop;
          }
          stat[id].sigma[ilev]=sqrt(stat[id].sigma[ilev])/5.;
          if(gp[id].mol[0].pops[ilev] > 1e-12) c++;

          if(gp[id].mol[0].pops[ilev] > 1e-12 && stat[id].sigma[ilev] > 0.){
            snr+=gp[id].mol[0].pops[ilev]/stat[id].sigma[ilev];
            n++;
          }
        }
        if(n>0) snr=snr/n;
        else if(n==0) snr=1e6;
        if(snr > 3.) gp[id].conv=2;
        if(snr <= 3 && gp[id].conv==2) gp[id].conv=1;
      }

      median=malloc(sizeof(*median)*gsl_max(c,1));
      c=0;
      for(id=0;id<par->pIntensity;id++){
        for(ilev=0;ilev<md[0].nlev;ilev++){
          if(gp[id].mol[0].pops[ilev] > 1e-12) median[c++]=gp[id].mol[0].pops[ilev]/stat[id].sigma[ilev];
        }
      }

      gsl_sort(median, 1, c);
      if(nItersDone>1){
        result1=median[0];
        result2 =gsl_stats_median_from_sorted_data(median, 1, c);
      }
      free(median);

      if(!silent) progressbar2(par, 1, nItersDone, percent, result1, result2);
      if(par->outputfile != NULL) popsout(par,gp,md);
      nItersDone++;
    }
    gsl_set_error_handler(defaultErrorHandler);

    freeMolsWithBlends(blends.mols, blends.numMolsWithBlends);

    for (i=0;i<par->nThreads;i++){
      gsl_rng_free(threadRans[i]);
    }
    free(threadRans);
    gsl_rng_free(ran);

    for(id=0;id<par->pIntensity;id++){
      free(stat[id].pop);
      free(stat[id].ave);
      free(stat[id].sigma);
    }
    free(stat);
  }

  par->dataFlags |= (1 << DS_bit_populations);

  if(par->binoutputfile != NULL) binpopsout(par,gp,md);

  *popsdone=1;
}

_Bool allBitsSet(const int flags, const int mask){
  /* Returns true only if all the masked bits of flags are set. */

  if(~flags & mask)
    return 0;
  else
    return 1;
}

_Bool anyBitSet(const int flags, const int mask){
  /* Returns true if any of the masked bits of flags are set. */

  if(flags & mask)
    return 1;
  else
    return 0;
}

_Bool bitIsSet(const int flags, const int bitI){
  /* Returns true if the designated bit of flags is set. */

  if(flags & (1 << bitI))
    return 1;
  else
    return 0;
}

_Bool onlyBitsSet(const int flags, const int mask){
  /* Returns true if flags has no bits set apart from those which are true in mask. */

  if(flags & ~mask)
    return 0;
  else
    return 1;
}

