/*
 *  main.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2016 The LIME development team
 *
 */
#include <locale.h>

#include "lime.h"


/* Forward declaration of functions only used in this file */
int initParImg(inputPars *par, image **img);
int main ();



#ifdef FASTEXP
double EXP_TABLE_2D[128][10];
double EXP_TABLE_3D[256][2][10];
/* I've hard-wired the dimensions of these arrays, but it would be better perhaps to declare them as pointers, and calculate the dimensions with the help of the function call:
  calcFastExpRange(FAST_EXP_MAX_TAYLOR, FAST_EXP_NUM_BITS, &numMantissaFields, &lowestExponent, &numExponentsUsed)
*/
#else
double EXP_TABLE_2D[1][1]; // nominal definitions so the fastexp.c module will compile.
double EXP_TABLE_3D[1][1][1];
#endif

int
initParImg(inputPars *par, image **img)
{
  /* Initialize par with default values, allocate space for the
     output fits images, initialize the images with default values,
     and finally call the input() routine from model.c to set both the
     par and image values.
  */

  int i,id,nImages;
  const double defaultAngle=-999.0;

  /* Set 'impossible' default values for mandatory parameters */
  par->radius    = 0;
  par->minScale  = 0;
  par->pIntensity= 0;
  par->sinkPoints= 0;

  /* Set default values for optional parameters */
  par->dust  	    = NULL;
  par->outputfile   = NULL;
  par->binoutputfile= NULL;
  par->gridfile     = NULL;
  par->pregrid      = NULL;
  par->restart      = NULL;

  par->collPartIds  = malloc(sizeof(int)*MAX_N_COLL_PART);
  for(i=0;i<MAX_N_COLL_PART;i++) par->collPartIds[i] = 0;
  par->nMolWeights  = malloc(sizeof(double)*MAX_N_COLL_PART);
  for(i=0;i<MAX_N_COLL_PART;i++) par->nMolWeights[i] = -1.0;
  par->dustWeights  = malloc(sizeof(double)*MAX_N_COLL_PART);
  for(i=0;i<MAX_N_COLL_PART;i++) par->dustWeights[i] = -1.0;

  par->tcmb = 2.728;
  par->lte_only=0;
  par->init_lte=0;
  par->sampling=2;
  par->blend=0;
  par->antialias=1;
  par->polarization=0;
  par->nThreads = NTHREADS;
  par->traceRayAlgorithm=0;

  /* Allocate initial space for molecular data file names */
  par->moldatfile=malloc(sizeof(char *)*MAX_NSPECIES);
  for(id=0;id<MAX_NSPECIES;id++){
    par->moldatfile[id]=NULL;
  }

  /* Allocate initial space for output fits images */
  (*img)=malloc(sizeof(image)*MAX_NIMAGES);
  for(i=0;i<MAX_NIMAGES;i++){
    (*img)[i].filename=NULL;
    (*img)[i].imgunits= NULL;
  }

  /* First call to the user function which sets par, img values. Note that, as far as img is concerned, here we just want to find out how many images the user wants, so we can malloc the array properly. We call input() a second time then to get the actual per-image parameter values.
  */
  input(par, *img);

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

  /* If the user has provided a list of image filenames, the corresponding elements of (*img).filename will be non-NULL. Thus we can deduce the number of images from the number of non-NULL elements. */
  nImages=0;
  while((*img)[nImages].filename!=NULL && nImages<MAX_NIMAGES)
    nImages++;

  if(nImages==0) {
    if(!silent) bail_out("No images defined (or you haven't set the 1st filename).");
    exit(1);
  }

  /* Set img defaults. */
  for(i=0;i<nImages;i++) {
    (*img)[i].source_vel=0.0;
    (*img)[i].phi=0.0;
    (*img)[i].nchan=0;
    (*img)[i].velres=-1.;
    (*img)[i].trans=-1;
    (*img)[i].molI=-1;
    (*img)[i].freq=-1.;
    (*img)[i].bandwidth=-1.;
    (*img)[i].incl    = defaultAngle;
    (*img)[i].azimuth = defaultAngle;
    (*img)[i].posang  = defaultAngle;
    (*img)[i].units = "0";
  }

  /* Second-pass reading of the user-set parameters (this time just to read the par->moldatfile and img stuff). */
  input(par,*img);

  return nImages;
}


void
run(inputPars inpars, image *img)
{
  /* Run LIME with inpars and the output fits files specified.

     This routine may be used as an interface to LIME from external
     programs. In this case, inpars and img must be specified by the
     external program.
  */
  int i,j,nLineImages;
  int initime=time(0);
  int popsdone=0;
  molData*     m = NULL;
  configInfo  par;
  struct grid* g = NULL;
  char fits_filename[100] = "";

  /*Set locale to avoid trouble when reading files*/
  setlocale(LC_ALL, "C");

  if(!silent) greetings();
  if(!silent) screenInfo();

#ifdef FASTEXP
  calcTableEntries(FAST_EXP_MAX_TAYLOR, FAST_EXP_NUM_BITS);
#endif

  parseInput(inpars, &par, &img, &m); /* Sets par.numDensities for !(par.doPregrid || par.restart) */

  if(par.doPregrid)
    {
      gridAlloc(&par,&g);
      predefinedGrid(&par,g); /* Sets par.numDensities */
      checkUserDensWeights(&par); /* Needs par.numDensities */
    }
  else if(par.restart)
    {
      popsin(&par,&g,&m,&popsdone);
    }
  else
    {
      checkUserDensWeights(&par); /* Needs par.numDensities */
      gridAlloc(&par,&g);
      buildGrid(&par,g);
    }

  /* Make all the continuum images, and count the non-continuum images at the same time:
  */
  nLineImages = 0;
  for(i=0;i<par.nImages;i++){
    if(img[i].doline)
      nLineImages++;
    else{
      if(par.restart)
        img[i].trans=0;
      else
        continuumSetup(i,img,m,&par,g);
      raytrace(i,&par,g,m,img);
      for(j=0;j<img[i].numunits;j++) {
        fitsFilename(fits_filename, &par, img, i, j);
        writeFits(fits_filename, &par, m, img, i, j);
        if(!silent) output(fits_filename);
      }
    }
  }

  if(nLineImages>0 && !popsdone)
    levelPops(m,&par,g,&popsdone);

  /* Now make the line images.
  */
  for(i=0;i<par.nImages;i++){
    if(img[i].doline){
      raytrace(i,&par,g,m,img);
      for(j=0;j<img[i].numunits;j++){
        fitsFilename(fits_filename, &par, img, i, j);
        writeFits(fits_filename, &par, m, img, i, j);
        if(!silent) output(fits_filename);
      }
    }
  }
  
  if(!silent) goodnight(initime);

  freeGrid( &par, m, g);
  freeMolData(par.nSpecies, m);
  free(par.moldatfile);
  free(par.collPartIds);
  free(par.nMolWeights);
  free(par.dustWeights);
}

void writeFits(const char *fits_filename, configInfo *par, molData *m, image *img, const int im, const int unit){
  if(img[im].doline == 1 || (img[im].doline==0 && par->polarization)){
    write3Dfits(fits_filename, par, m, img, im, unit);
  }
  else if(img[im].doline == 0){
    write2Dfits(fits_filename, par, m, img, im, unit);
  }
  else{
    if(!silent) bail_out("FITS output unclear");
  }
}

int main () {
  /* Main program for stand-alone LIME */

  inputPars par;
  image	*img = NULL;
  int nImages;

  silent = 0;

  nImages = initParImg(&par, &img);

  run(par, img);

  freeParImg(nImages, &par, img);

  return 0;

}

