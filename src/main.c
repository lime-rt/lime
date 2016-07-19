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
void freeParImg(const int nImages, inputPars *par, image *img);
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
  FILE *fp;

  /* Set 'impossible' default values for mandatory parameters */
  par->radius    =-1;
  par->minScale  =-1;
  par->pIntensity=-1;
  par->sinkPoints=-1;

  /* Set default values for optional parameters */
  par->dust  	    = NULL;
  par->outputfile   = NULL;
  par->binoutputfile= NULL;
  par->gridfile     = NULL;
  par->pregrid      = NULL;
  par->restart      = NULL;

  par->tcmb = 2.728;
  par->lte_only=0;
  par->init_lte=0;
  par->sampling=2;
  par->blend=0;
  par->antialias=1;
  par->polarization=0;
  par->nThreads = NTHREADS;

  /* Allocate initial space for molecular data file names */
  par->moldatfile=malloc(sizeof(char *)*MAX_NSPECIES);
  for(id=0;id<MAX_NSPECIES;id++){
    par->moldatfile[id]=NULL;
  }

  /* Allocate initial space for output fits images */
  (*img)=malloc(sizeof(image)*MAX_NIMAGES);
  for(i=0;i<MAX_NIMAGES;i++){
    (*img)[i].filename=NULL;
  }

  /* First-pass reading of the user-set parameters */
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
    (*img)[i].freq=-1.;
    (*img)[i].bandwidth=-1.;
  }

  /* Second-pass reading of the user-set parameters (this time just to read the par->moldatfile and img stuff). */
  input(par,*img);

  return nImages;
}


void
freeParImg(const int nImages, inputPars *par, image *img)
{
  /* Release memory allocated for the output fits images
     and for par->moldatfile
  */
  int i, id;
  for(i=0;i<nImages;i++){
    for(id=0;id<(img[i].pxls*img[i].pxls);id++){
      free( img[i].pixel[id].intense );
      free( img[i].pixel[id].tau );
    }
    free(img[i].pixel);
  }
  if( img != NULL )
    {
      free(img);
    }
  if( par->moldatfile != NULL )
    {
      free(par->moldatfile);
    }
}


void
run(inputPars inpars, image *img)
{
  /* Run LIME with inpars and the output fits files specified.

     This routine may be used as an interface to LIME from external
     programs. In this case, inpars and img must be specified by the
     external program.
  */
  int i,nLineImages;
  int initime=time(0);
  int popsdone=0;
  molData*     m = NULL;
  configInfo  par;
  struct grid* g = NULL;

  /*Set locale to avoid trouble when reading files*/
  setlocale(LC_ALL, "C");

  if(!silent) greetings();
  if(!silent) screenInfo();

#ifdef FASTEXP
  calcTableEntries(FAST_EXP_MAX_TAYLOR, FAST_EXP_NUM_BITS);
#endif

  parseInput(inpars, &par, &img, &m);

  if(par.doPregrid)
    {
      gridAlloc(&par,&g);
      predefinedGrid(&par,g);
    }
  else if(par.restart)
    {
      popsin(&par,&g,&m,&popsdone);
    }
  else
    {
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
      continuumSetup(i,img,m,&par,g);
      raytrace(i,&par,g,m,img);
      writefits(i,&par,m,img);
    }
  }

  if(nLineImages>0 && !popsdone)
    levelPops(m,&par,g,&popsdone);

  /* Now make the line images.
  */
  for(i=0;i<par.nImages;i++){
    if(img[i].doline){
      raytrace(i,&par,g,m,img);
      writefits(i,&par,m,img);
    }
  }
  
  if(!silent) goodnight(initime,img[0].filename);

  freeGrid( &par, m, g);
  freeMoldata(par.nSpecies, m);
  free(par.moldatfile);
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
