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

void
initParImg(inputPars *par, image **img)
{
  /* Initialize par with default values, allocate space for the
     output fits images, initialize the images with default values,
     and finally call the input routine from model.c to set both the
     par and image values.
  */

  int i, id;
  FILE *fp;

  /* Set default values */
  par->dust  	    = NULL;
  par->inputfile    = NULL;
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
  par->pIntensity=0;
  par->sinkPoints=0;
  par->doPregrid=0;

  par->nThreads = NTHREADS;

  /* Allocate space for output fits images */
  (*img)=malloc(sizeof(image)*MAX_NSPECIES);
  par->moldatfile=malloc(sizeof(char *) * MAX_NSPECIES);
  for(id=0;id<MAX_NSPECIES;id++){
    (*img)[id].filename=NULL;
    par->moldatfile[id]=NULL;
  }
  input(par, *img);
  id=-1;
  while((*img)[++id].filename!=NULL);
  par->nImages=id;
  if(par->nImages==0) {
    if(!silent) bail_out("Error: no images defined");
    exit(1);
  }

  *img=realloc(*img, sizeof(image)*par->nImages);

  id=-1;
  while(par->moldatfile[++id]!=NULL);
  par->nSpecies=id;
  if( par->nSpecies == 0 )
    {
      par->nSpecies = 1;
      free(par->moldatfile);
      par->moldatfile = NULL;
    }
  else
    {
      par->moldatfile=realloc(par->moldatfile, sizeof(char *)*par->nSpecies);
      /* Check if files exists */
      for(id=0;id<par->nSpecies;id++){
        if((fp=fopen(par->moldatfile[id], "r"))==NULL) {
          openSocket(par, id);
        }
        else {
          fclose(fp);
        }
      }
    }

  if(par->dust != NULL){
    if((fp=fopen(par->dust, "r"))==NULL){
      if(!silent) bail_out("Error opening dust opacity data file!");
      exit(1);
    }
    else  {
      fclose(fp);
    }
  }

  /* Set defaults and read inputPars and img[] */
  for(i=0;i<par->nImages;i++) {
    (*img)[i].source_vel=0.0;
    (*img)[i].phi=0.0;
    (*img)[i].nchan=0;
    (*img)[i].velres=-1.;
    (*img)[i].trans=-1;
    (*img)[i].freq=-1.;
    (*img)[i].bandwidth=-1.;
  }
  input(par,*img);
}


void
freeParImg(inputPars *par, image *img)
{
  /* Release memory allocated for the output fits images
     and for par->moldatfile
  */
  int i, id;
  for(i=0;i<par->nImages;i++){
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
run(inputPars *par, image *img)
{
  /* Run LIME with par and the output fits files specified.

     This routine may be used as an interface to LIME from external
     programs. In this case, par and img must be specified by the
     external program.
  */
  int i;
  int initime=time(0);
  int popsdone=0;
  molData*     m = NULL;
  struct grid* g = NULL;

  /*Set locale to avoid trouble when reading files*/
  setlocale(LC_ALL, "C");

  if(!silent) greetings();
  if(!silent) screenInfo();

#ifdef FASTEXP
  calcTableEntries(FAST_EXP_MAX_TAYLOR, FAST_EXP_NUM_BITS);
#endif

  parseInput(par, &img, &m);

  if(par->doPregrid)
    {
      gridAlloc(par,&g);
      predefinedGrid(par,g);
    }
  else if(par->restart)
    {
      popsin(par,&g,&m,&popsdone);
    }
  else
    {
      gridAlloc(par,&g);
      buildGrid(par,g);
    }

  for(i=0;i<par->nImages;i++){
    if(img[i].doline==1 && popsdone==0) {
      levelPops(m,par,g,&popsdone);
    }
    if(img[i].doline==0) {
      continuumSetup(i,img,m,par,g);
    }

    raytrace(i,par,g,m,img);
    writefits(i,par,m,img);
  }

  if(!silent) goodnight(initime,img[0].filename);

  freeGrid( par, m, g);
  freeMoldata(par, m);
}


int main () {
  /* Main program for stand-alone LIME */

  inputPars	par;
  image		*img = NULL;

  silent = 0;

  initParImg(&par, &img);

  run(&par, img);

  freeParImg(&par, img);

  return 0;
}
