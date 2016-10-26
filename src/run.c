/*
 *  run.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2016 The LIME development team
 *
 */

#include <locale.h>

#include "lime.h"

int
run(inputPars inpars, image *inimg, const int nImages)
{
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
  double *lamtab=NULL, *kaptab=NULL; 

  /*Set locale to avoid trouble when reading files*/
  setlocale(LC_ALL, "C");

  if(!silent) greetings();
  if(!silent) screenInfo();

#ifdef FASTEXP
  calcTableEntries(FAST_EXP_MAX_TAYLOR, FAST_EXP_NUM_BITS);
#endif

  par.nImages = nImages;
  parseInput(inpars, inimg, &par, &img, &md); /* Sets par.numDensities for !(par.doPregrid || par.restart) */

  if(!silent && par.nThreads>1){
    sprintf(message, "Number of threads used: %d", par.nThreads);
    printMessage(message);
  }

  if(par.doPregrid)
    {
      mallocAndSetDefaultGrid(&gp, (unsigned int)par.ncell);
      predefinedGrid(&par,gp); /* Sets par.numDensities */
      checkUserDensWeights(&par); /* Needs par.numDensities */
    }
  else if(par.restart)
    {
      popsin(&par,&gp,&md,&popsdone);
    }
  else
    {
      checkUserDensWeights(&par); /* Needs par.numDensities */
      mallocAndSetDefaultGrid(&gp, (unsigned int)par.ncell);
      buildGrid(&par,gp);
    }

  if(par.dust != NULL)
    readDustFile(par.dust, &lamtab, &kaptab, &nEntries);

  /* Make all the continuum images:
  */
  for(i=0;i<par.nImages;i++){
    if(!img[i].doline){
      raytrace(i, &par, gp, md, img, lamtab, kaptab, nEntries);
      writeFits(i,&par,md,img);
    }
  }

  if(par.nLineImages>0){
    molInit(&par, md);

    if(!popsdone){
      for(gi=0;gi<par.ncell;gi++){
        gp[gi].mol = malloc(sizeof(*(gp[gi].mol))*par.nSpecies);
        for(si=0;si<par.nSpecies;si++){
          gp[gi].mol[si].pops    = NULL;
          gp[gi].mol[si].partner = NULL;
          gp[gi].mol[si].cont    = NULL;
        }
      }
    }

    for(gi=0;gi<par.ncell;gi++){
      for(si=0;si<par.nSpecies;si++)
        gp[gi].mol[si].specNumDens = malloc(sizeof(double)*md[si].nlev);
    }
    calcGridMolDoppler(&par, md, gp);
    calcGridMolDensities(&par,gp);

    if(!popsdone)
      levelPops(md, &par, gp, &popsdone, lamtab, kaptab, nEntries);

    calcGridMolSpecNumDens(&par,md,gp);
  }

  freeSomeGridFields((unsigned int)par.ncell, (unsigned short)par.nSpecies, gp);

  /* Now make the line images.
  */
  for(i=0;i<par.nImages;i++){
    if(img[i].doline){
      raytrace(i, &par, gp, md, img, lamtab, kaptab, nEntries);
      writeFits(i,&par,md,img);
    }
  }
  
  if(!silent) goodnight(initime,img[0].filename);

  freeGrid((unsigned int)par.ncell, (unsigned short)par.nSpecies, gp);
  freeMolData(par.nSpecies, md);
  freeImg(par.nImages, img);
  freeConfigInfo(par);

  if(par.dust != NULL){
    free(kaptab);
    free(lamtab);
  }

  return status; /* This is a bit of a placeholder for now. Ideally we would like all the functions called to return status values rather than exiting. This would allow python-calling versions of Lime to exit 'nicely' at the top level. */
}


