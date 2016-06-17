/*
 *  main.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015 The LIME development team
 *
 */

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

int main () {
  int i,nLineImages;
  int initime=time(0);
  int popsdone=0;
  molData*     md = NULL;
  inputPars    par;
  struct grid* gp = NULL;
  image*       img = NULL;
  char message[80];

  if(!silent) greetings();
  if(!silent) screenInfo();

#ifdef FASTEXP
  calcTableEntries(FAST_EXP_MAX_TAYLOR, FAST_EXP_NUM_BITS);
#endif

  parseInput(&par,&img,&md);

  if(!silent && par.nThreads>1){
    sprintf(message, "Number of threads used: %d", par.nThreads);
    printMessage(message);
  }

  if(par.doPregrid)
    {
      gridAlloc(&par,&gp);
      predefinedGrid(&par,gp);
      par.dataStageI = 3; // Sort of.
    }
  else if(par.restart)
    {
      popsin(&par,&gp,&md,&popsdone);
      par.dataStageI = 4; // Sort of.
    }
  else
    {
      readOrBuildGrid(&par,&gp);
    }

  /* Make all the continuum images, and count the non-continuum images at the same time:
  */
  nLineImages = 0;
  for(i=0;i<par.nImages;i++){
    if(img[i].doline)
      nLineImages++;
    else{
      continuumSetup(i,img,md,&par,gp);
      raytrace(i,&par,gp,md,img);
      writefits(i,&par,md,img);
    }
  }

  if(nLineImages>0 && !popsdone){ // eventually, replace !popdone by (!popsdone || dataStageI<4)? *Really* eventually we want to get rid of popsdone.
    levelPops(md,&par,gp,&popsdone);
    par.dataStageI = 4;
/* Disable the next lines for now, since we have not tested dataStageI<4 in the 'if' of this block, because we can't use an input grid file at dataStageI==4 yet: we have to disentangle all the functionality of molinit() before we can contemplate doing that. 
  }else if(par.dataStageI==4 && par->nSolveIters>0 && par.writeGridAtStage[par.dataStageI-1]){
    sprintf(message, "You just read a grid file at data stage %d, now you want to write it again?", par.dataStageI);
    if(!silent) warning(message);
*/
  }

  if(par.dataStageI==4)
    writeGridIfRequired(&par, gp, md, lime_FITS);

  /* Now make the line images.
  */
  for(i=0;i<par.nImages;i++){
    if(img[i].doline){
      raytrace(i,&par,gp,md,img);
      writefits(i,&par,md,img);
    }
  }

  if(!silent) goodnight(initime,img[0].filename);

  freeGrid( &par, md, gp);
  freeInput(&par, img, md);
  return 0;
}
