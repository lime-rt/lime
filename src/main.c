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
  molData*     m = NULL;
  inputPars    par;
  struct grid* g = NULL;
  image*       img = NULL;

  if(!silent) greetings();
  if(!silent) screenInfo();

#ifdef FASTEXP
  calcTableEntries(FAST_EXP_MAX_TAYLOR, FAST_EXP_NUM_BITS);
#endif

  parseInput(&par,&img,&m); /* Sets par.numDensities for !(par.doPregrid || par.restart) */

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
  freeMolData(&par, m);
  freeInput(&par, img);
  return 0;
}
