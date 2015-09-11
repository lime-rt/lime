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
  int i;
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

  parseInput(&par,&img,&m);

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

  for(i=0;i<par.nImages;i++){
    if(img[i].doline==1 && popsdone==0) {
      levelPops(m,&par,g,&popsdone);
    }
    if(img[i].doline==0) {
      continuumSetup(i,img,m,&par,g);
    }

    raytrace(i,&par,g,m,img);

    if(img[i].unit<5)
      write3Dfits(i,&par,m,img);
    else if(img[i].unit==5)
      write2Dfits(i,&par,m,img);
    else{
      if(!silent) bail_out("Image unit number invalid");
      exit(0);
    }
  }

  if(!silent) goodnight(initime,img[0].filename);

  freeGrid( &par, m, g);
  freeInput(&par, img, m);
  return 0;
}
