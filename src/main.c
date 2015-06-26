/*
 *  main.c
 *  LIME, The versatile 3D line modeling environment
 *
 *  Created by Christian Brinch on 16/11/06.
 *  Copyright 2006-2014, Christian Brinch,
 *  <brinch@nbi.dk>
 *  Niels Bohr institutet
 *  University of Copenhagen
 *	All rights reserved.
 *
 *  LIME is derived from RATRAN by Michiel Hogerheijde and Floris van der Tak,
 *  Copyright 2000, Hogerheijde and van der Tak, A&A, 362, 697, 2000.
 *
 *  DISCLAIMER: LIME is provided as is and the author does not accept any 
 *  responsibility for erroneous results due to bugs in the code.
 *
 *  Any publication with results obtain using LIME should refer to
 *  Brinch & Hogerheijde, A&A, 523, A25, 2010
 *
 */

#include "lime.h"

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
    writefits(i,&par,m,img);
  }

  if(!silent) goodnight(initime,img[0].filename);

  freeGrid( &par, m, g);
  freeInput(&par, img, m);
  return 0;
}
