/*
 *  main.c
 *  LIME, The versatile 3D line modeling environment 
 *
 *  Created by Christian Brinch on 16/11/06.
 *  Copyright 2006-2012, Christian Brinch, 
 *  <brinch@nbi.dk>
 *  Niels Bohr institutet
 *  University of Copenhagen
 *	All rights reserved.
 *
 *  LIME is derived from RATRAN by Michiel Hogerheijde and Floris van der Tak,
 *  Copyright 2000, Hogerheijde and van der Tak, A&A, 362, 697, 2000.
 *
 *  DISCLAIMER: LIME is provided as is and the author does not resume any 
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
  molData	  *m;			
  inputPars	  par;
  struct grid *g;
  image		  *img;

  par.dust  	  = NULL;
  par.inputfile = NULL;
  par.outputfile= NULL;
  par.gridfile  = NULL;
  par.pregrid	  = NULL;

  if(!silent) greetings();
  if(!silent) screenInfo();
	
  parseInput(&par,&img,&m);
  gridAlloc(&par,&g);

  if(par.pregrid) predefinedGrid(&par,g); 
  else buildGrid(&par,g);

  for(i=0;i<par.nImages;i++){
    if(img[i].doline==1 && popsdone==0) {
      levelPops(m,&par,g);
      popsdone=1;
    }
    if(img[i].doline==0) {
      continuumSetup(i,img,m,&par,g);
    }
		
    raytrace(i,&par,g,m,img);			  	
    writefits(i,&par,m,img);
  }
  if(!silent) goodnight(initime,img[0].filename);
  return 0;
}
