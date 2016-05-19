/*
 *  LTEsolution.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015 The LIME development team
 *
 */

#include "lime.h"

void
LTE(inputPars *par, struct grid *g, molData *m){
  int id,ilev,ispec;
  double z;

  for(ispec=0;ispec<par->nSpecies;ispec++){
    for(id=0;id<par->pIntensity;id++){
      z=0;
      for(ilev=0;ilev<m[ispec].nlev;ilev++){
        z+=m[ispec].gstat[ilev]*exp(-100*CLIGHT*HPLANCK*m[ispec].eterm[ilev]/(KBOLTZ*g[id].t[0]));
      }
      for(ilev=0;ilev<m[ispec].nlev;ilev++){
        g[id].mol[ispec].pops[ilev]=m[ispec].gstat[ilev]*exp(-100*CLIGHT*HPLANCK*m[ispec].eterm[ilev]/(KBOLTZ*g[id].t[0]))/z;
      }
    }
  }
  if(par->outputfile) popsout(par,g,m);
}

