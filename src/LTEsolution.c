/*
 *  LTEsolution.c
 *  LIME, The versatile 3D line modeling environment
 *
 *  Created by Christian Brinch on 18/09/09.
 *  Copyright 2006-2014, Christian Brinch,
 *  <brinch@nbi.dk>
 *  Niels Bohr institutet
 *  University of Copenhagen
 *	All rights reserved.
 *
 */
#include "lime.h"

void
LTE(inputPars *par, struct grid *g, molData *m){
  int id,ilev,ispec;
  double z;

  for(ispec=0;ispec<par->nSpecies;ispec++){
    for(id=0;id<par->pIntensity;id++){
      g[id].nmol[ispec]=g[id].abun[ispec]*g[id].dens[0];
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

