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
  int id,ispec;

  for(id=0;id<par->pIntensity;id++){
    for(ispec=0;ispec<par->nSpecies;ispec++){
      g[id].nmol[ispec]=g[id].abun[ispec]*g[id].dens[0];
      lteOnePoint(par, m, ispec, g[id].t[0], g[id].mol[ispec].pops);
    }
  }
  if(par->outputfile) popsout(par,g,m);
}

void lteOnePoint(inputPars *par, molData *m, const int ispec, const double temp, double *pops){
  int ilev;
  double sum;

  sum = 0.0;
  for(ilev=0;ilev<m[ispec].nlev;ilev++){
    pops[ilev] = m[ispec].gstat[ilev]*exp(-100*CLIGHT*HPLANCK*m[ispec].eterm[ilev]/(KBOLTZ*temp));
    sum += pops[ilev];
  }
  for(ilev=0;ilev<m[ispec].nlev;ilev++)
    pops[ilev] /= sum;
}

