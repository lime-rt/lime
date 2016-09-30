/*
 *  LTEsolution.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2016 The LIME development team
 *
 */

#include "lime.h"

void
LTE(configInfo *par, struct grid *g, molData *m){
  int id,ispec;

  for(id=0;id<par->pIntensity;id++){
    for(ispec=0;ispec<par->nSpecies;ispec++){
      g[id].mol[ispec].nmol = g[id].abun[ispec]*g[id].dens[0];
      lteOnePoint(m, ispec, g[id].t[0], g[id].mol[ispec].pops);
    }
  }
  if(par->outputfile) popsout(par,g,m);
}

void lteOnePoint(molData *m, const int ispec, const double temp, double *pops){
  int ilev;
  double sum;

  sum = 0.0;
  for(ilev=0;ilev<m[ispec].nlev;ilev++){
    pops[ilev] = m[ispec].gstat[ilev]*exp(-HCKB*m[ispec].eterm[ilev]/temp);
    sum += pops[ilev];
  }
  for(ilev=0;ilev<m[ispec].nlev;ilev++)
    pops[ilev] /= sum;
}

