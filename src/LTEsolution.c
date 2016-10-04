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
LTE(configInfo *par, struct grid *gp, molData *md){
  int id,ispec;

  for(id=0;id<par->pIntensity;id++){
    for(ispec=0;ispec<par->nSpecies;ispec++){
      lteOnePoint(md, ispec, gp[id].t[0], gp[id].mol[ispec].pops);
    }
  }
  if(par->outputfile) popsout(par,gp,md);
}

void lteOnePoint(molData *md, const int ispec, const double temp, double *pops){
  int ilev;
  double sum;

  sum = 0.0;
  for(ilev=0;ilev<md[ispec].nlev;ilev++){
    pops[ilev] = md[ispec].gstat[ilev]*exp(-HCKB*md[ispec].eterm[ilev]/temp);
    sum += pops[ilev];
  }
  for(ilev=0;ilev<md[ispec].nlev;ilev++)
    pops[ilev] /= sum;
}

