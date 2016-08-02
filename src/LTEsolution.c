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
LTE(configInfo *par, struct grid *gp, molData *m){
  int id,ispec;

  for(id=0;id<par->pIntensity;id++){
    for(ispec=0;ispec<par->nSpecies;ispec++){
      gp[id].mol[ispec].nmol = gp[id].abun[ispec]*gp[id].dens[0];
      lteOnePoint(par, m, ispec, gp[id].t[0], gp[id].mol[ispec].pops);
    }
  }
  if(par->outputfile) popsout(par,gp,m);
}

void lteOnePoint(configInfo *par, molData *m, const int ispec, const double temp, double *pops){
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

