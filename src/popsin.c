/*
 *  popsin.c
 *  LIME, The versatile 3D line modeling environment
 *
 *  Created by Christian Brinch on 08/26/10.
 *  Copyright 2006-2013, Christian Brinch,
 *  <brinch@nbi.dk>
 *  Niels Bohr institutet
 *  University of Copenhagen
 *	All rights reserved.
 *
 */

#include "lime.h"


void
popsin(inputPars *par, struct grid **g, molData **m, int *popsdone){
  FILE *fp;
  int i,j,k;
  
  if((fp=fopen(par->restart, "rb"))==NULL){
    if(!silent) bail_out("Error writing binary output populations file!");
    exit(1);
  }
  
  fread(&par->radius,   sizeof(double), 1, fp);
  fread(&par->ncell,    sizeof(int), 1, fp);
  fread(&par->nSpecies, sizeof(int), 1, fp);
  
  free(*m);
  *m=malloc(sizeof(molData)*par->nSpecies);
  
  for(i=0;i<par->nSpecies;i++){
    fread(&(*m)[i].nlev,  sizeof(int),        1,fp);
    fread(&(*m)[i].nline, sizeof(int),        1,fp);
    fread(&(*m)[i].npart, sizeof(int),        1,fp);
    (*m)[i].ntrans=malloc(sizeof(int)*(*m)[i].npart);
    for(j=0;j<(*m)[i].npart;j++) fread(&(*m)[i].ntrans[j], sizeof(int), 1,fp);
    (*m)[i].lal=malloc(sizeof(int)*(*m)[i].nline);
    for(j=0;j<(*m)[i].nline;j++) fread(&(*m)[i].lal[j],    sizeof(int), 1,fp);
    (*m)[i].lau=malloc(sizeof(int)*(*m)[i].nline);
    for(j=0;j<(*m)[i].nline;j++) fread(&(*m)[i].lau[j],    sizeof(int), 1,fp);
    (*m)[i].aeinst=malloc(sizeof(double)*(*m)[i].nline);
    for(j=0;j<(*m)[i].nline;j++) fread(&(*m)[i].aeinst[j], sizeof(double), 1,fp);
    (*m)[i].freq=malloc(sizeof(double)*(*m)[i].nline);
    for(j=0;j<(*m)[i].nline;j++) fread(&(*m)[i].freq[j],   sizeof(double), 1,fp);
    (*m)[i].beinstl=malloc(sizeof(double)*(*m)[i].nline);
    for(j=0;j<(*m)[i].nline;j++) fread(&(*m)[i].beinstl[j],sizeof(double), 1,fp);
    (*m)[i].beinstu=malloc(sizeof(double)*(*m)[i].nline);
    for(j=0;j<(*m)[i].nline;j++) fread(&(*m)[i].beinstu[j],sizeof(double), 1,fp);
    (*m)[i].local_cmb = malloc(sizeof(double)*(*m)[i].nline);
    for(j=0;j<(*m)[i].nline;j++) fread(&(*m)[i].local_cmb[j],sizeof(double), 1,fp);
  }
  
  free(*g);
  *g=malloc(sizeof(struct grid)*par->ncell);
  
  for(i=0;i<par->ncell;i++){
    fread(&(*g)[i].id, sizeof (*g)[i].id, 1, fp);
    fread(&(*g)[i].x, sizeof (*g)[i].x, 1, fp);
    fread(&(*g)[i].vel, sizeof (*g)[i].vel, 1, fp);
    fread(&(*g)[i].sink, sizeof (*g)[i].sink, 1, fp);
    (*g)[i].nmol=malloc(par->nSpecies*sizeof(double));
    for(j=0;j<par->nSpecies;j++) fread(&(*g)[i].nmol[j], sizeof(double), 1, fp);
    fread(&(*g)[i].dopb, sizeof (*g)[i].dopb, 1, fp);
    (*g)[i].mol=malloc(par->nSpecies*sizeof(struct populations));
    for(j=0;j<par->nSpecies;j++){
      (*g)[i].mol[j].pops=malloc(sizeof(double)*(*m)[j].nlev);
      for(k=0;k<(*m)[j].nlev;k++) fread(&(*g)[i].mol[j].pops[k], sizeof(double), 1, fp);
      (*g)[i].mol[j].knu=malloc(sizeof(double)*(*m)[j].nline);
      for(k=0;k<(*m)[j].nline;k++) fread(&(*g)[i].mol[j].knu[k], sizeof(double), 1, fp);
      (*g)[i].mol[j].dust=malloc(sizeof(double)*(*m)[j].nline);
      for(k=0;k<(*m)[j].nline;k++) fread(&(*g)[i].mol[j].dust[k],sizeof(double), 1, fp);
      fread(&(*g)[i].mol[j].dopb,sizeof(double), 1, fp);
      fread(&(*g)[i].mol[j].binv,sizeof(double), 1, fp);
    }
  }
  fclose(fp);

  qhull(par, *g);
  distCalc(par, *g);
  getVelosplines(par,*g);
  *popsdone=1;
}

