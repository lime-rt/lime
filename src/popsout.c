/*
 *  popsout.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015 The LIME development team
 *
 */

#include "lime.h"

void
popsout(inputPars *par, struct grid *gp, molData *md){
  FILE *fp;
  int j,k,l;
  double dens;
  /* int i,mi,c,q=0,best; */
  /* double vel[3],ra[100],rb[100],za[100],zb[100],min; */
  
  if((fp=fopen(par->outputfile, "w"))==NULL){
    if(!silent) bail_out("Error writing output populations file!");
    exit(1);
  }
  fprintf(fp,"# Column definition: x, y, z, H2 density, kinetic gas temperature, molecular abundance, convergence flag, pops_0...pops_n\n");
  for(j=0;j<par->pIntensity;j++){
    dens=0.;
    for(l=0;l<par->collPart;l++) dens+=gp[j].dens[l];
    fprintf(fp,"%e %e %e %e %e %e %d ", gp[j].x[0], gp[j].x[1], gp[j].x[2], dens, gp[j].t[0], gp[j].nmol[0]/dens, gp[j].conv);
    for(k=0;k<md[0].nlev;k++) fprintf(fp,"%e ",gp[j].mol[0].pops[k]);
    fprintf(fp,"\n");
    //fprintf(fp,"%i %lf %lf %lf %lf %lf %lf %lf %lf\n", gp[j].id, gp[j].x[0], gp[j].x[1], gp[j].x[2],  gp[j].dens[0], gp[j].t[0], gp[j].vel[0], gp[j].vel[1], gp[j].vel[2]);
  }
  fclose(fp);
}


void
binpopsout(inputPars *par, struct grid *g, molData *m){
  FILE *fp;
  int i,j;
  
  if((fp=fopen(par->binoutputfile, "wb"))==NULL){
    if(!silent) bail_out("Error writing binary output populations file!");
    exit(1);
  }

  fwrite(&par->radius,   sizeof(double), 1, fp);
  fwrite(&par->ncell,    sizeof(int), 1, fp);
  fwrite(&par->nSpecies, sizeof(int), 1, fp);
  
  for(i=0;i<par->nSpecies;i++){
    fwrite(&m[i].nlev,  sizeof(int),               1,fp);
    fwrite(&m[i].nline, sizeof(int),               1,fp);
    fwrite(&m[i].npart, sizeof(int),               1,fp);
    fwrite(m[i].ntrans, sizeof(int)*m[i].npart,    1,fp);
    fwrite(m[i].lal,    sizeof(int)*m[i].nline,    1,fp);
    fwrite(m[i].lau,    sizeof(int)*m[i].nline,    1,fp);
    fwrite(m[i].aeinst, sizeof(double)*m[i].nline, 1,fp);
    fwrite(m[i].freq,   sizeof(double)*m[i].nline, 1,fp);
    fwrite(m[i].beinstl,sizeof(double)*m[i].nline, 1,fp);
    fwrite(m[i].beinstu,sizeof(double)*m[i].nline, 1,fp);
    fwrite(m[i].local_cmb, sizeof(double)*m[i].nline,1,fp);
    fwrite(&m[i].norm,  sizeof(double),1,fp);
    fwrite(&m[i].norminv,sizeof(double),1,fp);
  }
  
  for(i=0;i<par->ncell;i++){
    fwrite(&g[i].id,   sizeof(int),      1, fp);
    fwrite(&g[i].x,  3*sizeof(double),   1, fp);
    fwrite(&g[i].vel,3*sizeof(double),   1, fp);
    fwrite(&g[i].sink, sizeof(int),      1, fp);
    fwrite(g[i].nmol,  sizeof(double)*par->nSpecies,1, fp);
    fwrite(&g[i].dopb, sizeof g[i].dopb, 1, fp);
    for(j=0;j<par->nSpecies;j++){
      fwrite(g[i].mol[j].pops,  sizeof(double)*m[j].nlev, 1, fp);
      fwrite(g[i].mol[j].knu,   sizeof(double)*m[j].nline,1, fp);
      fwrite(g[i].mol[j].dust,  sizeof(double)*m[j].nline,1, fp);
      fwrite(&g[i].mol[j].dopb, sizeof(double),           1, fp);
      fwrite(&g[i].mol[j].binv, sizeof(double),           1, fp);
    }
    fwrite(&g[i].dens[0], sizeof(double), 1, fp);
    fwrite(&g[i].t[0],    sizeof(double), 1, fp);
    fwrite(&g[i].abun[0], sizeof(double), 1, fp);
 }
  
  
  fclose(fp);

}

  
  
