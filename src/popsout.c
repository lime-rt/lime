/*
 *  popsout.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
***TODO:
	- Change the definition of the file format so that nmol is now read with the other mol[] scalars.
 */

#include "lime.h"

void
popsout(configInfo *par, struct grid *gp, molData *md){
  FILE *fp;
  int j,k,l;
  double dens;
  /* int i,mi,c,q=0,best; */
  /* double vel[3],ra[100],rb[100],za[100],zb[100],min; */

  if((fp=fopen(par->outputfile, "w"))==NULL){
    if(!silent) bail_out("Error writing output populations file!");
    exit(1);
  }
  fprintf(fp,"# x y z H2_density kinetic_gas_temperature molecular_abundance convergence_flag");
  for(k=0;k<md[0].nlev;k++) fprintf(fp," pops_%d",k);
  fprintf(fp,"\n");
  for(j=0;j<par->pIntensity;j++){
    dens=0.;
    for(l=0;l<par->numDensities;l++) dens+=gp[j].dens[l]*par->nMolWeights[l];
    fprintf(fp,"%e %e %e %e %e %e %d ", gp[j].x[0], gp[j].x[1], gp[j].x[2], dens, gp[j].t[0], gp[j].mol[0].nmol/dens, gp[j].conv);
    for(k=0;k<md[0].nlev;k++) fprintf(fp,"%e ",gp[j].mol[0].pops[k]);
    fprintf(fp,"\n");
    //fprintf(fp,"%i %lf %lf %lf %lf %lf %lf %lf %lf\n", gp[j].id, gp[j].x[0], gp[j].x[1], gp[j].x[2],  gp[j].dens[0], gp[j].t[0], gp[j].vel[0], gp[j].vel[1], gp[j].vel[2]);
  }
  fclose(fp);
}

void
binpopsout(configInfo *par, struct grid *gp, molData *md){
  FILE *fp;
  int i,j,k,*nTrans=NULL;
  double dummy=-1.0;
  struct oldPop {
    double *dust, *knu;
  } *dummyMol=NULL;
  double **dummy2=NULL;

  nTrans = malloc(sizeof(int)*1);

  dummyMol = malloc(sizeof(*dummyMol)*par->nSpecies);
  dummy2 = malloc(sizeof(*dummy2)*par->nSpecies);
  for(j=0;j<par->nSpecies;j++){
    dummyMol[j].dust = malloc(sizeof(double)*md[j].nline);
    dummyMol[j].knu  = malloc(sizeof(double)*md[j].nline);
    dummy2[  j]      = malloc(sizeof(double)*md[j].nline);
    for(k=0;k<md[j].nline;k++){
      dummyMol[j].dust[k] = 0.0;
      dummyMol[j].knu[ k] = 0.0;
      dummy2[  j][     k] = 0.0;
    } 
  }

  if((fp=fopen(par->binoutputfile, "wb"))==NULL){
    if(!silent) bail_out("Error writing binary output populations file!");
    exit(1);
  }

  checkFwrite(fwrite(&par->radius,   sizeof(double), 1, fp), 1, "radius");
  checkFwrite(fwrite(&par->ncell,    sizeof(int),    1, fp), 1, "ncell");
  checkFwrite(fwrite(&par->nSpecies, sizeof(int),    1, fp), 1, "nSpecies");

  for(i=0;i<par->nSpecies;i++){
    if(md[i].part==NULL)
      nTrans[0] = 1;
    else
      nTrans[0] = md[i].part[0].ntrans;

    checkFwrite(fwrite(&md[i].nlev,  sizeof(int),                1,fp), 1, "nlev");
    checkFwrite(fwrite(&md[i].nline, sizeof(int),                1,fp), 1, "nline");
    checkFwrite(fwrite(&md[i].npart, sizeof(int),                1,fp), 1, "npart");
    checkFwrite(fwrite(nTrans,       sizeof(int),                1,fp), 1, "nTrans");
    checkFwrite(fwrite(md[i].lal,    sizeof(int)   *md[i].nline, 1,fp), 1, "lal");
    checkFwrite(fwrite(md[i].lau,    sizeof(int)   *md[i].nline, 1,fp), 1, "lau");
    checkFwrite(fwrite(md[i].aeinst, sizeof(double)*md[i].nline, 1,fp), 1, "aeinst");
    checkFwrite(fwrite(md[i].freq,   sizeof(double)*md[i].nline, 1,fp), 1, "freq");
    checkFwrite(fwrite(md[i].beinstl,sizeof(double)*md[i].nline, 1,fp), 1, "beinstl");
    checkFwrite(fwrite(md[i].beinstu,sizeof(double)*md[i].nline, 1,fp), 1, "beinstu");
    checkFwrite(fwrite(dummy2[i],    sizeof(double)*md[i].nline, 1,fp), 1, "dummy2");
    checkFwrite(fwrite(&dummy,       sizeof(double),             1,fp), 1, "dummy");
    checkFwrite(fwrite(&dummy,       sizeof(double),             1,fp), 1, "dummy");
  }

  for(i=0;i<par->ncell;i++){
    checkFwrite(fwrite(&gp[i].id,     sizeof(int),      1, fp), 1, "id");
    checkFwrite(fwrite(&gp[i].x,  DIM*sizeof(double),   1, fp), 1, "x");
    checkFwrite(fwrite(&gp[i].vel,DIM*sizeof(double),   1, fp), 1, "vel");
    checkFwrite(fwrite(&gp[i].sink,   sizeof(int),      1, fp), 1, "sink");
    for(j=0;j<par->nSpecies;j++)
      checkFwrite(fwrite(&gp[i].mol[j].nmol,  sizeof(double), 1, fp), 1, "nmol");
    checkFwrite(fwrite(&gp[i].dopb_turb, sizeof(double), 1, fp), 1, "dopb_turb");
    for(j=0;j<par->nSpecies;j++){
      checkFwrite(fwrite(gp[i].mol[j].pops,  sizeof(double)*md[j].nlev, 1, fp), 1, "pops");
      checkFwrite(fwrite(dummyMol[j].knu,    sizeof(double)*md[j].nline,1, fp), 1, "knu");
      checkFwrite(fwrite(dummyMol[j].dust,   sizeof(double)*md[j].nline,1, fp), 1, "dust");
      checkFwrite(fwrite(&gp[i].mol[j].dopb, sizeof(double),            1, fp), 1, "dopb");
      checkFwrite(fwrite(&gp[i].mol[j].binv, sizeof(double),            1, fp), 1, "binv");
    }
    checkFwrite(fwrite(&gp[i].dens[0],     sizeof(double), 1, fp), 1, "dens");
    checkFwrite(fwrite(&gp[i].t[0],        sizeof(double), 1, fp), 1, "t");
    checkFwrite(fwrite(&gp[i].mol[0].abun, sizeof(double), 1, fp), 1, "abun");
  }

  fclose(fp);

  for(j=0;j<par->nSpecies;j++){
    free(dummyMol[j].dust);
    free(dummyMol[j].knu);
    free(dummy2[j]);
  }
  free(dummyMol);
  free(dummy2);

  free(nTrans);
}

