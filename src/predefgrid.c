/*
 *  predefgrid.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2016 The LIME development team
 *
 */

#include "lime.h"

void
predefinedGrid(configInfo *par, struct grid *g){
  FILE *fp;
  int i,j;
  double x,y,z,scale;
  struct cell *dc=NULL; /* Not used at present. */
  unsigned long numCells;

  gsl_rng *ran = gsl_rng_alloc(gsl_rng_ranlxs2);
#ifdef TEST
  gsl_rng_set(ran,6611304);
#else
  gsl_rng_set(ran,time(0));
#endif

  fp=fopen(par->pregrid,"r");
  par->ncell=par->pIntensity+par->sinkPoints;

  for(i=0;i<par->pIntensity;i++){
    //    fscanf(fp,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &g[i].id, &g[i].x[0], &g[i].x[1], &g[i].x[2],  &g[i].dens[0], &g[i].t[0], &abun, &g[i].dopb, &g[i].vel[0], &g[i].vel[1], &g[i].vel[2]);
    //    fscanf(fp,"%d %lf %lf %lf %lf %lf %lf %lf\n", &g[i].id, &g[i].x[0], &g[i].x[1], &g[i].x[2],  &g[i].dens[0], &g[i].t[0], &abun, &g[i].dopb);
    int nRead = fscanf(fp,"%d %lf %lf %lf %lf %lf %lf %lf %lf\n", &g[i].id, &g[i].x[0], &g[i].x[1], &g[i].x[2],  &g[i].dens[0], &g[i].t[0], &g[i].vel[0], &g[i].vel[1], &g[i].vel[2]);
    if( nRead != 9 || g[i].id < 0 || g[i].id > par->ncell)
      {
        if(!silent) bail_out("Reading Grid File error");
        exit(0);
      }

    g[i].dopb=200;
    g[i].abun[0]=1e-9;

    g[i].sink=0;
    g[i].t[1]=g[i].t[0];
    g[i].mol[0].nmol=g[i].abun[0]*g[i].dens[0];
    g[i].B[0]=0.0;
    g[i].B[1]=0.0;
    g[i].B[2]=0.0;

    /* This next step needs to be done, even though it looks stupid */
    g[i].dir=malloc(sizeof(point)*1);
    g[i].ds =malloc(sizeof(double)*1);
    g[i].neigh =malloc(sizeof(struct grid *)*1);
    if(!silent) progressbar((double) i/((double)par->pIntensity-1), 4);	
  }

  checkGridDensities(par, g);

  for(i=par->pIntensity;i<par->ncell;i++){
    x=2*gsl_rng_uniform(ran)-1.;
    y=2*gsl_rng_uniform(ran)-1.;
    z=2*gsl_rng_uniform(ran)-1.;
    if(x*x+y*y+z*z<1){
      scale=par->radius*sqrt(1/(x*x+y*y+z*z));
      g[i].id=i;
      g[i].x[0]=scale*x;
      g[i].x[1]=scale*y;
      g[i].x[2]=scale*z;
      g[i].sink=1;
      g[i].abun[0]=0;
      g[i].dens[0]=1e-30;
      g[i].mol[0].nmol=0.0; /* Just to give it a value. */
      g[i].t[0]=par->tcmb;
      g[i].t[1]=par->tcmb;
      g[i].B[0]=0.0;
      g[i].B[1]=0.0;
      g[i].B[2]=0.0;
      g[i].dopb=0.;
      for(j=0;j<DIM;j++) g[i].vel[j]=0.;
    } else i--;
  }
  fclose(fp);

  delaunay(DIM, g, (unsigned long)par->ncell, 0, &dc, &numCells);
  distCalc(par,g);
  //  getArea(par,g, ran);
  //  getMass(par,g, ran);
  calcInterpCoeffs_lin(par,g);

  if(par->gridfile) write_VTK_unstructured_Points(par, g);
  gsl_rng_free(ran);
  free(dc);

  par->numDensities = 1;
}

