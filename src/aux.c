/*
 *  aux.c
 *  LIME, The versatile 3D line modeling environment
 *
 *  Created by Christian Brinch on 16/11/06.
 *  Copyright 2006-2014, Christian Brinch,
 *  <brinch@nbi.dk>
 *  Niels Bohr institutet
 *  University of Copenhagen
 *	All rights reserved.
 *
 */

#include "lime.h"
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>


void
parseInput( char* input_file, inputPars *par, image **img, molData **m){
  FILE *fp;
  int i,id;
  double BB[3];

  /* Set default values */
  par->restart[0]       = '\0';;
  par->pregrid[0]       = '\0';
  par->gridfile[0]      = '\0';
  par->binoutputfile[0] = '\0';
  par->outputfile[0]    = '\0';
  par->dust[0]          = '\0';

  par->tcmb = 2.728;
  par->lte_only=0;
  par->sampling=0;
  par->blend=0;
  par->antialias=1;
  par->polarization=0;
  par->pIntensity=0;
  par->sinkPoints=0;

  strcpy(par->python_module_name, "model");
  strcpy(par->python_module_path, "./");
  strcpy(par->density_func_name, "density");
  strcpy(par->velocity_func_name, "velocity");
  strcpy(par->temperature_func_name, "temperature");
  strcpy(par->doppler_func_name, "doppler");
  strcpy(par->abundance_func_name, "abundance");

  /* Allocate space for output fits images */
  (*img)=malloc(sizeof(image)*MAX_NSPECIES);
  par->moldatfile=malloc(sizeof( filename_t ) * MAX_NSPECIES);
  for(id=0;id<MAX_NSPECIES;id++){
    (*img)[id].filename[0]='\0';
    par->moldatfile[id][0]='\0';

    (*img)[id].source_vel=0.0;
    (*img)[id].phi=0.0;
    (*img)[id].nchan=0;
    (*img)[id].velres=-1.;
    (*img)[id].trans=-1;
    (*img)[id].freq=-1.;
    (*img)[id].bandwidth=-1.;
  }

  if( input( input_file, par, *img) != EXIT_SUCCESS )
    {
      if(!silent) bail_out("Error: Cannot Read input file");
      exit(1);
    }

  if( python_call_initialize( par ) != EXIT_SUCCESS )
    {
      if(!silent) bail_out("Error: Cannot initialize python");
      exit(1);
    }

  id=0;
  while(strlen( (*img)[++id].filename ) != 0 );
  par->nImages=id;
  if(par->nImages==0) {
    if(!silent) bail_out("Error: no images defined");
    exit(1);
  }

  *img=realloc(*img, sizeof(image)*par->nImages);

  id=-1;
  while( strlen( par->moldatfile[++id] ) != 0);
  par->nSpecies=id;
  if( par->nSpecies <= 0 )
    {
      if(!silent) bail_out("Error: no moldatfile provided");
      exit(1);
    }

  par->moldatfile=realloc(par->moldatfile, sizeof( filename_t )*par->nSpecies);


  par->ncell=par->pIntensity+par->sinkPoints;

  /* Check if files exists */
  for(id=0;id<par->nSpecies;id++){
    if((fp=fopen(par->moldatfile[id], "r"))==NULL) {
      openSocket(par, id);
    }
    else {
      fclose(fp);
    }
  }
  if(par->dust != NULL){
    if((fp=fopen(par->dust, "r"))==NULL){
      if(!silent) bail_out("Error opening dust opacity data file!");
      exit(1);
    }
    else  {
      fclose(fp);
    }
  }





  /* Allocate pixel space and parse image information */
  for(i=0;i<par->nImages;i++){
    if((*img)[i].nchan == 0 && (*img)[i].velres<0 ){
      /* Assume continuum image */

      /* Check for polarization */
      BB[0]=0.;
      magfield(par->minScale,par->minScale,par->minScale,BB);
      if(fabs(BB[0]) > 0.) par->polarization=1;

      if(par->polarization) (*img)[i].nchan=3;
      else (*img)[i].nchan=1;
      if((*img)[i].trans>-1 || (*img)[i].bandwidth>-1. || (*img)[i].freq==0 || par->dust==NULL){
        if(!silent) bail_out("Error: Image keywords are ambiguous");
        exit(1);
      }
      (*img)[i].doline=0;
    } else if (((*img)[i].nchan>0 || (*img)[i].velres > 0)){
      /* Assume line image */
      par->polarization=0;
      if(par->moldatfile==NULL){
        if(!silent) bail_out("Error: No data file is specified for line image.");
        exit(1);
      }
      if(((*img)[i].trans>-1 && (*img)[i].freq>-1) || ((*img)[i].trans<0 && (*img)[i].freq<0)){
        if(!silent) bail_out("Error: Specify either frequency or transition ");
        exit(1);
      }
      if(((*img)[i].nchan==0 && (*img)[i].bandwidth<0) || ((*img)[i].bandwidth<0 && (*img)[i].velres<0)){
        if(!silent) bail_out("Error: Image keywords are not set properly");
        exit(1);
      }
      (*img)[i].doline=1;
    }
    (*img)[i].imgres=(*img)[i].imgres/206264.806;
    (*img)[i].pixel = malloc(sizeof(spec)*(*img)[i].pxls*(*img)[i].pxls);
    for(id=0;id<((*img)[i].pxls*(*img)[i].pxls);id++){
      (*img)[i].pixel[id].intense = malloc(sizeof(double)*(*img)[i].nchan);
      (*img)[i].pixel[id].tau = malloc(sizeof(double)*(*img)[i].nchan);
    }
  }

  /* Allocate moldata array */
  (*m)=malloc(sizeof(molData)*par->nSpecies);
  for( i=0; i<par->nSpecies; i++ )
    {
      (*m)[i].ntrans = NULL;
      (*m)[i].lal = NULL;
      (*m)[i].lau = NULL;
      (*m)[i].lcl = NULL;
      (*m)[i].lcu = NULL;
      (*m)[i].aeinst = NULL;
      (*m)[i].freq = NULL;
      (*m)[i].beinstu = NULL;
      (*m)[i].beinstl = NULL;
      (*m)[i].up = NULL;
      (*m)[i].down = NULL;
      (*m)[i].eterm = NULL;
      (*m)[i].gstat = NULL;
      (*m)[i].jbar = NULL;
      (*m)[i].cmb = NULL;
      (*m)[i].local_cmb = NULL;
      (*m)[i].phot = NULL;
      (*m)[i].ds = NULL;
      (*m)[i].vfac = NULL;
      (*m)[i].weight = NULL;
    }
}

void
freeInput( inputPars *par, image* img, molData* mol )
{
  python_call_finalize();

  int i,id;
  if( mol!= 0 )
    {
      for( i=0; i<par->nSpecies; i++ )
        {
          if( mol[i].ntrans != NULL )
            {
              free(mol[i].ntrans);
            }
          if( mol[i].lal != NULL )
            {
              free(mol[i].lal);
            }
          if( mol[i].lau != NULL )
            {
              free(mol[i].lau);
            }
          if( mol[i].lcl != NULL )
            {
              free(mol[i].lcl);
            }
          if( mol[i].lcu != NULL )
            {
              free(mol[i].lcu);
            }
          if( mol[i].aeinst != NULL )
            {
              free(mol[i].aeinst);
            }
          if( mol[i].freq != NULL )
            {
              free(mol[i].freq);
            }
          if( mol[i].beinstu != NULL )
            {
              free(mol[i].beinstu);
            }
          if( mol[i].beinstl != NULL )
            {
              free(mol[i].beinstl);
            }
          if( mol[i].up != NULL )
            {
              free(mol[i].up);
            }
          if( mol[i].down != NULL )
            {
              free(mol[i].down);
            }
          if( mol[i].eterm != NULL )
            {
              free(mol[i].eterm);
            }
          if( mol[i].gstat != NULL )
            {
              free(mol[i].gstat);
            }
          if( mol[i].jbar != NULL )
            {
              free(mol[i].jbar);
            }
          if( mol[i].cmb != NULL )
            {
              free(mol[i].cmb);
            }
          if( mol[i].local_cmb != NULL )
            {
              free(mol[i].local_cmb);
            }
          if( mol[i].phot != NULL )
            {
              free(mol[i].phot);
            }
          if( mol[i].ds != NULL )
            {
              free(mol[i].ds);
            }
          if( mol[i].vfac != NULL )
            {
              free(mol[i].vfac);
            }
          if( mol[i].weight != NULL )
            {
              free(mol[i].weight);
            }
        }
      free(mol);
    }
  for(i=0;i<par->nImages;i++){
    for(id=0;id<(img[i].pxls*img[i].pxls);id++){
      free( img[i].pixel[id].intense );
      free( img[i].pixel[id].tau );
    }
    free(img[i].pixel);
  }
  free(img);
  free(par->moldatfile);
}



float
invSqrt(float x){
  /* The magic Quake(TM) fast inverse square root algorithm   */
  /* Can _only_ be used on 32-bit machine architectures       */
  float xhalf = 0.5f*x;
  int i = *(int*)&x;
  i = 0x5f3759df - (i>>1);
  x = *(float*)&i;
  x = x*(1.5f - xhalf*x*x);
  return x;
}

void
continuumSetup(int im, image *img, molData *m, inputPars *par, struct grid *g){
  int id;
  img[im].trans=0;
  m[0].nline=1;
  m[0].freq= malloc(sizeof(double));
  m[0].freq[0]=img[im].freq;
  for(id=0;id<par->ncell;id++) {
    freePopulation( par, m, g[id].mol );
    g[id].mol=malloc(sizeof(struct populations)*1);
    g[id].mol[0].dust = malloc(sizeof(double)*m[0].nline);
    g[id].mol[0].knu  = malloc(sizeof(double)*m[0].nline);
    g[id].mol[0].pops = NULL;
    g[id].mol[0].partner = NULL;
  }
  if(par->outputfile) popsout(par,g,m);
  kappa(m,g,par,0);
}

void
lineCount(int n,molData *m,int **counta,int **countb,int *nlinetot){
  int ispec,iline,count;

  *nlinetot=0;
  for(ispec=0;ispec<n;ispec++) *nlinetot+=m[ispec].nline;
  *counta=malloc(sizeof(int)* *nlinetot);
  *countb=malloc(sizeof(int)* *nlinetot);
  count=0;
  for(ispec=0;ispec<n;ispec++) {
    for(iline=0;iline<m[ispec].nline;iline++){
      (*counta)[count]=ispec;
      (*countb)[count++]=iline;
    }
  }
}

void
lineBlend(molData *m, inputPars *par, blend **matrix){
  int iline, jline, nlinetot=0,c;
  int *counta,*countb;

  lineCount(par->nSpecies, m, &counta, &countb, &nlinetot);

  c=0;
  for(iline=0;iline<nlinetot;iline++){
    for(jline=0;jline<nlinetot;jline++){
      if(fabs((m[counta[jline]].freq[countb[jline]]-m[counta[iline]].freq[countb[iline]])/m[counta[iline]].freq[countb[iline]]*CLIGHT) < blendmask
         && iline !=jline) c++;
    }
  }
  if(c>0){
    if(par->blend){
      if(!silent) warning("There are blended lines (Line blending is switched on)");
    } else {
      if(!silent) warning("There are blended lines (Line blending is switched off)");
    }

    (*matrix)=malloc(sizeof(blend)*c);

    c=0;
    for(iline=0;iline<nlinetot;iline++){
      for(jline=0;jline<nlinetot;jline++){
        if(fabs((m[counta[jline]].freq[countb[jline]]-m[counta[iline]].freq[countb[iline]])/m[counta[iline]].freq[countb[iline]]*CLIGHT) < blendmask
           && iline != jline){
          (*matrix)[c].line1=iline;
          (*matrix)[c].line2=jline;
          (*matrix)[c++].deltav=-(m[counta[jline]].freq[countb[jline]]-m[counta[iline]].freq[countb[iline]])/m[counta[iline]].freq[countb[iline]]*CLIGHT;
        }
      }
    }
  }
  free(counta);
  free(countb);

}

void
levelPops(molData *m, inputPars *par, struct grid *g, int *popsdone){
  int id,conv=0,iter,ilev,prog=0,ispec,c=0,n;
  double percent=0.,pstate,*median,result1=0,result2=0,snr;
  blend *matrix;
  struct statistics { double *pop, *ave, *sigma; } *stat;

  stat=malloc(sizeof(struct statistics)*par->pIntensity);

  for(id=0;id<par->ncell;id++) {
    freePopulation( par, m, g[id].mol );
    g[id].mol=malloc(sizeof(struct populations)*par->nSpecies);
    int i;
    for( i=0; i<par->nSpecies; i++ )
      {
        g[id].mol[i].dust = NULL;
        g[id].mol[i].knu  = NULL;
        g[id].mol[i].pops = NULL;
        g[id].mol[i].partner = NULL;
      }
  }

  /* Random number generator */
  gsl_rng *ran = gsl_rng_alloc(gsl_rng_ranlxs2);
  gsl_rng_set(ran,time(0));

  /* Read in all molecular data */
  for(id=0;id<par->nSpecies;id++) molinit(m,par,g,id);

  /* Check for blended lines */
  lineBlend(m,par,&matrix);

  if(par->lte_only) LTE(par,g,m);

  for(id=0;id<par->pIntensity;id++){
    stat[id].pop=malloc(sizeof(double)*m[0].nlev*5);
    stat[id].ave=malloc(sizeof(double)*m[0].nlev);
    stat[id].sigma=malloc(sizeof(double)*m[0].nlev);
    for(ilev=0;ilev<m[0].nlev;ilev++) {
      for(iter=0;iter<5;iter++) stat[id].pop[ilev+m[0].nlev*iter]=g[id].mol[0].pops[ilev];
    }
  }

  if(par->outputfile) popsout(par,g,m);


  /* Initialize convergence flag */
  for(id=0;id<par->ncell;id++){
    g[id].conv=0;
  }

  if(par->lte_only==0){
    do{
      if(!silent) progressbar2(prog++, 0, result1, result2);
      pstate=0.;

      for(id=0;id<par->ncell && !g[id].sink;id++){
        if(!silent) progressbar((double)id/par->pIntensity,10);
        for(ilev=0;ilev<m[0].nlev;ilev++) {
          for(iter=0;iter<4;iter++) stat[id].pop[ilev+m[0].nlev*iter]=stat[id].pop[ilev+m[0].nlev*(iter+1)];
          stat[id].pop[ilev+m[0].nlev*4]=g[id].mol[0].pops[ilev];
        }
        if(g[id].dens[0] > 0 && g[id].t[0] > 0){
          photon(id,g,m,0,ran,par,matrix);
          for(ispec=0;ispec<par->nSpecies;ispec++) stateq(id,g,m,&pstate,ispec,par);
        }
        if(!silent) warning("");

        snr=0;
        n=0;
        for(ilev=0;ilev<m[0].nlev;ilev++) {
          stat[id].ave[ilev]=0;
          for(iter=0;iter<5;iter++) stat[id].ave[ilev]+=stat[id].pop[ilev+m[0].nlev*iter];
          stat[id].ave[ilev]=stat[id].ave[ilev]/5.;
          stat[id].sigma[ilev]=0;
          for(iter=0;iter<5;iter++) stat[id].sigma[ilev]+=pow(stat[id].pop[ilev+m[0].nlev*iter]-stat[id].ave[ilev],2);
          stat[id].sigma[ilev]=sqrt(stat[id].sigma[ilev])/5.;
          if(g[id].mol[0].pops[ilev] > 1e-12) c++;

          if(g[id].mol[0].pops[ilev] > 1e-12 && stat[id].sigma[ilev] > 0.){
            snr+=g[id].mol[0].pops[ilev]/stat[id].sigma[ilev];
            n++;
          }
        }
        if(n>0) snr=snr/n;
        else if(n==0) snr=1e6;
        if(snr > 3.) g[id].conv=2;
        if(snr <= 3 && g[id].conv==2) g[id].conv=1;
      }

      median=malloc(sizeof(double)*gsl_max(c,1));
      c=0;
      for(id=0;id<par->pIntensity;id++){
        for(ilev=0;ilev<m[0].nlev;ilev++){
          if(g[id].mol[0].pops[ilev] > 1e-12) median[c++]=g[id].mol[0].pops[ilev]/stat[id].sigma[ilev];
        }
      }

      gsl_sort(median, 1, c);
      if(conv>1){
        result1=median[0];
        result2 =gsl_stats_median_from_sorted_data(median, 1, c);
      }
      free(median);

      if(!silent) progressbar2(prog, percent, result1, result2);
      if(par->outputfile) popsout(par,g,m);
    } while(conv++<NITERATIONS);
    if(par->binoutputfile) binpopsout(par,g,m);
  }
  gsl_rng_free(ran);
  for(id=0;id<par->pIntensity;id++){
    free(stat[id].pop);
    free(stat[id].ave);
    free(stat[id].sigma);
  }
  free(stat);
  *popsdone=1;
}



