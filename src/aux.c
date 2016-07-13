/*
 *  aux.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2016 The LIME development team
 *
 */

#include "lime.h"
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <float.h>


void
parseInput(inputPars *par, image **img, molData **m){
  int i,id, ispec;
  double BB[3],normBSquared;
  double cosPhi,sinPhi,cosTheta,sinTheta,dummyVel[DIM];

  par->ncell=par->pIntensity+par->sinkPoints;
  par->radiusSqu=par->radius*par->radius;
  par->minScaleSqu=par->minScale*par->minScale;
  if(par->pregrid!=NULL) par->doPregrid=1;

  /* Check that the user has supplied this function (needed unless par->pregrid):
  */
  if(!par->pregrid)
    velocity(0.0,0.0,0.0, dummyVel);

  /*
Now we need to calculate the cutoff value used in calcSourceFn(). The issue is to decide between

  y = (1 - exp[-x])/x

or the approximation using the Taylor expansion of exp(-x), which to 3rd order is

  y ~ 1 - x/2 + x^2/6.

The cutoff will be the value of abs(x) for which the error in the exact expression equals the next unused Taylor term, which is x^3/24. This error can be shown to be given for small |x| by epsilon/|x|, where epsilon is the floating-point precision of the computer. Hence the cutoff evaluates to

  |x|_cutoff = (24*epsilon)^{1/4}.

  */
  par->taylorCutoff = pow(24.*DBL_EPSILON, 0.25);

  /* Allocate pixel space and parse image information */
  for(i=0;i<par->nImages;i++){
    if((*img)[i].nchan == 0 && (*img)[i].velres<0 ){
      /* Assume continuum image */

      if(par->polarization){
        (*img)[i].nchan=3;

        if(!silent){
          /* Do a sketchy check which might indicate if the user has forgotten to supply a magfield function, and warn if this comes up positive. Note: there is no really robust way at present to distinguish the default magfield function (which, if called, indicates that the user forgot to supply their own) from one the user has supplied but which happens to set the B field to 0 at the origin.
          */
          magfield(par->minScale,par->minScale,par->minScale,BB);
          normBSquared = BB[0]*BB[0] + BB[1]*BB[1] + BB[2]*BB[2];
          if(normBSquared <= 0.) warning("Zero B field - did you remember to supply a magfield function?");
        }
      }else
        (*img)[i].nchan=1;

      if((*img)[i].trans>-1 || (*img)[i].bandwidth>-1. || (*img)[i].freq==0 || par->dust==NULL){
        if(!silent) bail_out("Error: Image keywords are ambiguous");
        exit(1);
      }
      (*img)[i].doline=0;
    } else if (((*img)[i].nchan>0 || (*img)[i].velres > 0)){
      /* Assume line image. */
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

    /* Rotation matrix

            |1          0           0   |
     R_x(a)=|0        cos(a)      sin(a)|
            |0       -sin(a)      cos(a)|

            |cos(b)     0       -sin(b)|
     R_y(b)=|  0        1          0   |
            |sin(b)     0        cos(b)|

            |      cos(b)       0          -sin(b)|
     Rot =  |sin(a)sin(b)     cos(a)  sin(a)cos(b)|
            |cos(a)sin(b)    -sin(a)  cos(a)cos(b)|

    */

    cosPhi   = cos((*img)[i].phi);
    sinPhi   = sin((*img)[i].phi);
    cosTheta = cos((*img)[i].theta);
    sinTheta = sin((*img)[i].theta);
    (*img)[i].rotMat[0][0] =           cosPhi;
    (*img)[i].rotMat[0][1] =  0.0;
    (*img)[i].rotMat[0][2] =          -sinPhi;
    (*img)[i].rotMat[1][0] =  sinTheta*sinPhi;
    (*img)[i].rotMat[1][1] =  cosTheta;
    (*img)[i].rotMat[1][2] =  sinTheta*cosPhi;
    (*img)[i].rotMat[2][0] =  cosTheta*sinPhi;
    (*img)[i].rotMat[2][1] = -sinTheta;
    (*img)[i].rotMat[2][2] =  cosTheta*cosPhi;
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
      (*m)[i].down = NULL;
      (*m)[i].ntemp = NULL;
      (*m)[i].eterm = NULL;
      (*m)[i].gstat = NULL;
      (*m)[i].cmb = NULL;
      (*m)[i].local_cmb = NULL;
    }
}

void
freeMoldata( inputPars *par, molData* mol )
{
  int i;
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
          if( mol[i].eterm != NULL )
            {
              free(mol[i].eterm);
            }
          if( mol[i].gstat != NULL )
            {
              free(mol[i].gstat);
            }
          if( mol[i].cmb != NULL )
            {
              free(mol[i].cmb);
            }
          if( mol[i].local_cmb != NULL )
            {
              free(mol[i].local_cmb);
            }

          if( mol[i].down != NULL )
            {
              int j=0;
              for (j=0;j<mol[i].npart;j++) free(mol[i].down[j]);
              free(mol[i].down);
            }
	 free(mol[i].ntemp);

        }
      free(mol);
    }
}

void
freeGridPointData(inputPars *par, gridPointData *mol){
  int i;
  if (mol!= 0){
    for (i=0;i<par->nSpecies;i++){
      if (mol[i].jbar != NULL){
        free(mol[i].jbar);
      }
      if (mol[i].phot != NULL){
        free(mol[i].phot);
      }
      if (mol[i].vfac != NULL){
        free(mol[i].vfac);
      }
    }
    free(mol);
  }
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

void checkGridDensities(inputPars *par, struct grid *g){
  int i;
  static _Bool warningAlreadyIssued=0;
  char errStr[80];

  if(!silent){ /* Warn if any densities too low. */
    i = 0;
    while(i<par->pIntensity && !warningAlreadyIssued){
      if(g[i].dens[0]<TYPICAL_ISM_DENS){
        warningAlreadyIssued = 1;
        sprintf(errStr, "g[%d].dens[0] at %.1e is below typical values for the ISM (~%.1e).", i, g[i].dens[0], TYPICAL_ISM_DENS);
        warning(errStr);
        warning("This could give you convergence problems. NOTE: no further warnings will be issued.");
      }
      i++;
    }
  }
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
  if(*nlinetot > 0){
  *counta=malloc(sizeof(*counta)* *nlinetot);
  *countb=malloc(sizeof(*countb)* *nlinetot);
  } else {
    if(!silent) bail_out("Error: Line count finds no lines");
    exit(0);
  }
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
  int id,conv=0,iter,ilev,prog=0,ispec,c=0,n,i,threadI,nVerticesDone;
  double percent=0.,*median,result1=0,result2=0,snr,delta_pop;
  blend *matrix;
  struct statistics { double *pop, *ave, *sigma; } *stat;
  const gsl_rng_type *ranNumGenType = gsl_rng_ranlxs2;
  _Bool luWarningGiven=0;
  gsl_error_handler_t *defaultErrorHandler=NULL;

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
  gsl_rng *ran = gsl_rng_alloc(ranNumGenType);
#ifdef TEST
  gsl_rng_set(ran, 1237106) ;
#else 
  gsl_rng_set(ran,time(0));
#endif

  gsl_rng **threadRans;
  threadRans = malloc(sizeof(gsl_rng *)*par->nThreads);

  for (i=0;i<par->nThreads;i++){
    threadRans[i] = gsl_rng_alloc(ranNumGenType);
    gsl_rng_set(threadRans[i],(int)(gsl_rng_uniform(ran)*1e6));
  }

  /* Read in all molecular data */
  for(id=0;id<par->nSpecies;id++) molinit(m,par,g,id);

  /* Check for blended lines */
  lineBlend(m,par,&matrix);

  if(par->lte_only || par->init_lte) LTE(par,g,m);

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
    defaultErrorHandler = gsl_set_error_handler_off();
    /*
This is done to allow proper handling of errors which may arise in the LU solver within stateq(). It is done here because the GSL documentation does not recommend leaving the error handler at the default within multi-threaded code.

While this is off however, other gsl_* etc calls will not exit if they encounter a problem. We may need to pay some attention to trapping their errors.
    */

    do{
      if(!silent) progressbar2(0, prog++, 0, result1, result2);

      for(id=0;id<par->ncell && !g[id].sink;id++){
        for(ilev=0;ilev<m[0].nlev;ilev++) {
          for(iter=0;iter<4;iter++) stat[id].pop[ilev+m[0].nlev*iter]=stat[id].pop[ilev+m[0].nlev*(iter+1)];
          stat[id].pop[ilev+m[0].nlev*4]=g[id].mol[0].pops[ilev];
        }
      }

      nVerticesDone=0;
      omp_set_dynamic(0);
#pragma omp parallel private(i,id,ispec,threadI) num_threads(par->nThreads)
      {
        threadI = omp_get_thread_num();

        /* Declare and allocate thread-private variables */
        gridPointData *mp;	// Could have declared them earlier
        double *halfFirstDs;	// and included them in private() I guess.
        mp=malloc(sizeof(gridPointData)*par->nSpecies);
        for (i=0;i<par->nSpecies;i++){
          mp[i].phot = malloc(sizeof(double)*m[i].nline*max_phot);
          mp[i].vfac = malloc(sizeof(double)*           max_phot);
          mp[i].jbar = malloc(sizeof(double)*m[i].nline);
        }
        halfFirstDs = malloc(sizeof(*halfFirstDs)*max_phot);

#pragma omp for
        for(id=0;id<par->pIntensity;id++){
#pragma omp atomic
          ++nVerticesDone;

          if (threadI == 0){ // i.e., is master thread
            if(!silent) progressbar((double)nVerticesDone/par->pIntensity,10);
          }
          if(g[id].dens[0] > 0 && g[id].t[0] > 0){
            photon(id,g,m,0,threadRans[threadI],par,matrix,mp,halfFirstDs);
            for(ispec=0;ispec<par->nSpecies;ispec++)
              stateq(id,g,m,ispec,par,mp,halfFirstDs,&luWarningGiven);
          }
          if (threadI == 0){ // i.e., is master thread
            if(!silent) warning("");
          }
        }

        freeGridPointData(par, mp);
        free(halfFirstDs);
      } // end parallel block.

      for(id=0;id<par->ncell && !g[id].sink;id++){
        snr=0;
        n=0;
        for(ilev=0;ilev<m[0].nlev;ilev++) {
          stat[id].ave[ilev]=0;
          for(iter=0;iter<5;iter++) stat[id].ave[ilev]+=stat[id].pop[ilev+m[0].nlev*iter];
          stat[id].ave[ilev]=stat[id].ave[ilev]/5.;
          stat[id].sigma[ilev]=0;
          for(iter=0;iter<5;iter++) {
            delta_pop = stat[id].pop[ilev+m[0].nlev*iter]-stat[id].ave[ilev];
            stat[id].sigma[ilev]+=delta_pop*delta_pop;
          }
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

      median=malloc(sizeof(*median)*gsl_max(c,1));
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

      if(!silent) progressbar2(1, prog, percent, result1, result2);
      if(par->outputfile) popsout(par,g,m);
    } while(conv++<NITERATIONS);
    gsl_set_error_handler(defaultErrorHandler);
    if(par->binoutputfile) binpopsout(par,g,m);
  }

  for (i=0;i<par->nThreads;i++){
    gsl_rng_free(threadRans[i]);
  }
  free(threadRans);
  gsl_rng_free(ran);
  for(id=0;id<par->pIntensity;id++){
    free(stat[id].pop);
    free(stat[id].ave);
    free(stat[id].sigma);
  }
  free(stat);
  *popsdone=1;
}



