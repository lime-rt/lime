/*
 *  dust.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#include "lime.h"

/*....................................................................*/
double planckfunc(const double freq, const double tKelvin){
  double bb=10.,wn;
  if(tKelvin<EPS) bb = 0.0;
  else {
    wn=freq/CLIGHT;
    if (HPLANCK*freq>100.*KBOLTZ*tKelvin) 
      bb=2.*HPLANCK*wn*wn*freq*exp(-HPLANCK*freq/KBOLTZ/tKelvin);
    else 
      bb=2.*HPLANCK*wn*wn*freq/(exp(HPLANCK*freq/KBOLTZ/tKelvin)-1);
  }
  return bb;
}

/*....................................................................*/
void readDustFile(char *dustFileName, double **lamtab, double **kaptab\
  , int *nEntries){

  /* NOTE! The calling routine must free lamtab and kaptab after use.
  */
  FILE *fp;
  int i=0,k;
  char string[80];

  /* Open the file and count the values it contains.
  */
  if((fp=fopen(dustFileName, "r"))==NULL){
    if(!silent) bail_out("Error opening dust opacity data file!");
    exit(1);
  }


  while(fgets(string,80,fp) != NULL){
    i++;
  }

  rewind(fp);

  /* Now read the values.
  */
  if(i>0){
    *lamtab=malloc(sizeof(**lamtab)*i);
    *kaptab=malloc(sizeof(**kaptab)*i);
  } else {
    if(!silent) bail_out("No opacities read");
    exit(1);
  }
  for(k=0;k<i;k++){
    checkFscanf(fscanf(fp,"%lf %lf\n", &(*lamtab)[k], &(*kaptab)[k]), 2, "dust opacities.");
    (*lamtab)[k]=log10((*lamtab)[k]/1e6);
    (*kaptab)[k]=log10((*kaptab)[k]);
  }
  fclose(fp);

  *nEntries = i;
}

/*....................................................................*/
double interpolateKappa(const double freq, double *lamtab, double *kaptab\
  , const int nEntries, gsl_spline *spline, gsl_interp_accel *acc){
  /* Note that the multiplications by 0.1 below are to convert cm^2/g to m^2/kg. */

  double loglam, kappa;

  loglam=log10(CLIGHT/freq);
  if(loglam < lamtab[0])
    kappa = 0.1*pow(10.,kaptab[0] + (loglam-lamtab[0])\
          *(kaptab[1]-kaptab[0])/(lamtab[1]-lamtab[0]));
  else if(loglam > lamtab[nEntries-1])
    kappa = 0.1*pow(10.,kaptab[nEntries-2] + (loglam-lamtab[nEntries-2])\
          *(kaptab[nEntries-1]-kaptab[nEntries-2])\
          /(lamtab[nEntries-1]-lamtab[nEntries-2]));
  else
    kappa = 0.1*pow(10.,gsl_spline_eval(spline,loglam,acc));

  return kappa;
}

/*....................................................................*/
void calcDustData(configInfo *par, double *dens, double *freqs\
  , const double gtd, double *kappatab, const int numLines, const double tsKelvin[]\
  , double *knus, double *dusts){

  double tKelvin,gasMassDensityAMUs,dustToGas;
  int di,iline;

  /* Check if input model supplies a dust temperature. Otherwise use the kinetic temperature. */
  if(tsKelvin[1]<=0.0) { /* Flags that the user has not set it. */
    tKelvin = tsKelvin[0];
  } else {
    tKelvin = tsKelvin[1];
  }

  if(par->collPartUserSetFlags==0){ /* this means the user did not set any of the collision-partner-related parameters. Use the old formula. */
    dustToGas = AMU*2.4*dens[0]/gtd;
  }else{
    gasMassDensityAMUs = 0.0;
    for(di=0;di<par->numDensities;di++)
      gasMassDensityAMUs += dens[di]*par->collPartMolWeights[di];

    dustToGas = AMU*gasMassDensityAMUs/gtd;
  }

  for(iline=0;iline<numLines;iline++){
    knus[iline] = kappatab[iline]*dustToGas;
    dusts[iline] = planckfunc(freqs[iline],tKelvin);
  }
}

