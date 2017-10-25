/*
 *  aux.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2017 The LIME development team
 *
TODO:
 */

#include "lime.h"

/*....................................................................*/
void sigintHandler(int sigI){
#ifdef IS_PYTHON
  Py_Finalize();
#endif

//*** write output file?

exit(1);
}

/*....................................................................*/
void checkFgets(char *fgetsResult, char *message){
  const size_t buflen=80;
  char string[buflen];

  if(fgetsResult==NULL){
    if(!silent){
      snprintf(string, buflen, "fgets() failed to read %s", message);
      bail_out(string);
    }
exit(1);
  }
}

/*....................................................................*/
void checkFscanf(const int fscanfResult, const int expectedNum, char *message){
  const size_t buflen=80;
  char string[buflen];

  if(fscanfResult!=expectedNum){
    if(!silent){
      snprintf(string, buflen, "fscanf() failed to read %s - read %d bytes when %d expected.", message, fscanfResult, expectedNum);
      bail_out(string);
    }
exit(1);
  }
}

/*....................................................................*/
void checkFread(const size_t freadResult, const size_t expectedNum, char *message){
  const size_t buflen=80;
  char string[buflen];

  if(freadResult!=expectedNum){
    if(!silent){
      snprintf(string, buflen, "fread() failed to read %s. Expected %d got %d", message, (int)expectedNum, (int)freadResult);
      bail_out(string);
    }
exit(1);
  }
}

/*....................................................................*/
void checkFwrite(const size_t fwriteResult, const size_t expectedNum, char *message){
  const size_t buflen=80;
  char string[buflen];

  if(fwriteResult!=expectedNum){
    if(!silent){
      snprintf(string, buflen, "fwrite() failed to write %s. Expected %d got %d", message, (int)expectedNum, (int)fwriteResult);
      bail_out(string);
    }
exit(1);
  }
}

/*....................................................................*/
void mallocAndSetDefaultMolData(const int nSpecies, molData **md){
  int i;

  (*md)=malloc(sizeof(molData)*nSpecies);
  for(i=0;i<nSpecies;i++){
    (*md)[i].nlev  = -1;
    (*md)[i].nline = -1;
    (*md)[i].npart = -1;
    (*md)[i].amass = -1.0;
    (*md)[i].part    = NULL;
    (*md)[i].lal     = NULL;
    (*md)[i].lau     = NULL;
    (*md)[i].aeinst  = NULL;
    (*md)[i].gir     = NULL;
    (*md)[i].freq    = NULL;
    (*md)[i].beinstu = NULL;
    (*md)[i].beinstl = NULL;
    (*md)[i].eterm   = NULL;
    (*md)[i].gstat   = NULL;
    (*md)[i].cmb     = NULL;
  }
}

/*....................................................................*/
void calcSourceFn(double dTau, const configInfo *par, double *remnantSnu, double *expDTau){
  /*
  The source function S is defined as j_nu/alpha, which is clearly not
  defined for alpha==0. However S is used in the algorithm only in the
  term (1-exp[-alpha*ds])*S, which is defined for all values of alpha.
  The present function calculates this term and returns it in the
  argument remnantSnu. For values of abs(alpha*ds) less than a pre-
  calculated cutoff supplied in configInfo, a Taylor approximation is
  used.

  Note that the same cutoff condition holds for replacement of
  exp(-dTau) by its Taylor expansion to 3rd order.

  Note that this is called from within the multi-threaded block.
  */

#ifdef FASTEXP
  *expDTau = FastExp(dTau);
  if (fabs(dTau)<par->taylorCutoff){
    *remnantSnu = 1. - dTau*(1. - dTau*(1./3.))*(1./2.);
  } else {
    *remnantSnu = (1.-(*expDTau))/dTau;
  }
#else
  if (fabs(dTau)<par->taylorCutoff){
    *remnantSnu = 1. - dTau*(1. - dTau*(1./3.))*(1./2.);
    *expDTau = 1. - dTau*(*remnantSnu);
  } else {
    *expDTau = exp(-dTau);
    *remnantSnu = (1.-(*expDTau))/dTau;
  }
#endif
}

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
double
gaussline(const double v, const double oneOnSigma){
  double val;
  val = v*v*oneOnSigma*oneOnSigma;
#ifdef FASTEXP
  return FastExp(val);
#else
  return exp(-val);
#endif
}

/*....................................................................*/
double
dotProduct3D(const double *vA, const double *vB){
  return vA[0]*vB[0] + vA[1]*vB[1] + vA[2]*vB[2];
}

/*....................................................................*/
void
copyInparStr(const char *inStr, char **outStr){
  if(inStr==NULL || strlen(inStr)<=0 || strlen(inStr)>STR_LEN_0){
    *outStr = NULL;
  }else{
    *outStr = malloc(sizeof(**outStr)*(STR_LEN_0+1));
    strcpy(*outStr, inStr);
  }
}

/*....................................................................*/
_Bool
charPtrIsNullOrEmpty(const char *inStr){
  if(inStr==NULL || strlen(inStr)<=0)
    return TRUE;
  else
    return FALSE;
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

/*....................................................................*/
_Bool allBitsSet(const int flags, const int mask){
  /* Returns true only if all the masked bits of flags are set. */

  if(~flags & mask)
    return 0;
  else
    return 1;
}

/*....................................................................*/
_Bool anyBitSet(const int flags, const int mask){
  /* Returns true if any of the masked bits of flags are set. */

  if(flags & mask)
    return 1;
  else
    return 0;
}

/*....................................................................*/
_Bool bitIsSet(const int flags, const int bitI){
  /* Returns true if the designated bit of flags is set. */

  if(flags & (1 << bitI))
    return 1;
  else
    return 0;
}

/*....................................................................*/
_Bool onlyBitsSet(const int flags, const int mask){
  /* Returns true if flags has no bits set apart from those which are true in mask. */

  if(flags & ~mask)
    return 0;
  else
    return 1;
}

