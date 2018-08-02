/*
 *  run.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#include "lime.h"
#include <locale.h>
#include "gridio.h" /* For countDensityCols() */

/*....................................................................*/
void
_parseInput_rump(configInfo *par, _Bool checkForSingularities){
  int numFuncDensities;

  if(par->pregrid==NULL && par->restart){
    par->nSpecies=0; /* This will get set during popsin(). */
    par->girdatfile = NULL;
  }

  numFuncDensities = checkUserFunctions(par, checkForSingularities);
  par->numDensities = numFuncDensities;

  if(!(par->doPregrid || par->restart || par->gridInFile!=NULL)){
    calcGridDensGlobalMax(par);
  }

#ifndef KLUDGE_FOR_BAD_DENS
  if(par->samplingAlgorithm==0)
    defaultDensyPower = DENSITY_POWER;
  else
    defaultDensyPower = TREE_POWER;
#endif

}

/*....................................................................*/
void
_furtherParChecks(configInfo *par){
  FILE *fp;
  char message[STR_LEN_1+1];

  if(par->nContImages>0){
    if(par->dust==NULL){
      if(!silent) bail_out("You must point par->dust to a dust opacity file for a continuum image.");
exit(1);
    }else{
      if((fp=fopen(par->dust, "r"))==NULL){
        if(!silent){
          snprintf(message, STR_LEN_1, "Couldn't open dust opacity data file %s", par->dust);
          bail_out(message);
        }
exit(1);
      }
      fclose(fp);
    }
  }

  par->useVelFuncInRaytrace = FALSE; /* the default */
  if(par->nLineImages>0 && par->traceRayAlgorithm==0 && !par->doPregrid){
    if(bitIsSet(defaultFuncFlags, USERFUNC_velocity)){
      if(!silent) warning("No velocity function supplied - raytracing will have lower precision.");
    }else
      par->useVelFuncInRaytrace = TRUE;
  }

#ifdef IS_PYTHON
  if(par->nThreads>1 && par->useVelFuncInRaytrace){
    par->useVelFuncInRaytrace = FALSE;
    if(!silent)
      warning("You cannot call a python velocity function when multi-threaded.");
  }
#endif

  par->edgeVelsAvailable=0; /* default value, this is set within getEdgeVelocities(). */

  if(!silent){
    if(par->lte_only && par->nSolveIters>0)
      warning("Requesting par->nSolveIters>0 will have no effect if LTE calculation is also requested.");

    else if(par->nSolveIters<=par->nSolveItersDone && !allBitsSet(par->dataFlags, DS_mask_populations))
      warning("No supplied pops values, and par->nSolveIters <= par->nSolveItersDone.");
  }

  if(par->restart){ /* pops will be read from file */
    par->doSolveRTE = FALSE;
    par->doMolCalcs = (par->nLineImages>0);
    par->popsHasBeenInit = TRUE;

  }else{
    if(par->nSolveIters>par->nSolveItersDone || par->lte_only) /* To save the user having to set par->doSolveRTE as well as par->nSolveIters>0 or par->lte_only. */
      par->doSolveRTE = TRUE;

    par->popsHasBeenInit = FALSE;
  }

  if(!par->popsHasBeenInit && par->doSolveRTE)
    par->needToInitPops = TRUE;
  else
    par->needToInitPops = FALSE;

  if((!par->lte_only && par->doSolveRTE) || par->nLineImages>0)
    par->needToInitSND = TRUE;
  else
    par->needToInitSND = FALSE;

  par->SNDhasBeenInit = FALSE; /* default */

  par->doMolCalcs = par->doSolveRTE || par->nLineImages>0;
  if(par->doMolCalcs && par->moldatfile==NULL){
    if(!silent) bail_out("You must point par->moldatfile to a data file.");
exit(1);
  }

  if(!silent && !par->doMolCalcs && par->init_lte)
    warning("Your choice of par->init_lte will have no effect.");

  if(par->nSpecies>0 && !par->doMolCalcs){
    if(!silent) bail_out("If you want only continuum calculations you must supply zero moldatfiles.");
exit(1);
  }
}

/*....................................................................*/
int
run(inputPars inpars, image *inimg, const int nImages){
  /* Run LIME with inpars and the output fits files specified.

     This routine may be used as an interface to LIME from external
     programs. In this case, inpars and img must be specified by the
     external program.
  */
  int nSpecies,i,status=0,sigactionStatus=0;
  int initime=time(0);
  int popsdone=0,nExtraSolverIters=0;
  molData *md=NULL;
  configInfo par;
  imageInfo *img=NULL;
  struct grid *gp=NULL;
  char message[STR_LEN_1+1];
  int nEntries=0;
  double *lamtab=NULL,*kaptab=NULL;

  if(inpars.pregrid==NULL && inpars.restart==NULL){
    return run_new(inpars, inimg, nImages);
  }

  struct sigaction sigact = {.sa_handler = sigintHandler};
  sigactionStatus = sigaction(SIGINT, &sigact, NULL);
  if(sigactionStatus){
    if(!silent){
      snprintf(message, STR_LEN_1, "Call to sigaction() returned with status %d", sigactionStatus);
      bail_out(message);
    }
exit(1);
  }

  /*Set locale to avoid trouble when reading files*/
  setlocale(LC_ALL, "C");

  if(!silent) greetings(VERSION);
  if(!silent) screenInfo();

#ifdef FASTEXP
  calcExpTableEntries(FAST_EXP_MAX_TAYLOR, FAST_EXP_NUM_BITS);
#endif
  fillErfTable();

  nSpecies = copyInpars(inpars, inimg, nImages, &par, &img); /* In init.c */
  par.nSpecies = nSpecies;
  setOtherEasyConfigValues(nImages, &par, &img); /* In init.c */

  //*** read files here instead of later? That would set par.dataFlags and maybe allow use of fewer bespoke par check routines.

  parChecks(&par);
  _parseInput_rump(&par, TRUE); /* Sets par.numDensities for !(par.doPregrid || par.restart) */
  parseImagePars(&par, &img); /* In init.c */
  _furtherParChecks(&par);

  /* Allocate moldata array.
  */
  if(par.nSpecies>0){
    mallocAndSetDefaultMolData(par.nSpecies, &md);
  } /* otherwise leave it at NULL - we will not be using it. */

  if(!silent && par.nThreads>1){
    snprintf(message, STR_LEN_0, "Number of threads used: %d", par.nThreads);
    printMessage(message);
  }

  if(par.doPregrid){
    mallocAndSetDefaultGrid(&gp, (size_t)par.ncell, (size_t)par.nSpecies);
    predefinedGrid(&par,gp); /* Sets par.numDensities. */
    checkUserDensWeights(&par); /* In collparts.c. Needs par.numDensities. */

  }else{ /* must mean par.restart */
    popsin(&par,&gp,&md,&popsdone);
  }

  if(par.dust != NULL)
    readDustFile(par.dust, &lamtab, &kaptab, &nEntries);

  /* Make all the continuum images:
  */
  if(par.nContImages>0){
    for(i=0;i<par.nImages;i++){
      if(!img[i].doline){
        raytrace(i, &par, gp, md, img, lamtab, kaptab, nEntries);
        writeFitsAllUnits(i, &par, img);
      }
    }
  }

  if(par.doMolCalcs){
    if(!popsdone){
      molInit(&par, md);
      calcGridMolDoppler(&par, md, gp);
    }
    if(par.useAbun)
      calcGridMolDensities(&par, &gp);
  }

  if(par.needToInitPops)
    gridPopsInit(&par,md,gp);

  if(par.needToInitSND)
    specNumDensInit(&par,md,gp);

  if(par.doSolveRTE){
    nExtraSolverIters = levelPops(md, &par, gp, &popsdone, lamtab, kaptab, nEntries);
    par.nSolveItersDone += nExtraSolverIters;
  }

  if(onlyBitsSet(par.dataFlags & DS_mask_all_but_mag, DS_mask_3))
    writeGridIfRequired(&par, gp, NULL, 3);
  else if(onlyBitsSet(par.dataFlags & DS_mask_all_but_mag, DS_mask_5)){
    writeGridIfRequired(&par, gp, md, 5);
  }else if(!silent){
    snprintf(message, STR_LEN_0, "Data flags %x match neither mask 3 %x (cont.) or 5 %x (line).", par.dataFlags, DS_mask_3, DS_mask_5);
    warning(message);
  }

  if(par.nLineImages>0)
    calcGridMolSpecNumDens(&par,md,gp); /* In grid_aux.c */

  freeSomeGridFields((unsigned int)par.ncell, (unsigned short)par.nSpecies, gp);

  /* Now make the line images.
  */
  if(par.nLineImages>0){
    for(i=0;i<par.nImages;i++){
      if(img[i].doline){
        raytrace(i, &par, gp, md, img, lamtab, kaptab, nEntries);
        writeFitsAllUnits(i, &par, img);
      }
    }
  }

  if(!silent){
    if(par.nImages>0) reportOutput(img[0].filename);
    goodnight(initime);
  }

  freeGrid((unsigned int)par.ncell, (unsigned short)par.nSpecies, gp);
  freeMolData(par.nSpecies, md);
  freeImgInfo(par.nImages, img);
  freeConfigInfo(&par);

  if(par.dust != NULL){
    free(kaptab);
    free(lamtab);
  }

  return status; /* This is a bit of a placeholder for now. Ideally we would like all the functions called to return status values rather than exiting. This would allow python-calling versions of Lime to exit 'nicely' at the top level. */
}

