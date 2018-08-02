/*
 *  run_new.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#include "lime.h"
#include <locale.h>

#ifdef TEST
_Bool fixRandomSeeds = TRUE;
#else
_Bool fixRandomSeeds = FALSE;
#endif

double defaultDensyPower = DENSITY_POWER;

/*....................................................................*/
int
run_new(inputPars inpars, image *inimg, const int nImages){
  /* Run LIME with inpars and the output fits files specified.

     This routine may be used as an interface to LIME from external
     programs. In this case, inpars and img must be specified by the
     external program.
  */
  int nSpecies,i,status=0,sigactionStatus=0,numCollPartRead=0;
  int initime=time(0);
  int dummyPopsdone=0,nExtraSolverIters=0;
  molData *md=NULL;
  configInfo par;
  imageInfo *img=NULL;
  struct grid *gp=NULL;
  char message[STR_LEN_1+1];
  int nEntries=0;
  double *lamtab=NULL,*kaptab=NULL;
  char **collPartNamesRead=NULL;

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

  if(par.gridInFile!=NULL)
    readGridWrapper(&par, &gp, &collPartNamesRead, &numCollPartRead); /* In gridio.c. Sets values in par.dataFlags. */

  parChecks(&par);
  parseInputWhenIncompleteGridFile(&par, TRUE); /* In init.c. Sets par.numDensities */
  parseImagePars(&par, &img); /* In init.c */
  furtherParChecks(&par); /* In init.c */

  /* Allocate moldata array.
  */
  if(par.nSpecies>0){
    mallocAndSetDefaultMolData(par.nSpecies, &md);
  } /* otherwise leave it at NULL - we will not be using it. */

  if(!silent && par.nThreads>1){
    snprintf(message, STR_LEN_1, "Number of threads used: %d", par.nThreads);
    printMessage(message);
  }

  checkUserDensWeights(&par); /* In collparts.c. Needs par.numDensities. */
  buildGrid(&par,&gp);

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
    molInit(&par, md);
    calcGridMolDoppler(&par, md, gp);

    if(par.useAbun)
      calcGridMolDensities(&par, &gp);
  }

  if(par.needToInitPops)
    gridPopsInit(&par,md,gp);

  if(par.needToInitSND)
    specNumDensInit(&par,md,gp);

  if(par.doSolveRTE){
    nExtraSolverIters = levelPops(md, &par, gp, &dummyPopsdone, lamtab, kaptab, nEntries); /* In solver.c */
    par.nSolveItersDone += nExtraSolverIters;
  }

  if(onlyBitsSet(par.dataFlags & DS_mask_all_but_mag, DS_mask_3))
    writeGridIfRequired(&par, gp, NULL, 3);
  else if(onlyBitsSet(par.dataFlags & DS_mask_all_but_mag, DS_mask_5)){
    writeGridIfRequired(&par, gp, md, 5);
  }else if(!silent){
    snprintf(message, STR_LEN_1, "Data flags %x match neither mask 3 %x (cont.) or 5 %x (line).", par.dataFlags, DS_mask_3, DS_mask_5);
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

  freeArrayOfStrings(collPartNamesRead, numCollPartRead);

  return status; /* This is a bit of a placeholder for now. Ideally we would like all the functions called to return status values rather than exiting. This would allow python-calling versions of Lime to exit 'nicely' at the top level. */
}
