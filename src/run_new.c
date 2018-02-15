/*
 *  run_new.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2017 The LIME development team
 *
TODO:
 */

#include "lime.h"
#include <locale.h>
#include "defaults.h"

#ifdef NO_STDOUT
#include "pyshared_io.h"
#endif

#ifdef TEST
_Bool fixRandomSeeds = TRUE;
#else
_Bool fixRandomSeeds = FALSE;
#endif

double defaultDensyPower = DENSITY_POWER;

/*....................................................................*/
void _gridPopsInit(configInfo *par, molData *md, struct grid *gp){
  int i,id,ilev;

  for(i=0;i<par->nSpecies;i++){
    /* Allocate space for populations etc */
    for(id=0;id<par->ncell; id++){
      free(gp[id].mol[i].pops);
      gp[id].mol[i].pops = malloc(sizeof(double)*md[i].nlev);
      for(ilev=0;ilev<md[i].nlev;ilev++)
        gp[id].mol[i].pops[ilev] = 0.0;
    }
  }
}

/*....................................................................*/
int
run_new(inputPars inpars, image *inimg, const int nImages){
  /* Run LIME with inpars and the output fits files specified.

     This routine may be used as an interface to LIME from external
     programs. In this case, inpars and img must be specified by the
     external program.
  */
  int i,gi,si,ei,status=0,sigactionStatus=0,numCollPartRead=0;
  int initime=time(0);
  int dummyPopsdone=0,nExtraSolverIters=0;
  molData *md=NULL;
  configInfo par;
  imageInfo *img=NULL;
  struct grid *gp=NULL;
  char message[STR_LEN_1];
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

  par.nSpecies = copyInpars(inpars, inimg, nImages, &par, &img);
  setOtherEasyConfigValues(nImages, &par, &img);

  if(par.gridInFile!=NULL)
    readGridWrapper(&par, &gp, &collPartNamesRead, &numCollPartRead); /* In gridio.c */

  parChecks(&par);
  parseInputWhenIncompleteGridFile(&par, TRUE); /* Sets par.numDensities */
  parseImagePars(&par, &img);

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

    for(gi=0;gi<par.ncell;gi++){
      for(si=0;si<par.nSpecies;si++){
        gp[gi].mol[si].specNumDens = malloc(sizeof(double)*md[si].nlev);
        for(ei=0;ei<md[si].nlev;ei++){
          gp[gi].mol[si].specNumDens[ei] = 0.0;
        }
      }
    }

    if(par.doSolveRTE){
      _gridPopsInit(&par,md,gp);
      nExtraSolverIters = levelPops(md, &par, gp, &dummyPopsdone, lamtab, kaptab, nEntries); /* In solver.c */
    }

    calcGridMolSpecNumDens(&par,md,gp); /* In grid_aux.c */

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
