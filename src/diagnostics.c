/*
 *  diagnostics.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#include "lime.h"

/*....................................................................*/
void
diagPrintInput(const inputPars inpars, image *inimg, const int nImages){
  /*
This is intended purely as a diagnostic program to compare the parameter values passed in by the various wrappers to LIME.
  */
  int i,j;

  printf("              radius = %e\n", inpars.radius);
  printf("            minScale = %e\n", inpars.minScale);
  printf("                tcmb = %e\n", inpars.tcmb);

  for(i=0;i<MAX_N_COLL_PART;i++){
    if(inpars.nMolWeights[i]>=0.0)
      printf("       nMolWeights[%d] = %e\n", i, inpars.nMolWeights[i]);
    if(inpars.dustWeights[i]>=0.0)
      printf("       dustWeights[%d] = %e\n", i, inpars.dustWeights[i]);
    if(inpars.collPartMolWeights[i]>=0.0)
      printf("collPartMolWeights[%d] = %e\n", i, inpars.collPartMolWeights[i]);
    if(inpars.collPartIds[i]>0)
      printf("       collPartIds[%d] = %d\n", i, inpars.collPartIds[i]);
  }

  for(i=0;i<MAX_N_HIGH;i++){
    if(inpars.gridDensMaxValues[i]<0.0)
      break;

    printf(" gridDensMaxValues[%d] = %e\n", i, inpars.gridDensMaxValues[i]);
    printf("    gridDensMaxLoc[%d] = [", i);
    for(j=0;j<DIM;j++){
      printf("%e,", inpars.gridDensMaxLoc[i][j]);
    }
    printf("]\n");
  }

  if(i<=0) printf(" gridDensMaxValues unset.\n");

  printf("          sinkPoints = %d\n", inpars.sinkPoints);
  printf("          pIntensity = %d\n", inpars.pIntensity);
  printf("               blend = %d\n", inpars.blend);
  printf("   traceRayAlgorithm = %d\n", inpars.traceRayAlgorithm);
  printf("   samplingAlgorithm = %d\n", inpars.samplingAlgorithm);
  printf("            sampling = %d\n", inpars.sampling);
  printf("            lte_only = %d\n", inpars.lte_only);
  printf("            init_lte = %d\n", inpars.init_lte);
  printf("           antialias = %d\n", inpars.antialias);
  printf("        polarization = %d\n", inpars.polarization);
  printf("            nThreads = %d\n", inpars.nThreads);
  printf("         nSolveIters = %d\n", inpars.nSolveIters);

  if(inpars.moldatfile!=NULL && inpars.girdatfile!=NULL){
    for(i=0;i<MAX_NSPECIES;i++){
      if(!charPtrIsNullOrEmpty(inpars.moldatfile[i]))
        printf("        moldatfile[%d] = %s\n", i, inpars.moldatfile[i]);

      if (!charPtrIsNullOrEmpty(inpars.girdatfile[i]))
        printf("        girdatfile[%d] = %s\n", i, inpars.girdatfile[i]);
    }
  }else{
    if(inpars.moldatfile==NULL) printf("            moldatfile = NULL\n");
    if(inpars.girdatfile==NULL) printf("            girdatfile = NULL\n");
  }

  if(inpars.collPartNames!=NULL){
    for(i=0;i<MAX_N_COLL_PART;i++){
      if (!charPtrIsNullOrEmpty(inpars.collPartNames[i]))
        printf("     collPartNames[%d] = %s\n", i, inpars.collPartNames[i]);
    }
  }else
    printf("         collPartNames = NULL\n");

  if (   inpars.outputfile!=NULL)
    printf("          outputfile = %s\n", inpars.outputfile);
  else
    printf("          outputfile = NULL\n");

  if (inpars.binoutputfile!=NULL)
    printf("       binoutputfile = %s\n", inpars.binoutputfile);
  else
    printf("       binoutputfile = NULL\n");

  if (     inpars.gridfile!=NULL)
    printf("            gridfile = %s\n", inpars.gridfile);
  else
    printf("            gridfile = NULL\n");

  if (      inpars.pregrid!=NULL)
    printf("             pregrid = %s\n", inpars.pregrid);
  else
    printf("             pregrid = NULL\n");

  if (      inpars.restart!=NULL)
    printf("             restart = %s\n", inpars.restart);
  else
    printf("             restart = NULL\n");

  if (         inpars.dust!=NULL)
    printf("                dust = %s\n", inpars.dust);
  else
    printf("                dust = NULL\n");

  if (   inpars.gridInFile!=NULL)
    printf("          gridInFile = %s\n", inpars.gridInFile);
  else
    printf("          gridInFile = NULL\n");

  if(inpars.gridOutFiles!=NULL){
    for(i=0;i<NUM_GRID_STAGES;i++){
      if (inpars.gridOutFiles[i]!=NULL)
        printf("    gridOutFiles[%d] = %s\n", i, inpars.gridOutFiles[i]);
      else
        printf("    gridOutFiles[%d] = NULL\n", i);
    }
  }else
    printf("        gridOutFiles = NULL\n");

  if(inpars.resetRNG)
    printf("            resetRNG = TRUE\n");
  else
    printf("            resetRNG = FALSE\n");

  if(inpars.doSolveRTE)
    printf("          doSolveRTE = TRUE\n");
  else
    printf("          doSolveRTE = FALSE\n");

  for(i=0;i<nImages;i++){
    printf("\n");
    printf("Image %d\n", i);
    printf("\n");

    printf("              velres = %e\n", inimg[i].velres);
    printf("              imgres = %e\n", inimg[i].imgres);
    printf("                freq = %e\n", inimg[i].freq);
    printf("           bandwidth = %e\n", inimg[i].bandwidth);
    printf("          source_vel = %e\n", inimg[i].source_vel);
    printf("               theta = %e\n", inimg[i].theta);
    printf("                 phi = %e\n", inimg[i].phi);
    printf("                incl = %e\n", inimg[i].incl);
    printf("              posang = %e\n", inimg[i].posang);
    printf("             azimuth = %e\n", inimg[i].azimuth);

    printf("               nchan = %d\n", inimg[i].nchan);
    printf("               trans = %d\n", inimg[i].trans);
    printf("                molI = %d\n", inimg[i].molI);
    printf("                pxls = %d\n", inimg[i].pxls);
    printf("                unit = %d\n", inimg[i].unit);

    if (   inimg[i].units!=NULL)
      printf("               units = %s\n", inimg[i].units);
    else
      printf("               units = NULL\n");

    if (inimg[i].filename!=NULL)
      printf("            filename = %s\n", inimg[i].filename);
    else
      printf("            filename = NULL\n");

    if(inimg[i].doInterpolateVels)
      printf("   doInterpolateVels = TRUE\n");
    else
      printf("   doInterpolateVels = FALSE\n");
  }
}

/*....................................................................*/
void diagPrintTestGridOutput(struct grid *gp, const int id, const int i\
  , const int ilev){

  if(gp[id].mol==NULL){
    printf("gp[%d].mol==NULL!\n", id);
  }else{
    printf("gp[%d].mol[%d].binv=%e\n", id, i, gp[id].mol[i].binv);
    printf("gp[%d].mol[%d].nmol=%e\n", id, i, gp[id].mol[i].nmol);

    if(gp[id].mol[i].pops==NULL)
      printf("gp[%d].mol[%d].pops==NULL\n", id, i);
    else
      printf("gp[%d].mol[%d].pops[%d]=%e\n", id, i, ilev, gp[id].mol[i].pops[ilev]);

    if(gp[id].mol[i].specNumDens==NULL)
      printf("gp[%d].mol[%d].specNumDens==NULL\n", id, i);
    else
      printf("gp[%d].mol[%d].specNumDens[%d]=%e\n", id, i, ilev, gp[id].mol[i].specNumDens[ilev]);

    if(gp[id].mol[i].cont==NULL)
      printf("gp[%d].mol[%d].cont==NULL\n", id, i);
    else{
      printf("gp[%d].mol[%d].cont[%d].dust=%e\n", id, i, ilev, gp[id].mol[i].cont[ilev].dust);
      printf("gp[%d].mol[%d].cont[%d].knu=%e\n",  id, i, ilev, gp[id].mol[i].cont[ilev].knu);
    }
  }
  printf("\n");
}

