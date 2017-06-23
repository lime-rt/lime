/*
 *  frees.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2017 The LIME development team
 *
 */

#include "lime.h"

void
freeConfigInfo(configInfo *par){
  int i;

  free(par->nMolWeights);
  free(par->collPartMolWeights);
  free(par->collPartIds);
  free(par->outputfile);
  free(par->binoutputfile);
  free(par->gridfile);
  free(par->pregrid);
  free(par->restart);
  free(par->dust);
  free(par->gridDensMaxValues);
  free(par->gridDensMaxLoc);
  free(par->gridInFile);

  if(par->collPartNames!= NULL){
    for(i=0;i<par->numDensities;i++)
      free(par->collPartNames[i]);
    free(par->collPartNames);
  }
  if(par->moldatfile!= NULL){
    for(i=0;i<par->nSpecies;i++)
      free(par->moldatfile[i]);
    free(par->moldatfile);
  }
  if(par->girdatfile!= NULL){
    for(i=0;i<par->nSpecies;i++)
      free(par->girdatfile[i]);
    free(par->girdatfile);
  }
  if(par->gridOutFiles!= NULL){
    for(i=0;i<NUM_GRID_STAGES;i++)
      free(par->gridOutFiles[i]);
    free(par->gridOutFiles);
  }
}

void
freeInputPars(inputPars *par){
  free(par->collPartIds);
  free(par->nMolWeights);
  free(par->dustWeights);
  free(par->collPartMolWeights);
  free(par->moldatfile);
  free(par->girdatfile);
  free(par->collPartNames);
  free(par->gridOutFiles);
  free(par->gridDensMaxValues);
  free(par->gridDensMaxLoc);
}

void
freeGrid(const unsigned int numPoints, const unsigned short numSpecies\
  , struct grid *gp){

  unsigned int i_u;

  if(gp != NULL){
    for(i_u=0;i_u<numPoints;i_u++){
      free(gp[i_u].v1);
      free(gp[i_u].v2);
      free(gp[i_u].v3);
      free(gp[i_u].dir);
      free(gp[i_u].neigh);
      free(gp[i_u].w);
      free(gp[i_u].dens);
      free(gp[i_u].ds);
      freePopulation(numSpecies, gp[i_u].mol);
    }
    free(gp);
  }
}

void
freeGridPointData(const int nSpecies, gridPointData *mol){
  int i;
  if(mol!= NULL){
    for(i=0;i<nSpecies;i++){
      free(mol[i].jbar);
      free(mol[i].phot);
      free(mol[i].vfac);
      free(mol[i].vfac_loc);
    }
  }
}

void
freeImgInfo(const int nImages, imageInfo *img){
  int i,id;
  for(i=0;i<nImages;i++){
    for(id=0;id<(img[i].pxls*img[i].pxls);id++){
      free( img[i].pixel[id].intense );
      free( img[i].pixel[id].tau );
    }
    free(img[i].pixel);
    free(img[i].filename);
    free(img[i].imgunits);
  }
  free(img);
}

void
freeMolData(const int nSpecies, molData *md){
  int i,j;
  if(md!= NULL){
    for(i=0;i<nSpecies;i++){
      if(md[i].part != NULL){
        for(j=0; j<md[i].npart; j++){
          free(md[i].part[j].down);
          free(md[i].part[j].temp);
          free(md[i].part[j].lcl);
          free(md[i].part[j].lcu);
          free(md[i].part[j].name);
        }
        free(md[i].part);
      }
      free(md[i].lal);
      free(md[i].lau);
      free(md[i].aeinst);
      free(md[i].gir);
      free(md[i].freq);
      free(md[i].beinstu);
      free(md[i].beinstl);
      free(md[i].eterm);
      free(md[i].gstat);
      free(md[i].cmb);
    }
    free(md);
  }
}

void
freeMolsWithBlends(struct molWithBlends *mols, const int numMolsWithBlends){
  int mi, li;

  if(mols != NULL){
    for(mi=0;mi<numMolsWithBlends;mi++){
      if(mols[mi].lines != NULL){
        for(li=0;li<mols[mi].numLinesWithBlends;li++)
          free(mols[mi].lines[li].blends);
        free(mols[mi].lines);
      }
    }
    free(mols);
  }
}

void
freePopulation(const unsigned short numSpecies, struct populations *pop){
  if(pop != NULL){
    unsigned short i_s;
    for(i_s=0;i_s<numSpecies;i_s++){
      free(pop[i_s].pops);
      free(pop[i_s].partner);
      free(pop[i_s].specNumDens);
      free(pop[i_s].cont);
    }
    free(pop);
  }
}

void
freeSomeGridFields(const unsigned int numPoints, const unsigned short numSpecies\
  , struct grid *gp){

  unsigned int i_u;
  unsigned short i_s;

  if(gp != NULL){
    for(i_u=0;i_u<numPoints;i_u++){
      free(gp[i_u].w);
      gp[i_u].w    = NULL;
      free(gp[i_u].ds);
      gp[i_u].ds   = NULL;

      if(gp[i_u].mol != NULL){
        for(i_s=0;i_s<numSpecies;i_s++){
          free(gp[i_u].mol[i_s].pops);
          gp[i_u].mol[i_s].pops = NULL;
          free(gp[i_u].mol[i_s].partner);
          gp[i_u].mol[i_s].partner = NULL;
          free(gp[i_u].mol[i_s].cont);
          gp[i_u].mol[i_s].cont = NULL;
        }
      }
    }
  }
}

