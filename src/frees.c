/*
 *  frees.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2016 The LIME development team
 *
 */

#include "lime.h"

void freeGrid(const unsigned int numPoints, const unsigned short numSpecies\
  , struct grid *gp){

  unsigned int i_u;

  if(gp != NULL){
    for(i_u=0;i_u<numPoints;i_u++){
      free(gp[i_u].a0);
      free(gp[i_u].a1);
      free(gp[i_u].a2);
      free(gp[i_u].a3);
      free(gp[i_u].a4);
      free(gp[i_u].dir);
      free(gp[i_u].neigh);
      free(gp[i_u].w);
      free(gp[i_u].dens);
      free(gp[i_u].abun);
      free(gp[i_u].ds);
      freePopulation(numSpecies, gp[i_u].mol);
    }
    free(gp);
  }
}

void
freeGridPointData(configInfo *par, gridPointData *mol){
  int i;
  if(mol!= NULL){
    for(i=0;i<par->nSpecies;i++){
      free(mol[i].jbar);
      free(mol[i].phot);
      free(mol[i].vfac);
    }
    free(mol);
  }
}

void
freeMolData(const int nSpecies, molData *mol){
  int i,j;
  if(mol!= NULL){
    for(i=0;i<nSpecies;i++){
      if(mol[i].part != NULL){
        for(j=0; j<mol[i].npart; j++){
          free(mol[i].part[j].down);
          free(mol[i].part[j].temp);
          free(mol[i].part[j].lcl);
          free(mol[i].part[j].lcu);
        }
        free(mol[i].part);
      }
      free(mol[i].lal);
      free(mol[i].lau);
      free(mol[i].aeinst);
      free(mol[i].freq);
      free(mol[i].beinstu);
      free(mol[i].beinstl);
      free(mol[i].eterm);
      free(mol[i].gstat);
      free(mol[i].cmb);
    }
    free(mol);
  }
}

void freeMolsWithBlends(struct molWithBlends *mols, const int numMolsWithBlends){
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
freeParImg(const int nImages, inputPars *par, image *img){
  /* Release memory allocated for the output fits images
     and for inputPars pointers.
  */
  int i,id;
  for(i=0;i<nImages;i++){
    for(id=0;id<(img[i].pxls*img[i].pxls);id++){
      free( img[i].pixel[id].intense );
      free( img[i].pixel[id].tau );
    }
    free(img[i].pixel);
  }
  free(img);

  free(par->moldatfile);
  free(par->collPartIds);
  free(par->nMolWeights);
  free(par->dustWeights);
}

void freePopulation(const unsigned short numSpecies, struct populations *pop){
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

void freeSomeGridFields(const unsigned int numPoints, const unsigned short numSpecies\
  , struct grid *gp){

  unsigned int i_u;
  unsigned short i_s;

  if(gp != NULL){
    for(i_u=0;i_u<numPoints;i_u++){
      free(gp[i_u].a0);
      gp[i_u].a0 = NULL;
      free(gp[i_u].a1);
      gp[i_u].a1 = NULL;
      free(gp[i_u].a2);
      gp[i_u].a2 = NULL;
      free(gp[i_u].a3);
      gp[i_u].a3 = NULL;
      free(gp[i_u].a4);
      gp[i_u].a4 = NULL;
      free(gp[i_u].w);
      gp[i_u].w    = NULL;
      free(gp[i_u].abun);
      gp[i_u].abun = NULL;
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

