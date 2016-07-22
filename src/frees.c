/*
 *  frees.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015 The LIME development team
 *
 */

#include "lime.h"

void
freeGAux(const unsigned long numPoints, const int numSpecies, struct gAuxType *gAux) {
  unsigned long ppi;

  if(gAux !=NULL){
    for(ppi=0;ppi<numPoints;ppi++)
      freePop2(numSpecies, gAux[ppi].mol);
    free(gAux);
  }
}

void
freeGrid(configInfo *par, molData *m , struct grid *g){
  int i;
  if(g != NULL){
    for(i=0;i<(par->pIntensity+par->sinkPoints); i++){
      free(g[i].a0);
      free(g[i].a1);
      free(g[i].a2);
      free(g[i].a3);
      free(g[i].a4);
      free(g[i].dir);
      free(g[i].neigh);
      free(g[i].w);
      free(g[i].dens);
      free(g[i].abun);
      free(g[i].ds);
      freePopulation( par, m, g[i].mol );
    }
    free(g);
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
      free(mol[i].local_cmb);
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

void
freePop2(const int numSpecies, struct pop2 *mol){
  int i;

  if(mol !=NULL){
    for(i=0;i<numSpecies;i++){
      free(mol[i].specNumDens);
      free(mol[i].knu);
      free(mol[i].dust);
    }
    free(mol);
  }
}

void
freePopulation(configInfo *par, molData *m, struct populations *pop ){
  if(pop != NULL){
    int j,k;
    for(j=0; j<par->nSpecies; j++){
      free( pop[j].pops );
      free( pop[j].knu );
      free( pop[j].dust );
      free( pop[j].partner );
    }
    free(pop);
  }
}


