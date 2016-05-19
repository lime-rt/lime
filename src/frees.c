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
freeInput( inputPars *par, image* img )
{
  int i,id;
  for(i=0;i<par->nImages;i++){
    for(id=0;id<(img[i].pxls*img[i].pxls);id++){
      free( img[i].pixel[id].intense );
      free( img[i].pixel[id].tau );
    }
    free(img[i].pixel);
  }
  if( img != NULL )
    {
      free(img);
    }
  if( par->moldatfile != NULL )
    {
      free(par->moldatfile);
    }
}

void
freeMolData( inputPars *par, molData* mol )
{
  int i,j;
  if( mol!= 0 )
    {
      for( i=0; i<par->nSpecies; i++ )
        {
          if( mol[i].part != NULL )
            {

              for( j=0; j<mol[i].npart; j++ )
                {
                  if( mol[i].part[j].lcl != NULL )
                    {
                      free(mol[i].part[j].lcl);
                    }
                  if( mol[i].part[j].lcu != NULL )
                    {
                      free(mol[i].part[j].lcu);
                    }
                }
              free(mol[i].part);
            }
          if( mol[i].lal != NULL )
            {
              free(mol[i].lal);
            }
          if( mol[i].lau != NULL )
            {
              free(mol[i].lau);
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


void
freePopulation(const inputPars *par, const molData* m, struct populations* pop ) {
  if( pop !=NULL )
    {
      int j,k;
      for( j=0; j<par->nSpecies; j++ )
        {
          if( pop[j].pops != NULL )
            {
              free( pop[j].pops );
            }
          if( pop[j].knu != NULL )
            {
              free( pop[j].knu );
            }
          if( pop[j].dust != NULL )
            {
              free( pop[j].dust );
            }
          if( pop[j].partner != NULL )
            {
              if( m != NULL )
                {
                  for(k=0; k<m[j].npart; k++)
                    {
                      if( pop[j].partner[k].up != NULL )
                        {
                          free(pop[j].partner[k].up);
                        }
                      if( pop[j].partner[k].down != NULL )
                        {
                          free(pop[j].partner[k].down);
                        }
                    }
                }
              free( pop[j].partner );
            }
        }
      free(pop);
    }
}
void
freeGrid(const inputPars *par, const molData* m ,struct grid* g){
  int i;
  if( g != NULL )
    {
      for(i=0;i<(par->pIntensity+par->sinkPoints); i++){
        if(g[i].a0 != NULL)
          {
            free(g[i].a0);
          }
        if(g[i].a1 != NULL)
          {
            free(g[i].a1);
          }
        if(g[i].a2 != NULL)
          {
            free(g[i].a2);
          }
        if(g[i].a3 != NULL)
          {
            free(g[i].a3);
          }
        if(g[i].a4 != NULL)
          {
            free(g[i].a4);
          }
        if(g[i].dir != NULL)
          {
            free(g[i].dir);
          }
        if(g[i].neigh != NULL)
          {
            free(g[i].neigh);
          }
        if(g[i].w != NULL)
          {
            free(g[i].w);
          }
        if(g[i].dens != NULL)
          {
            free(g[i].dens);
          }
        if(g[i].nmol != NULL)
          {
            free(g[i].nmol);
          }
        if(g[i].abun != NULL)
          {
            free(g[i].abun);
          }
        if(g[i].ds != NULL)
          {
            free(g[i].ds);
          }
        if(g[i].mol != NULL)
          {
            freePopulation( par, m, g[i].mol );
          }
      }
      free(g);
    }
}


