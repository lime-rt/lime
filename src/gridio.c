/*
 *  gridio.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015 The LIME development team
 *
 */

#include "lime.h"

/*
This module contains generic routines for writing grid data to, and reading it from, a file on disk. The problem with doing this is that the grid struct contains different amounts of information at different times in the running of the code. In an attempt to regulate this, the following five states of completeness (encoded in values of the variable 'dataStageI') have been defined:

	dataStageI==0: No useable data in grid.

	dataStageI==1. At this stage the vector of grid objects has been malloc'd and values have been generated for the following struct elements:
		id
		x
		sink

	dataStageI==2. This stage is entered after the Delaunay neighbours of each grid point have been determined. The following further struct elements are expected to have been malloc'd (in the case of pointers) and given values:
		numNeigh
		neigh

	dataStageI==3. This is entered after sampling the user-supplied functions for density, velocity etc. The following further struct elements are expected to have been malloc'd (in the case of pointers) and given values:
		vel
		a0, a1 etc.
		dens
		t
		dopb
		abun
Note that it is necessary to make this stage dependent on the previous because we need information about the near-neighbours to calculate the a0, a1 etc coefficients.

	dataStageI==4. After 1 or more iterations of populating the levels, we enter the final stage, in which all the grid struct elements have been malloc'd and given at least preliminary values; specifically now the element
		mol

	*NOTE* that LIME may run to completion without ever reaching stage 4 - if all the images required were continuum ones, for example.


A note about the vectors 'links', 'nnLinks' and 'firstNearNeigh':
-----------------------------------------------------------------
The 'neigh' element of struct grid is a clunky way to store the information about the nearest neighbours of each grid point, because it leads to a host of separate little mallocs, thus memory for the grid is spread all over the place. There is also some duplication of quantities which refer to the link between 2 grid points rather than to the points themselves. This arrangement does not lend itself easily to storage in a file.

To meet the needs of file storage we divide the information contained in 'neigh' between 3 new vectors 'links', 'nnLinks' and 'firstNearNeigh'. The 1st of these is a list of all Delaunay edges between grid points. The 2nd contains indices to the 1st, and is arranged such that indices to all the links connected to a single grid point occur in a connected sequence. The 3rd vector has the same number of entries as the total number of grid points. It contains the index of the first entry in the 2nd vector which corresponds to that grid point. Thus in order, for example, to iterate over all links connected to the ith grid point, we need only do something as follows:

  for(j=0;j<g[i].numNeigh;j++){
    k = j + firstNearNeigh[i];
    link = links[nnLinks[k]];
  }

In addition, any link-specific stuff (at present just the 'a' coefficients used to interpolate the velocity along the length of the link) is stored in the first vector 'links'.

Since 'firstNearNeigh' has a single entry for each grid point, it can be stored with other like quantities in the GRID block. The 'links' and 'nnLinks' vectors have lengths which differ from this and thus must be stored in separate blocks.
*/

/*....................................................................*/
void writeGridIfRequired(inputPars *par, struct grid *gp, molData *md, const int fileFormatI){
  int status = 0;
  char **collPartNames=NULL; /*** this is a placeholder until we start reading these. */
  char message[80];

  if(par->writeGridAtStage[par->dataStageI-1]){
    status = writeGrid(par->gridOutFiles[par->dataStageI-1], fileFormatI\
      , *par, DIM, NUM_VEL_COEFFS, gp, md, collPartNames, par->dataStageI);

    if(status){
      sprintf(message, "writeGrid at data stage %d returned with status %d", par->dataStageI, status);
      if(!silent) bail_out(message);
      exit(1);
    }
  }

  free(collPartNames);
}

/*....................................................................*/
lime_fptr *openFileForRead(char *inFileName, const int fileFormatI, int *dataStageI){
  lime_fptr *fptr=NULL;

  if(fileFormatI==lime_FITS){
    fptr = openFITSFileForRead(inFileName, dataStageI);
  }else{
    return NULL;
  }

  return fptr;
}

/*....................................................................*/
lime_fptr *openFileForWrite(char *outFileName, const int fileFormatI, const int dataStageI){
  lime_fptr *fptr=NULL;

  if(fileFormatI==lime_FITS){
    fptr = openFITSFileForWrite(outFileName, dataStageI);
  }else{
    return NULL;
  }

  return fptr;
}

/*....................................................................*/
void closeAndFree(lime_fptr *fptr, const int fileFormatI\
  , unsigned int *firstNearNeigh, struct linkType **nnLinks\
  , struct linkType *links, const unsigned int totalNumLinks){

  unsigned int li;

  closeFile(fptr, fileFormatI);
  free(firstNearNeigh);
  free(nnLinks);
  if(links!=NULL){
    for(li=0;li<totalNumLinks;li++)
      free(links[li].aCoeffs);
   free(links);
  }
}

/*....................................................................*/
void freeReadGrid(struct grid **gp, struct gridInfoType gridInfoRead){
  unsigned int i_u;
  unsigned short i_s;

  if(*gp != NULL){
    for(i_u=0;i_u<(gridInfoRead.nInternalPoints+gridInfoRead.nSinkPoints);i_u++){
      free((*gp)[i_u].a0);
      free((*gp)[i_u].a1);
      free((*gp)[i_u].a2);
      free((*gp)[i_u].a3);
      free((*gp)[i_u].a4);
      free((*gp)[i_u].dir);
      free((*gp)[i_u].neigh);
      free((*gp)[i_u].w);
      free((*gp)[i_u].dens);
      free((*gp)[i_u].nmol);
      free((*gp)[i_u].abun);
      free((*gp)[i_u].ds);
      if((*gp)[i_u].mol != NULL){
        for(i_s=0;i_s<gridInfoRead.nSpecies;i_s++){
          free((*gp)[i_u].mol[i_s].pops );
          free((*gp)[i_u].mol[i_s].knu );
          free((*gp)[i_u].mol[i_s].dust );
        }
        free((*gp)[i_u].mol);
      }
    }
    free(*gp);
  }
}

/*....................................................................*/
void closeFile(lime_fptr *fptr, const int fileFormatI){

  if(fileFormatI==lime_FITS){
    closeFITSFile(fptr);
  }
}

/*....................................................................*/
void constructLinkArrays(const unsigned int numGridPoints, struct grid *g\
  , struct linkType **links, unsigned int *totalNumLinks, struct linkType ***nnLinks\
  , unsigned int **firstNearNeigh, unsigned int *totalNumNeigh, const int dataStageI){
  /*
See the comment at the beginning of the present module for a description of how the pointers 'links', 'nnLinks' and 'firstNearNeigh' relate to the grid struct.
  */

  unsigned int li, ni, ci, idA, idB, trialIdA, linkId;
  int i, jA, jB, k, nearI;
  _Bool *pointIsDone=NULL, linkNotFound;
  struct grid *gAPtr=NULL, *gBPtr=NULL;
  int *nnLinkIs=NULL;
  char message[80];
  double evenCoeffSign;

  pointIsDone     = malloc(sizeof(*pointIsDone)    *numGridPoints);
  *firstNearNeigh = malloc(sizeof(**firstNearNeigh)*numGridPoints);

  /* First add up all the numNeigh.
  */
  ni = 0;
  for(i=0;i<numGridPoints;i++){
    ni += g[i].numNeigh;
    pointIsDone[g[i].id] = 0;
  }
  *totalNumNeigh = ni;

  *links   = malloc(sizeof(**links)*(*totalNumNeigh)); /* Bigger than needed, we'll reduce it later. */
  nnLinkIs = malloc(sizeof(*nnLinkIs)*(*totalNumNeigh));

  li = 0;
  ni = 0;
  for(i=0;i<numGridPoints;i++){
    gAPtr = g+i;
    idA = g[i].id;
    (*firstNearNeigh)[idA] = ni;

    for(jA=0;jA<gAPtr->numNeigh;jA++){
      /* Check to see if the NN has been done.
      */
      gBPtr = gAPtr->neigh[jA];
      idB = gBPtr->id;
      if(pointIsDone[idB]){
        /* The link exists; we need to find a pointer to it for the next entry of nnLinks. To do that, we need to go through the list of NN links of the point whose id is idB and identify that link which connects back to idA.
        */
        linkNotFound = 1; /* default */
        for(jB=0;jB<gBPtr->numNeigh;jB++){
          trialIdA = gBPtr->neigh[jB]->id;
          if(trialIdA==idA){
            linkNotFound = 0;
            k = (*firstNearNeigh)[idB] + jB;
            linkId = nnLinkIs[k];
          }
        }
        if(linkNotFound){
          if(!silent) bail_out("Link not found.");
          exit(1);
        }

      }else{
        /* The link does not yet exist; we must create it, load it into links, and find a pointer to it for the next entry of nnLinks.
        */
        linkId = li;
        (*links)[li].id = linkId;
        (*links)[li].g[0] = gAPtr;
        (*links)[li].g[1] = gBPtr;

        li++;
      }

      nnLinkIs[ni] = linkId;
      ni++;
    }

    pointIsDone[idA] = 1;
  }
  *totalNumLinks = li;

  *links = realloc(*links, sizeof(struct linkType)*(*totalNumLinks));

  if(dataStageI>2){
    for(li=0;li<*totalNumLinks;li++)
      (*links)[li].aCoeffs = malloc(sizeof(double)*NUM_VEL_COEFFS);

    for(li=0;li<*totalNumLinks;li++){
      if((*links)[li].g[0]->sink && (*links)[li].g[1]->sink){
        for(ci=0;ci<NUM_VEL_COEFFS;ci++)
          (*links)[li].aCoeffs[ci] = 0.0;
      }else{

        /* If g[0] is a sink point then try g[1], because sink point a*'s are set to zero. Remember to invert the sign of the even coefficients if we are reading them from g[1]. (See the header comments to loadNnIntoGrid() for an explanation of the reason for the sign inversion.)
        */
        if((*links)[li].g[0]->sink){ /* If we get to here, this ensures that g[1] is not a sink point. */
          nearI = 1;
          evenCoeffSign = -1.0;
        }else{
          nearI = 0;
          evenCoeffSign = 1.0;
        }

        gAPtr = (*links)[li].g[nearI];
        /* Find which neighbour of gAPtr corresponds to the link: */
        linkNotFound = 1;
        for(jA=0;jA<gAPtr->numNeigh;jA++){
          idB = gAPtr->neigh[jA]->id;
          if(linkNotFound && idB==(*links)[li].g[1-nearI]->id){
            (*links)[li].aCoeffs[0] = evenCoeffSign*gAPtr->a0[jA];
            (*links)[li].aCoeffs[1] =               gAPtr->a1[jA];
            (*links)[li].aCoeffs[2] = evenCoeffSign*gAPtr->a2[jA];
            (*links)[li].aCoeffs[3] =               gAPtr->a3[jA];
            (*links)[li].aCoeffs[4] = evenCoeffSign*gAPtr->a4[jA];
            linkNotFound = 0;
          }
        }
        if(linkNotFound){
          sprintf(message, "SNAFU with link %d grid point %d", li, gAPtr->id);
          if(!silent) bail_out(message);
          exit(1);
        }
      } /* end if not both sink */
    } /* end loop over links */
  }else{ /* dataStageI<=2 */
    for(li=0;li<*totalNumLinks;li++)
      (*links)[li].aCoeffs = NULL;
  }

  *nnLinks = malloc(sizeof(**nnLinks)*(*totalNumNeigh));
  for(ni=0;ni<(*totalNumNeigh);ni++)
    (*nnLinks)[ni] = &(*links)[nnLinkIs[ni]];

  free(nnLinkIs);
  free(pointIsDone);
  /* The calling routine must free firstNearNeigh, nnLinks, links. */
}


/*....................................................................*/
int writeGrid(char *outFileName, const int fileFormatI, inputPars par\
  , unsigned short numDims, unsigned short numACoeffs, struct grid *g\
  , molData *md, char **collPartNames, const int dataStageI){

  lime_fptr *fptr=NULL;
  int status = 0;
  unsigned short speciesI;
  char message[80];
  struct linkType *links=NULL, **nnLinks=NULL;
  unsigned int totalNumLinks, totalNumNeigh, *firstNearNeigh=NULL;

  if (outFileName==""){
    if(!silent) warning("Cannot write grid list to file, filename is blank.");
    return 1;
  }

  if (g==NULL || dataStageI<1){
    if(!silent) warning("Cannot write grid list to file, there are no entries in it.");
    return 2;
  }

  sprintf(message, "Writing grid-point list to file %s", outFileName);
  if(!silent) printMessage(message);

  if (dataStageI>1){
    constructLinkArrays((unsigned int)par.ncell, g, &links, &totalNumLinks\
      , &nnLinks, &firstNearNeigh, &totalNumNeigh, dataStageI);
  }

  fptr = openFileForWrite(outFileName, fileFormatI, dataStageI);
  if(fptr==NULL){
    closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, totalNumLinks);
    return 3;
  }

  status = writeGridBlock(fptr, fileFormatI, par, numDims, g, firstNearNeigh, collPartNames, dataStageI);
  if(status){
    closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, totalNumLinks);
    return 4;
  }

  if (dataStageI>1){
    status = writeNnIndicesBlock(fptr, fileFormatI, totalNumNeigh, nnLinks, links);
    if(status){
      closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, totalNumLinks);
      return 5;
    }

    status = writeLinksBlock(fptr, fileFormatI, totalNumLinks, numACoeffs, links, dataStageI);
    if(status){
      closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, totalNumLinks);
      return 6;
    }

    if (dataStageI>3){
      for(speciesI=0;speciesI<par.nSpecies;speciesI++){
        status = writePopsBlock(fptr, fileFormatI, (unsigned int)par.ncell, md, speciesI, g);
        if(status){
          closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, totalNumLinks);
          return 7;
        }
      }
    }
  }

  closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, totalNumLinks);
  return 0;
}

/*....................................................................*/
int writeGridBlock(lime_fptr *fptr, const int fileFormatI, inputPars par\
  , unsigned short numDims, struct grid *g, unsigned int *firstNearNeigh\
  , char **collPartNames, const int dataStageI){

  int status=0;

  if(fileFormatI==lime_FITS){
    writeGridExtToFits(fptr, par, numDims, g, firstNearNeigh, collPartNames, dataStageI);
  }else{
    status = 1;
  }

  return status;
}

/*....................................................................*/
int writeNnIndicesBlock(lime_fptr *fptr, const int fileFormatI\
  , const unsigned int totalNumNeigh, struct linkType **nnLinks\
  , struct linkType *links){

  int status=0;

  if(fileFormatI==lime_FITS){
    writeNnIndicesExtToFits(fptr, totalNumNeigh, nnLinks, links);
  }else{
    status = 1;
  }

  return status;
}

/*....................................................................*/
int writeLinksBlock(lime_fptr *fptr, const int fileFormatI\
  , const unsigned int totalNumLinks, const unsigned short numACoeffs\
  , struct linkType *links, const int dataStageI){

  int status=0;

  if(fileFormatI==lime_FITS){
    writeLinksExtToFits(fptr, totalNumLinks, numACoeffs, links, dataStageI);
  }else{
    status = 1;
  }

  return status;
}

/*....................................................................*/
int writePopsBlock(lime_fptr *fptr, const int fileFormatI\
  , unsigned int numGridPoints, molData *md, unsigned short speciesI\
  , struct grid *g){

  int status=0;

  if(fileFormatI==lime_FITS){
    writePopsExtToFits(fptr, numGridPoints, md, speciesI, g);
  }else{
    status = 1;
  }

  return status;
}

/*....................................................................*/
int readGrid(char *inFileName, const int fileFormatI\
  , struct gridInfoType *gridInfoRead, struct grid **gp\
  , char ***collPartNames, int *numCollPartRead, int *dataStageI){
//********** free grid before call, and set gridInfoRead elements to defaults.

  lime_fptr *fptr;
  int status=0;
  unsigned short i_s, numBlocks;
  unsigned int *firstNearNeigh=NULL, totalNumGridPoints, i_u;
  struct linkType *links=NULL, **nnLinks=NULL;
  char message[80];

  sprintf(message, "Reading grid-point list from file %s", inFileName);
  if(!silent) printMessage(message);

  /* Open the file and also return the data stage. */
  fptr = openFileForRead(inFileName, fileFormatI, dataStageI);

  if(*dataStageI<1 || *dataStageI>NUM_GRID_STAGES){
    if(!silent) bail_out("Data stage out of range.");
    exit(1);
  }

  sprintf(message, "Grid-point data stage determined as %d", *dataStageI);
  if(!silent) printMessage(message);

  /* Read the values which should be in grid for every stage.
  */
  status = readGridBlock(fptr, fileFormatI, gridInfoRead, gp, &firstNearNeigh\
    , collPartNames, numCollPartRead, *dataStageI);

  if(status || gridInfoRead->nSinkPoints<=0\
  || gridInfoRead->nInternalPoints<=0 || gridInfoRead->nDims<=0){
    closeAndFree(fptr, fileFormatI, firstNearNeigh, NULL, NULL, 0);
    freeReadGrid(gp, *gridInfoRead);

    if(status)
      return 1;
    else if(gridInfoRead->nSinkPoints<=0)
      return 2;
    else if(gridInfoRead->nInternalPoints<=0)
      return 3;
    else if(gridInfoRead->nDims<=0)
      return 4;
    else{
      sprintf(message, "This indicates a programming error. Please contact the developer.");
      if(!silent) bail_out(message);
      exit(1);
    }
  }

  if((*dataStageI)>2 && (gridInfoRead->nDensities<=0 || gridInfoRead->nSpecies<=0)){
    closeAndFree(fptr, fileFormatI, firstNearNeigh, NULL, NULL, 0);
    freeReadGrid(gp, *gridInfoRead);

    if(gridInfoRead->nDensities<=0)
      return 5;
    else if(gridInfoRead->nSpecies<=0) //***** what if all continuum images??
      return 6;
    else{
      sprintf(message, "This indicates a programming error. Please contact the developer.");
      if(!silent) bail_out(message);
      exit(1);
    }
  }

  totalNumGridPoints = gridInfoRead->nInternalPoints+gridInfoRead->nSinkPoints;

  if((*dataStageI)>1){ /* there should be blocks for the nnIndices and links. */
    status = readLinksBlock(fptr, fileFormatI, gridInfoRead, *gp, &links, *dataStageI);
    if(status){
      closeAndFree(fptr, fileFormatI, firstNearNeigh, NULL, links, gridInfoRead->nLinks);
      freeReadGrid(gp, *gridInfoRead);
      return 7;
    }

    if ((*dataStageI)>2 && gridInfoRead->nACoeffs<=0){
      closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, gridInfoRead->nLinks);
      freeReadGrid(gp, *gridInfoRead);
      return 8;
    }

    status = readNnIndicesBlock(fptr, fileFormatI, links, &nnLinks, gridInfoRead);
    if(status){
      closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, gridInfoRead->nLinks);
      freeReadGrid(gp, *gridInfoRead);
      return 9;
    }

    /* Convert the NN information back to the standard LIME grid struct format.
    */
    loadNnIntoGrid(firstNearNeigh, nnLinks, links, *gridInfoRead, *gp, *dataStageI);

    if ((*dataStageI)>3){ /* there should be pops blocks. */
      status = getNumPopsBlocks(fptr, fileFormatI, &numBlocks);

      if(status || numBlocks<=0 || numBlocks!=gridInfoRead->nSpecies){
        closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, gridInfoRead->nLinks);
        freeReadGrid(gp, *gridInfoRead);
        if(status)
          return 10;
        else if(numBlocks<=0)
          return 11;
        else if(numBlocks!=gridInfoRead->nSpecies)
          return 12;
        else{
          sprintf(message, "This indicates a programming error. Please contact the developer.");
          if(!silent) bail_out(message);
          exit(1);
        }
      }

      gridInfoRead->mols = malloc(sizeof(struct molInfoType)*gridInfoRead->nSpecies);

      for(i_u=0;i_u<totalNumGridPoints;i_u++){
        (*gp)[i_u].mol = malloc(sizeof(struct populations)*gridInfoRead->nSpecies);
        for(i_s=0;i_s<gridInfoRead->nSpecies;i_s++){
          (*gp)[i_u].mol[i_s].dust    = NULL;
          (*gp)[i_u].mol[i_s].knu     = NULL;
          (*gp)[i_u].mol[i_s].pops    = NULL;
          (*gp)[i_u].mol[i_s].partner = NULL;
        }
      }

      for(i_s=0;i_s<gridInfoRead->nSpecies;i_s++){
        gridInfoRead->mols[i_s].molName = NULL;
        gridInfoRead->mols[i_s].nLevels = -1;
        gridInfoRead->mols[i_s].nLines = -1;

        status = readPopsBlock(fptr, fileFormatI, i_s, *gp, gridInfoRead);

        if(status || gridInfoRead->mols[i_s].nLevels<=0){
          closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, gridInfoRead->nLinks);
          freeReadGrid(gp, *gridInfoRead);
          if(status)
            return 13;
          else if(gridInfoRead->mols[i_s].nLevels<=0)
            return 14;
          else{
            sprintf(message, "This indicates a programming error. Please contact the developer.");
            if(!silent) bail_out(message);
            exit(1);
          }
        }
      }
    }
  }

  /* Set unread pointers to NULL:
  */
  if((*dataStageI)<4){
    for(i_u=0;i_u<totalNumGridPoints;i_u++)
      (*gp)[i_u].mol = NULL;

    if((*dataStageI)<3){
      for(i_u=0;i_u<totalNumGridPoints;i_u++){
        (*gp)[i_u].a0   = NULL;
        (*gp)[i_u].a1   = NULL;
        (*gp)[i_u].a2   = NULL;
        (*gp)[i_u].a3   = NULL;
        (*gp)[i_u].a4   = NULL;
        (*gp)[i_u].dens = NULL;
        (*gp)[i_u].abun = NULL;
      }

      if((*dataStageI)<2){
        for(i_u=0;i_u<totalNumGridPoints;i_u++)
          (*gp)[i_u].neigh = NULL;
      }
    }
  }

  /* The following pointers are for secondary information which is either calculated later from the read values, or used later to store temporary values:
  */
  for(i_u=0;i_u<totalNumGridPoints;i_u++){
    (*gp)[i_u].dir  = NULL;
    (*gp)[i_u].w    = NULL;
    (*gp)[i_u].ds   = NULL;
  }

  closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, gridInfoRead->nLinks);
  return 0;
}

/*....................................................................*/
int readGridBlock(lime_fptr *fptr, const int fileFormatI\
  , struct gridInfoType *gridInfoRead, struct grid **gp\
  , unsigned int **firstNearNeigh, char ***collPartNames\
  , int *numCollPartRead, const int dataStageI){

  int status=0;

  if(fileFormatI==lime_FITS){
    readGridExtFromFits(fptr, gridInfoRead, gp, firstNearNeigh\
      , collPartNames, numCollPartRead, dataStageI);
  }else{
    status = 1;
  }

  return status;
}

/*....................................................................*/
int readLinksBlock(lime_fptr *fptr, const int fileFormatI\
  , struct gridInfoType *gridInfoRead, struct grid *g\
  , struct linkType **links, const int dataStageI){

  int status=0;

  if(fileFormatI==lime_FITS){
    readLinksExtFromFits(fptr, gridInfoRead, g, links, dataStageI);
  }else{
    status = 1;
  }

  return status;
}

/*....................................................................*/
int readNnIndicesBlock(lime_fptr *fptr, const int fileFormatI, struct linkType *links\
  , struct linkType ***nnLinks, struct gridInfoType *gridInfoRead){

  int status=0;

  if(fileFormatI==lime_FITS){
    readNnIndicesExtFromFits(fptr, links, nnLinks, gridInfoRead);
  }else{
    status = 1;
  }

  return status;
}

/*....................................................................*/
void loadNnIntoGrid(unsigned int *firstNearNeigh, struct linkType **nnLinks\
  , struct linkType *links, struct gridInfoType gridInfoRead, struct grid *gp\
  , const int dataStageI){
  /*
See the comment at the beginning of the present module for a description of how the pointers 'links', 'nnLinks' and 'firstNearNeigh' relate to the grid struct.

The function mallocs the following extensions of struct g for each grid point:
	neigh
if dataStageI>2:
	a0
	a1
	a2
	a3
	a4

The business here with the evenCoeffSign is due to the fact that, for each grid point A, the 'a' coefficients are calculated for a polynomial with independent coordinate being the distance from that grid point towards the appropriate neighbour B. Point B will also store its own 'a' coefficients, which produce identically the same polynomial, but which in general will have different values, because in this case the coefficients are calculated appropriate to the distance from B towards A. Now that the definition of the independent coordinate has been defined as the fractional distance minus 0.5, thus symmetrical about the midpoint between A and B, the coefficients for B can be obtained from those of A by inverting the sign of the even-power coefficients. Note that we invert the even coefficients, not the odd, because not only has the sign of the independent coordinate (i.e the fractional distance along the link minus 0.5) reversed, but also the that of the function being interpolated (the component of velocity in the direction of the neighbour point).

Further: the 'a' coefficients for the link or edge between A and B are stored in a struct of type linkType. They are stored appropriate to point g[0] of that struct. We cycle through all the grid points in the present function and, for each point/neighbour pair, we look first to see which end of the matching link element 'point' is - i.e., whether it is g[0] or g[1] of the link. If it is g[0], then the 'a' coefficients stored in link are appropriate for it; if g[1], then the even-power ones have to be inverted in sign.
  */

  const unsigned int totalNumGridPoints = gridInfoRead.nInternalPoints+gridInfoRead.nSinkPoints;
  unsigned int i_u;
  int j;
  struct linkType *linkPtr;
  char message[80];
  double evenCoeffSign;

  for(i_u=0;i_u<totalNumGridPoints;i_u++){
    gp[i_u].neigh = malloc(sizeof(struct grid *)*gp[i_u].numNeigh);

    for(j=0;j<gp[i_u].numNeigh;j++){
      linkPtr = nnLinks[firstNearNeigh[i_u]+j];
      if(linkPtr->g[0]->id==(int)i_u){
        gp[i_u].neigh[j] = linkPtr->g[1];
      }else{
        gp[i_u].neigh[j] = linkPtr->g[0];
      }
    }
  }

  if(dataStageI>2){
    /* Just for the time being: */
    if(gridInfoRead.nACoeffs!=5){
      if(!silent){
        sprintf(message, "There should be %d ACOEFF_n columns, but %d were read.", NUM_VEL_COEFFS, (int)gridInfoRead.nACoeffs);
        bail_out(message);
      }
      exit(1);
    }

    for(i_u=0;i_u<totalNumGridPoints;i_u++){
      gp[i_u].a0 = malloc(sizeof(double)*gp[i_u].numNeigh);
      gp[i_u].a1 = malloc(sizeof(double)*gp[i_u].numNeigh);
      gp[i_u].a2 = malloc(sizeof(double)*gp[i_u].numNeigh);
      gp[i_u].a3 = malloc(sizeof(double)*gp[i_u].numNeigh);
      gp[i_u].a4 = malloc(sizeof(double)*gp[i_u].numNeigh);
    }
    for(i_u=gridInfoRead.nInternalPoints;i_u<gridInfoRead.nSinkPoints;i_u++){
      for(j=0;j<gp[i_u].numNeigh;j++){
        gp[i_u].a0[j] = 0.0;
        gp[i_u].a1[j] = 0.0;
        gp[i_u].a2[j] = 0.0;
        gp[i_u].a3[j] = 0.0;
        gp[i_u].a4[j] = 0.0;
      }
    }
    for(i_u=0;i_u<gridInfoRead.nInternalPoints;i_u++){
      for(j=0;j<gp[i_u].numNeigh;j++){
        linkPtr = nnLinks[firstNearNeigh[i_u]+j];
        if(linkPtr->g[0]->id==(int)i_u){
          evenCoeffSign = 1.0;
        }else{
          evenCoeffSign = -1.0;
        }

        gp[i_u].a0[j] = evenCoeffSign*linkPtr->aCoeffs[0];
        gp[i_u].a1[j] =               linkPtr->aCoeffs[1];
        gp[i_u].a2[j] = evenCoeffSign*linkPtr->aCoeffs[2];
        gp[i_u].a3[j] =               linkPtr->aCoeffs[3];
        gp[i_u].a4[j] = evenCoeffSign*linkPtr->aCoeffs[4];
      }
    }
  }
}

/*....................................................................*/
int getNumPopsBlocks(lime_fptr *fptr, const int fileFormatI, unsigned short *numBlocks){
  int status = 0;
  _Bool blockFound = 1;

  *numBlocks = 0;
  while(blockFound && !status){
    status = checkPopsBlockExists(fptr, fileFormatI, *numBlocks, &blockFound);
    (*numBlocks)++;
  }

  (*numBlocks)--;

  return status;
}

/*....................................................................*/
int checkPopsBlockExists(lime_fptr *fptr, const int fileFormatI\
  , const unsigned short speciesI, _Bool *blockFound){
  int status=0;

  *blockFound = 0;

  if(fileFormatI==lime_FITS){
    *blockFound = checkPopsFitsExtExists(fptr, speciesI);
  }else{
    status = 1;
  }

  return status;
}

/*....................................................................*/
int readPopsBlock(lime_fptr *fptr, const int fileFormatI\
  , const unsigned short speciesI, struct grid *g\
  , struct gridInfoType *gridInfoRead){

  int status=0;

  if(fileFormatI==lime_FITS){
    readPopsExtFromFits(fptr, speciesI, g, gridInfoRead);
  }else{
    status = 1;
  }

  return status;
}


