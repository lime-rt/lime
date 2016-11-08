/*
 *  gridio.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2016 The LIME development team
 *
TODO:
 */

#include "lime.h"
#include "gridio.h"

/*
This module contains generic routines for writing grid data to, and reading it from, a file on disk. The problem with doing this is that the grid struct contains different amounts of information at different times in the running of the code. In order to quantify and regulate this, a dataFlags integer is used to record the presence or absence (as indicated by the value of the appropriate bit in the mask) of particular types of information. The bits associated with certain fields of struct grid are given in the lime.h header.

Notes:
  - If several fields of the struct are listed for a given bit, a value of 1 for that bit indicates that all fields are expected to be present, and a value of 0 indicates that all fields of that bit will be ignored.
  - Only a few combinations are not allowed, as follows:
    * No other bit may be set if DS_bit_x is not.
    * DS_bit_ACOEFF may not be set if either DS_bit_neighbours or DS_bit_velocity is not.
    * DS_bit_populations may not be set unless all the others are set as well (except DS_bit_magfield may be set or not).

Writing the grid file:
----------------------
To make things simpler, four stages have been defined at which the user may write the grid data to file. These are described in the following table:

	Data mask bits set:	    stage 0	    stage 1	    stage 2	    stage 3	    stage 4
	.....................................................................................................
	DS_bit_x             		0		1		1		1		1
	DS_bit_neighbours    		0		0		1		1		1
	DS_bit_velocity      		0		0		0		1		1
	DS_bit_density       		0		0		0		1		1
	DS_bit_abundance     		0		0		0		1		1
	DS_bit_turb_doppler  		0		0		0		1		1
	DS_bit_temperatures  		0		0		0		1		1
	DS_bit_magfield  		0		0		0		x		x
	DS_bit_ACOEFF        		0		0		0		1		1
	DS_bit_populations   		0		0		0		0		1
	.....................................................................................................

Notes:
  - dataStageI==0 has been included for completeness/robustness but the user may not write a file with nothing in it.
  - LIME may run to completion without ever reaching stage 4 - if all the images required were continuum ones, for example.

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
lime_fptr *
openFileForWrite(char *outFileName, const int fileFormatI){
  lime_fptr *fptr=NULL;

  if(fileFormatI==lime_FITS){
    fptr = openFITSFileForWrite(outFileName);
  }else{
    return NULL;
  }

  return fptr;
}

/*....................................................................*/
void
constructLinkArrays(struct gridInfoType gridInfo, struct grid *gp, struct linkType **links\
  , unsigned int *totalNumLinks, struct linkType ***nnLinks\
  , unsigned int **firstNearNeigh, unsigned int *totalNumNeigh, const int dataFlags){
  /*
See the comment at the beginning of the present module for a description of how the pointers 'links', 'nnLinks' and 'firstNearNeigh' relate to the grid struct.
  */

  const unsigned int totalNumGridPoints = gridInfo.nInternalPoints+gridInfo.nSinkPoints;
  unsigned int li,ni,idA,idB,trialIdA,linkId,i_ui,gi0,gi1;
  unsigned short i_us,j_us;
  int jA,jB,k,nearI;
  _Bool *pointIsDone=NULL, linkNotFound;
  struct grid *gAPtr=NULL, *gBPtr=NULL;
  int *nnLinkIs=NULL;
  char message[80];

  if(!allBitsSet(dataFlags, DS_mask_neighbours)) /* Grid fields id, sink, numNeigh and neigh must be available. */
    return;

  pointIsDone     = malloc(sizeof(*pointIsDone)    *totalNumGridPoints);
  *firstNearNeigh = malloc(sizeof(**firstNearNeigh)*totalNumGridPoints);

  /* First add up all the numNeigh.
  */
  ni = 0;
  for(i_ui=0;i_ui<totalNumGridPoints;i_ui++){
    ni += gp[i_ui].numNeigh;
    pointIsDone[gp[i_ui].id] = 0;
  }
  *totalNumNeigh = ni;

  *links   = malloc(sizeof(**links)*(*totalNumNeigh)); /* Bigger than needed, we'll reduce it later. */
  nnLinkIs = malloc(sizeof(*nnLinkIs)*(*totalNumNeigh));

  li = 0;
  ni = 0;
  for(i_ui=0;i_ui<totalNumGridPoints;i_ui++){
    gAPtr = gp+i_ui;
    idA = gAPtr->id;
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
        (*links)[li].gis[0] = idA;
        (*links)[li].gis[1] = idB;

        li++;
      }

      nnLinkIs[ni] = linkId;
      ni++;
    }

    pointIsDone[idA] = 1;
  }
  *totalNumLinks = li;

  *links = realloc(*links, sizeof(struct linkType)*(*totalNumLinks));

  if(allBitsSet(dataFlags, DS_mask_ACOEFF)){
    for(i_ui=0;i_ui<*totalNumLinks;i_ui++)
      (*links)[i_ui].vels = malloc(sizeof(double)*gridInfo.nLinkVels*gridInfo.nDims);

    for(i_ui=0;i_ui<*totalNumLinks;i_ui++){
      gi0 = (*links)[i_ui].gis[0];
      gi1 = (*links)[i_ui].gis[1];
      if(gp[gi0].sink && gp[gi1].sink){
        for(i_us=0;i_us<gridInfo.nLinkVels;i_us++){
          for(j_us=0;j_us<gridInfo.nDims;j_us++)
            (*links)[i_ui].vels[gridInfo.nDims*i_us + j_us] = 0.0;
        }
      }else{ /* One of the two vertices must not be a sink point. */
        /*
We want to read in the intra-edge velocity samples which are currently still being stored in each grid point. This is x2 redundant storage, because if grid points A and B are at the ends of a given edge, gp[A].v1 == gp[B].v3 and so forth. We want to define the order of velocities stored in the .vels list of each links element as that of the grid point that comes first in the .gis list, i.e. .gis[0]. The only time this does not work is if .gis[0] is a sink point, because we are not storing intra-edge velocity samples for sink points. In this case we get the .vels values from .gis[1] and reverse their order.
        */

        if(gp[gi0].sink){
          /* If we get to here, we can be certain that .gis[1] is not a sink point. */
          nearI = 1;
        }else{
          nearI = 0;
        }

        gAPtr = &gp[(*links)[i_ui].gis[nearI]];
        /* Find which neighbour of gAPtr corresponds to the link: */
        linkNotFound = 1;
        for(jA=0;jA<gAPtr->numNeigh;jA++){
          idB = gAPtr->neigh[jA]->id;
          if(linkNotFound && idB==(*links)[i_ui].gis[1-nearI]){
            if(nearI==0){
              for(j_us=0;j_us<gridInfo.nDims;j_us++){
                (*links)[i_ui].vels[gridInfo.nDims*0 + j_us] = gAPtr->v1[gridInfo.nDims*jA + j_us];
                (*links)[i_ui].vels[gridInfo.nDims*1 + j_us] = gAPtr->v2[gridInfo.nDims*jA + j_us];
                (*links)[i_ui].vels[gridInfo.nDims*2 + j_us] = gAPtr->v3[gridInfo.nDims*jA + j_us];
              }
            }else{
              for(j_us=0;j_us<gridInfo.nDims;j_us++){
                (*links)[i_ui].vels[gridInfo.nDims*0 + j_us] = gAPtr->v3[gridInfo.nDims*jA + j_us];
                (*links)[i_ui].vels[gridInfo.nDims*1 + j_us] = gAPtr->v2[gridInfo.nDims*jA + j_us];
                (*links)[i_ui].vels[gridInfo.nDims*2 + j_us] = gAPtr->v1[gridInfo.nDims*jA + j_us];
              }
            }
            linkNotFound = 0;
          }
        }
        if(linkNotFound){
          sprintf(message, "SNAFU with link %d grid point %d", (int)i_ui, gAPtr->id);
          if(!silent) bail_out(message);
          exit(1);
        }
      } /* end if not both sink */
    } /* end loop over links */
  }else{ /* a coeffs not available */
    for(i_ui=0;i_ui<*totalNumLinks;i_ui++)
      (*links)[i_ui].vels = NULL;
  }

  *nnLinks = malloc(sizeof(**nnLinks)*(*totalNumNeigh));
  for(i_ui=0;i_ui<(*totalNumNeigh);i_ui++)
    (*nnLinks)[i_ui] = &(*links)[nnLinkIs[i_ui]];

  free(nnLinkIs);
  free(pointIsDone);
  /* The calling routine must free firstNearNeigh, nnLinks, links. */
}

/*....................................................................*/
void
closeFile(lime_fptr *fptr, const int fileFormatI){
  if(fileFormatI==lime_FITS)
    closeFITSFile(fptr);
  //**** error if not?
}

/*....................................................................*/
void
closeAndFree(lime_fptr *fptr, const int fileFormatI\
  , unsigned int *firstNearNeigh, struct linkType **nnLinks\
  , struct linkType *links, const unsigned int totalNumLinks){

  unsigned int i_ui;

  closeFile(fptr, fileFormatI);

  free(firstNearNeigh);
  free(nnLinks);
  if(links!=NULL){
    for(i_ui=0;i_ui<totalNumLinks;i_ui++)
      free(links[i_ui].vels);
   free(links);
  }
}

/*....................................................................*/
int
writeKeywords(lime_fptr *fptr, const int fileFormatI\
  , struct keywordType *kwds, const int numKeywords){

  int status=0;

  if(fileFormatI==lime_FITS){
    writeKeywordsToFits(fptr, kwds, numKeywords);
  }else{
    status = 1;
  }

  return status;
}

/*....................................................................*/
int
writeGridTable(lime_fptr *fptr, const int fileFormatI\
  , struct gridInfoType gridInfo, struct grid *gp, unsigned int *firstNearNeigh\
  , char **collPartNames, const int dataFlags){

  int status=0;

  if(fileFormatI==lime_FITS){
    writeGridExtToFits(fptr, gridInfo, gp, firstNearNeigh, collPartNames, dataFlags);
  }else{
    status = 1;
  }

  return status;
}

/*....................................................................*/
int
writeNnIndicesTable(lime_fptr *fptr, const int fileFormatI\
  , struct gridInfoType gridInfo, struct linkType **nnLinks){\

  int status=0;

  if(fileFormatI==lime_FITS){
    writeNnIndicesExtToFits(fptr, gridInfo, nnLinks);
  }else{
    status = 1;
  }

  return status;
}

/*....................................................................*/
int
writeLinksTable(lime_fptr *fptr, const int fileFormatI\
  , struct gridInfoType gridInfo, struct linkType *links){

  int status=0;

  if(fileFormatI==lime_FITS){
    writeLinksExtToFits(fptr, gridInfo, links);
  }else{
    status = 1;
  }

  return status;
}

/*....................................................................*/
int
writePopsTable(lime_fptr *fptr, const int fileFormatI\
  , struct gridInfoType gridInfo, unsigned short speciesI\
  , struct grid *gp){


  int status=0;

  if(fileFormatI==lime_FITS){
    writePopsExtToFits(fptr, gridInfo, speciesI, gp);
  }else{
    status = 1;
  }

  return status;
}

/*....................................................................*/
int
writeGrid(char *outFileName, const int fileFormatI, struct gridInfoType gridInfo\
  , struct keywordType *primaryKwds, const int numKeywords, struct grid *gp, char **collPartNames, const int dataFlags){

  /*
This is designed to be a generic function to write the grid data (in any of its accepted degrees of completeness) to file. It is assumed that several tables of different size will need to be written, corresponding to the different dimensionalities of the elements of the 'grid' struct. These are described in the following list.

	Table:		Contains:
	---------------------------------------------------------------
	grid		All the scalar elements of struct grid.

	nnIndices	Basically stores the information in the element 'neigh'.

	links		Stores the elements a0, a1 etc.

	pops		Stores the element 'pops' of each element 'mol'. There will be a separate table for each radiating species.
	---------------------------------------------------------------
  */

  lime_fptr *fptr=NULL;
  int status = 0;
  unsigned short i_us;
  char message[80];
  struct linkType *links=NULL, **nnLinks=NULL;
  unsigned int *firstNearNeigh=NULL;

  if (outFileName==""){
    if(!silent) warning("Cannot write grid list to file, filename is blank.");
    return 1;
  }

  if (gp==NULL || !allBitsSet(dataFlags, DS_mask_x)){
    if(!silent) warning("Cannot write grid list to file, there are no entries in it.");
    return 2;
  }

  sprintf(message, "Writing grid-point list to file %s", outFileName);
  if(!silent) printMessage(message);

  constructLinkArrays(gridInfo, gp, &links, &gridInfo.nLinks\
    , &nnLinks, &firstNearNeigh, &gridInfo.nNNIndices, dataFlags);

  fptr = openFileForWrite(outFileName, fileFormatI);
  if(fptr==NULL){
    closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, gridInfo.nLinks);
    return 3;
  }

  /* Write header keywords:
  */
  status = writeKeywords(fptr, fileFormatI, primaryKwds, numKeywords);
  if(status){
    closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, gridInfo.nLinks);
    return 4;
  }

  status = writeGridTable(fptr, fileFormatI, gridInfo, gp, firstNearNeigh, collPartNames, dataFlags);
  if(status){
    closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, gridInfo.nLinks);
    return 5;
  }

  if (links!=NULL && nnLinks!=NULL){
    status = writeNnIndicesTable(fptr, fileFormatI, gridInfo, nnLinks);
    if(status){
      closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, gridInfo.nLinks);
      return 6;
    }

    status = writeLinksTable(fptr, fileFormatI, gridInfo, links);
    if(status){
      closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, gridInfo.nLinks);
      return 7;
    }
  }

  if (allBitsSet(dataFlags, DS_mask_populations)){
    for(i_us=0;i_us<gridInfo.nSpecies;i_us++){
      status = writePopsTable(fptr, fileFormatI, gridInfo, i_us, gp);
      if(status){
        closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, gridInfo.nLinks);
        return 8;
      }
    }
  }

  closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, gridInfo.nLinks);
  return 0;
}


/*....................................................................*/
lime_fptr *
openFileForRead(char *inFileName, const int fileFormatI){
  lime_fptr *fptr=NULL;

  if(fileFormatI==lime_FITS){
    fptr = openFITSFileForRead(inFileName);
  }else{
    return NULL;
  }

  return fptr;
}

/*....................................................................*/
int
readKeywords(lime_fptr *fptr, const int fileFormatI\
  , struct keywordType *kwds, const int numKeywords){

  int status=0;

  if(fileFormatI==lime_FITS){
    readKeywordsFromFits(fptr, kwds, numKeywords);
  }else{
    status = 1;
  }

  return status;
}

/*....................................................................*/
int
readGridTable(lime_fptr *fptr, const int fileFormatI\
  , struct gridInfoType *gridInfoRead, struct grid **gp\
  , unsigned int **firstNearNeigh, char ***collPartNames\
  , int *numCollPartRead, int *dataFlags){
  /*
Individual routines called should set the appropriate bits of dataFlags; also malloc gp and set all its defaults. (Note there is a bespoke routine grid.c:mallocAndSetDefaultGrid() to do the latter.)
  */

  int status=0;

  if(fileFormatI==lime_FITS){
    readGridExtFromFits(fptr, gridInfoRead, gp, firstNearNeigh\
      , collPartNames, numCollPartRead, dataFlags);
  }else{
    status = 1;
  }

  return status;
}

/*....................................................................*/
int
readLinksTable(lime_fptr *fptr, const int fileFormatI\
  , struct gridInfoType *gridInfoRead, struct grid *gp\
  , struct linkType **links, int *dataFlags){
  /*
Individual routines called should set the appropriate bits of dataFlags.
  */

  int status=0;

  if(fileFormatI==lime_FITS){
    readLinksExtFromFits(fptr, gridInfoRead, gp, links, dataFlags);
  }else{
    status = 1;
  }

  return status;
}

/*....................................................................*/
int
readNnIndicesTable(lime_fptr *fptr, const int fileFormatI, struct linkType *links\
  , struct linkType ***nnLinks, struct gridInfoType *gridInfoRead, int *dataFlags){
  /*
Individual routines called should set the appropriate bits of dataFlags.
  */

  int status=0;

  if(fileFormatI==lime_FITS){
    readNnIndicesExtFromFits(fptr, links, nnLinks, gridInfoRead, dataFlags);
  }else{
    status = 1;
  }

  return status;
}

/*....................................................................*/
void
loadNnIntoGrid(unsigned int *firstNearNeigh, struct linkType **nnLinks\
  , struct gridInfoType gridInfoRead, struct grid *gp){
  /*
See the comment at the beginning of the present module for a description of how the pointers 'links', 'nnLinks' and 'firstNearNeigh' relate to the grid struct.

The function mallocs extension 'neigh' of struct g for each grid point.
  */

  const unsigned int totalNumGridPoints = gridInfoRead.nInternalPoints+gridInfoRead.nSinkPoints;
  unsigned int i_ui;
  int j;
  struct linkType *linkPtr;

  for(i_ui=0;i_ui<totalNumGridPoints;i_ui++){
    gp[i_ui].neigh = malloc(sizeof(struct grid *)*gp[i_ui].numNeigh);

    for(j=0;j<gp[i_ui].numNeigh;j++){
      linkPtr = nnLinks[firstNearNeigh[i_ui]+j];
      if(linkPtr->gis[0]==i_ui){
        gp[i_ui].neigh[j] = &gp[linkPtr->gis[1]];
      }else{
        gp[i_ui].neigh[j] = &gp[linkPtr->gis[0]];
      }
    }
  }
}

/*....................................................................*/
void
loadLinkVelsIntoGrid(unsigned int *firstNearNeigh, struct linkType **nnLinks\
  , struct gridInfoType gridInfoRead, struct grid *gp){
  /*
See the comment at the beginning of the present module for a description of how the pointers 'links', 'nnLinks' and 'firstNearNeigh' relate to the grid struct.

The function mallocs (then writes to) the following extensions of struct g for each grid point:
	v1
	v2
	v3

The values are obtained from those stored in the .vels extension of elements of a list of edges (aka links) which are pointed to by elements of the input argument nnLinks. Storing these values in the grid points is x2 redundant, because if grid points A and B are at the ends of a given edge, the same set of values is stored in each.

There is a subtle point about ordering to grasp. The order of the .vels entries in each edge element is prescribed to be the same as the order as stored in the vertex which is listed first in the edge's .gis attribute, i.e. .gis[0]. As we work through the list of grid points, if a given point corresponds to .gis[0] of the edge in question, obvisously we just copy over the velocities; but if it corresponds to .gis[1], we have to reverse the velocity order.
  */

  const unsigned int totalNumGridPoints = gridInfoRead.nInternalPoints+gridInfoRead.nSinkPoints;
  unsigned int i_ui;
  unsigned short i_us;
  int j;
  struct linkType *linkPtr;
  char message[80];

  /* Just for the time being: */
  if(gridInfoRead.nLinkVels!=NUM_VEL_COEFFS){
    if(!silent){
      sprintf(message, "There should be %d VEL_n columns, but %d were read.", NUM_VEL_COEFFS, (int)gridInfoRead.nLinkVels);
      bail_out(message);
    }
    exit(1);
  }

  for(i_ui=0;i_ui<totalNumGridPoints;i_ui++){
    gp[i_ui].v1 = malloc(gridInfoRead.nDims*gp[i_ui].numNeigh*sizeof(double));
    gp[i_ui].v2 = malloc(gridInfoRead.nDims*gp[i_ui].numNeigh*sizeof(double));
    gp[i_ui].v3 = malloc(gridInfoRead.nDims*gp[i_ui].numNeigh*sizeof(double));
  }
  for(i_ui=gridInfoRead.nInternalPoints;i_ui<gridInfoRead.nSinkPoints;i_ui++){
    for(j=0;j<gp[i_ui].numNeigh;j++){
      for(i_us=0;i_us<gridInfoRead.nDims;i_us++){
        gp[i_ui].v1[gridInfoRead.nDims*j+i_us] = 0.0;
        gp[i_ui].v2[gridInfoRead.nDims*j+i_us] = 0.0;
        gp[i_ui].v3[gridInfoRead.nDims*j+i_us] = 0.0;
      }
    }
  }
  for(i_ui=0;i_ui<gridInfoRead.nInternalPoints;i_ui++){
    for(j=0;j<gp[i_ui].numNeigh;j++){
      linkPtr = nnLinks[firstNearNeigh[i_ui]+j];
      if(linkPtr->gis[0]==i_ui){
        for(i_us=0;i_us<gridInfoRead.nDims;i_us++){
          gp[i_ui].v1[gridInfoRead.nDims*j+i_us] = (linkPtr->vels)[gridInfoRead.nDims*0 + i_us];
          gp[i_ui].v2[gridInfoRead.nDims*j+i_us] = (linkPtr->vels)[gridInfoRead.nDims*1 + i_us];
          gp[i_ui].v3[gridInfoRead.nDims*j+i_us] = (linkPtr->vels)[gridInfoRead.nDims*2 + i_us];
        }
      }else{
        for(i_us=0;i_us<gridInfoRead.nDims;i_us++){
          gp[i_ui].v1[gridInfoRead.nDims*j+i_us] = (linkPtr->vels)[gridInfoRead.nDims*2 + i_us];
          gp[i_ui].v2[gridInfoRead.nDims*j+i_us] = (linkPtr->vels)[gridInfoRead.nDims*1 + i_us];
          gp[i_ui].v3[gridInfoRead.nDims*j+i_us] = (linkPtr->vels)[gridInfoRead.nDims*0 + i_us];
        }
      }
    }
  }
}

/*....................................................................*/
int
checkPopsTableExists(lime_fptr *fptr, const int fileFormatI\
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
int
getNumPopsTables(lime_fptr *fptr, const int fileFormatI, unsigned short *numTables){
  int status = 0;
  _Bool blockFound = 1;

  *numTables = 0;
  while(blockFound && !status){
    status = checkPopsTableExists(fptr, fileFormatI, *numTables, &blockFound);
    (*numTables)++;
  }

  (*numTables)--;

  return status;
}

/*....................................................................*/
int
readPopsTable(lime_fptr *fptr, const int fileFormatI\
  , const unsigned short speciesI, struct grid *gp\
  , struct gridInfoType *gridInfoRead){

  int status=0;
  unsigned int totalNumGridPoints, i_ui;

  totalNumGridPoints = gridInfoRead->nInternalPoints + gridInfoRead->nSinkPoints;
  for(i_ui=0;i_ui<totalNumGridPoints;i_ui++){
    gp[i_ui].mol[speciesI].cont    = NULL;
    gp[i_ui].mol[speciesI].pops    = NULL;
    gp[i_ui].mol[speciesI].partner = NULL;
  }

  if(fileFormatI==lime_FITS){
    readPopsExtFromFits(fptr, speciesI, gp, gridInfoRead);
  }else{
    status = 1;
  }

  return status;
}

/*....................................................................*/
int
readGrid(char *inFileName, const int fileFormatI\
  , struct gridInfoType *gridInfoRead, struct keywordType *primaryKwds\
  , const int numKeywords, struct grid **gp, char ***collPartNames, int *numCollPartRead\
  , int *dataFlags){

  /*
This is designed to be a generic function to read the grid data from file. It is assumed that the data will be stored in several tables of different size, corresponding to the different dimensionalities of the elements of the 'grid' struct. See 'writeGrid' for a description.

Some sanity checks are performed here and also in the deeper functions, but any check of the validity of the state of completeness of the grid data (as encoded in the returned argument dataFlags) is left to the calling routine.

NOTE that gp should not be allocated before this routine is called.
  */

  lime_fptr *fptr;
  int status=0;
  unsigned short i_us, numTables;
  unsigned int *firstNearNeigh=NULL, totalNumGridPoints, i_ui;
  struct linkType *links=NULL, **nnLinks=NULL;
  char message[80];

  sprintf(message, "Reading grid-point list from file %s", inFileName);
  if(!silent) printMessage(message);

  *dataFlags = 0;

  /* Set defaults for *gridInfoRead:
  */
  gridInfoRead->nInternalPoints = 0;
  gridInfoRead->nSinkPoints = 0;
  gridInfoRead->nLinks = 0;
  gridInfoRead->nNNIndices = 0;
  gridInfoRead->nDims = 0;
  gridInfoRead->nSpecies = 0;
  gridInfoRead->nDensities = 0;
  gridInfoRead->nLinkVels = 0;
  gridInfoRead->mols = NULL;

  /* Open the file and also return the data stage. */
  fptr = openFileForRead(inFileName, fileFormatI);

  /* Read header keywords:
  */
  status = readKeywords(fptr, fileFormatI, primaryKwds, numKeywords);
  if(status){
    closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, 0);
    return 1;
  }

  /* Read the values which should be in grid for every stage.
  */
  status = readGridTable(fptr, fileFormatI, gridInfoRead, gp, &firstNearNeigh\
    , collPartNames, numCollPartRead, dataFlags); /* Sets appropriate bits of dataFlags; also mallocs gp and sets all its defaults. */
  totalNumGridPoints = gridInfoRead->nSinkPoints + gridInfoRead->nInternalPoints;
  if(status){
    closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, 0);
    freeGrid(totalNumGridPoints, gridInfoRead->nSpecies, *gp);
    return 2;
  }

  /* Some sanity checks:
  */
  if((*dataFlags)!=0 && gridInfoRead->nSinkPoints<=0 || gridInfoRead->nInternalPoints<=0 || gridInfoRead->nDims<=0){
    closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, 0);
    freeGrid(totalNumGridPoints, gridInfoRead->nSpecies, *gp);

    if(gridInfoRead->nSinkPoints<=0)
      return 3;
    else if(gridInfoRead->nInternalPoints<=0)
      return 4;
    else if(gridInfoRead->nDims<=0)
      return 5;
    else{
      sprintf(message, "This indicates a programming error. Please contact the developer.");
      if(!silent) bail_out(message);
      exit(1);
    }
  }

  /* Some more sanity checks:
  */
  if((allBitsSet(*dataFlags, DS_mask_density)   && gridInfoRead->nDensities<=0)\
  || (allBitsSet(*dataFlags, DS_mask_abundance) && gridInfoRead->nSpecies<=0)){
    closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, 0);
    freeGrid(totalNumGridPoints, gridInfoRead->nSpecies, *gp);

    if(gridInfoRead->nDensities<=0)
      return 6;
    else if(gridInfoRead->nSpecies<=0) //***** what if all continuum images??
      return 7;
    else{
      sprintf(message, "This indicates a programming error. Please contact the developer.");
      if(!silent) bail_out(message);
      exit(1);
    }
  }

  status = readLinksTable(fptr, fileFormatI, gridInfoRead, *gp, &links, dataFlags); /* Sets appropriate bits of dataFlags. */
  if(status){
    closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, gridInfoRead->nLinks);
    freeGrid(totalNumGridPoints, gridInfoRead->nSpecies, *gp);
    return 8;
  }

  /* Sanity check:
  */
  if (allBitsSet(*dataFlags, DS_mask_ACOEFF) && gridInfoRead->nLinkVels<=0){
    closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, gridInfoRead->nLinks);
    freeGrid(totalNumGridPoints, gridInfoRead->nSpecies, *gp);
    return 9;
  }

  status = readNnIndicesTable(fptr, fileFormatI, links, &nnLinks, gridInfoRead, dataFlags); /* Sets appropriate bits of dataFlags. */
  if(status){
    closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, gridInfoRead->nLinks);
    freeGrid(totalNumGridPoints, gridInfoRead->nSpecies, *gp);
    return 10;
  }

  if(allBitsSet(*dataFlags, DS_mask_neighbours)){
    /* Convert the NN information back to the standard LIME grid struct format.
    */
    loadNnIntoGrid(firstNearNeigh, nnLinks, *gridInfoRead, *gp); /* mallocs extension 'neigh' of struct g for each grid point. */

    if(allBitsSet(*dataFlags, DS_mask_ACOEFF))
      loadLinkVelsIntoGrid(firstNearNeigh, nnLinks, *gridInfoRead, *gp);
  }

  status = getNumPopsTables(fptr, fileFormatI, &numTables);
  if(status){
    closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, gridInfoRead->nLinks);
    freeGrid(totalNumGridPoints, gridInfoRead->nSpecies, *gp);
    return 11;
  }

  /* Sanity check:
  */
  if(numTables>0 && numTables!=gridInfoRead->nSpecies){
    closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, gridInfoRead->nLinks);
    freeGrid(totalNumGridPoints, gridInfoRead->nSpecies, *gp);
    return 12;
  }

  if(numTables>0){
    (*dataFlags) |= (1 << DS_bit_populations);

    for(i_ui=0;i_ui<totalNumGridPoints;i_ui++)
      (*gp)[i_ui].mol = malloc(sizeof(struct populations)*gridInfoRead->nSpecies);

    gridInfoRead->mols = malloc(sizeof(struct molInfoType)*gridInfoRead->nSpecies);
    for(i_us=0;i_us<gridInfoRead->nSpecies;i_us++){
      gridInfoRead->mols[i_us].molName = NULL;
      gridInfoRead->mols[i_us].nLevels = -1;
      gridInfoRead->mols[i_us].nLines = -1;

      status = readPopsTable(fptr, fileFormatI, i_us, *gp, gridInfoRead); /* Sets defaults for all the fields under grid.mol. */
      if(status){
        closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, gridInfoRead->nLinks);
        freeGrid(totalNumGridPoints, gridInfoRead->nSpecies, *gp);
        return 13;
      }

      /* Sanity check:
      */
      if(gridInfoRead->mols[i_us].nLevels<=0){
        closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, gridInfoRead->nLinks);
        freeGrid(totalNumGridPoints, gridInfoRead->nSpecies, *gp);
        return 14;
      }
    }
  }

  closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, gridInfoRead->nLinks);
  return 0;
}



