/*
 *  gridio.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
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

	Data mask bits set:	    stage 0	    stage 1	    stage 2	    stage 3	    stage 4	    stage 5
	.....................................................................................................................
	DS_bit_x             		0		1		1		1		1		1
	DS_bit_neighbours    		0		0		1		1		1		1
	DS_bit_velocity      		0		0		0		0		1		1
	DS_bit_density       		0		0		0		1		1		1
	DS_bit_abundance     		0		0		0		0		1		1
	DS_bit_turb_doppler  		0		0		0		0		1		1
	DS_bit_temperatures  		0		0		0		1		1		1
	DS_bit_magfield  		0		0		0		x		x		x
	DS_bit_ACOEFF        		0		0		0		0		1		1
	DS_bit_populations   		0		0		0		0		0		1
	.....................................................................................................................

Notes:
  - dataStageI==0 has been included for completeness/robustness but the user may not write a file with nothing in it.
  - Stage 3 is as far as LIME gets when only continuum images/calculations are required.

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
void
initializeKeyword(struct keywordType *kwd){
  (*kwd).datatype = 0;
  (*kwd).keyname = malloc(sizeof(char)*STRLEN_KNAME);
  (*kwd).comment = malloc(sizeof(char)*STRLEN_KCOMM);
  (*kwd).intValue = 0;
  (*kwd).floatValue = 0.0;
  (*kwd).doubleValue = 0.0;
  (*kwd).charValue = malloc(sizeof(char)*STRLEN_KCHAR);
}

/*....................................................................*/
void
freeKeyword(struct keywordType kwd){
  free(kwd.keyname);
  free(kwd.comment);
  free(kwd.charValue);
}

/*....................................................................*/
void
freeKeywords(struct keywordType *kwds, const int numKwds){
  int i;

  if(kwds!=NULL){
    for(i=0;i<numKwds;i++)
      freeKeyword(kwds[i]);
    free(kwds);
  }
}

/*....................................................................*/
void
freeGridInfo(struct gridInfoType *gridInfo){
  unsigned short i_us;

  if(gridInfo->mols!=NULL){
    for(i_us=0;i_us<gridInfo->nSpecies;i_us++)
      free(gridInfo->mols[i_us].molName);
    free(gridInfo->mols);
  }
}

/*....................................................................*/
lime_fptr
openFileForWrite(char *outFileName){
  lime_fptr fptr=lime_init;

#if defined(lime_IO) && lime_IO==lime_HDF5
  fptr = openHDF5FileForWrite(outFileName);
#else
  fptr = openFITSFileForWrite(outFileName);
#endif

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
  unsigned int li,ni,idA,idB,trialIdA,linkId=0,i_ui,gi0,gi1;
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
closeFile(lime_fptr fptr){
#if defined(lime_IO) && lime_IO==lime_HDF5
  closeHDF5File(fptr);
#else
  closeFITSFile(fptr);
#endif
}

/*....................................................................*/
void
closeAndFree(lime_fptr fptr\
  , unsigned int *firstNearNeigh, struct linkType **nnLinks\
  , struct linkType *links, const unsigned int totalNumLinks){

  unsigned int i_ui;

  closeFile(fptr);

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
writeKeywords(lime_fptr fptr\
  , struct keywordType *kwds, const int numKeywords){

  int status=0;

#if defined(lime_IO) && lime_IO==lime_HDF5
  writeKeywordsToHDF5(fptr, kwds, numKeywords);
#else
  writeKeywordsToFITS(fptr, kwds, numKeywords);
#endif

  return status;
}

/*....................................................................*/
int
writeGridTable(lime_fptr fptr\
  , struct gridInfoType gridInfo, struct grid *gp, unsigned int *firstNearNeigh\
  , char **collPartNames, const int dataFlags){

  int status=0;

#if defined(lime_IO) && lime_IO==lime_HDF5
  writeGridExtToHDF5(fptr, gridInfo, gp, firstNearNeigh, collPartNames, dataFlags);
#else
  writeGridExtToFITS(fptr, gridInfo, gp, firstNearNeigh, collPartNames, dataFlags);
#endif

  return status;
}

/*....................................................................*/
int
writeNnIndicesTable(lime_fptr fptr\
  , struct gridInfoType gridInfo, struct linkType **nnLinks){\

  int status=0;

#if defined(lime_IO) && lime_IO==lime_HDF5
  writeNnIndicesExtToHDF5(fptr, gridInfo, nnLinks);
#else
  writeNnIndicesExtToFITS(fptr, gridInfo, nnLinks);
#endif

  return status;
}

/*....................................................................*/
int
writeLinksTable(lime_fptr fptr\
  , struct gridInfoType gridInfo, struct linkType *links){

  int status=0;

#if defined(lime_IO) && lime_IO==lime_HDF5
  writeLinksExtToHDF5(fptr, gridInfo, links);
#else
  writeLinksExtToFITS(fptr, gridInfo, links);
#endif

  return status;
}

/*....................................................................*/
int
writePopsTable(lime_fptr fptr\
  , struct gridInfoType gridInfo, unsigned short speciesI\
  , struct grid *gp){


  int status=0;

#if defined(lime_IO) && lime_IO==lime_HDF5
  writePopsGroupToHDF5(fptr, gridInfo, speciesI, gp);
#else
  writePopsExtToFITS(fptr, gridInfo, speciesI, gp);
#endif

  return status;
}

/*....................................................................*/
int
writeGrid(char *outFileName, struct gridInfoType gridInfo\
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

  lime_fptr fptr=lime_init;
  int status = 0;
  unsigned short i_us;
  char message[80];
  struct linkType *links=NULL, **nnLinks=NULL;
  unsigned int *firstNearNeigh=NULL;

  if (outFileName==NULL || strlen(outFileName)<=0){
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

  fptr = openFileForWrite(outFileName);
  if(fptr _FAILED_TO_OPEN){
    closeAndFree(fptr, firstNearNeigh, nnLinks, links, gridInfo.nLinks);
    return 3;
  }

  /* Write header keywords:
  */
  status = writeKeywords(fptr, primaryKwds, numKeywords);
  if(status){
    closeAndFree(fptr, firstNearNeigh, nnLinks, links, gridInfo.nLinks);
    return 4;
  }

  status = writeGridTable(fptr, gridInfo, gp, firstNearNeigh, collPartNames, dataFlags);
  if(status){
    closeAndFree(fptr, firstNearNeigh, nnLinks, links, gridInfo.nLinks);
    return 5;
  }

  if (links!=NULL && nnLinks!=NULL){
    status = writeNnIndicesTable(fptr, gridInfo, nnLinks);
    if(status){
      closeAndFree(fptr, firstNearNeigh, nnLinks, links, gridInfo.nLinks);
      return 6;
    }

    status = writeLinksTable(fptr, gridInfo, links);
    if(status){
      closeAndFree(fptr, firstNearNeigh, nnLinks, links, gridInfo.nLinks);
      return 7;
    }
  }

  if (allBitsSet(dataFlags, DS_mask_populations)){
    for(i_us=0;i_us<gridInfo.nSpecies;i_us++){
      status = writePopsTable(fptr, gridInfo, i_us, gp);
      if(status){
        closeAndFree(fptr, firstNearNeigh, nnLinks, links, gridInfo.nLinks);
        return 8;
      }
    }
  }

  closeAndFree(fptr, firstNearNeigh, nnLinks, links, gridInfo.nLinks);
  return 0;
}


/*....................................................................*/
lime_fptr
openFileForRead(char *inFileName){
  lime_fptr fptr=lime_init;

#if defined(lime_IO) && lime_IO==lime_HDF5
  fptr = openHDF5FileForRead(inFileName);
#else
  fptr = openFITSFileForRead(inFileName);
#endif

  return fptr;
}

/*....................................................................*/
int
readKeywords(lime_fptr fptr\
  , struct keywordType *kwds, const int numKeywords){

  int status=0;

#if defined(lime_IO) && lime_IO==lime_HDF5
  readKeywordsFromHDF5(fptr, kwds, numKeywords);
#else
  readKeywordsFromFITS(fptr, kwds, numKeywords);
#endif

  return status;
}

/*....................................................................*/
int
_readGridTable(lime_fptr fptr, struct gridInfoType *gridInfoRead\
  , struct grid **gp, unsigned int **firstNearNeigh, char ***collPartNames\
  , int *numCollPartRead, int *dataFlags, _Bool *densMolColsExists){
  /*
Individual routines called should set the appropriate bits of dataFlags; also malloc gp and set all its defaults. (Note there is a bespoke routine grid.c:mallocAndSetDefaultGrid() to do the latter.)
  */

  int status=0;

#if defined(lime_IO) && lime_IO==lime_HDF5
  readGridExtFromHDF5(fptr, gridInfoRead, gp, firstNearNeigh\
    , collPartNames, numCollPartRead, dataFlags, densMolColsExists);
#else
  readGridExtFromFITS(fptr, gridInfoRead, gp, firstNearNeigh\
    , collPartNames, numCollPartRead, dataFlags, densMolColsExists);
#endif

  return status;
}

/*....................................................................*/
int
readLinksTable(lime_fptr fptr\
  , struct gridInfoType *gridInfoRead, struct grid *gp\
  , struct linkType **links, int *dataFlags){
  /*
Individual routines called should set the appropriate bits of dataFlags.
  */

  int status=0;

#if defined(lime_IO) && lime_IO==lime_HDF5
  readLinksExtFromHDF5(fptr, gridInfoRead, gp, links, dataFlags);
#else
  readLinksExtFromFITS(fptr, gridInfoRead, gp, links, dataFlags);
#endif

  return status;
}

/*....................................................................*/
int
readNnIndicesTable(lime_fptr fptr, struct linkType *links\
  , struct linkType ***nnLinks, struct gridInfoType *gridInfoRead, int *dataFlags){
  /*
Individual routines called should set the appropriate bits of dataFlags.
  */

  int status=0;

#if defined(lime_IO) && lime_IO==lime_HDF5
  readNnIndicesExtFromHDF5(fptr, links, nnLinks, gridInfoRead, dataFlags);
#else
  readNnIndicesExtFromFITS(fptr, links, nnLinks, gridInfoRead, dataFlags);
#endif

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
  for(i_ui=0;i_ui<totalNumGridPoints;i_ui++){
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
checkPopsTableExists(lime_fptr fptr\
  , const unsigned short speciesI, _Bool *blockFound){
  int status=0;

  *blockFound = 0;

#if defined(lime_IO) && lime_IO==lime_HDF5
  *blockFound = checkPopsHDF5GroupExists(fptr, speciesI);
#else
  *blockFound = checkPopsFITSExtExists(fptr, speciesI);
#endif

  return status;
}

/*....................................................................*/
int
getNumPopsTables(lime_fptr fptr, unsigned short *numTables){
  int status = 0;
  _Bool blockFound = 1;

  *numTables = 0;
  while(blockFound && !status){
    status = checkPopsTableExists(fptr, *numTables, &blockFound);
    (*numTables)++;
  }

  (*numTables)--;

  return status;
}

/*....................................................................*/
int
readPopsTable(lime_fptr fptr\
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

#if defined(lime_IO) && lime_IO==lime_HDF5
  readPopsGroupFromHDF5(fptr, speciesI, gp, gridInfoRead);
#else
  readPopsExtFromFITS(  fptr, speciesI, gp, gridInfoRead);
#endif

  return status;
}

/*....................................................................*/
int
readGrid(char *inFileName, struct gridInfoType *gridInfoRead\
  , struct keywordType *primaryKwds, const int numKeywords\
  , struct grid **gp, char ***collPartNames, int *numCollPartRead\
  , int *dataFlags, _Bool *densMolColsExists){

  /*
This is designed to be a generic function to read the grid data from file. It is assumed that the data will be stored in several tables of different size, corresponding to the different dimensionalities of the elements of the 'grid' struct. See 'writeGrid' for a description.

Some sanity checks are performed here and also in the deeper functions, but any check of the validity of the state of completeness of the grid data (as encoded in the returned argument dataFlags) is left to the calling routine.

NOTE that gp should not be allocated before this routine is called, but it must be freed by the calling routine after use.

NOTE that collPartNames and its components must be freed after use.
  */

  lime_fptr fptr=lime_init;
  int status=0;
  unsigned short i_us, numTables;
  unsigned int *firstNearNeigh=NULL, totalNumGridPoints;
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
  fptr = openFileForRead(inFileName);

  /* Read header keywords:
  */
  status = readKeywords(fptr, primaryKwds, numKeywords);
  if(status){
    closeAndFree(fptr, firstNearNeigh, nnLinks, links, 0);
    return 1;
  }

  /* Read the values which should be in grid for every stage.
  */
  status = _readGridTable(fptr, gridInfoRead, gp, &firstNearNeigh\
    , collPartNames, numCollPartRead, dataFlags, densMolColsExists); /* Sets appropriate bits of dataFlags; also mallocs gp and sets all its defaults. */
  totalNumGridPoints = gridInfoRead->nSinkPoints + gridInfoRead->nInternalPoints;
  if(status){
    closeAndFree(fptr, firstNearNeigh, nnLinks, links, 0);
    freeGrid(totalNumGridPoints, gridInfoRead->nSpecies, *gp);
    return 2;
  }

  /* Some sanity checks:
  */
  if((*dataFlags)!=0 && (gridInfoRead->nSinkPoints<=0 || gridInfoRead->nInternalPoints<=0 || gridInfoRead->nDims<=0)){
    closeAndFree(fptr, firstNearNeigh, nnLinks, links, 0);
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
    closeAndFree(fptr, firstNearNeigh, nnLinks, links, 0);
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

  status = readLinksTable(fptr, gridInfoRead, *gp, &links, dataFlags); /* Sets appropriate bits of dataFlags. */
  if(status){
    closeAndFree(fptr, firstNearNeigh, nnLinks, links, gridInfoRead->nLinks);
    freeGrid(totalNumGridPoints, gridInfoRead->nSpecies, *gp);
    return 8;
  }

  /* Sanity check:
  */
  if (allBitsSet(*dataFlags, DS_mask_ACOEFF) && gridInfoRead->nLinkVels<=0){
    closeAndFree(fptr, firstNearNeigh, nnLinks, links, gridInfoRead->nLinks);
    freeGrid(totalNumGridPoints, gridInfoRead->nSpecies, *gp);
    return 9;
  }

  status = readNnIndicesTable(fptr, links, &nnLinks, gridInfoRead, dataFlags); /* Sets appropriate bits of dataFlags. */
  if(status){
    closeAndFree(fptr, firstNearNeigh, nnLinks, links, gridInfoRead->nLinks);
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

  status = getNumPopsTables(fptr, &numTables);
  if(status){
    closeAndFree(fptr, firstNearNeigh, nnLinks, links, gridInfoRead->nLinks);
    freeGrid(totalNumGridPoints, gridInfoRead->nSpecies, *gp);
    return 11;
  }

  /* Sanity check:
  */
  if(numTables>0 && numTables!=gridInfoRead->nSpecies){
    closeAndFree(fptr, firstNearNeigh, nnLinks, links, gridInfoRead->nLinks);
    freeGrid(totalNumGridPoints, gridInfoRead->nSpecies, *gp);
    return 12;
  }

  if(numTables>0){
    (*dataFlags) |= (1 << DS_bit_populations);

    gridInfoRead->mols = malloc(sizeof(struct molInfoType)*gridInfoRead->nSpecies);
    for(i_us=0;i_us<gridInfoRead->nSpecies;i_us++){
      gridInfoRead->mols[i_us].molName = NULL;
      gridInfoRead->mols[i_us].nLevels = -1;
      gridInfoRead->mols[i_us].nLines = -1;

      status = readPopsTable(fptr, i_us, *gp, gridInfoRead); /* Sets defaults for all the fields under grid.mol. */
      if(status){
        closeAndFree(fptr, firstNearNeigh, nnLinks, links, gridInfoRead->nLinks);
        freeGrid(totalNumGridPoints, gridInfoRead->nSpecies, *gp);
        return 13;
      }

      /* Sanity check:
      */
      if(gridInfoRead->mols[i_us].nLevels<=0){
        closeAndFree(fptr, firstNearNeigh, nnLinks, links, gridInfoRead->nLinks);
        freeGrid(totalNumGridPoints, gridInfoRead->nSpecies, *gp);
        return 14;
      }
    }
  }

  closeAndFree(fptr, firstNearNeigh, nnLinks, links, gridInfoRead->nLinks);

  return 0;
}

/*....................................................................*/
int
countDensityCols(char *inFileName, int *numDensities){
  int status=0;

#if defined(lime_IO) && lime_IO==lime_HDF5
  *numDensities = countDensityColsHDF5(inFileName);
#else
  *numDensities = countDensityColsFITS(inFileName);
#endif

  return status;
}

/*....................................................................*/
void
sanityCheckOfRead(const int status, configInfo *par, struct gridInfoType gridInfoRead){
  char message[STR_LEN_0];

  if(status){
    if(!silent){
      snprintf(message, STR_LEN_0, "Read of grid file failed with status return %d", status);
      bail_out(message);
    }
exit(1);
  }

  /* Test that dataFlags obeys the rules. */
  /* No other bit may be set if DS_bit_x is not: */
  if(anyBitSet(par->dataFlags, (DS_mask_all & ~(1 << DS_bit_x))) && !bitIsSet(par->dataFlags, DS_bit_x)){
    if(!silent) bail_out("You may not read a grid file without X, ID or IS_SINK data.");
exit(1);
  }

  /* DS_bit_ACOEFF may not be set if either DS_bit_neighbours or DS_bit_velocity is not: */
  if(bitIsSet(par->dataFlags, DS_bit_ACOEFF)\
  && !(bitIsSet(par->dataFlags, DS_bit_neighbours) && bitIsSet(par->dataFlags, DS_bit_velocity))){
    if(!silent) bail_out("You may not read a grid file with ACOEFF but no VEL or neighbour data.");
exit(1);
  }

  /* DS_bit_populations may not be set unless all the others (except DS_bit_magfield) are set as well: */
  if(bitIsSet(par->dataFlags, DS_bit_populations)\
  && !allBitsSet(par->dataFlags & DS_mask_all_but_mag, DS_mask_populations)){
    if(!silent) bail_out("You may not read a grid file with pop data unless all other data is present.");
exit(1);
  }

  /* Test gridInfoRead values against par values and overwrite the latter, with a warning, if necessary.
  */
  if(gridInfoRead.nSinkPoints>0 && par->sinkPoints>0){
    if((int)gridInfoRead.nSinkPoints!=par->sinkPoints){
      if(!silent) warning("par->sinkPoints will be overwritten");
    }
    if((int)gridInfoRead.nInternalPoints!=par->pIntensity){
      if(!silent) warning("par->pIntensity will be overwritten");
    }
  }
  par->sinkPoints = (int)gridInfoRead.nSinkPoints;
  par->pIntensity = (int)gridInfoRead.nInternalPoints;
  par->ncell = par->sinkPoints + par->pIntensity;

  if(gridInfoRead.nDims!=DIM){ /* At present this situation is already detected and handled inside readGridExtFromFits(), but it may not be in future. The test here has no present functionality but saves trouble later if we change grid.x from an array to a pointer. */
    if(!silent){
      snprintf(message, STR_LEN_0, "Grid file had %d dimensions but there should be %d.", (int)gridInfoRead.nDims, DIM);
      bail_out(message);
    }
exit(1);
  }
  if(gridInfoRead.nSpecies > 0){
    if((int)gridInfoRead.nSpecies!=par->nSpecies && par->doMolCalcs){
      if(!silent){
        snprintf(message, STR_LEN_0, "Grid file had %d species but you have provided moldata files for %d."\
          , (int)gridInfoRead.nSpecies, par->nSpecies);
        bail_out(message);
      }
exit(1);
/**** should compare name to name - at some later time after we have read these from the moldata files? */
    }
  }
  if(gridInfoRead.nDensities>0 && par->numDensities>0 && (int)gridInfoRead.nDensities!=par->numDensities){
    if(!silent){
      snprintf(message, STR_LEN_0, "Grid file had %d densities but you have provided %d."\
        , (int)gridInfoRead.nDensities, par->numDensities);
      bail_out(message);
    }
exit(1);
  }
  par->numDensities = gridInfoRead.nDensities;

  if(par->radius>0){ /* means the user has set it to something in the model file. */
    if(!silent)
      warning("Your value of par->radius will be overwritten by that read from file.");
  }

  if(par->minScale>0){ /* means the user has set it to something in the model file. */
    if(!silent)
      warning("Your value of par->minScale will be overwritten by that read from file.");
  }

  if(par->nSolveItersDone>0 && (par->init_lte || par->lte_only)){
    if(!silent)
      warning("Your choice of LTE calculation will erase the RTE solution you read from file.");
  }

  if(allBitsSet(par->dataFlags, DS_mask_populations) && par->nSolveItersDone<=0){
    if(!silent)
      bail_out("Populations were read but par->nSolveItersDone<=0.");
exit(1);
  }
}

/*....................................................................*/
void
readGridWrapper(configInfo *par, struct grid **gp, char ***collPartNames\
  , int *numCollPartRead){

  const int numDesiredKwds=3;
  struct keywordType *desiredKwds=malloc(sizeof(struct keywordType)*numDesiredKwds);
  int i,status=0;
  struct gridInfoType gridInfoRead;
  _Bool densMolColsExists;

  i = 0;
  initializeKeyword(&desiredKwds[i]);
  desiredKwds[i].datatype = lime_DOUBLE;
  snprintf(desiredKwds[i].keyname, STRLEN_KNAME, "RADIUS  ");

  i++;
  initializeKeyword(&desiredKwds[i]);
  desiredKwds[i].datatype = lime_DOUBLE;
  snprintf(desiredKwds[i].keyname, STRLEN_KNAME, "MINSCALE");

  i++;
  initializeKeyword(&desiredKwds[i]);
  desiredKwds[i].datatype = lime_INT;
  snprintf(desiredKwds[i].keyname, STRLEN_KNAME, "NSOLITER");

  status = readGrid(par->gridInFile, &gridInfoRead, desiredKwds\
    , numDesiredKwds, gp, collPartNames, numCollPartRead\
    , &(par->dataFlags), &densMolColsExists);

  par->nSolveItersDone = desiredKwds[2].intValue;

  sanityCheckOfRead(status, par, gridInfoRead);

  if(gridInfoRead.nSpecies>0)
    par->useAbun = !densMolColsExists;

  par->radius          = desiredKwds[0].doubleValue;
  par->minScale        = desiredKwds[1].doubleValue;

  par->radiusSqu   = par->radius*par->radius;
  par->minScaleSqu = par->minScale*par->minScale;

  freeKeywords(desiredKwds, numDesiredKwds);
  freeGridInfo(&gridInfoRead);

/*
**** Ideally we should also have a test on nACoeffs.

**** Ideally we should also have a test on the mols entries - at some later time after we have read the corresponding values from the moldata files?
*/
}

