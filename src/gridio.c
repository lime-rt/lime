/*
 *  gridio.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2016 The LIME development team
 *
 */

#include "lime.h"

/*
This module contains generic routines for writing grid data to, and reading it from, a file on disk. The problem with doing this is that the grid struct contains different amounts of information at different times in the running of the code. In order to quantify and regulate this, a dataFlags integer is used to record the presence or absence (as indicated by the value of the appropriate bit in the mask) of particular types of information. The bits associated with certain fields of struct grid are given in the lime.h header.

Notes:
  - If several fields of the struct are listed for a given bit, a value of 1 for that bit indicates that all fields are expected to be present, and a value of 0 indicates that all fields of that bit will be ignored.
  - Only a few combinations are not allowed, as follows:
    * No other bit may be set if DS_bit_x is not.
    * DS_bit_ACOEFF may not be set if either DS_bit_neighbours or DS_bit_velocity is not.
    * DS_bit_populations may not be set unless all the others are set as well.

Writing the grid file:
----------------------
To make things simpler, four stages have been defined at which the user may write the grid data to file. These are described in the following table:

	Data mask bits set:	dataStageI=0	dataStageI=1	dataStageI=2	dataStageI=3	dataStageI=4
	.....................................................................................................
	DS_bit_x             		0		1		1		1		1
	DS_bit_neighbours    		0		0		1		1		1
	DS_bit_velocity      		0		0		0		1		1
	DS_bit_density       		0		0		0		1		1
	DS_bit_abundance     		0		0		0		1		1
	DS_bit_turb_doppler  		0		0		0		1		1
	DS_bit_temperatures  		0		0		0		1		1
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
void writeGridIfRequired(configInfo *par, struct grid *gp, molData *md, const int fileFormatI){
  int status = 0;
  char **collPartNames=NULL; /*** this is a placeholder until we start reading these. */
  char message[80];
  int dataStageI=0;

  /* Work out the data stage:
  */
  if(!allBitsSet(par->dataFlags, DS_mask_1)){
    if(!silent) warning("Trying to write at data stage 0.");
    return;
  }

  if(      allBitsSet(par->dataFlags, DS_mask_4)){
    dataStageI = 4;
  }else if(allBitsSet(par->dataFlags, DS_mask_3)){
    dataStageI = 3;
  }else if(allBitsSet(par->dataFlags, DS_mask_2)){
    dataStageI = 2;
  }else{
    dataStageI = 1;
  }

  if(par->writeGridAtStage[dataStageI-1]){
    status = writeGrid(par->gridOutFiles[dataStageI-1], fileFormatI\
      , *par, DIM, NUM_VEL_COEFFS, gp, md, collPartNames, par->dataFlags);

    if(status){
      sprintf(message, "writeGrid at data stage %d returned with status %d", dataStageI, status);
      if(!silent) bail_out(message);
      exit(1);
    }
  }

  free(collPartNames);
}

/*....................................................................*/
lime_fptr *openFileForWrite(char *outFileName, const int fileFormatI){
  lime_fptr *fptr=NULL;

  if(fileFormatI==lime_FITS){
    fptr = openFITSFileForWrite(outFileName);
  }else{
    return NULL;
  }

  return fptr;
}

/*....................................................................*/
void constructLinkArrays(const unsigned int numGridPoints, struct grid *gp\
  , struct linkType **links, unsigned int *totalNumLinks, struct linkType ***nnLinks\
  , unsigned int **firstNearNeigh, unsigned int *totalNumNeigh, const int dataFlags){
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

  if(!allBitsSet(dataFlags, DS_mask_neighbours)) /* Grid fields id, sink, numNeigh and neigh must be available. */
    return;

  pointIsDone     = malloc(sizeof(*pointIsDone)    *numGridPoints);
  *firstNearNeigh = malloc(sizeof(**firstNearNeigh)*numGridPoints);

  /* First add up all the numNeigh.
  */
  ni = 0;
  for(i=0;i<numGridPoints;i++){
    ni += gp[i].numNeigh;
    pointIsDone[gp[i].id] = 0;
  }
  *totalNumNeigh = ni;

  *links   = malloc(sizeof(**links)*(*totalNumNeigh)); /* Bigger than needed, we'll reduce it later. */
  nnLinkIs = malloc(sizeof(*nnLinkIs)*(*totalNumNeigh));

  li = 0;
  ni = 0;
  for(i=0;i<numGridPoints;i++){
    gAPtr = gp+i;
    idA = gp[i].id;
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
        (*links)[li].gp[0] = gAPtr;
        (*links)[li].gp[1] = gBPtr;

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
    for(li=0;li<*totalNumLinks;li++)
      (*links)[li].aCoeffs = malloc(sizeof(double)*NUM_VEL_COEFFS);

    for(li=0;li<*totalNumLinks;li++){
      if((*links)[li].gp[0]->sink && (*links)[li].gp[1]->sink){
        for(ci=0;ci<NUM_VEL_COEFFS;ci++)
          (*links)[li].aCoeffs[ci] = 0.0;
      }else{

        /* If gp[0] is a sink point then try gp[1], because sink point a*'s are set to zero. Remember to invert the sign of the even coefficients if we are reading them from gp[1]. (See the header comments to loadNnIntoGrid() for an explanation of the reason for the sign inversion.)
        */
        if((*links)[li].gp[0]->sink){ /* If we get to here, this ensures that gp[1] is not a sink point. */
          nearI = 1;
          evenCoeffSign = -1.0;
        }else{
          nearI = 0;
          evenCoeffSign = 1.0;
        }

        gAPtr = (*links)[li].gp[nearI];
        /* Find which neighbour of gAPtr corresponds to the link: */
        linkNotFound = 1;
        for(jA=0;jA<gAPtr->numNeigh;jA++){
          idB = gAPtr->neigh[jA]->id;
          if(linkNotFound && idB==(*links)[li].gp[1-nearI]->id){
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
  }else{ /* a coeffs not available */
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
void closeFile(lime_fptr *fptr, const int fileFormatI){
  if(fileFormatI==lime_FITS){
    closeFITSFile(fptr);
  }//**** error if not?
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
int writeGridTable(lime_fptr *fptr, const int fileFormatI, configInfo par\
  , unsigned short numDims, struct grid *gp, unsigned int *firstNearNeigh\
  , char **collPartNames, const int dataFlags){

  int status=0;

  if(fileFormatI==lime_FITS){
    writeGridExtToFits(fptr, par, numDims, gp, firstNearNeigh, collPartNames, dataFlags);
  }else{
    status = 1;
  }

  return status;
}

/*....................................................................*/
int writeNnIndicesTable(lime_fptr *fptr, const int fileFormatI\
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
int writeLinksTable(lime_fptr *fptr, const int fileFormatI\
  , const unsigned int totalNumLinks, const unsigned short numACoeffs\
  , struct linkType *links){

  int status=0;

  if(fileFormatI==lime_FITS){
    writeLinksExtToFits(fptr, totalNumLinks, numACoeffs, links);
  }else{
    status = 1;
  }

  return status;
}

/*....................................................................*/
int writePopsTable(lime_fptr *fptr, const int fileFormatI\
  , unsigned int numGridPoints, molData *md, unsigned short speciesI\
  , struct grid *gp){

  int status=0;

  if(fileFormatI==lime_FITS){
    writePopsExtToFits(fptr, numGridPoints, md, speciesI, gp);
  }else{
    status = 1;
  }

  return status;
}

/*....................................................................*/
int writeGrid(char *outFileName, const int fileFormatI, configInfo par\
  , unsigned short numDims, unsigned short numACoeffs, struct grid *gp\
  , molData *md, char **collPartNames, const int dataFlags){

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
  unsigned short speciesI;
  char message[80];
  struct linkType *links=NULL, **nnLinks=NULL;
  unsigned int totalNumLinks=-1, totalNumNeigh=-1, *firstNearNeigh=NULL;

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

  constructLinkArrays((unsigned int)par.ncell, gp, &links, &totalNumLinks\
    , &nnLinks, &firstNearNeigh, &totalNumNeigh, dataFlags);

  fptr = openFileForWrite(outFileName, fileFormatI);
  if(fptr==NULL){
    closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, totalNumLinks);
    return 3;
  }

  status = writeGridTable(fptr, fileFormatI, par, numDims, gp, firstNearNeigh, collPartNames, dataFlags);
  if(status){
    closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, totalNumLinks);
    return 4;
  }

  if (links!=NULL && nnLinks!=NULL){
    status = writeNnIndicesTable(fptr, fileFormatI, totalNumNeigh, nnLinks, links);
    if(status){
      closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, totalNumLinks);
      return 5;
    }

    status = writeLinksTable(fptr, fileFormatI, totalNumLinks, numACoeffs, links);
    if(status){
      closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, totalNumLinks);
      return 6;
    }
  }

  if (allBitsSet(dataFlags, DS_mask_populations)){
    for(speciesI=0;speciesI<par.nSpecies;speciesI++){
      status = writePopsTable(fptr, fileFormatI, (unsigned int)par.ncell, md, speciesI, gp);
      if(status){
        closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, totalNumLinks);
        return 7;
      }
    }
  }

  closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, totalNumLinks);
  return 0;
}


/*....................................................................*/
lime_fptr *openFileForRead(char *inFileName, const int fileFormatI){
  lime_fptr *fptr=NULL;

  if(fileFormatI==lime_FITS){
    fptr = openFITSFileForRead(inFileName);
  }else{
    return NULL;
  }

  return fptr;
}

/*....................................................................*/
int readGridTable(lime_fptr *fptr, const int fileFormatI\
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
int readLinksTable(lime_fptr *fptr, const int fileFormatI\
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
int readNnIndicesTable(lime_fptr *fptr, const int fileFormatI, struct linkType *links\
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
void loadNnIntoGrid(unsigned int *firstNearNeigh, struct linkType **nnLinks\
  , struct gridInfoType gridInfoRead, struct grid *gp){
  /*
See the comment at the beginning of the present module for a description of how the pointers 'links', 'nnLinks' and 'firstNearNeigh' relate to the grid struct.

The function mallocs extension 'neigh' of struct g for each grid point.
  */

  const unsigned int totalNumGridPoints = gridInfoRead.nInternalPoints+gridInfoRead.nSinkPoints;
  unsigned int i_u;
  int j;
  struct linkType *linkPtr;

  for(i_u=0;i_u<totalNumGridPoints;i_u++){
    gp[i_u].neigh = malloc(sizeof(struct grid *)*gp[i_u].numNeigh);

    for(j=0;j<gp[i_u].numNeigh;j++){
      linkPtr = nnLinks[firstNearNeigh[i_u]+j];
      if(linkPtr->gp[0]->id==(int)i_u){
        gp[i_u].neigh[j] = linkPtr->gp[1];
      }else{
        gp[i_u].neigh[j] = linkPtr->gp[0];
      }
    }
  }
}

/*....................................................................*/
void loadACoeffsIntoGrid(unsigned int *firstNearNeigh, struct linkType **nnLinks\
  , struct gridInfoType gridInfoRead, struct grid *gp){
  /*
See the comment at the beginning of the present module for a description of how the pointers 'links', 'nnLinks' and 'firstNearNeigh' relate to the grid struct.

The function mallocs the following extensions of struct g for each grid point:
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
      if(linkPtr->gp[0]->id==(int)i_u){
        evenCoeffSign =  1.0;
      }else{
        evenCoeffSign = -1.0;
      }

      gp[i_u].a0[j] = evenCoeffSign*(linkPtr->aCoeffs)[0];
      gp[i_u].a1[j] =               (linkPtr->aCoeffs)[1];
      gp[i_u].a2[j] = evenCoeffSign*(linkPtr->aCoeffs)[2];
      gp[i_u].a3[j] =               (linkPtr->aCoeffs)[3];
      gp[i_u].a4[j] = evenCoeffSign*(linkPtr->aCoeffs)[4];
    }
  }
}

/*....................................................................*/
int checkPopsTableExists(lime_fptr *fptr, const int fileFormatI\
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
int getNumPopsTables(lime_fptr *fptr, const int fileFormatI, unsigned short *numTables){
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
int readPopsTable(lime_fptr *fptr, const int fileFormatI\
  , const unsigned short speciesI, struct grid *gp\
  , struct gridInfoType *gridInfoRead){

  int status=0;
  unsigned int totalNumGridPoints, i_u;

  totalNumGridPoints = gridInfoRead->nInternalPoints + gridInfoRead->nSinkPoints;
  for(i_u=0;i_u<totalNumGridPoints;i_u++){
    gp[i_u].mol[speciesI].dust    = NULL;
    gp[i_u].mol[speciesI].knu     = NULL;
    gp[i_u].mol[speciesI].pops    = NULL;
    gp[i_u].mol[speciesI].partner = NULL;
  }

  if(fileFormatI==lime_FITS){
    readPopsExtFromFits(fptr, speciesI, gp, gridInfoRead);
  }else{
    status = 1;
  }

  return status;
}

/*....................................................................*/
int readGrid(char *inFileName, const int fileFormatI\
  , struct gridInfoType *gridInfoRead, struct grid **gp\
  , char ***collPartNames, int *numCollPartRead, int *dataFlags){

  /*
This is designed to be a generic function to read the grid data from file. It is assumed that the data will be stored in several tables of different size, corresponding to the different dimensionalities of the elements of the 'grid' struct. See 'writeGrid' for a description.

Some sanity checks are performed here and also in the deeper functions, but any check of the validity of the state of completeness of the grid data (as encoded in the returned argument dataFlags) is left to the calling routine.

NOTE that gp should not be allocated before this routine is called.
  */

  lime_fptr *fptr;
  int status=0;
  unsigned short i_s, numTables;
  unsigned int *firstNearNeigh=NULL, totalNumGridPoints, i_u;
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
  gridInfoRead->nACoeffs = 0;
  gridInfoRead->mols = NULL;

  /* Open the file and also return the data stage. */
  fptr = openFileForRead(inFileName, fileFormatI);

  /* Read the values which should be in grid for every stage.
  */
  status = readGridTable(fptr, fileFormatI, gridInfoRead, gp, &firstNearNeigh\
    , collPartNames, numCollPartRead, dataFlags); /* Sets appropriate bits of dataFlags; also mallocs gp and sets all its defaults. */
  totalNumGridPoints = gridInfoRead->nSinkPoints + gridInfoRead->nInternalPoints;
  if(status){
    closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, 0);
    freeGrid(totalNumGridPoints, gridInfoRead->nSpecies, *gp);
    return 1;
  }

  /* Some sanity checks:
  */
  if((*dataFlags)!=0 && gridInfoRead->nSinkPoints<=0 || gridInfoRead->nInternalPoints<=0 || gridInfoRead->nDims<=0){
    closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, 0);
    freeGrid(totalNumGridPoints, gridInfoRead->nSpecies, *gp);

    if(gridInfoRead->nSinkPoints<=0)
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

  /* Some more sanity checks:
  */
  if((allBitsSet(*dataFlags, DS_mask_density)   && gridInfoRead->nDensities<=0)\
  || (allBitsSet(*dataFlags, DS_mask_abundance) && gridInfoRead->nSpecies<=0)){
    closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, 0);
    freeGrid(totalNumGridPoints, gridInfoRead->nSpecies, *gp);

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

  status = readLinksTable(fptr, fileFormatI, gridInfoRead, *gp, &links, dataFlags); /* Sets appropriate bits of dataFlags. */
  if(status){
    closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, gridInfoRead->nLinks);
    freeGrid(totalNumGridPoints, gridInfoRead->nSpecies, *gp);
    return 7;
  }

  /* Sanity check:
  */
  if (allBitsSet(*dataFlags, DS_mask_ACOEFF) && gridInfoRead->nACoeffs<=0){
    closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, gridInfoRead->nLinks);
    freeGrid(totalNumGridPoints, gridInfoRead->nSpecies, *gp);
    return 8;
  }

  status = readNnIndicesTable(fptr, fileFormatI, links, &nnLinks, gridInfoRead, dataFlags); /* Sets appropriate bits of dataFlags. */
  if(status){
    closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, gridInfoRead->nLinks);
    freeGrid(totalNumGridPoints, gridInfoRead->nSpecies, *gp);
    return 9;
  }

  if(allBitsSet(*dataFlags, DS_mask_neighbours)){
    /* Convert the NN information back to the standard LIME grid struct format.
    */
    loadNnIntoGrid(firstNearNeigh, nnLinks, *gridInfoRead, *gp); /* mallocs extension 'neigh' of struct g for each grid point. */

    if(allBitsSet(*dataFlags, DS_mask_ACOEFF))
      loadACoeffsIntoGrid(firstNearNeigh, nnLinks, *gridInfoRead, *gp);
  }

  status = getNumPopsTables(fptr, fileFormatI, &numTables);
  if(status){
    closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, gridInfoRead->nLinks);
    freeGrid(totalNumGridPoints, gridInfoRead->nSpecies, *gp);
    return 10;
  }

  /* Sanity check:
  */
  if(numTables>0 && numTables!=gridInfoRead->nSpecies){
    closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, gridInfoRead->nLinks);
    freeGrid(totalNumGridPoints, gridInfoRead->nSpecies, *gp);
    return 11;
  }

  if(numTables>0){
    (*dataFlags) |= (1 << DS_bit_populations);

    for(i_u=0;i_u<totalNumGridPoints;i_u++)
      (*gp)[i_u].mol = malloc(sizeof(struct populations)*gridInfoRead->nSpecies);

    gridInfoRead->mols = malloc(sizeof(struct molInfoType)*gridInfoRead->nSpecies);
    for(i_s=0;i_s<gridInfoRead->nSpecies;i_s++){
      gridInfoRead->mols[i_s].molName = NULL;
      gridInfoRead->mols[i_s].nLevels = -1;
      gridInfoRead->mols[i_s].nLines = -1;

      status = readPopsTable(fptr, fileFormatI, i_s, *gp, gridInfoRead); /* Sets defaults for all the fields under grid.mol. */
      if(status){
        closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, gridInfoRead->nLinks);
        freeGrid(totalNumGridPoints, gridInfoRead->nSpecies, *gp);
        return 12;
      }

      /* Sanity check:
      */
      if(gridInfoRead->mols[i_s].nLevels<=0){
        closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, gridInfoRead->nLinks);
        freeGrid(totalNumGridPoints, gridInfoRead->nSpecies, *gp);
        return 13;
      }
    }
  }

  closeAndFree(fptr, fileFormatI, firstNearNeigh, nnLinks, links, gridInfoRead->nLinks);
  return 0;
}



