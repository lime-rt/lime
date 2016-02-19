/*
 *  grid2fits.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015 The LIME development team
 *
TODOS (some day):
	- Change the type of g[i].id to unsigned int.
	- Change the type of g[i].numNeigh to unsigned short.
	- Change g[i].t, g[i].dopb and g[i].abun to float.
 */

#include "lime.h"


/*
The present module contains routines for transferring the LIME grid point data to or from a FITS format file. The purpose of the present comment block is to describe the FITS file format.

1) Different amounts of information in grid.
--------------------------------------------
Describing the format of the FITS file is complicated by the fact that the amount of information stored in the grid data structure tends to increase during the run of the LIME program. One can distingish 4 stages of completion which are at least conceptually different, even if they are not cleanly separated in the present implementation of the code. These are summarized as follows.

	A) At this stage the vector of grid objects has been malloc'd and values have been generated for the following struct elements:
		id
		x
		sink
	This stage is flagged to the software by the element 'neigh' of the first grid point being set to NULL.

	B) This stage is entered after the Delaunay neighbours of each grid point have been determined. The following further struct elements are expected to have been malloc'd (in the case of pointers) and given values:
		numNeigh
		neigh
	This stage is flagged by the element 'dens' of the first grid point being set to NULL.

	C) This stage is entered after sampling the user-supplied functions for density, velocity etc. The following further struct elements are expected to have been malloc'd (in the case of pointers) and given values:
		vel
		a0, a1 etc.
		dens
		t
		dopb
		abun
	This stage is flagged by the element 'mol' of the first grid point being set to NULL.

	D) After 1 or more iterations of populating the levels, we enter stage D, in which all the grid struct elements have been malloc'd and given at least preliminary values; specifically now the element
		mol

(Other grid struct values are deducible from the ones given or are used for temporary storage only.)

The dependency relationship tree for the stages may be diagrammed as follows:

	A -> B -> C -> D.

2) The FITS file format for stage D
-----------------------------------
This is defined in the following list of extensions. All extensions are binary table except where indicated. The letter in the second row for column descriptions gives the FITS data type. See eg

  https://heasarc.gsfc.nasa.gov/docs/software/fitsio/c/c_user/node20.html

for a key to these.

	1) GRID
	Number of rows = number of grid points.

	Keywords:
		RADIUS
		COLLPARn		# 1 for each nth collision partner.

	Columns:
		ID		V
		Xj		D	# Cartesian components of the point location, 1 col per jth dimension.
		IS_SINK		L	# =True iff the point lies on the edge of the model.
		NUMNEIGH	U
		FIRST_NN	V	# See explanation in section 3 below.
		VELj		D	# 1 col per jth dimension.
		TURBDPLR	E	# Given Gaussian lineshape exp(-v^2/[B^2 + 2*k*T/m]), this is B.
		DENSITYn	D	# 1 per nth collision partner.
		TEMPKNTC	E	# From t[0].
		TEMPDUST	E	# From t[1].
		ABUNMOLm	E	# 1 per mth molecular species.

	2) NN_INDICES (see explanation in section 3 below)
	Number of rows = number of grid points * average number of Delaunay links per point.

	Columns:
		LINK_I		V

	3) LINKS (see explanation in section 3 below)
	Number of rows = number of Delaunay links.

	Columns:
		GRID_I_1	V
		GRID_I_2	V
		ACOEFF_p	D	# 1 per pth order of the velocity polynomial.

	4 etc) LEVEL_POPS_m (1 per mth molecular species)
	This is an image extension of size (number of grid cells)*(number of energy levels this species).

	Keywords:
		MOL_NAME

The subsets of these FITS structures which are written under the other 3 stages are described in the header comment to the function writeGridToFits().

3) The NN_INDICES and LINKS extensions and the FIRST_NN column.
---------------------------------------------------------------
The 'neigh' element of struct grid is a clunky way to store the information about the nearest neighbours of each grid point, because it leads to a host of separate little mallocs, thus memory for the grid is spread all over the place. There is also some duplication of quantities which refer to the link between 2 grid points rather than to the points themselves. This arrangement does not lend itself easily to storage in a file.

To meet the needs of FITS storage we divide the information contained in 'neigh' between 3 new vectors 'links', 'nnLinks' and 'firstNearNeigh'. The 1st of these is a list of all Delaunay edges between grid points. The 2nd contains indices to the 1st, and is arranged such that indices to all the links connected to a single grid point occur in a connected sequence. The 3rd vector has the same number of entries as the total number of grid points. It contains the index of the first entry in the 2nd vector which corresponds to that grid point. Thus in order, for example, to iterate over all links connected to the ith grid point, we need only do something as follows:

  for(j=0;j<g[i].numNeigh;j++){
    k = j + firstNearNeigh[i];
    link = links[nnLinks[k]];
  }

In addition, any link-specific stuff (at present just the 'a' coefficients used to interpolate the velocity along the length of the link) is stored in the first vector 'links'.

Since 'firstNearNeigh' has a single entry for each grid point, it can be stored as a column FIRST_NN in the GRID extension. The 'links' and 'nnLinks' vectors must be stored in separate extensions LINKS and NN_INDICES.
*/

/*....................................................................*/
int writeGridToFits(char *outFileName, inputPars par, unsigned short numDims\
  , unsigned short numACoeffs, struct grid *g, molData *md, char **collPartNames){

  /*
The present function writes information contained in the grid variable to FITS file. The extensions/keywords/columns written depend on the state of the data in grid, as described in the header comment to the present module. This is described below via pseudo-code.

The shorthand used below for the FITS info is
#		Ext.
#			Keyword
#			--------
#			Column

The pseudo-code is as follows:

  if (g!=NULL){ # = at least stage A, write:
#		GRID
#			RADIUS
#			--------
#			ID
#			Xj
#			IS_SINK

    if (g[0].neigh!=NULL){ # = at least stage B, write:
#		GRID
#			--------
#			NUMNEIGH
#			FIRST_NN
#
#		NN_INDICES
#			--------
#			LINK_I
#
#		LINKS
#			--------
#			GRID_I_0
#			GRID_I_1

      if (g[0].dens!=NULL){ # = at least stage C, write:
#		GRID
#			COLLPARn
#			--------
#			VELj
#			TURBDPLR
#			DENSITYn
#			TEMPKNTC
#			TEMPDUST
#			ABUNMOLm
#
#		LINKS
#			--------
#			ACOEFF_p

        if (g[0].mol!=NULL){ # = stage D, write:
#		LEVEL_POPS_m
#			MOL_NAME
#			--------
        }
      }
    }
  }
  */

  fitsfile *fptr;
  int status = 0;
  unsigned short speciesI;
  char negfile[100]="! ", message[80];
  struct linkType *links=NULL, **nnLinks=NULL;
  unsigned int totalNumLinks, totalNumNeigh, *firstNearNeigh=NULL;
  _Bool stageC=0;

  if (outFileName==""){
    if(!silent) warning("Cannot write grid list to file, filename is blank.");
    return 1;
  }

  if (g==NULL){
    if(!silent) warning("Cannot write grid list to file, there are no entries in it.");
    return 2;
  }

  sprintf(message, "Writing grid list to file %s", outFileName);
  if(!silent) printMessage(message);

  /* If we get to here, we have at least stage A info in the grid. */

  if (g[0].neigh!=NULL){ /* => at least stage B. */
    constructLinkArrays((unsigned int)par.ncell, g, &links, &totalNumLinks\
      , &nnLinks, &firstNearNeigh, &totalNumNeigh);
  }

  fits_create_file(&fptr, outFileName, &status);
  if(status!=0){
    if(!silent) warning("Overwriting existing fits file");
    status=0;
    strcat(negfile,outFileName);
    fits_create_file(&fptr, negfile, &status);
    processFitsError(status);
  }

  /* GRID
  */
  writeGridExtToFits(fptr, par, numDims, g, firstNearNeigh, collPartNames); /* Deals internally with stages B and C. */

  if (g[0].neigh!=NULL){ /* => at least stage B. */
    /* NN_INDICES
    */
    writeNnIndicesExtToFits(fptr, totalNumNeigh, nnLinks, links);

    /* LINKS
    */
    if(g[0].dens!=NULL) stageC = 1;
    writeLinksExtToFits(fptr, stageC, totalNumLinks, numACoeffs, links);

    if (g[0].mol!=NULL){ /* => stage D. */
      for(speciesI=0;speciesI<par.nSpecies;speciesI++){
        /* LEVEL_POPS_m
        */
        writePopsExtToFits(fptr, (unsigned int)par.ncell, md, speciesI, g);
      }
    }
  }

  fits_close_file(fptr, &status);
  processFitsError(status);

  free(firstNearNeigh);
  free(nnLinks);
  free(links);

  return 0;
}

/*....................................................................*/
void constructLinkArrays(unsigned int numGridPoints, struct grid *g\
  , struct linkType **links, unsigned int *totalNumLinks, struct linkType ***nnLinks\
  , unsigned int **firstNearNeigh, unsigned int *totalNumNeigh){
  /*
See the comment at the beginning of the present module for a description of how the pointers 'links', 'nnLinks' and 'firstNearNeigh' relate to the grid struct.
  */

  unsigned int li, ni, idA, idB, trialIdA, linkId;
  int i, jA, jB, k;
  _Bool *pointIsDone=NULL, linkNotFound;
  struct grid *gAPtr=NULL, *gBPtr=NULL;
  int *nnLinkIs=NULL;
  char message[80];

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

  if(g[0].dens!=NULL){ /* => at least stage C. */
    for(li=0;li<*totalNumLinks;li++){
      gAPtr = (*links)[li].g[0];
      /* Find which neighbour of gAPtr corresponds to the link: */
      linkNotFound = 1;
      for(jA=0;jA<gAPtr->numNeigh;jA++){
        idB = gAPtr->neigh[jA]->id;
        if(linkNotFound && idB==(*links)[li].g[1]->id){
          (*links)[li].aCoeffs[0] = gAPtr->a0[jA];
          (*links)[li].aCoeffs[1] = gAPtr->a1[jA];
          (*links)[li].aCoeffs[2] = gAPtr->a2[jA];
          (*links)[li].aCoeffs[3] = gAPtr->a3[jA];
          (*links)[li].aCoeffs[4] = gAPtr->a4[jA];
          linkNotFound = 0;
        }
      }
      if(linkNotFound){
        sprintf(message, "SNAFU with link %d grid point %d", li, gAPtr->id);
        if(!silent) bail_out(message);
        exit(1);
      }
    }
  }

  *nnLinks = malloc(sizeof(**nnLinks)*(*totalNumNeigh));
  for(ni=0;ni<(*totalNumNeigh);ni++)
    (*nnLinks)[ni] = &(*links)[nnLinkIs[ni]];

  free(nnLinkIs);
  free(pointIsDone);
  /* The calling routine must free firstNearNeigh, nnLinks, links. */
}

/*....................................................................*/
void defineGridExtColumns(unsigned short numKwdChars, unsigned short numDims\
  , unsigned short numCollPart, unsigned short numSpecies, char *ttype[]\
  , char *tform[], char *tunit[], int dataTypes[]){

  /* Note that data types in all capitals are defined in fitsio.h. */

  const int maxNumDims=9, maxNumSpecies=9, maxNumCollPart=9;
  int colI=0, i;
  char message[80];

  if(numDims>maxNumDims){
    if(!silent){
      sprintf(message, "Caller asked for %d dims but colnames can only be written for %d.", numDims, maxNumDims);
      bail_out(message);
    }
    exit(1);
  }

  if(numSpecies>maxNumSpecies){
    if(!silent){
      sprintf(message, "Caller asked for %d species but colnames can only be written for %d.", numSpecies, maxNumSpecies);
      bail_out(message);
    }
    exit(1);
  }

  if(numCollPart>maxNumCollPart){
    if(!silent){
      sprintf(message, "Caller asked for %d coll. part. but colnames can only be written for %d.", numCollPart, maxNumCollPart);
      bail_out(message);
    }
    exit(1);
  }

  colI = 0;
  ttype[colI] = malloc(numKwdChars);
  sprintf(ttype[colI], "ID");
  tform[colI] = "V";
  tunit[colI] = "\0";
  dataTypes[colI] = TUINT;

  /* should rather have a vector column? */
  for(i=0;i<numDims;i++){
    colI++;
    ttype[colI] = malloc(numKwdChars);
    sprintf(ttype[colI], "X%d", i+1);
    tform[colI] = "D";
    tunit[colI] = "m";
    dataTypes[colI] = TDOUBLE;
  }

  colI++;
  ttype[colI] = malloc(numKwdChars);
  sprintf(ttype[colI], "IS_SINK");
  tform[colI] = "L";
  tunit[colI] = "\0";
  dataTypes[colI] = TLOGICAL;

  colI++;
  ttype[colI] = malloc(numKwdChars);
  sprintf(ttype[colI], "NUMNEIGH");
  tform[colI] = "U";
  tunit[colI] = "\0";
  dataTypes[colI] = TUSHORT;

  colI++;
  ttype[colI] = malloc(numKwdChars);
  sprintf(ttype[colI], "FIRST_NN");
  tform[colI] = "V";
  tunit[colI] = "\0";
  dataTypes[colI] = TUINT;

  /* should rather have a vector column? */
  for(i=0;i<numDims;i++){
    colI++;
    ttype[colI] = malloc(numKwdChars);
    sprintf(ttype[colI], "VEL%d", i+1);
    tform[colI] = "D";
    tunit[colI] = "m/s";
    dataTypes[colI] = TDOUBLE;
  }

  /* should rather have a vector column? */
  for(i=0;i<numCollPart;i++){
    colI++;
    ttype[colI] = malloc(numKwdChars);
    sprintf(ttype[colI], "DENSITY%d", i+1);
    tform[colI] = "D";
    tunit[colI] = "kg/m^3";
    dataTypes[colI] = TDOUBLE;
  }

  colI++;
  ttype[colI] = malloc(numKwdChars);
  sprintf(ttype[colI], "TURBDPLR");
  tform[colI] = "E";
  tunit[colI] = "m/s";
  dataTypes[colI] = TFLOAT;

  colI++;
  ttype[colI] = malloc(numKwdChars);
  sprintf(ttype[colI], "TEMPKNTC");
  tform[colI] = "E";
  tunit[colI] = "K";
  dataTypes[colI] = TFLOAT;

  colI++;
  ttype[colI] = malloc(numKwdChars);
  sprintf(ttype[colI], "TEMPDUST");
  tform[colI] = "E";
  tunit[colI] = "K";
  dataTypes[colI] = TFLOAT;

  /* should rather have a vector column? */
  for(i=0;i<numSpecies;i++){
    colI++;
    ttype[colI] = malloc(numKwdChars);
    sprintf(ttype[colI], "ABUNMOL%d", i+1);
    tform[colI] = "E";
    tunit[colI] = "\0";
    dataTypes[colI] = TFLOAT;
  }
}

/*....................................................................*/
void writeGridExtToFits(fitsfile *fptr, inputPars par, unsigned short numDims\
  , struct grid *g, unsigned int *firstNearNeigh, char **collPartNames){

  /* Note that data types in all capitals are defined in fitsio.h. */

  unsigned int *ids=NULL;
  double *xj=NULL;
  _Bool *sink=NULL;/* Probably should be char* but this seems to work. */
  unsigned short *numNeigh=NULL;
  double *velj=NULL, *densn=NULL;
  float *dopb=NULL, *t=NULL, *abunm=NULL;
  const unsigned short numKwdChars = 9; /* 8 characters + \0. */
  const int maxNumCollPart = 9;
  int status=0, colI=0, i, j, m, n, localNumCollPart;
  LONGLONG firstRow=1, firstElem=1;
  int numCols = 7 + numDims*2 + par.collPart + par.nSpecies;
  char extname[] = "GRID";
  char genericComment[80];
  char genericKwd[numKwdChars], message[80];

  /* Define the name, datatype, and physical units for the columns.
  */
  char *ttype[numCols];
  char *tform[numCols];
  char *tunit[numCols];
  int dataTypes[numCols];

  defineGridExtColumns(numKwdChars, numDims, (unsigned short)par.collPart\
    ,(unsigned short)par.nSpecies, ttype, tform, tunit, dataTypes);

  /* Append a new empty binary table onto the FITS file.
  */
  fits_create_tbl( fptr, BINARY_TBL, 0, numCols, ttype, tform, tunit, extname, &status);
  processFitsError(status);

  for(colI=0;colI<numCols;colI++)
    free(ttype[colI]);

  /* Write the columns:
  */
  colI = 0;
  ids = malloc(sizeof(*ids)*par.ncell);
  for(i=0;i<par.ncell;i++) {
    ids[i] = (unsigned int)g[i].id;
  }
  fits_write_col(fptr, dataTypes[colI], colI+1, firstRow, firstElem, (LONGLONG)par.ncell, ids, &status);
  processFitsError(status);
  free(ids);

  xj = malloc(sizeof(*xj)*par.ncell);
  for(j=0;j<numDims;j++){
    colI++;
    for(i=0;i<par.ncell;i++) xj[i] = g[i].x[j];
    fits_write_col(fptr, dataTypes[colI], colI+1, firstRow, firstElem, (LONGLONG)par.ncell, xj, &status);
    processFitsError(status);
  }
  free(xj);

  colI++;
  sink = malloc(sizeof(*sink)*par.ncell);
  for(i=0;i<par.ncell;i++) sink[i] = (_Bool)g[i].sink;
  fits_write_col(fptr, dataTypes[colI], colI+1, firstRow, firstElem, (LONGLONG)par.ncell, sink, &status);
  processFitsError(status);
  free(sink);

  if (firstNearNeigh!=NULL){ /* => at least stage B. */
    colI++;
    numNeigh = malloc(sizeof(*numNeigh)*par.ncell);
    for(i=0;i<par.ncell;i++) numNeigh[i] = (unsigned short)g[i].numNeigh;
    fits_write_col(fptr, dataTypes[colI], colI+1, firstRow, firstElem, (LONGLONG)par.ncell, numNeigh, &status);
    processFitsError(status);
    free(numNeigh);

    colI++;
    fits_write_col(fptr, dataTypes[colI], colI+1, firstRow, firstElem, (LONGLONG)par.ncell, firstNearNeigh, &status);
    processFitsError(status);

    if (g[0].dens!=NULL){ /* => at least stage C. */
      velj = malloc(sizeof(*velj)*par.ncell);
      for(j=0;j<numDims;j++){
        colI++;
        for(i=0;i<par.ncell;i++) velj[i] = g[i].vel[j];
        fits_write_col(fptr, dataTypes[colI], colI+1, firstRow, firstElem, (LONGLONG)par.ncell, velj, &status);
        processFitsError(status);
      }
      free(velj);

      densn = malloc(sizeof(*densn)*par.ncell);
      for(n=0;n<par.collPart;n++){
        colI++;
        for(i=0;i<par.ncell;i++) densn[i] = g[i].dens[n];
        fits_write_col(fptr, dataTypes[colI], colI+1, firstRow, firstElem, (LONGLONG)par.ncell, densn, &status);
        processFitsError(status);
      }
      free(densn);

      colI++;
      dopb = malloc(sizeof(*dopb)*par.ncell);
      for(i=0;i<par.ncell;i++) dopb[i] = (float)g[i].dopb;
      fits_write_col(fptr, dataTypes[colI], colI+1, firstRow, firstElem, (LONGLONG)par.ncell, dopb, &status);
      processFitsError(status);
      free(dopb);

      t = malloc(sizeof(*t)*par.ncell);
      colI++;
      for(i=0;i<par.ncell;i++) t[i] = (float)g[i].t[0];
      fits_write_col(fptr, dataTypes[colI], colI+1, firstRow, firstElem, (LONGLONG)par.ncell, t, &status);
      processFitsError(status);

      colI++;
      for(i=0;i<par.ncell;i++) t[i] = (float)g[i].t[1];
      fits_write_col(fptr, dataTypes[colI], colI+1, firstRow, firstElem, (LONGLONG)par.ncell, t, &status);
      processFitsError(status);
      free(t);

      abunm = malloc(sizeof(*abunm)*par.ncell);
      for(m=0;m<par.nSpecies;m++){
        colI++;
        for(i=0;i<par.ncell;i++) abunm[i] = (float)g[i].abun[m];
        fits_write_col(fptr, dataTypes[colI], colI+1, firstRow, firstElem, (LONGLONG)par.ncell, abunm, &status);
        processFitsError(status);
      }
      free(abunm);
    }
  }

  /* Write keywords.
  */
  fits_write_key(fptr, TDOUBLE, "RADIUS  ", &(par.radius), "[m] Model radius.", &status);
  processFitsError(status);

  if (collPartNames!=NULL){
    if(par.collPart>maxNumCollPart){
      if(!silent){
        sprintf(message, "There seem to be %d collision partners but keywords can only be written for %d.", par.collPart, maxNumCollPart);
        warning(message);
      }
      localNumCollPart = maxNumCollPart;
    }else{
      localNumCollPart = par.collPart;
    }

    for(i=0;i<localNumCollPart;i++){
      sprintf(genericKwd, "COLLPAR%d", i+1);
      sprintf(genericComment, "Collision partner %d", i+1);
      fits_write_key(fptr, TSTRING, genericKwd, collPartNames[i], genericComment, &status);
      processFitsError(status);
    }
  }
}

/*....................................................................*/
void writeNnIndicesExtToFits(fitsfile *fptr, unsigned int totalNumNeigh\
  , struct linkType **nnLinks, struct linkType *links){
  /*
See the comment at the beginning of the present module for a description of how the NN_INDICES extension relates to the grid struct.

	Extension name: NN_INDICES
	Number of rows = number of grid points * average number of Delaunay links per point.

	Columns:
		LINK_I		V

Note that data types in all capitals are defined in fitsio.h.
  */

  unsigned int *linkIs=NULL;
  int status=0, i;
  LONGLONG firstRow=1, firstElem=1;
  int numCols = 1;
  char extname[] = "NN_INDICES";

  /* Define the name, datatype, and physical units for the columns.
  */
  char *ttype[] = { "LINK_I" };
  char *tform[] = { "V"      };
  char *tunit[] = { "\0"     };
  int dataType = TUINT;

  /* Append a new empty binary table onto the FITS file.
  */
  fits_create_tbl( fptr, BINARY_TBL, 0, numCols, ttype, tform, tunit, extname, &status);
  processFitsError(status);

  linkIs = malloc(sizeof(*linkIs)*totalNumNeigh);
  for(i=0;i<totalNumNeigh;i++)
    linkIs[i] = nnLinks[i]->id;
  fits_write_col(fptr, dataType, 1, firstRow, firstElem, (LONGLONG)totalNumNeigh, linkIs, &status);
  processFitsError(status);
  free(linkIs);
}

/*....................................................................*/
void writeLinksExtToFits(fitsfile *fptr, _Bool stageC, unsigned int totalNumLinks\
  , unsigned short numACoeffs, struct linkType *links){
  /*
See the comment at the beginning of the present module for a description of how the LINKS extension relates to the grid struct.

Note that data types in all capitals are defined in fitsio.h.
  */

  unsigned int *ids=NULL;
  double *aCoeffs=NULL;
  int status=0, colI=0, i, n;
  LONGLONG firstRow=1, firstElem=1;
  int numCols = 2 + numACoeffs;
  char extname[] = "LINKS";
  const int numKwdChars = 9; /* 8 characters + \0. */

  /* Define the name, datatype, and physical units for the columns.
  */
  char *ttype[numCols];
  char *tform[numCols];
  char *tunit[numCols];
  int dataTypes[numCols];

  colI = 0;
  ttype[colI] = malloc(numKwdChars);
  sprintf(ttype[colI], "GRID_I_1");
  tform[colI] = "V";
  tunit[colI] = "\0";
  dataTypes[colI] = TUINT;

  colI++;
  ttype[colI] = malloc(numKwdChars);
  sprintf(ttype[colI], "GRID_I_2");
  tform[colI] = "V";
  tunit[colI] = "\0";
  dataTypes[colI] = TUINT;

  /* should rather have a vector column? */
  for(n=0;n<numACoeffs;n++){
    colI++;
    ttype[colI] = malloc(numKwdChars);
    sprintf(ttype[colI], "ACOEFF_%d", n+1);
    tform[colI] = "D";
    tunit[colI] = "\0";
    dataTypes[colI] = TDOUBLE;
  }

  /* Append a new empty binary table onto the FITS file.
  */
  fits_create_tbl(fptr, BINARY_TBL, 0, numCols, ttype, tform, tunit, extname, &status);
  processFitsError(status);

  for(colI=0;colI<numCols;colI++)
    free(ttype[colI]);

  /* Write columns.
  */
  colI = 0;
  ids = malloc(sizeof(*ids)*totalNumLinks);
  for(i=0;i<totalNumLinks;i++)
    ids[i] = (unsigned int)links[i].g[0]->id;
  fits_write_col(fptr, dataTypes[colI], colI+1, firstRow, firstElem, (LONGLONG)totalNumLinks, ids, &status);
  processFitsError(status);

  colI++;
  for(i=0;i<totalNumLinks;i++)
    ids[i] = (unsigned int)links[i].g[1]->id;
  fits_write_col(fptr, dataTypes[colI], colI+1, firstRow, firstElem, (LONGLONG)totalNumLinks, ids, &status);
  processFitsError(status);
  free(ids);

  if (stageC){
    aCoeffs = malloc(sizeof(*aCoeffs)*totalNumLinks);
    for(n=0;n<numACoeffs;n++){
      colI++;
      for(i=0;i<totalNumLinks;i++) aCoeffs[i] = links[i].aCoeffs[n];
      fits_write_col(fptr, dataTypes[colI], colI+1, firstRow, firstElem, (LONGLONG)totalNumLinks, aCoeffs, &status);
      processFitsError(status);
    }
    free(aCoeffs);
  }
}

/*....................................................................*/
void writePopsExtToFits(fitsfile *fptr, unsigned int numGridPoints\
  , molData *md, unsigned short speciesI, struct grid *g){
  /*
See the comment at the beginning of the present module for a description of how the LEVEL_POPS_m extensions relate to the grid struct.

	Extension name: LEVEL_POPS_m (1 per molecular species)
	This is an image extension of size (number of grid cells)*(number of energy levels this species)

	Keywords:
		MOL_NAME

Note that data types in all capitals are defined in fitsio.h.
  */

  int status=0, xi, yi;
  char extname[13];
  float *row=NULL;
  int bitpix = FLOAT_IMG;
  const long naxis = 2;  /* i.e. 2-dimensional image */    
  unsigned int numEnergyLevels = md[speciesI].nlev;
  long naxes[] = { numEnergyLevels, numGridPoints };
  long fpixels[naxis],lpixels[naxis];

  sprintf(extname, "LEVEL_POPS_%d", (int)speciesI+1);

  fits_create_img(fptr, bitpix, naxis, naxes, &status);
  processFitsError(status);

  row = malloc(sizeof(*row)*numEnergyLevels);

  /* Write FITS data.
  */
  for(yi=0;yi<numGridPoints;yi++){
    for(xi=0;xi<numEnergyLevels;xi++)
      row[xi] = (float)g[yi].mol[speciesI].pops[xi]; 

    fpixels[0]=1;
    fpixels[1]=yi+1;
    lpixels[0]=numEnergyLevels;
    lpixels[1]=yi+1;

    fits_write_subset(fptr, TFLOAT, fpixels, lpixels, row, &status);
    processFitsError(status);
  }

  free(row);

  /* write keywords:
  */
  fits_write_key(fptr, TSTRING, "MOL_NAME ", md[speciesI].molName, "\0", &status);
  processFitsError(status);

  fits_write_key(fptr, TSTRING, "EXTNAME ", extname, "\0", &status);
  processFitsError(status);
}

/*....................................................................*/
void readGridFromFits(char *inFileName, inputPars par, unsigned short numDims\
  , unsigned short numACoeffs, struct grid **g, _Bool *userValuesFound\
  , _Bool *nnValuesFound, _Bool *popsValuesFound, molData *md\
  , char ***collPartNames, int *numCollPartRead){
  /*
See the comment at the beginning of the present module for information about the expected structure of the FITS file and the categorization of the amount of information present in 4 stages from A to D.

Briefly, the stage appropriate to the file is assessed according to the following criteria:

	- If the routine is called at all, it is assumed that the stage is at least A.

	- If a NUMNEIGH column is found in the GRID extension, the stage is at least B.

	- If a VEL1 column is found in the GRID extension, the stage is at least C.

	- If a LEVEL_POPS_1 extension is found, the stage is D.
  */

  fitsfile *fptr;
  int status=0, i;
  char extname[13];
  unsigned short si=0;
  unsigned int *firstNearNeigh=NULL;
  struct linkType *links=NULL, **nnLinks=NULL;

  fits_open_file(&fptr, inFileName, READONLY, &status);
  processFitsError(status);

  /* Go to the GRID extension, then read it.
  */
  fits_movnam_hdu(fptr, BINARY_TBL, "GRID", 0, &status);
  processFitsError(status);

  readGridExtFromFits(fptr, par, numDims, g, &firstNearNeigh\
    , userValuesFound, collPartNames, numCollPartRead);

  if(firstNearNeigh==NULL){ /* In this case we are only at stage A. */
    *nnValuesFound = 0;

    for(i=0;i<par.ncell;i++){
      (*g)[i].neigh = NULL;
      (*g)[i].a0    = NULL;
      (*g)[i].a1    = NULL;
      (*g)[i].a2    = NULL;
      (*g)[i].a3    = NULL;
      (*g)[i].a4    = NULL;
    }

  }else{ /* we are at least at stage B: there should be NN_INDICES and LINKS extensions. */
    /* Go to the LINKS extension, then read it.
    */
    fits_movnam_hdu(fptr, BINARY_TBL, "LINKS", 0, &status);
    processFitsError(status);

    readLinksExtFromFits(fptr, *userValuesFound, numACoeffs, *g, &links);

    /* Go to the NN_INDICES extension, then read it.
    */
    fits_movnam_hdu(fptr, BINARY_TBL, "NN_INDICES", 0, &status);
    processFitsError(status);

    readNnIndicesExtFromFits(fptr, links, &nnLinks);
    *nnValuesFound = 1;

    /* Convert the NN information back to the standard LIME grid struct format.
    */
    loadNnIntoGrid(firstNearNeigh, nnLinks, links, (unsigned int)par.ncell, g);

    free(nnLinks);
    free(links);
    free(firstNearNeigh);
  }

  *popsValuesFound = 0; /* Default value. */
  if (userValuesFound){ /* we are at least at stage C. */
    /*
See if there are LEVEL_POPS_m extensions. The number of such extensions should either equal zero or equal to par.nSpecies.

Try to move to the first LEVEL_POPS extension:
    */
    si = 0;
    sprintf(extname, "LEVEL_POPS_%d", (int)si+1);
    fits_movnam_hdu(fptr, IMAGE_HDU, extname, 0, &status);
    if(status!=BAD_HDU_NUM){ /* Success, at least partially. That is, ok we found the extension, but now we check for other errors: */
      processFitsError(status);
      *popsValuesFound = 1;

      for(i=0;i<par.ncell;i++)
        (*g)[i].mol = malloc(sizeof(struct populations)*par.nSpecies);

      for(si==0;si<par.nSpecies;si++){
        sprintf(extname, "LEVEL_POPS_%d", (int)si+1);

        /* Try to move to extension LEVEL_POPS_<si>:
        */
        fits_movnam_hdu(fptr, IMAGE_HDU, extname, 0, &status);
        processFitsError(status);

        readPopsExtFromFits(fptr, (unsigned int)par.ncell, md, si, g);
      }
    }
  }

  if(*popsValuesFound==0){
    for(i=0;i<par.ncell;i++)
      (*g)[i].mol = NULL;
  }

  fits_close_file(fptr, &status);
  processFitsError(status);
}

/*....................................................................*/
void readGridExtFromFits(fitsfile *fptr, inputPars par, unsigned short numDims\
  , struct grid **g, unsigned int **firstNearNeigh, _Bool *userValuesFound\
  , char ***collPartNames, int *numCollPartRead){
  /*
The present function mallocs 'g'. If (nnValuesFound) then it also mallocs 'firstNearNeigh'.

If (*userValuesFound) then the function also mallocs the following elements of struct grid for all grid points:
	dens
	abun
Otherwise these are set to NULL.

If a COLLPARn keywords are found in the GRID extension header then collPartNames is malloc'd to the number of these.
  */

  LONGLONG numGridCells, firstRow=1, firstElem=1;
  int status=0, colNum, anynul=0, i, j, n;//, localNumCollPart;
  char colName[20];
  _Bool nnValuesFound;
  double modelRadius;
  char *genericComment;
  char genericKwd[9], message[80];
  const int maxNumCollPart = 9;
  char dummyCollPartName[80];
  unsigned int *ids=NULL;
  double *xj=NULL;
  _Bool *sink=NULL;/* Probably should be char* but this seems to work. */
  unsigned short *numNeigh=NULL;
  double *velj=NULL, *densn=NULL;
  float *dopb=NULL, *t=NULL, *abunm=NULL;

  /* Find out how many rows there are, then malloc the array.
  */
  fits_get_num_rowsll(fptr, &numGridCells, &status);
  processFitsError(status);

  if(numGridCells!=par.ncell && !silent){
    sprintf(message, "The inputPars struct says there should be %d grid points, but %u were found.\n", par.ncell, (unsigned int)numGridCells);
    warning(message);
  }

  *g = malloc(sizeof(**g)*numGridCells);

  /* Read the columns.
  */
  fits_get_colnum(fptr, CASEINSEN, "ID", &colNum, &status);
  processFitsError(status);

  ids = malloc(sizeof(*ids)*numGridCells);
  fits_read_col(fptr, TUINT, colNum, firstRow, firstElem, numGridCells, 0, ids, &anynul, &status);
  processFitsError(status);

  for(i=0;i<numGridCells;i++) {
    (*g)[i].id = (int)ids[i];
  }
  free(ids);

  xj = malloc(sizeof(*xj)*numGridCells);
  for(i=0;i<numDims;i++){
    sprintf(colName, "X%d", i+1);
    fits_get_colnum(fptr, CASEINSEN, colName, &colNum, &status);
    processFitsError(status);

    fits_read_col(fptr, TDOUBLE, colNum, firstRow, firstElem, numGridCells, 0, xj, &anynul, &status);
    processFitsError(status);

    for(j=0;j<numGridCells;j++) {
      (*g)[j].x[i] = xj[j];
    }
  }
  free(xj);

  fits_get_colnum(fptr, CASEINSEN, "IS_SINK", &colNum, &status);
  processFitsError(status);

  sink = malloc(sizeof(*sink)*numGridCells);
  fits_read_col(fptr, TLOGICAL, colNum, firstRow, firstElem, numGridCells, 0, sink, &anynul, &status);
  processFitsError(status);

  for(i=0;i<numGridCells;i++) {
    (*g)[i].sink = (int)sink[i];
  }
  free(sink);

  /* Check for the existence of a NUMNEIGH column:
  */
  nnValuesFound = 0; /* Default value. */
  fits_get_colnum(fptr, CASEINSEN, "NUMNEIGH", &colNum, &status);
  if(status!=COL_NOT_FOUND){
    processFitsError(status);

    nnValuesFound = 1;
    numNeigh = malloc(sizeof(*numNeigh)*numGridCells);
    fits_read_col(fptr, TUSHORT, colNum, firstRow, firstElem, numGridCells, 0, numNeigh, &anynul, &status);
    processFitsError(status);

    for(i=0;i<numGridCells;i++) {
      (*g)[i].numNeigh = (int)numNeigh[i];
    }
    free(numNeigh);

    *firstNearNeigh = malloc(sizeof(**firstNearNeigh)*numGridCells); /* This not being NULL will signal that these columns were found. */

    fits_get_colnum(fptr, CASEINSEN, "FIRST_NN", &colNum, &status);
    processFitsError(status);

    fits_read_col(fptr, TUINT, colNum, firstRow, firstElem, numGridCells, 0, *firstNearNeigh, &anynul, &status);
    processFitsError(status);

    /* Check for the existence of a VEL1 column:
    */
    *userValuesFound = 0; /* Default value. */
    i = 0;
    sprintf(colName, "VEL%d", i+1);
    fits_get_colnum(fptr, CASEINSEN, colName, &colNum, &status);
    if(status!=COL_NOT_FOUND){
      processFitsError(status);

      *userValuesFound = 1;

      /* Dimension the appropriate elements of g:
      */
      for(i=0;i<numGridCells;i++) {
        (*g)[i].dens = malloc(sizeof(double)*par.collPart);
        (*g)[i].abun = malloc(sizeof(double)*par.nSpecies);
      }

      /* Read the VEL columns:
      */
      velj = malloc(sizeof(*velj)*numGridCells);
      for(i=0;i<numDims;i++){
        sprintf(colName, "VEL%d", i+1);
        fits_get_colnum(fptr, CASEINSEN, colName, &colNum, &status);
        processFitsError(status);

        fits_read_col(fptr, TDOUBLE, colNum, firstRow, firstElem, numGridCells, 0, velj, &anynul, &status);
        processFitsError(status);

        for(j=0;j<numGridCells;j++) {
          (*g)[j].vel[i] = velj[j];
        }
      }
      free(velj);

      /* Read the TURBDPLR column:
      */
      fits_get_colnum(fptr, CASEINSEN, "TURBDPLR", &colNum, &status);
      processFitsError(status);

      dopb = malloc(sizeof(*dopb)*numGridCells);
      fits_read_col(fptr, TFLOAT, colNum, firstRow, firstElem, numGridCells, 0, dopb, &anynul, &status);
      processFitsError(status);

      for(i=0;i<numGridCells;i++) {
        (*g)[i].dopb = (double)dopb[i];
      }
      free(dopb);

      /* Read the DENSITY columns:
      */
      densn = malloc(sizeof(*densn)*numGridCells);
      for(n=0;n<par.collPart;n++){
        sprintf(colName, "DENSITY%d", n+1);
        fits_get_colnum(fptr, CASEINSEN, colName, &colNum, &status);
        processFitsError(status);

        fits_read_col(fptr, TDOUBLE, colNum, firstRow, firstElem, numGridCells, 0, densn, &anynul, &status);
        processFitsError(status);

        for(j=0;j<numGridCells;j++) {
          (*g)[j].dens[n] = densn[j];
        }
      }
      free(densn);

      /* Read the TEMPKNTC column:
      */
      fits_get_colnum(fptr, CASEINSEN, "TEMPKNTC", &colNum, &status);
      processFitsError(status);

      t = malloc(sizeof(*t)*numGridCells);
      fits_read_col(fptr, TFLOAT, colNum, firstRow, firstElem, numGridCells, 0, t, &anynul, &status);
      processFitsError(status);

      for(i=0;i<numGridCells;i++) {
        (*g)[i].t[0] = (double)t[i];
      }

      /* Read the TEMPDUST column:
      */
      fits_get_colnum(fptr, CASEINSEN, "TEMPDUST", &colNum, &status);
      processFitsError(status);

      fits_read_col(fptr, TFLOAT, colNum, firstRow, firstElem, numGridCells, 0, t, &anynul, &status);
      processFitsError(status);

      for(i=0;i<numGridCells;i++) {
        (*g)[i].t[1] = (double)t[i];
      }
      free(t);

      /* Read the ABUNMOL columns:
      */
      abunm = malloc(sizeof(*abunm)*numGridCells);
      for(i=0;i<par.nSpecies;i++){
        sprintf(colName, "ABUNMOL%d", i+1);
        fits_get_colnum(fptr, CASEINSEN, colName, &colNum, &status);
        processFitsError(status);

        fits_read_col(fptr, TFLOAT, colNum, firstRow, firstElem, numGridCells, 0, abunm, &anynul, &status);
        processFitsError(status);

        for(j=0;j<numGridCells;j++) {
          (*g)[j].abun[i] = (double)abunm[j];
        }
      }
      free(abunm);
    }else{ /* user values not found */
      /* Null the appropriate elements of g:
      */
      for(i=0;i<numGridCells;i++) {
        (*g)[i].dens = NULL;
        (*g)[i].abun = NULL;
      }
    } /* end if(*userValuesFound) block. */
  } /* end if(nnValuesFound) block. */


  /* Read kwds:
  */
  fits_read_key(fptr, TDOUBLE, "RADIUS  ", &modelRadius, NULL, &status);
  processFitsError(status);
  /************compare against par.radius??*/

  /* Check if there are any COLLPAR keywords.
  */
  *numCollPartRead = 0; /* Default. */
  i = 0;
  sprintf(genericKwd, "COLLPAR%d", i+1);
  fits_read_key(fptr, TSTRING, genericKwd, dummyCollPartName, NULL, &status);
  if (status!=KEY_NO_EXIST){
    processFitsError(status);

    if(par.collPart>maxNumCollPart){
      if(!silent){
        sprintf(message, "There seem to be %d collision partners but keywords can only be written for %d.", par.collPart, maxNumCollPart);
        warning(message);
      }
      *numCollPartRead = maxNumCollPart;
    }else{
      *numCollPartRead = par.collPart;
    }

    *collPartNames = malloc(sizeof(**collPartNames)*(*numCollPartRead));

    for(i=0;i<(*numCollPartRead);i++){
      sprintf(genericKwd, "COLLPAR%d", i+1);
      (*collPartNames)[i] = malloc(sizeof(char)*100);
      fits_read_key(fptr, TSTRING, genericKwd, (*collPartNames)[i], NULL, &status);
      processFitsError(status);
    }
  }
}

/*....................................................................*/
void readLinksExtFromFits(fitsfile *fptr, _Bool userValuesFound\
  , unsigned short numACoeffs, struct grid *g, struct linkType **links){
  /*
See the comment at the beginning of the present module for a description of how the LINKS extension relates to the grid struct.

The present function mallocs *links.
  */

  LONGLONG totalNumLinks, firstRow=1, firstElem=1;
  int status=0, colNum, anynul=0, li, n;
  char colName[21];
  unsigned int *ids=NULL;
  double *aCoeffs=NULL;

  /* Find out how many rows there are, then malloc the array.
  */
  fits_get_num_rowsll(fptr, &totalNumLinks, &status);
  processFitsError(status);

  *links = malloc(sizeof(**links)*totalNumLinks);

  /* Read GRID_I_1 column.
  */
  fits_get_colnum(fptr, CASEINSEN, "GRID_I_1", &colNum, &status);
  processFitsError(status);

  ids = malloc(sizeof(*ids)*totalNumLinks);
  fits_read_col(fptr, TUINT, colNum, firstRow, firstElem, totalNumLinks, 0, ids, &anynul, &status);
  processFitsError(status);

  for(li=0;li<totalNumLinks;li++) {
    (*links)[li].g[0] = &g[ids[li]];
  }

  /* Read GRID_I_2 column.
  */
  fits_get_colnum(fptr, CASEINSEN, "GRID_I_2", &colNum, &status);
  processFitsError(status);

  fits_read_col(fptr, TUINT, colNum, firstRow, firstElem, totalNumLinks, 0, ids, &anynul, &status);
  processFitsError(status);

  for(li=0;li<totalNumLinks;li++) {
    (*links)[li].g[1] = &g[ids[li]];
  }
  free(ids);

  if(userValuesFound>0){ /* Means we are in at least stage C. */
    aCoeffs = malloc(sizeof(*aCoeffs)*totalNumLinks);
    for(n=0;n<numACoeffs;n++){
      /* Read the ACOEFF_n columns.
      */
      sprintf(colName, "ACOEFF_%d", n+1);
      fits_get_colnum(fptr, CASEINSEN, colName, &colNum, &status);
      processFitsError(status);

      fits_read_col(fptr, TDOUBLE, colNum, firstRow, firstElem, totalNumLinks, 0, aCoeffs, &anynul, &status);
      processFitsError(status);

      for(li=0;li<totalNumLinks;li++) {
        (*links)[li].aCoeffs[n] = aCoeffs[li];
      }
    }
    free(aCoeffs);

  }else{
    for(li=0;li<totalNumLinks;li++) {
      for(n=0;n<numACoeffs;n++){
        (*links)[li].aCoeffs[n] = 0.0;
      }
    }
  } 
}

/*....................................................................*/
void readNnIndicesExtFromFits(fitsfile *fptr, struct linkType *links\
  , struct linkType ***nnLinks){
  /*
See the comment at the beginning of the present module for a description of how the NN_INDICES extension relates to the grid struct.

The function mallocs *nnLinks.
  */

  LONGLONG totalNumNeigh, firstRow=1, firstElem=1;
  int status=0, colNum, anynul=0, i;
  unsigned int *linkIs=NULL;

  /* Find out how many rows there are, then malloc the array.
  */
  fits_get_num_rowsll(fptr, &totalNumNeigh, &status);
  processFitsError(status);

  *nnLinks = malloc(sizeof(**nnLinks)*totalNumNeigh);

  /* Read LINK_I column.
  */
  fits_get_colnum(fptr, CASEINSEN, "LINK_I", &colNum, &status);
  processFitsError(status);

  linkIs = malloc(sizeof(*linkIs)*totalNumNeigh);
  fits_read_col(fptr, TUINT, colNum, firstRow, firstElem, totalNumNeigh, 0, linkIs, &anynul, &status);
  processFitsError(status);

  for(i=0;i<totalNumNeigh;i++) {
    (*nnLinks)[i] = &links[linkIs[i]];
  }
  free(linkIs);
}

/*....................................................................*/
void loadNnIntoGrid(unsigned int *firstNearNeigh, struct linkType **nnLinks\
  , struct linkType *links, unsigned int numGridCells, struct grid **g){
  /*
See the comment at the beginning of the present module for a description of how the pointers 'links', 'nnLinks' and 'firstNearNeigh' relate to the grid struct.

The function mallocs the following extensions of struct g for each grid point:
	neigh
	a0
	a1
	a2
	a3
	a4
  */

  unsigned int i;
  unsigned short j;
  struct linkType *linkPtr;

  for(i=0;i<numGridCells;i++){
    (*g)[i].neigh = malloc(sizeof(struct grid *)*(*g)[i].numNeigh);
    (*g)[i].a0    = malloc(sizeof(double)       *(*g)[i].numNeigh);
    (*g)[i].a1    = malloc(sizeof(double)       *(*g)[i].numNeigh);
    (*g)[i].a2    = malloc(sizeof(double)       *(*g)[i].numNeigh);
    (*g)[i].a3    = malloc(sizeof(double)       *(*g)[i].numNeigh);
    (*g)[i].a4    = malloc(sizeof(double)       *(*g)[i].numNeigh);
    for(j=0;j<((*g)[i].numNeigh);j++){
      linkPtr = nnLinks[firstNearNeigh[i]+j];
      if(linkPtr->g[0]->id==i){
        (*g)[i].neigh[j] = linkPtr->g[1];
      }else{
        (*g)[i].neigh[j] = linkPtr->g[0];
      }

      (*g)[i].a0[j] = linkPtr->aCoeffs[0];
      (*g)[i].a1[j] = linkPtr->aCoeffs[1];
      (*g)[i].a2[j] = linkPtr->aCoeffs[2];
      (*g)[i].a3[j] = linkPtr->aCoeffs[3];
      (*g)[i].a4[j] = linkPtr->aCoeffs[4];
    }
  }
}

/*....................................................................*/
void readPopsExtFromFits(fitsfile *fptr, unsigned int numGridPoints\
  , molData *md, unsigned short speciesI, struct grid **g){
  /*
See the comment at the beginning of the present module for a description of how the LEVEL_POPS_m extensions relate to the grid struct.

The function mallocs (*g)[yi].mol[speciesI].pops for each grid point yi and species.
  */

  const int maxLenMolName = 8;
  int bitpix, naxis, status=0, xi, yi, anynul=0;
  long *naxes=NULL;
  float *row=NULL;
  long fpixels[2],lpixels[2];
  long inc[2] = {1,1};
  char molNameRead[maxLenMolName+1];
  unsigned int numEnergyLevels = md[speciesI].nlev;

  fits_get_img_param(fptr, 2, &bitpix, &naxis, naxes, &status);
  processFitsError(status);
  /****** error if naxes != {numEnergyLevels, numGridPoints} or bitpix!=FLOAT_IMG.*/
  free(naxes);

  row = malloc(sizeof(*row)*numEnergyLevels);

  /* Read FITS data.
  */
  for(yi=0;yi<numGridPoints;yi++){
    fpixels[0]=1;
    fpixels[1]=yi+1;
    lpixels[0]=numEnergyLevels;
    lpixels[1]=yi+1;

    fits_read_subset(fptr, TFLOAT, fpixels, lpixels, inc, 0, row, &anynul, &status);
    processFitsError(status);

    (*g)[yi].mol[speciesI].pops = malloc(sizeof(double)*numEnergyLevels);
    for(xi=0;xi<numEnergyLevels;xi++)
      (*g)[yi].mol[speciesI].pops[xi] = (double)row[xi];
  }

  free(row);

  /* Read kwds:
  */
  fits_read_key(fptr, TSTRING, "MOL_NAME", molNameRead, NULL, &status);
  processFitsError(status);
  /****** raise exception if it is not the same as md[speciesI].molName. */
}



