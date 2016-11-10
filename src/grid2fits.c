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
	- Vector columns in defineGridExtColumns()?
*/

#include "lime.h"
#include "gridio.h"

/*
The present module contains routines for transferring the LIME grid point data to or from a FITS format file. The purpose of the present comment block is to describe the FITS file format. Note that the amount and type of information stored depends on the 'data stage' of the grid struct, as described in the header remarks to module gridio.c.

In the description below, the data stage bit associated with the presence of a particular extension, column or keyword is given on the leftmost place of each line. When the file is read, all the objects associated with a bit must be present for the bit to be set.

Note that all extensions are binary table except where indicated. The letter in the second row for column descriptions gives the FITS data type. See eg

  https://heasarc.gsfc.nasa.gov/docs/software/fitsio/c/c_user/node20.html

for a key to these.

Where column names contain a lower-case letter, this is a placeholder for a digit as explained in the repective comment.

0	0) The primary HDU

	Keywords:
0		RADIUS		D	# The model radius in metres.

0	1) GRID
	Number of rows = number of grid points.

	Keywords:
0		COLLPARn	A	# 1 for each nth collision partner.

	Columns:
0		ID		V
0		Xj		D	# Cartesian components of the point location, 1 col per jth dimension.
0		IS_SINK		L	# =True iff the point lies on the edge of the model.
1		NUMNEIGH	U
1		FIRST_NN	V	# See explanation in section 3 below.
2		VELj		D	# 1 col per jth dimension.
3		DENSITYn	D	# 1 per nth collision partner.
4		ABUNMOLm	E	# 1 per mth molecular species.
5		TURBDPLR	E	# Given Gaussian lineshape exp(-v^2/[B^2 + 2*k*T/m]), this is B.
6		TEMPKNTC	E	# From t[0].
6		TEMPDUST	E	# From t[1].

1	2) NN_INDICES (see explanation in the header to gridio.c)
	Number of rows = number of grid points * average number of Delaunay links per point.

	Columns:
1		LINK_I		V

1	3) LINKS (see explanation in the header to gridio.c)
	Number of rows = number of Delaunay links.

	Columns:
1		GRID_I_1	V
1		GRID_I_2	V
7		V_p_j		D	# 1 per pth velocity sample per jth dimension.

8	4 etc) LEVEL_POPS_m (1 per mth molecular species)
	This is an image extension of size (number of grid cells)*(number of energy levels this species).

	Keywords:
8		MOL_NAME

Note that at present the data in the 'partner' element of grid.mol is *NOT* being stored.
*/

/*....................................................................*/
fitsfile *
openFITSFileForWrite(char *outFileName){
  fitsfile *fptr=NULL;
  int status=0;
  char negfile[100]="! ";

  fits_create_file(&fptr, outFileName, &status);
  if(status!=0){
    if(!silent) warning("Overwriting existing fits file");
    status=0;
    strcat(negfile,outFileName);
    fits_create_file(&fptr, negfile, &status);
    processFitsError(status);
  }

  fits_create_img(fptr, 8, 0, NULL, &status);
  processFitsError(status);

  return fptr;
}

/*....................................................................*/
void
initializeKeyword(struct keywordType *kwd){
  (*kwd).datatype = 0;
  (*kwd).keyname = NULL;
  (*kwd).comment = NULL;
  (*kwd).intValue = 0;
  (*kwd).floatValue = 0.0;
  (*kwd).doubleValue = 0.0;
  (*kwd).charValue = NULL;
}

/*....................................................................*/
void
closeFITSFile(fitsfile *fptr){
  int status=0;

  fits_close_file(fptr, &status);
  processFitsError(status);
}

/*....................................................................*/
void
writeKeywordsToFits(lime_fptr *fptr, struct keywordType *kwds\
  , const int numKeywords){

  int i, status;
  char message[80];

  for(i=0;i<numKeywords;i++){
    status = 0;

    if(     kwds[i].datatype==TSTRING)
      fits_write_key(fptr, TSTRING, kwds[i].keyname, kwds[i].charValue, kwds[i].comment, &status);
    else if(kwds[i].datatype==TINT)
      fits_write_key(fptr, TINT, kwds[i].keyname, &kwds[i].intValue, kwds[i].comment, &status);
    else if(kwds[i].datatype==TFLOAT)
      fits_write_key(fptr, TFLOAT, kwds[i].keyname, &kwds[i].floatValue, kwds[i].comment, &status);
    else if(kwds[i].datatype==TDOUBLE)
      fits_write_key(fptr, TDOUBLE, kwds[i].keyname, &kwds[i].doubleValue, kwds[i].comment, &status);
    else{
      if(!silent){
        sprintf(message, "Keyword %d dataype %d is not currently accepted.", i, kwds[i].datatype);
        bail_out(message);
      }
      exit(1);
    }
    processFitsError(status);
  }
}

/*....................................................................*/
void
defineAndLoadColumns(fitsfile *fptr, struct gridInfoType gridInfo\
  , const int dataFlags, const unsigned short numColNameChars, char ***allColNames, int **allColNumbers\
  , int *maxNumCols, int **colDataTypes){
  /*
To understand what is going on in the present function one needs to be aware that there is a set of all possible columns, corresponding to the set of quasi-scalar elements of the grid struct; and also a subset of that (the 'columns to write'), which is the columns for which data currently exists in the grid struct. The latter set is defined by the bits set in dataFlags.

The function does two things:
  - Sets up the binary table extension GRID with dimensions corresponding to the available data in the grid struct, as encoded in dataFlags.
  - Returns 3 vectors: allColNames, allColNumbers and colDataTypes, described as follows:

      allColNames: pretty self-explanatory.

      allColNumbers: this is a vector of size equal to the total number of possible columns, but its values are indices in the range {1,...,numColsToWrite}. If one of the possible columns is not scheduled to be written (because the appropriate bit of dataFlags was not set), then the matching value of allColNumbers will be 0. The number of columns scheduled to be written is therefore the number of non-zero elements of allColNumbers.

      colDataTypes: this contains data types, but NOT for all columns, just for those to be written.

NOTES:
  - The calling routine needs to free allColNames, allColNumbers and colDataTypes after it is finshed with them.
  - Data types in all capitals are defined in fitsio.h.
  */

  const unsigned short maxNumDims=9, maxNumSpecies=9, maxNumDensities=9;
  unsigned short i_us;
  int colI,i,colToWriteI,status=0;
  char message[80];
  char **tformAllCols=NULL;
  char **tunitAllCols=NULL;
  int *dataTypeAllCols=NULL;

  if(gridInfo.nDims>maxNumDims){
    if(!silent){
      sprintf(message, "Caller asked for %d dims but colnames can only be written for %d.", (int)gridInfo.nDims, (int)maxNumDims);
      bail_out(message);
    }
    exit(1);
  }

  if(gridInfo.nSpecies>maxNumSpecies){
    if(!silent){
      sprintf(message, "Caller asked for %d species but colnames can only be written for %d.", (int)gridInfo.nSpecies, (int)maxNumSpecies);
      bail_out(message);
    }
    exit(1);
  }

  if(gridInfo.nDensities>maxNumDensities){
    if(!silent){
      sprintf(message, "Caller asked for %d coll. part. but colnames can only be written for %d.", (int)gridInfo.nDensities, (int)maxNumDensities);
      bail_out(message);
    }
    exit(1);
  }

  *maxNumCols = 10 + gridInfo.nDims*2 + gridInfo.nSpecies + gridInfo.nDensities;

  *allColNames    = malloc(sizeof(**allColNames)   *(*maxNumCols));
  *allColNumbers  = malloc(sizeof(**allColNumbers) *(*maxNumCols));
  tformAllCols    = malloc(sizeof(*tformAllCols)   *(*maxNumCols));
  tunitAllCols    = malloc(sizeof(*tunitAllCols)   *(*maxNumCols));
  dataTypeAllCols = malloc(sizeof(*dataTypeAllCols)*(*maxNumCols));

  for(i=0;i<(*maxNumCols);i++){ /* Set default: */
    (*allColNames)[i] = malloc(sizeof(char)*numColNameChars);
    (*allColNumbers)[i] = 0;
  }

  colToWriteI = 0;
  colI = 0;

  sprintf((*allColNames)[colI], "ID");
  if(bitIsSet(dataFlags, DS_bit_x)){
    colToWriteI++;
    (*allColNumbers)[colI] = colToWriteI;
  }
  tformAllCols[colI] = "V";
  tunitAllCols[colI] = "\0";
  dataTypeAllCols[colI] = TUINT;

  /* should rather have a vector column? */
  for(i_us=0;i_us<gridInfo.nDims;i_us++){
    colI++;
    sprintf((*allColNames)[colI], "X%d", (int)i_us+1);
    if(bitIsSet(dataFlags, DS_bit_x)){
      colToWriteI++;
      (*allColNumbers)[colI] = colToWriteI;
    }
    tformAllCols[colI] = "D";
    tunitAllCols[colI] = "m";
    dataTypeAllCols[colI] = TDOUBLE;
  }

  colI++;
  sprintf((*allColNames)[colI], "IS_SINK");
  if(bitIsSet(dataFlags, DS_bit_x)){
    colToWriteI++;
    (*allColNumbers)[colI] = colToWriteI;
  }
  tformAllCols[colI] = "L";
  tunitAllCols[colI] = "\0";
  dataTypeAllCols[colI] = TLOGICAL;

  colI++;
  sprintf((*allColNames)[colI], "NUMNEIGH");
  if(bitIsSet(dataFlags, DS_bit_neighbours)){
    colToWriteI++;
    (*allColNumbers)[colI] = colToWriteI;
  }
  tformAllCols[colI] = "U";
  tunitAllCols[colI] = "\0";
  dataTypeAllCols[colI] = TUSHORT;

  colI++;
  sprintf((*allColNames)[colI], "FIRST_NN");
  if(bitIsSet(dataFlags, DS_bit_neighbours)){
    colToWriteI++;
    (*allColNumbers)[colI] = colToWriteI;
  }
  tformAllCols[colI] = "V";
  tunitAllCols[colI] = "\0";
  dataTypeAllCols[colI] = TUINT;

  /* should rather have a vector column? */
  for(i_us=0;i_us<gridInfo.nDims;i_us++){
    colI++;
    sprintf((*allColNames)[colI], "VEL%d", (int)i_us+1);
    if(bitIsSet(dataFlags, DS_bit_velocity)){
      colToWriteI++;
      (*allColNumbers)[colI] = colToWriteI;
    }
    tformAllCols[colI] = "D";
    tunitAllCols[colI] = "m/s";
    dataTypeAllCols[colI] = TDOUBLE;
  }

  /* should rather have a vector column? */
  for(i_us=0;i_us<gridInfo.nDensities;i_us++){
    colI++;
    sprintf((*allColNames)[colI], "DENSITY%d", (int)i_us+1);
    if(bitIsSet(dataFlags, DS_bit_density)){
      colToWriteI++;
      (*allColNumbers)[colI] = colToWriteI;
    }
    tformAllCols[colI] = "D";
    tunitAllCols[colI] = "kg/m^3";
    dataTypeAllCols[colI] = TDOUBLE;
  }

  /* should rather have a vector column? */
  for(i_us=0;i_us<gridInfo.nSpecies;i_us++){
    colI++;
    sprintf((*allColNames)[colI], "ABUNMOL%d", (int)i_us+1);
    if(bitIsSet(dataFlags, DS_bit_abundance)){
      colToWriteI++;
      (*allColNumbers)[colI] = colToWriteI;
    }
    tformAllCols[colI] = "E";
    tunitAllCols[colI] = "\0";
    dataTypeAllCols[colI] = TFLOAT;
  }

  colI++;
  sprintf((*allColNames)[colI], "TURBDPLR");
  if(bitIsSet(dataFlags, DS_bit_turb_doppler)){
    colToWriteI++;
    (*allColNumbers)[colI] = colToWriteI;
  }
  tformAllCols[colI] = "E";
  tunitAllCols[colI] = "m/s";
  dataTypeAllCols[colI] = TFLOAT;

  colI++;
  sprintf((*allColNames)[colI], "TEMPKNTC");
  if(bitIsSet(dataFlags, DS_bit_temperatures)){
    colToWriteI++;
    (*allColNumbers)[colI] = colToWriteI;
  }
  tformAllCols[colI] = "E";
  tunitAllCols[colI] = "K";
  dataTypeAllCols[colI] = TFLOAT;

  colI++;
  sprintf((*allColNames)[colI], "TEMPDUST");
  if(bitIsSet(dataFlags, DS_bit_temperatures)){
    colToWriteI++;
    (*allColNumbers)[colI] = colToWriteI;
  }
  tformAllCols[colI] = "E";
  tunitAllCols[colI] = "K";
  dataTypeAllCols[colI] = TFLOAT;

  /* should rather have a vector column? */
  for(i=0;i<3;i++){/* **** should rather loop to gridInfo.nDims but only entre here if it ==3?? */
    colI++;
    sprintf((*allColNames)[colI], "B_FIELD%d", i+1);
    if(bitIsSet(dataFlags, DS_bit_magfield)){
      colToWriteI++;
      (*allColNumbers)[colI] = colToWriteI;
    }
    tformAllCols[colI] = "E";
    tunitAllCols[colI] = "T";
    dataTypeAllCols[colI] = TFLOAT;
  }

  /* Define the name, format, datatype, and physical units for the columns which will actually be written.
  */
  char *ttype[colToWriteI];
  char *tform[colToWriteI];
  char *tunit[colToWriteI];
  *colDataTypes = malloc(sizeof(**colDataTypes)*colToWriteI);

  for(colI=0;colI<(*maxNumCols);colI++){
    i = (*allColNumbers)[colI];
    if(i>0){
      ttype[i-1] = (*allColNames)[colI];
      tform[i-1] = tformAllCols[colI];
      tunit[i-1] = tunitAllCols[colI];
      (*colDataTypes)[i-1] = dataTypeAllCols[colI];
    }
  }

  /* Append a new empty binary table onto the FITS file.
  */
  fits_create_tbl(fptr, BINARY_TBL, 0, colToWriteI, ttype, tform, tunit, "GRID", &status);
  processFitsError(status);

  free(tformAllCols);
  free(tunitAllCols);
  free(dataTypeAllCols);
}

/*....................................................................*/
int
getColIndex(char **allColNames, const int maxNumCols, char *colName){
  int i=0;
  int colFound=0; /* -> bool */

  while(i<maxNumCols && !colFound){
    if(!strcmp(allColNames[i],colName))
      colFound = 1;

    i++;
  }

  if(!colFound){
    if(!silent) bail_out("Column name not found - this should not happen.");
    exit(1);
  }

  return i-1;
}

/*....................................................................*/
void
writeGridExtToFits(fitsfile *fptr, struct gridInfoType gridInfo\
  , struct grid *gp, unsigned int *firstNearNeigh\
  , char **collPartNames, const int dataFlags){
  /*
This writes whatever information is in the grid struct (as specified by the dataFlags) and which also has a dimensionality which is some simple multiple of the number of grid points, to a single FITS binary table extension called GRID. The function tries to be fairly forgiving of screwy situations but it will exit if the minimum information is not present (defined as allBitsSet(dataFlags, DS_mask_x), which implies that elements .id, .x and .sink should all contain valid values).

Note that data types in all capitals are defined in fitsio.h.
  */

  const unsigned int totalNumGridPoints = gridInfo.nInternalPoints+gridInfo.nSinkPoints;
  const unsigned short numKwdChars=9; /* 8 characters + \0. */
  const unsigned short numColNameChars=21; /* 20 characters + \0. */
  const unsigned short maxNumCollPart = 9;
  unsigned int *ids=NULL,i_ui;
  double *xj=NULL;
  _Bool *sink=NULL;/* Probably should be char* but this seems to work. */
  unsigned short *numNeigh=NULL,i_us,localNumCollPart;
  double *velj=NULL,*densn=NULL;
  float *dopb=NULL, *t=NULL, *abunm=NULL, *bField=NULL;
  int status=0, colI=0, i, di, maxNumCols;
  LONGLONG firstRow=1, firstElem=1;
  char genericComment[80];
  char genericKwd[numKwdChars], message[80];
  char colName[numColNameChars];
  char **allColNames=NULL;
  int *allColNumbers=NULL, *colDataTypes=NULL;

  if(!allBitsSet(dataFlags, DS_mask_x)){
    if(!silent) bail_out("Data stage indicates no grid data!");
    exit(1);
  }

  /*
Ok we have a bit of a tricky situation here in that the number of columns we write is going to depend on the information available in gp, as encoded in the dataFlags. We need to work out which columns we are going to write ahead of time because we need the appropriate data on ALL the columns to set up the table size before we can start to write their individual data values. I also want to avoid checking dataFlags twice in two different contexts - that is how errors arise. So I've arranged that the following routine will do all the donkey work of setting up only those columns we can write and then using that information to define the table size. The routine also returns three more vectors:

	allColNames   - This is returned to remove what would otherwise be a hard-wired dependence that the column ordering was the same in the present routine as in defineAndLoadColumns(). With this vector, the present routine can search for a column name in it and then use the returned vector index to access (from the next vector) the number of the column in the (smaller) sequence of valid columns.

	allColNumbers - This contains the number of a column in the sequence (beginning at 1) of those for which gp has data. If gp contains no data for a given column name, its entry in allColNumbers will be 0.

	colDataTypes  - This just contains the data types for the valid columns.
  */
  defineAndLoadColumns(fptr, gridInfo\
    , dataFlags, numColNameChars, &allColNames, &allColNumbers, &maxNumCols, &colDataTypes);

  /* Write the columns:
  */
  colI = allColNumbers[getColIndex(allColNames, maxNumCols, "ID")];
  if(colI<=0){
    if(!silent) bail_out("This should not occur, it is some sort of bug.");
    exit(1);
  }
  ids = malloc(sizeof(*ids)*totalNumGridPoints);
  for(i_ui=0;i_ui<totalNumGridPoints;i_ui++)
    ids[i_ui] = (unsigned int)gp[i_ui].id;
  fits_write_col(fptr, colDataTypes[colI-1], colI, firstRow, firstElem, (LONGLONG)totalNumGridPoints, ids, &status);
  processFitsError(status);
  free(ids);

  xj = malloc(sizeof(*xj)*totalNumGridPoints);
  for(i_us=0;i_us<gridInfo.nDims;i_us++){
    sprintf(colName, "X%d", (int)i_us+1);
    colI = allColNumbers[getColIndex(allColNames, maxNumCols, colName)];
    if(colI<=0){
      if(!silent) bail_out("This should not occur, it is some sort of bug.");
      exit(1);
    }

    for(i_ui=0;i_ui<totalNumGridPoints;i_ui++)
      xj[i_ui] = gp[i_ui].x[i_us];
    fits_write_col(fptr, colDataTypes[colI-1], colI, firstRow, firstElem, (LONGLONG)totalNumGridPoints, xj, &status);
    processFitsError(status);
  }
  free(xj);

  colI = allColNumbers[getColIndex(allColNames, maxNumCols, "IS_SINK")];
  if(colI<=0){
    if(!silent) bail_out("This should not occur, it is some sort of bug.");
    exit(1);
  }
  sink = malloc(sizeof(*sink)*totalNumGridPoints);
  for(i_ui=0;i_ui<totalNumGridPoints;i_ui++)
    sink[i_ui] = (_Bool)gp[i_ui].sink;
  fits_write_col(fptr, colDataTypes[colI-1], colI, firstRow, firstElem, (LONGLONG)totalNumGridPoints, sink, &status);
  processFitsError(status);
  free(sink);

  colI = allColNumbers[getColIndex(allColNames, maxNumCols, "NUMNEIGH")];
  if(colI>0){
    numNeigh = malloc(sizeof(*numNeigh)*totalNumGridPoints);
    for(i_ui=0;i_ui<totalNumGridPoints;i_ui++)
      numNeigh[i_ui] = (unsigned short)gp[i_ui].numNeigh;
    fits_write_col(fptr, colDataTypes[colI-1], colI, firstRow, firstElem, (LONGLONG)totalNumGridPoints, numNeigh, &status);
    processFitsError(status);
    free(numNeigh);

    colI = allColNumbers[getColIndex(allColNames, maxNumCols, "FIRST_NN")];
    if(colI<=0){
      if(!silent) bail_out("This should not occur, it is some sort of bug.");
      exit(1);
    }
    fits_write_col(fptr, colDataTypes[colI-1], colI, firstRow, firstElem, (LONGLONG)totalNumGridPoints, firstNearNeigh, &status);
    processFitsError(status);
  }

  /* Check if first VEL column has info:
  */
  colI = allColNumbers[getColIndex(allColNames, maxNumCols, "VEL1")];
  if(colI>0){
    velj = malloc(sizeof(*velj)*totalNumGridPoints);
    for(i_us=0;i_us<gridInfo.nDims;i_us++){
      sprintf(colName, "VEL%d", (int)i_us+1);
      colI = allColNumbers[getColIndex(allColNames, maxNumCols, colName)];
      if(colI<=0){
        if(!silent) bail_out("This should not occur, it is some sort of bug.");
        exit(1);
      }

      for(i_ui=0;i_ui<totalNumGridPoints;i_ui++)
        velj[i_ui] = gp[i_ui].vel[i_us];
      fits_write_col(fptr, colDataTypes[colI-1], colI, firstRow, firstElem, (LONGLONG)totalNumGridPoints, velj, &status);
      processFitsError(status);
    }
    free(velj);
  }

  /* Check if first DENSITY column has info:
  */
  colI = allColNumbers[getColIndex(allColNames, maxNumCols, "DENSITY1")];
  if(colI>0){
    densn = malloc(sizeof(*densn)*totalNumGridPoints);
    for(i_us=0;i_us<gridInfo.nDensities;i_us++){
      sprintf(colName, "DENSITY%d", (int)i_us+1);
      colI = allColNumbers[getColIndex(allColNames, maxNumCols, colName)];
      if(colI<=0){
        if(!silent) bail_out("This should not occur, it is some sort of bug.");
        exit(1);
      }

      for(i_ui=0;i_ui<totalNumGridPoints;i_ui++)
        densn[i_ui] = gp[i_ui].dens[i_us];
      fits_write_col(fptr, colDataTypes[colI-1], colI, firstRow, firstElem, (LONGLONG)totalNumGridPoints, densn, &status);
      processFitsError(status);
    }
    free(densn);
  }

  /* Check if first ABUNMOL column has info:
  */
  colI = allColNumbers[getColIndex(allColNames, maxNumCols, "ABUNMOL1")];
  if(colI>0){
    abunm = malloc(sizeof(*abunm)*totalNumGridPoints);
    for(i_us=0;i_us<gridInfo.nSpecies;i_us++){
      sprintf(colName, "ABUNMOL%d", (int)i_us+1);
      colI = allColNumbers[getColIndex(allColNames, maxNumCols, colName)];
      if(colI<=0){
        if(!silent) bail_out("This should not occur, it is some sort of bug.");
        exit(1);
      }

      for(i_ui=0;i_ui<totalNumGridPoints;i_ui++)
        abunm[i_ui] = (float)gp[i_ui].abun[i_us];
      fits_write_col(fptr, colDataTypes[colI-1], colI, firstRow, firstElem, (LONGLONG)totalNumGridPoints, abunm, &status);
      processFitsError(status);
    }
    free(abunm);
  }

  colI = allColNumbers[getColIndex(allColNames, maxNumCols, "TURBDPLR")];
  if(colI>0){
    dopb = malloc(sizeof(*dopb)*totalNumGridPoints);
    for(i_ui=0;i_ui<totalNumGridPoints;i_ui++)
      dopb[i_ui] = (float)gp[i_ui].dopb_turb;
    fits_write_col(fptr, colDataTypes[colI-1], colI, firstRow, firstElem, (LONGLONG)totalNumGridPoints, dopb, &status);
    processFitsError(status);
    free(dopb);
  }

  colI = allColNumbers[getColIndex(allColNames, maxNumCols, "TEMPKNTC")];
  if(colI>0){
    t = malloc(sizeof(*t)*totalNumGridPoints);
    for(i_ui=0;i_ui<totalNumGridPoints;i_ui++)
      t[i_ui] = (float)gp[i_ui].t[0];
    fits_write_col(fptr, colDataTypes[colI-1], colI, firstRow, firstElem, (LONGLONG)totalNumGridPoints, t, &status);
    processFitsError(status);

    colI = allColNumbers[getColIndex(allColNames, maxNumCols, "TEMPDUST")];
    if(colI<=0){
      if(!silent) bail_out("This should not occur, it is some sort of bug.");
      exit(1);
    }

    for(i_ui=0;i_ui<totalNumGridPoints;i_ui++)
      t[i_ui] = (float)gp[i_ui].t[1];
    fits_write_col(fptr, colDataTypes[colI-1], colI, firstRow, firstElem, (LONGLONG)totalNumGridPoints, t, &status);
    processFitsError(status);
    free(t);
  }

  /* Check if first B_FIELD column has info:
  */
  colI = allColNumbers[getColIndex(allColNames, maxNumCols, "B_FIELD1")];
  if(colI>0){
    bField = malloc(sizeof(*bField)*totalNumGridPoints);
    for(di=0;di<3;di++){//**** keep this hard-wired or rather test that gridInfo.nDims==3??
      sprintf(colName, "B_FIELD%d", di+1);
      colI = allColNumbers[getColIndex(allColNames, maxNumCols, colName)];
      if(colI<=0){
        if(!silent) bail_out("This should not occur, it is some sort of bug.");
        exit(1);
      }

      for(i_ui=0;i_ui<totalNumGridPoints;i_ui++)
        bField[i_ui] = (float)gp[i_ui].B[di];
      fits_write_col(fptr, colDataTypes[colI-1], colI, firstRow, firstElem, (LONGLONG)totalNumGridPoints, bField, &status);
      processFitsError(status);
    }
    free(bField);
  }

  /* Write keywords.
  */
  if (collPartNames!=NULL){
    if(gridInfo.nDensities>maxNumCollPart){
      if(!silent){
        sprintf(message, "There seem to be %d collision partners but keywords can only be written for %d.", (int)gridInfo.nDensities, (int)maxNumCollPart);
        warning(message);
      }
      localNumCollPart = maxNumCollPart;
    }else{
      localNumCollPart = gridInfo.nDensities;
    }

    for(i_us=0;i_us<localNumCollPart;i_us++){
      sprintf(genericKwd, "COLLPAR%d", (int)i_us+1);
      sprintf(genericComment, "Collision partner %d", (int)i_us+1);
      fits_write_key(fptr, TSTRING, genericKwd, collPartNames[i_us], genericComment, &status);
      processFitsError(status);
    }
  }

  for(i=0;i<maxNumCols;i++) free(allColNames[i]);
  free(allColNames);
  free(allColNumbers);
  free(colDataTypes);
}

/*....................................................................*/
void
writeNnIndicesExtToFits(fitsfile *fptr, struct gridInfoType gridInfo\
  , struct linkType **nnLinks){
  /*
See the comment at the beginning of module gridio.c for a description of how the NN_INDICES extension relates to the grid struct.

	Extension name: NN_INDICES
	Number of rows = number of grid points * average number of Delaunay links per point.

	Columns:
		LINK_I		V

Note that data types in all capitals are defined in fitsio.h.
  */

  unsigned int *linkIs=NULL,i_ui;
  int status=0;
  LONGLONG firstRow=1, firstElem=1;
  int numCols = 1;
  char extname[] = "NN_INDICES";

  if (nnLinks==NULL){
    if(!silent) bail_out("No link or near-neighbour data!");
    exit(1);
  }

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

  linkIs = malloc(sizeof(*linkIs)*gridInfo.nNNIndices);
  for(i_ui=0;i_ui<gridInfo.nNNIndices;i_ui++)
    linkIs[i_ui] = nnLinks[i_ui]->id;
  fits_write_col(fptr, dataType, 1, firstRow, firstElem, (LONGLONG)gridInfo.nNNIndices, linkIs, &status);
  processFitsError(status);
  free(linkIs);
}

/*....................................................................*/
void
writeLinksExtToFits(fitsfile *fptr, struct gridInfoType gridInfo\
  , struct linkType *links){
  /*
See the comment at the beginning of module gridio.c for a description of how the LINKS extension relates to the grid struct.

Notes:
  - Data types in all capitals are defined in fitsio.h.

  - The business where the column names are first loaded into tempColNames, then the pointers to each one into ttype, is done because:
    . fits_create_tbl() demands a 5th argument of type char*[] and will not accept char**;

    . sprintf() cannot be used to write a string to an unallocated char*, thus something like the following is illegal:
        sprintf(ttype[colI], "ACOEFF_%d", n+1);

    . The following is legal, but cannot be used for all column names, because we have to generate some on the fly:
        ttype[colI] = <column name string literal>;

  */
  const int numColNameChars=21;
  unsigned int *ids=NULL,i_ui;
  unsigned short i_us,j_us;
  double *vels=NULL;
  int status=0, colI=0, i;
  LONGLONG firstRow=1, firstElem=1;
  int numCols;
  char extname[] = "LINKS";
  char **tempColNames=NULL;

  if(links==NULL){
    if(!silent) bail_out("No link data!");
    exit(1);
  }

  if(links[0].vels==NULL)
    numCols = 2;
  else
    numCols = 2 + gridInfo.nDims*gridInfo.nLinkVels;

  /* Define the name, datatype, and physical units for the columns.
  */
  char *ttype[numCols];
  char *tform[numCols];
  char *tunit[numCols];
  int dataTypes[numCols];

  tempColNames = malloc(sizeof(*tempColNames)*numCols);
  for(i=0;i<numCols;i++) tempColNames[i]=malloc(sizeof(char)*numColNameChars);

  colI = 0;
  sprintf(tempColNames[colI], "GRID_I_1");
  tform[colI] = "V";
  tunit[colI] = "\0";
  dataTypes[colI] = TUINT;

  colI++;
  sprintf(tempColNames[colI], "GRID_I_2");
  tform[colI] = "V";
  tunit[colI] = "\0";
  dataTypes[colI] = TUINT;

  if(links[0].vels!=NULL){
    /* Should rather have vector columns? */
    for(i_us=0;i_us<gridInfo.nLinkVels;i_us++){
      for(j_us=0;j_us<gridInfo.nDims;j_us++){
        colI++;
        sprintf(tempColNames[colI], "V_%d_%d", (int)i_us+1, (int)j_us+1);
        tform[colI] = "D";
        tunit[colI] = "\0";
        dataTypes[colI] = TDOUBLE;
      }
    }
  }

  for(i=0;i<numCols;i++)
    ttype[i] = tempColNames[i];

  /* Append a new empty binary table onto the FITS file.
  */
  fits_create_tbl(fptr, BINARY_TBL, 0, numCols, ttype, tform, tunit, extname, &status);
  processFitsError(status);

  for(i=0;i<numCols;i++)
    free(tempColNames[i]);
  free(tempColNames);

  /* Write columns.
  */
  colI = 0;
  ids = malloc(sizeof(*ids)*gridInfo.nLinks);
  for(i_ui=0;i_ui<gridInfo.nLinks;i_ui++)
    ids[i_ui] = links[i_ui].gis[0];
  fits_write_col(fptr, dataTypes[colI], colI+1, firstRow, firstElem, (LONGLONG)gridInfo.nLinks, ids, &status);
  processFitsError(status);

  colI++;
  for(i_ui=0;i_ui<gridInfo.nLinks;i_ui++)
    ids[i_ui] = links[i_ui].gis[1];
  fits_write_col(fptr, dataTypes[colI], colI+1, firstRow, firstElem, (LONGLONG)gridInfo.nLinks, ids, &status);
  processFitsError(status);
  free(ids);

  if(links[0].vels!=NULL){
    vels = malloc(sizeof(*vels)*gridInfo.nLinks);
    for(i_us=0;i_us<gridInfo.nLinkVels;i_us++){
      for(j_us=0;j_us<gridInfo.nDims;j_us++){
        colI++;
        for(i_ui=0;i_ui<gridInfo.nLinks;i_ui++)
          vels[i_ui] = links[i_ui].vels[gridInfo.nDims*i_us + j_us];
        fits_write_col(fptr, dataTypes[colI], colI+1, firstRow, firstElem, (LONGLONG)gridInfo.nLinks, vels, &status);
        processFitsError(status);
      }
    }
    free(vels);
  }
}

/*....................................................................*/
void
writePopsExtToFits(fitsfile *fptr, struct gridInfoType gridInfo\
  , const unsigned short speciesI, struct grid *gp){
  /*
	Extension name: LEVEL_POPS_m (1 per molecular species)
	This is an image extension of size (number of grid cells)*(number of energy levels this species)

	Keywords:
		MOL_NAME

Note that data types in all capitals are defined in fitsio.h.
  */

  const unsigned int totalNumGridPoints = gridInfo.nInternalPoints+gridInfo.nSinkPoints;
  unsigned int i_ui;
  int status=0, xi;
  char extname[13];
  float *row=NULL;
  int bitpix = FLOAT_IMG;
  const long naxis = 2;  /* i.e. 2-dimensional image */    
  int numEnergyLevels = (int)gridInfo.mols[speciesI].nLevels;
  long naxes[] = { (long)numEnergyLevels, (long)totalNumGridPoints };
  long fpixels[naxis],lpixels[naxis];

  sprintf(extname, "LEVEL_POPS_%d", (int)speciesI+1);

  fits_create_img(fptr, bitpix, naxis, naxes, &status);
  processFitsError(status);

  row = malloc(sizeof(*row)*numEnergyLevels);

  /* Write FITS data.
  */
  for(i_ui=0;i_ui<totalNumGridPoints;i_ui++){
    for(xi=0;xi<numEnergyLevels;xi++)
      row[xi] = (float)gp[i_ui].mol[speciesI].pops[xi]; 

    fpixels[0]=1;
    fpixels[1]=i_ui+1;
    lpixels[0]=numEnergyLevels;
    lpixels[1]=i_ui+1;

    fits_write_subset(fptr, TFLOAT, fpixels, lpixels, row, &status);
    processFitsError(status);
  }

  free(row);

  /* write keywords:
  */
  fits_write_key(fptr, TSTRING, "MOL_NAME ", gridInfo.mols[speciesI].molName, "\0", &status);
  processFitsError(status);

  fits_write_key(fptr, TSTRING, "EXTNAME ", extname, "\0", &status);
  processFitsError(status);
}

/*....................................................................*/
fitsfile *
openFITSFileForRead(char *inFileName){
  fitsfile *fptr=NULL;
  int status=0;

  fits_open_file(&fptr, inFileName, READONLY, &status);
  processFitsError(status);

  return fptr;
}

/*....................................................................*/
int
countColsBasePlusInt(fitsfile *fptr, char *baseName){
  char colName[20];
  int i, status, colNum;

  i = 0;
  status = 0;
  while(!status){
    sprintf(colName, "%s%d", baseName, i+1);
    fits_get_colnum(fptr, CASEINSEN, colName, &colNum, &status);
    i++;
  }
  if(status!=COL_NOT_FOUND)
    processFitsError(status);

  return i-1;
}

/*....................................................................*/
int
countKeywords(fitsfile *fptr, char *baseName){
  char kwdName[9];
  int i, status;
  char kwdValue[80];

  i = 0;
  status = 0;
  while(!status){
    sprintf(kwdName, "%s%d", baseName, i+1);
    fits_read_key(fptr, TSTRING, kwdName, kwdValue, NULL, &status);
    i++;
  }
  if(status!=KEY_NO_EXIST)
    processFitsError(status);

  return i-1;
}

/*....................................................................*/
void
readKeywordsFromFits(lime_fptr *fptr, struct keywordType *kwds\
  , const int numKeywords){

  int i, status;
  char message[80];

  for(i=0;i<numKeywords;i++){
    status = 0;

    if(     kwds[i].datatype==TSTRING)
      fits_read_key(fptr, TSTRING, kwds[i].keyname, kwds[i].charValue, kwds[i].comment, &status);
    else if(kwds[i].datatype==TINT)
      fits_read_key(fptr, TINT, kwds[i].keyname, &kwds[i].intValue, kwds[i].comment, &status);
    else if(kwds[i].datatype==TFLOAT)
      fits_read_key(fptr, TFLOAT, kwds[i].keyname, &kwds[i].floatValue, kwds[i].comment, &status);
    else if(kwds[i].datatype==TDOUBLE)
      fits_read_key(fptr, TDOUBLE, kwds[i].keyname, &kwds[i].doubleValue, kwds[i].comment, &status);
    else{
      if(!silent){
        sprintf(message, "Keyword %d dataype %d is not currently accepted.", i, kwds[i].datatype);
        bail_out(message);
      }
      exit(1);
    }
    processFitsError(status);
  }
}

/*....................................................................*/
void
readGridExtFromFits(fitsfile *fptr, struct gridInfoType *gridInfoRead\
  , struct grid **gp, unsigned int **firstNearNeigh, char ***collPartNames\
  , int *numCollPartRead, int *dataFlags){
  /*
The present function mallocs 'gp' and sets defaults for all the simple or first-level struct elements.

If a COLLPARn keywords are found in the GRID extension header then collPartNames is malloc'd to the number of these.
  */

  LONGLONG numGridCells, firstRow=1, firstElem=1, i_LL;
  int status=0, colNum, anynul=0, i;
  char colName[20];
  char genericKwd[9];
  char message[80];
  unsigned int *ids=NULL;
  double *xj=NULL;
  _Bool *sink=NULL;/* Probably should be char* but this seems to work. */
  unsigned short *numNeigh=NULL, i_us;
  double *velj=NULL, *densn=NULL;
  float *dopb=NULL, *t=NULL, *abunm=NULL, *bField=NULL;

  /* Go to the GRID extension.
  */
  fits_movnam_hdu(fptr, BINARY_TBL, "GRID", 0, &status);
  processFitsError(status);

  /* Find out how many rows there are, then malloc the array.
  */
  fits_get_num_rowsll(fptr, &numGridCells, &status);
  processFitsError(status);
  if(numGridCells<=0){
    if(!silent) warning("No rows found in grid dataset.");
    return; /* I.e. with dataFlags left unchanged. */
  }

  mallocAndSetDefaultGrid(gp, (unsigned int)numGridCells);

  /* Read the columns.
  */
  fits_get_colnum(fptr, CASEINSEN, "ID", &colNum, &status);
  if(status==COL_NOT_FOUND){
    if(!silent) warning("No ID column found in grid dataset.");
    return; /* I.e. with dataFlags left unchanged. */
  }
  processFitsError(status);

  ids = malloc(sizeof(*ids)*numGridCells);
  fits_read_col(fptr, TUINT, colNum, firstRow, firstElem, numGridCells, 0, ids, &anynul, &status);
  processFitsError(status);

  for(i_LL=0;i_LL<numGridCells;i_LL++) {
    (*gp)[i_LL].id = (int)ids[i_LL];
  }
  free(ids);

  gridInfoRead->nDims = (unsigned short)countColsBasePlusInt(fptr, "X");
  if(gridInfoRead->nDims<=0){
    if(!silent) warning("No X columns found in grid dataset.");
    return; /* I.e. with dataFlags left unchanged. */
  }

  /* We have to do this here (as well after the call to readGrid()) because grid.x is a pre-sized array rather than a pointer we can malloc. Later this should be changed to allow us to define the sizes of all arrays in grid purely from the data in the file.
  */
  if(gridInfoRead->nDims!=DIM){
    if(!silent){
      sprintf(message, "%d Xn columns read, but there should be %d.", (int)gridInfoRead->nDims, DIM);
      bail_out(message);
    }
    exit(1);
  }

  xj = malloc(sizeof(*xj)*numGridCells);
  for(i_us=0;i_us<gridInfoRead->nDims;i_us++){
    sprintf(colName, "X%d", (int)i_us+1);
    fits_get_colnum(fptr, CASEINSEN, colName, &colNum, &status);
    processFitsError(status);

    fits_read_col(fptr, TDOUBLE, colNum, firstRow, firstElem, numGridCells, 0, xj, &anynul, &status);
    processFitsError(status);

    for(i_LL=0;i_LL<numGridCells;i_LL++) {
      (*gp)[i_LL].x[i_us] = xj[i_LL];
    }
  }
  free(xj);

  fits_get_colnum(fptr, CASEINSEN, "IS_SINK", &colNum, &status);
  if(status==COL_NOT_FOUND){
    if(!silent) warning("No IS_SINK column found in grid dataset.");
    return; /* I.e. with dataFlags left unchanged. */
  }
  processFitsError(status);

  sink = malloc(sizeof(*sink)*numGridCells);
  fits_read_col(fptr, TLOGICAL, colNum, firstRow, firstElem, numGridCells, 0, sink, &anynul, &status);
  processFitsError(status);

  gridInfoRead->nSinkPoints = 0;
  for(i_LL=0;i_LL<numGridCells;i_LL++) {
    (*gp)[i_LL].sink = (int)sink[i_LL];
    if((*gp)[i_LL].sink)
      gridInfoRead->nSinkPoints++;
  }
  free(sink);

  gridInfoRead->nInternalPoints = numGridCells - gridInfoRead->nSinkPoints;

  /* If we have made it this far, we can set the first bit of dataFlags. Woot!
  */
  (*dataFlags) |= (1 << DS_bit_x);

  fits_get_colnum(fptr, CASEINSEN, "NUMNEIGH", &colNum, &status);
  if(status!=COL_NOT_FOUND){
    processFitsError(status);

    numNeigh = malloc(sizeof(*numNeigh)*numGridCells);
    fits_read_col(fptr, TUSHORT, colNum, firstRow, firstElem, numGridCells, 0, numNeigh, &anynul, &status);
    processFitsError(status);

    for(i_LL=0;i_LL<numGridCells;i_LL++) {
      (*gp)[i_LL].numNeigh = (int)numNeigh[i_LL];
    }
    free(numNeigh);

    fits_get_colnum(fptr, CASEINSEN, "FIRST_NN", &colNum, &status);
    if(status!=COL_NOT_FOUND){
      processFitsError(status);

      *firstNearNeigh = malloc(sizeof(**firstNearNeigh)*numGridCells);
      fits_read_col(fptr, TUINT, colNum, firstRow, firstElem, numGridCells, 0, *firstNearNeigh, &anynul, &status);
      processFitsError(status);

      /* If we made it to here, we can set the neighbour bit of dataFlags. Note however that this bit is only on trial til we check for the LINKS extension. So it had better behave itself.
      */
      (*dataFlags) |= (1 << DS_bit_neighbours);
    }
  }

  /* See if there are any VEL columns:
  */
  status = 0;
  sprintf(colName, "VEL1");
  fits_get_colnum(fptr, CASEINSEN, colName, &colNum, &status);
  if(status!=COL_NOT_FOUND){
    processFitsError(status);

    /* Read the VEL columns:
    */
    velj = malloc(sizeof(*velj)*numGridCells);
    for(i_us=0;i_us<gridInfoRead->nDims;i_us++){
      sprintf(colName, "VEL%d", (int)i_us+1);
      fits_get_colnum(fptr, CASEINSEN, colName, &colNum, &status);
      processFitsError(status);

      fits_read_col(fptr, TDOUBLE, colNum, firstRow, firstElem, numGridCells, 0, velj, &anynul, &status);
      processFitsError(status);

      for(i_LL=0;i_LL<numGridCells;i_LL++) {
        (*gp)[i_LL].vel[i_us] = velj[i_LL];
      }
    }
    free(velj);

    (*dataFlags) |= (1 << DS_bit_velocity);
  }

  /* Count the numbers of DENSITYn columns:
  */
  gridInfoRead->nDensities = (unsigned short)countColsBasePlusInt(fptr, "DENSITY");
  if(gridInfoRead->nDensities > 0){
    for(i_LL=0;i_LL<numGridCells;i_LL++)
      (*gp)[i_LL].dens = malloc(sizeof(double)*gridInfoRead->nDensities);

    /* Read the DENSITY columns:
    */
    densn = malloc(sizeof(*densn)*numGridCells);
    for(i_us=0;i_us<gridInfoRead->nDensities;i_us++){
      sprintf(colName, "DENSITY%d", (int)i_us+1);
      fits_get_colnum(fptr, CASEINSEN, colName, &colNum, &status);
      processFitsError(status);

      fits_read_col(fptr, TDOUBLE, colNum, firstRow, firstElem, numGridCells, 0, densn, &anynul, &status);
      processFitsError(status);

      for(i_LL=0;i_LL<numGridCells;i_LL++) {
        (*gp)[i_LL].dens[i_us] = densn[i_LL];
      }
    }
    free(densn);

    (*dataFlags) |= (1 << DS_bit_density);
  }

  /* Count the numbers of ABUNMOLm columns:
  */
  gridInfoRead->nSpecies = (unsigned short)countColsBasePlusInt(fptr, "ABUNMOL");
  if(gridInfoRead->nSpecies > 0){
    for(i_LL=0;i_LL<numGridCells;i_LL++) {
      (*gp)[i_LL].abun = malloc(sizeof(double)*gridInfoRead->nSpecies);
    }

    /* Read the ABUNMOL columns:
    */
    abunm = malloc(sizeof(*abunm)*numGridCells);
    for(i_us=0;i_us<gridInfoRead->nSpecies;i_us++){
      sprintf(colName, "ABUNMOL%d", (int)i_us+1);
      fits_get_colnum(fptr, CASEINSEN, colName, &colNum, &status);
      processFitsError(status);

      fits_read_col(fptr, TFLOAT, colNum, firstRow, firstElem, numGridCells, 0, abunm, &anynul, &status);
      processFitsError(status);

      for(i_LL=0;i_LL<numGridCells;i_LL++) {
        (*gp)[i_LL].abun[i_us] = (double)abunm[i_LL];
      }
    }
    free(abunm);

    (*dataFlags) |= (1 << DS_bit_abundance);
  }

  /* Read the TURBDPLR column:
  */
  fits_get_colnum(fptr, CASEINSEN, "TURBDPLR", &colNum, &status);
  if(status!=COL_NOT_FOUND){
    processFitsError(status);

    dopb = malloc(sizeof(*dopb)*numGridCells);
    fits_read_col(fptr, TFLOAT, colNum, firstRow, firstElem, numGridCells, 0, dopb, &anynul, &status);
    processFitsError(status);

    for(i_LL=0;i_LL<numGridCells;i_LL++) {
      (*gp)[i_LL].dopb_turb = (double)dopb[i_LL];
    }
    free(dopb);

    (*dataFlags) |= (1 << DS_bit_turb_doppler);
  }

//*** there is probably a neater way to do the temperatures. Resetting the t0 defaults is a bit ugly.
  /* Read the TEMPKNTC column:
  */
  status = 0;
  fits_get_colnum(fptr, CASEINSEN, "TEMPKNTC", &colNum, &status);
  if(status!=COL_NOT_FOUND){
    processFitsError(status);

    t = malloc(sizeof(*t)*numGridCells);
    fits_read_col(fptr, TFLOAT, colNum, firstRow, firstElem, numGridCells, 0, t, &anynul, &status);
    processFitsError(status);

    for(i_LL=0;i_LL<numGridCells;i_LL++) {
      (*gp)[i_LL].t[0] = (double)t[i_LL];
    }

    /* Read the TEMPDUST column:
    */
    fits_get_colnum(fptr, CASEINSEN, "TEMPDUST", &colNum, &status);
    if(status==COL_NOT_FOUND){ /* Set t0 back to defaults: */
      for(i_LL=0;i_LL<numGridCells;i_LL++) {
        (*gp)[i_LL].t[0] = -1;
      }
    }else{
      processFitsError(status);

      fits_read_col(fptr, TFLOAT, colNum, firstRow, firstElem, numGridCells, 0, t, &anynul, &status);
      processFitsError(status);

      for(i_LL=0;i_LL<numGridCells;i_LL++) {
        (*gp)[i_LL].t[1] = (double)t[i_LL];
      }

      (*dataFlags) |= (1 << DS_bit_temperatures);
    }
    free(t);
  }

  /* See if there are any B_FIELD columns:
  */
  status = 0;
  sprintf(colName, "B_FIELD1");
  fits_get_colnum(fptr, CASEINSEN, colName, &colNum, &status);
  if(status!=COL_NOT_FOUND){
    processFitsError(status);

    /* Read the B_FIELD columns:
    */
    bField = malloc(sizeof(*bField)*numGridCells);
    for(i=0;i<3;i++){//**** keep this hard-wired or rather test that gridInfo.nDims==3??
      sprintf(colName, "B_FIELD%d", i+1);
      fits_get_colnum(fptr, CASEINSEN, colName, &colNum, &status);
      processFitsError(status);

      fits_read_col(fptr, TFLOAT, colNum, firstRow, firstElem, numGridCells, 0, bField, &anynul, &status);
      processFitsError(status);

      for(i_LL=0;i_LL<numGridCells;i_LL++) {
        (*gp)[i_LL].B[i] = (double)bField[i_LL];
      }
    }
    free(bField);

    (*dataFlags) |= (1 << DS_bit_magfield);
  }

  /* Check if there are any COLLPAR keywords.
  */
  *numCollPartRead = countKeywords(fptr, "COLLPAR");
  if(*numCollPartRead <= 0){
    *collPartNames = NULL;
  }else{
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
void
readLinksExtFromFits(fitsfile *fptr, struct gridInfoType *gridInfoRead\
  , struct grid *gp, struct linkType **links, int *dataFlags){
  /*
See the comment at the beginning of gridio.c for a description of how the LINKS extension relates to the grid struct.

The present function mallocs the pointer *links.
  */

  LONGLONG totalNumLinks, firstRow=1, firstElem=1, i_LL;
  int status=0,colNum,anynul=0,i;
  char colName[21];
  unsigned int *ids=NULL, totalNumGridPoints, i_ui;
  double *vels=NULL;
  char message[80];
  unsigned short i_us,j_us;
  int colGrid1Found, colGrid2Found; //->bool

  if(!bitIsSet(*dataFlags, DS_bit_neighbours))
    return;

  totalNumGridPoints = gridInfoRead->nInternalPoints + gridInfoRead->nSinkPoints; /* Just for a bit more brevity. */

  /* Go to the LINKS extension.
  */
  fits_movnam_hdu(fptr, BINARY_TBL, "LINKS", 0, &status);
  if(status==BAD_HDU_NUM){
    if(!silent) warning("No LINKS extension found in grid dataset.");

    /* Unset the DS_mask_neighbours bit: */
    *dataFlags &= ~(1 << DS_bit_neighbours);
    return;
  }
  processFitsError(status);

  /* Find out how many rows there are, then malloc the array.
  */
  fits_get_num_rowsll(fptr, &totalNumLinks, &status);
  processFitsError(status);

  gridInfoRead->nLinks = (unsigned int)totalNumLinks;
  if(gridInfoRead->nLinks<=0){
    if(!silent) warning("No rows in LINKS extension of grid dataset.");

    /* Unset the neighbours bit: */
    *dataFlags &= ~(1 << DS_bit_neighbours);
    return;
  }
  *links = malloc(sizeof(**links)*totalNumLinks);

  for(i_LL=0;i_LL<totalNumLinks;i_LL++) {
    (*links)[i_LL].id = (unsigned int)i_LL;
  }

  /* Read GRID_I_1 column.
  */
  colGrid1Found = 0;
  colGrid2Found = 0;
  fits_get_colnum(fptr, CASEINSEN, "GRID_I_1", &colNum, &status);
  if(status!=COL_NOT_FOUND){
    processFitsError(status);
    colGrid1Found = 1;

    ids = malloc(sizeof(*ids)*totalNumLinks);
    fits_read_col(fptr, TUINT, colNum, firstRow, firstElem, totalNumLinks, 0, ids, &anynul, &status);
    processFitsError(status);

    for(i_LL=0;i_LL<totalNumLinks;i_LL++) {
      i_ui = ids[i_LL];
      if(i_ui<0 || i_ui>=totalNumGridPoints){
        if(!silent){
          sprintf(message, "GRID_I_1 %dth-row value %ud is outside range [0,%ud]", (int)i_LL, i_ui, totalNumGridPoints);
          bail_out(message);
        }
        exit(1);
      }
      (*links)[i_LL].gis[0] = gp[i_ui].id;
    }

    /* Read GRID_I_2 column.
    */
    fits_get_colnum(fptr, CASEINSEN, "GRID_I_2", &colNum, &status);
    if(status!=COL_NOT_FOUND){
      processFitsError(status);
      colGrid2Found = 1;

      fits_read_col(fptr, TUINT, colNum, firstRow, firstElem, totalNumLinks, 0, ids, &anynul, &status);
      processFitsError(status);

      for(i_LL=0;i_LL<totalNumLinks;i_LL++) {
        i_ui = ids[i_LL];
        if(i_ui<0 || i_ui>=totalNumGridPoints){
          if(!silent){
            sprintf(message, "GRID_I_2 %dth-row value %ud is outside range [0,%ud]", (int)i_LL, i_ui, totalNumGridPoints);
            bail_out(message);
          }
          exit(1);
        }
        (*links)[i_LL].gis[1] = gp[i_ui].id;
      }
    }
    free(ids);
  }

  if(!(colGrid1Found && colGrid2Found)){
    if(!silent) warning("Both GRID_I columns not found in LINKS extension of grid dataset.");
    free(*links);
    *links = NULL;
    /* Unset the neighbours bit: */
    *dataFlags &= ~(1 << DS_bit_neighbours);
    return;
  }

  /* Find out how many V_* columns there are.
  */
  i = 0;
  status = 0;
  while(!status){
    sprintf(colName, "V_%d_", i+1);
    if(countColsBasePlusInt(fptr, colName)!=gridInfoRead->nDims)
      status = 1;

    i++;
  }
  gridInfoRead->nLinkVels = i - 1;
  status = 0;

  if(gridInfoRead->nLinkVels<=0){
    for(i_LL=0;i_LL<totalNumLinks;i_LL++)
      (*links)[i_LL].vels = NULL;
    return;
  }

  for(i_LL=0;i_LL<totalNumLinks;i_LL++)
    (*links)[i_LL].vels = malloc(sizeof(double)*gridInfoRead->nLinkVels*gridInfoRead->nDims);

  vels = malloc(sizeof(*vels)*totalNumLinks);
  for(i_us=0;i_us<gridInfoRead->nLinkVels;i_us++){
    for(j_us=0;j_us<gridInfoRead->nDims;j_us++){
      /* Read the V_n_d columns.
      */
      sprintf(colName, "V_%d_%d", (int)i_us+1, (int)j_us+1);
      fits_get_colnum(fptr, CASEINSEN, colName, &colNum, &status);
      processFitsError(status);

      fits_read_col(fptr, TDOUBLE, colNum, firstRow, firstElem, totalNumLinks, 0, vels, &anynul, &status);
      processFitsError(status);

      for(i_LL=0;i_LL<totalNumLinks;i_LL++)
        (*links)[i_LL].vels[gridInfoRead->nDims*j_us + i_us] = vels[i_LL];
    }
  }
  free(vels);

  (*dataFlags) |= (1 << DS_bit_ACOEFF);

}

/*....................................................................*/
void
readNnIndicesExtFromFits(fitsfile *fptr, struct linkType *links\
  , struct linkType ***nnLinks, struct gridInfoType *gridInfoRead, int *dataFlags){
  /*
See the comment at the beginning of gridio.c for a description of how the NN_INDICES extension relates to the grid struct.

The function mallocs the pointer *nnLinks.
  */

  LONGLONG totalNumNeigh, firstRow=1, firstElem=1, i_LL;
  int status=0, colNum, anynul=0;
  unsigned int *linkIs=NULL;

  if(!bitIsSet(*dataFlags, DS_bit_neighbours))
    return;

  /* Go to the NN_INDICES extension.
  */
  fits_movnam_hdu(fptr, BINARY_TBL, "NN_INDICES", 0, &status);
  if(status==BAD_HDU_NUM){
    if(!silent) warning("No NN_INDICES extension found in grid dataset.");

    /* Unset the neighbours bit: */
    *dataFlags &= ~(1 << DS_bit_neighbours);
    return;
  }
  processFitsError(status);

  /* Find out how many rows there are, then malloc the array.
  */
  fits_get_num_rowsll(fptr, &totalNumNeigh, &status);
  processFitsError(status);

  gridInfoRead->nNNIndices = (unsigned int)totalNumNeigh;
  if(gridInfoRead->nNNIndices<=0){
    if(!silent) warning("No rows in NN_INDICES extension of grid dataset.");

    /* Unset the neighbours bit: */
    *dataFlags &= ~(1 << DS_bit_neighbours);
    return;
  }
  *nnLinks = malloc(sizeof(**nnLinks)*totalNumNeigh);

  /* Read LINK_I column.
  */
  fits_get_colnum(fptr, CASEINSEN, "LINK_I", &colNum, &status);
  if(status==COL_NOT_FOUND){
    if(!silent) warning("LINK_I column not found in NN_INDICES extension of grid dataset.");

    /* Unset the neighbours bit: */
    *dataFlags &= ~(1 << DS_bit_neighbours);
    return;
  }
  processFitsError(status);

  linkIs = malloc(sizeof(*linkIs)*totalNumNeigh);
  fits_read_col(fptr, TUINT, colNum, firstRow, firstElem, totalNumNeigh, 0, linkIs, &anynul, &status);
  processFitsError(status);

  for(i_LL=0;i_LL<totalNumNeigh;i_LL++)
    (*nnLinks)[i_LL] = &links[linkIs[i_LL]];

  free(linkIs);
}

/*....................................................................*/
_Bool
checkPopsFitsExtExists(fitsfile *fptr, const unsigned short speciesI){
  const unsigned short maxNumSpecies = 9;
  char message[80];
  char extname[13];
  int status=0;

  if(speciesI+1>maxNumSpecies){
    if(!silent){
      sprintf(message, "Species block %d is greater than the limit %d", (int)speciesI+1, (int)maxNumSpecies);
      bail_out(message);
    }
    exit(1);
  }

  sprintf(extname, "LEVEL_POPS_%d", (int)speciesI+1);

  /* Try to move to extension LEVEL_POPS_<speciesI>:
  */
  fits_movnam_hdu(fptr, IMAGE_HDU, extname, 0, &status);
  if(status==BAD_HDU_NUM)
    return 0;

  /* If we reached here it means we found the extension, but now we check for other errors:
  */
  processFitsError(status);

  return 1;
}

/*....................................................................*/
void
readPopsExtFromFits(fitsfile *fptr, const unsigned short speciesI\
  , struct grid *gp, struct gridInfoType *gridInfoRead){
  /*
See the comment at the beginning of the present module for a description of how the LEVEL_POPS_m extensions relate to the grid struct.

The function mallocs gp[yi].mol[speciesI].pops for each grid point yi and species.
  */

  const int maxLenMolName = 8;
  const unsigned short maxNumSpecies = 9;
  int naxis, status=0, xi, anynul=0, bitpix;
/*The interface to fits_get_img_param() says long* for argument 5 but I get a seg fault unless I use an array.
  long *naxes; */
long naxes[2];
  float *row=NULL;
  long fpixels[2],lpixels[2];
  long inc[2] = {1,1};
  char molNameRead[maxLenMolName+1];
  char message[80];
  char extname[13];
  unsigned int numGridPoints, i_ui;

  if(speciesI+1>maxNumSpecies){
    if(!silent){
      sprintf(message, "Species block %d is greater than the limit %d", (int)speciesI+1, (int)maxNumSpecies);
      bail_out(message);
    }
    exit(1);
  }
  sprintf(extname, "LEVEL_POPS_%d", (int)speciesI+1);

  /* Try to move to extension LEVEL_POPS_<speciesI>:
  */
  fits_movnam_hdu(fptr, IMAGE_HDU, extname, 0, &status);
  processFitsError(status);

  fits_get_img_param(fptr, 2, &bitpix, &naxis, naxes, &status);
  processFitsError(status);

  numGridPoints = gridInfoRead->nInternalPoints + gridInfoRead->nSinkPoints;
  if((long)numGridPoints != naxes[1]){
    if(!silent){
      sprintf(message, "Expected %ld grid points but extension %s has %ld"\
        , (long)(gridInfoRead->nInternalPoints + gridInfoRead->nSinkPoints), extname, naxes[1]);
      bail_out(message);
    }
    exit(1);
  }
  gridInfoRead->mols[speciesI].nLevels = (int)naxes[0];
  gridInfoRead->mols[speciesI].nLines = -1;
/*  free(naxes); */

  row = malloc(sizeof(*row)*gridInfoRead->mols[speciesI].nLevels);

  /* Read FITS data.
  */
  for(i_ui=0;i_ui<numGridPoints;i_ui++){
    fpixels[0]=1;
    fpixels[1]=(int)i_ui+1;
    lpixels[0]=gridInfoRead->mols[speciesI].nLevels;
    lpixels[1]=(int)i_ui+1;

    fits_read_subset(fptr, TFLOAT, fpixels, lpixels, inc, 0, row, &anynul, &status);
    processFitsError(status);

    gp[i_ui].mol[speciesI].pops = malloc(sizeof(double)*gridInfoRead->mols[speciesI].nLevels);
    for(xi=0;xi<gridInfoRead->mols[speciesI].nLevels;xi++)
      gp[i_ui].mol[speciesI].pops[xi] = (double)row[xi];

  }

  free(row);

  /* Read kwds:
  */
  fits_read_key(fptr, TSTRING, "MOL_NAME", molNameRead, NULL, &status);
  gridInfoRead->mols[speciesI].molName = molNameRead;//*****??
  processFitsError(status);
}


