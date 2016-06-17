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
	- If any of the stage 3 columns (VEL etc) are not present, get the values via user-supplied functions.
*/

#include "lime.h"


/*
The present module contains routines for transferring the LIME grid point data to or from a FITS format file. The purpose of the present comment block is to describe the FITS file format. Note that the amount and type of information stored depends on the 'data stage' of the grid struct, as described in the header remarks to module gridio.c. In the description below, the minimum stage integer associated with the presence of a particular extension, column or keyword is given on the leftmost place of each line.

Note that all extensions are binary table except where indicated. The letter in the second row for column descriptions gives the FITS data type. See eg

  https://heasarc.gsfc.nasa.gov/docs/software/fitsio/c/c_user/node20.html

for a key to these.

1	0) The primary HDU
	Keywords:
1		STAGE		I

1	1) GRID
	Number of rows = number of grid points.

	Keywords:
1		RADIUS		D
1		COLLPARn	A	# 1 for each nth collision partner.

	Columns:
1		ID		V
1		Xj		D	# Cartesian components of the point location, 1 col per jth dimension.
1		IS_SINK		L	# =True iff the point lies on the edge of the model.
2		NUMNEIGH	U
2		FIRST_NN	V	# See explanation in section 3 below.
3		VELj		D	# 1 col per jth dimension.
3		DENSITYn	D	# 1 per nth collision partner.
3		TEMPKNTC	E	# From t[0].
3		TEMPDUST	E	# From t[1].
3		TURBDPLR	E	# Given Gaussian lineshape exp(-v^2/[B^2 + 2*k*T/m]), this is B.
3		ABUNMOLm	E	# 1 per mth molecular species.

2	2) NN_INDICES (see explanation in the header to gridio.c)
	Number of rows = number of grid points * average number of Delaunay links per point.

	Columns:
2		LINK_I		V

2	3) LINKS (see explanation in the header to gridio.c)
	Number of rows = number of Delaunay links.

	Columns:
2		GRID_I_1	V
2		GRID_I_2	V
3		ACOEFF_p	D	# 1 per pth order of the velocity polynomial.

4	4 etc) LEVEL_POPS_m (1 per mth molecular species)
	This is an image extension of size (number of grid cells)*(number of energy levels this species).

	Keywords:
4		MOL_NAME

Note that at present the data in the 'partner' element of grid.mol is *NOT* being stored.
*/

/*....................................................................*/
fitsfile *openFITSFileForRead(char *inFileName, int *dataStageI){
  fitsfile *fptr=NULL;
  int status=0;
  short tempStageI;

  fits_open_file(&fptr, inFileName, READONLY, &status);
  processFitsError(status);

  fits_read_key(fptr, TSHORT, "STAGE   ", &tempStageI, NULL, &status);
  processFitsError(status);

  *dataStageI = (int)tempStageI;

  return fptr;
}

/*....................................................................*/
fitsfile *openFITSFileForWrite(char *outFileName, const int dataStageI){
  fitsfile *fptr=NULL;
  int status=0;
  char negfile[100]="! ";
  short dsiShort = (short)dataStageI, *ptrToShort=&dsiShort; /* Have to do this otherwise the compiler complains that the const int is not an 'lvalue'. */

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

  fits_write_key(fptr, TSHORT, "STAGE   ", ptrToShort, "Data stage", &status);
  processFitsError(status);

  return fptr;
}

/*....................................................................*/
void closeFITSFile(fitsfile *fptr){
  int status=0;

  fits_close_file(fptr, &status);
  processFitsError(status);
}

/*....................................................................*/
void defineGridExtColumns(const unsigned short numKwdChars, inputPars par\
  , const unsigned short numDims, const int dataStageI, char *ttype[]\
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

  if(par.nSpecies>maxNumSpecies){
    if(!silent){
      sprintf(message, "Caller asked for %d species but colnames can only be written for %d.", par.nSpecies, maxNumSpecies);
      bail_out(message);
    }
    exit(1);
  }

  if(par.collPart>maxNumCollPart){
    if(!silent){
      sprintf(message, "Caller asked for %d coll. part. but colnames can only be written for %d.", par.collPart, maxNumCollPart);
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

  if(dataStageI<2)
    return;

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

  if(dataStageI<3)
    return;

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
  for(i=0;i<par.collPart;i++){
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
  for(i=0;i<par.nSpecies;i++){
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
  , struct grid *g, unsigned int *firstNearNeigh, char **collPartNames\
  , const int dataStageI){

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
  int numCols;
  char extname[] = "GRID";
  char genericComment[80];
  char genericKwd[numKwdChars], message[80];

  if(dataStageI<1){
    if(!silent) bail_out("Data stage indicates no grid data!");
    exit(1);
  }

  if(dataStageI<2)
    numCols = 2 + numDims;
  else if(dataStageI<3)
    numCols = 4 + numDims;
  else
    numCols = 7 + numDims*2 + par.collPart + par.nSpecies;

  /* Define the name, datatype, and physical units for the columns.
  */
  char *ttype[numCols];
  char *tform[numCols];
  char *tunit[numCols];
  int dataTypes[numCols];

  defineGridExtColumns(numKwdChars, par, numDims, dataStageI, ttype, tform, tunit, dataTypes);

  /* Append a new empty binary table onto the FITS file.
  */
  fits_create_tbl( fptr, BINARY_TBL, 0, numCols, ttype, tform, tunit, extname, &status);
  processFitsError(status);

  for(colI=0;colI<numCols;colI++)
    free(ttype[colI]); /* This is freed because it is specifically malloc'd inside defineGridExtColumns. */

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

  if (dataStageI>1){
    colI++;
    numNeigh = malloc(sizeof(*numNeigh)*par.ncell);
    for(i=0;i<par.ncell;i++) numNeigh[i] = (unsigned short)g[i].numNeigh;
    fits_write_col(fptr, dataTypes[colI], colI+1, firstRow, firstElem, (LONGLONG)par.ncell, numNeigh, &status);
    processFitsError(status);
    free(numNeigh);

    colI++;
    fits_write_col(fptr, dataTypes[colI], colI+1, firstRow, firstElem, (LONGLONG)par.ncell, firstNearNeigh, &status);
    processFitsError(status);

    if (dataStageI>2){
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
void writeNnIndicesExtToFits(fitsfile *fptr, const unsigned int totalNumNeigh\
  , struct linkType **nnLinks, struct linkType *links){
  /*
See the comment at the beginning of module gridio.c for a description of how the NN_INDICES extension relates to the grid struct.

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
void writeLinksExtToFits(fitsfile *fptr, const unsigned int totalNumLinks\
  , const unsigned short numACoeffs, struct linkType *links, const int dataStageI){
  /*
See the comment at the beginning of module gridio.c for a description of how the LINKS extension relates to the grid struct.

Note that data types in all capitals are defined in fitsio.h.
  */

  unsigned int *ids=NULL;
  double *aCoeffs=NULL;
  int status=0, colI=0, i, n;
  LONGLONG firstRow=1, firstElem=1;
  int numCols;
  char extname[] = "LINKS";
  const int numKwdChars = 9; /* 8 characters + \0. */

  if(dataStageI<2){
    if(!silent) bail_out("Data stage indicates no link data!");
    exit(1);
  }

  if(dataStageI<3)
    numCols = 2;
  else
    numCols = 2 + numACoeffs;

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

  if(dataStageI>2){
    /* Should rather have a vector column? */
    for(n=0;n<numACoeffs;n++){
      colI++;
      ttype[colI] = malloc(numKwdChars);
      sprintf(ttype[colI], "ACOEFF_%d", n+1);
      tform[colI] = "D";
      tunit[colI] = "\0";
      dataTypes[colI] = TDOUBLE;
    }
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

  if(dataStageI>2){
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
void writePopsExtToFits(fitsfile *fptr, const unsigned int numGridPoints\
  , molData *md, const unsigned short speciesI, struct grid *g){
  /*
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
int countColumns(fitsfile *fptr, char *baseName){
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
int countKeywords(fitsfile *fptr, char *baseName){
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
void readGridExtFromFits(fitsfile *fptr, struct gridInfoType *gridInfoRead\
  , struct grid **gp, unsigned int **firstNearNeigh, char ***collPartNames\
  , int *numCollPartRead, const int dataStageI){
  /*
The present function mallocs 'g'. If dataStageI>1 then it also mallocs 'firstNearNeigh'. If dataStageI>2 then the function also mallocs the following elements of struct grid for all grid points:
	dens
	abun
Otherwise these are set to NULL.

If a COLLPARn keywords are found in the GRID extension header then collPartNames is malloc'd to the number of these.
  */

  LONGLONG numGridCells, firstRow=1, firstElem=1, i_LL;
  int status=0, colNum, anynul=0, i;
  char colName[20];
  double modelRadius;
  char genericKwd[9];
  char message[80];
  unsigned int *ids=NULL;
  double *xj=NULL;
  _Bool *sink=NULL;/* Probably should be char* but this seems to work. */
  unsigned short *numNeigh=NULL, i_us;
  double *velj=NULL, *densn=NULL;
  float *dopb=NULL, *t=NULL, *abunm=NULL;

  /* Go to the GRID extension.
  */
  fits_movnam_hdu(fptr, BINARY_TBL, "GRID", 0, &status);
  processFitsError(status);

  /* Find out how many rows there are, then malloc the array.
  */
  fits_get_num_rowsll(fptr, &numGridCells, &status);
  processFitsError(status);

  mallocAndSetDefaultGrid(gp, (unsigned int)numGridCells);

  /* Read the columns.
  */
  fits_get_colnum(fptr, CASEINSEN, "ID", &colNum, &status);
  processFitsError(status);

  ids = malloc(sizeof(*ids)*numGridCells);
  fits_read_col(fptr, TUINT, colNum, firstRow, firstElem, numGridCells, 0, ids, &anynul, &status);
  processFitsError(status);

  for(i_LL=0;i_LL<numGridCells;i_LL++) {
    (*gp)[i_LL].id = (int)ids[i_LL];
  }
  free(ids);

  gridInfoRead->nDims = (unsigned short)countColumns(fptr, "X");

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

  if(dataStageI>1){
    fits_get_colnum(fptr, CASEINSEN, "NUMNEIGH", &colNum, &status);
    processFitsError(status);

    numNeigh = malloc(sizeof(*numNeigh)*numGridCells);
    fits_read_col(fptr, TUSHORT, colNum, firstRow, firstElem, numGridCells, 0, numNeigh, &anynul, &status);
    processFitsError(status);

    for(i_LL=0;i_LL<numGridCells;i_LL++) {
      (*gp)[i_LL].numNeigh = (int)numNeigh[i_LL];
    }
    free(numNeigh);

    fits_get_colnum(fptr, CASEINSEN, "FIRST_NN", &colNum, &status);
    processFitsError(status);

    *firstNearNeigh = malloc(sizeof(**firstNearNeigh)*numGridCells);
    fits_read_col(fptr, TUINT, colNum, firstRow, firstElem, numGridCells, 0, *firstNearNeigh, &anynul, &status);
    processFitsError(status);
  }

  if(dataStageI>2){
    /* Count the numbers of DENSITYn and ABUNMOLn columns:
    */
    gridInfoRead->nDensities = (unsigned short)countColumns(fptr, "DENSITY");
    gridInfoRead->nSpecies   = (unsigned short)countColumns(fptr, "ABUNMOL");

    /* Dimension the appropriate elements of gp:
    */
    for(i_LL=0;i_LL<numGridCells;i_LL++) {
      (*gp)[i_LL].dens = malloc(sizeof(double)*gridInfoRead->nDensities);
      (*gp)[i_LL].abun = malloc(sizeof(double)*gridInfoRead->nSpecies);
      (*gp)[i_LL].nmol = malloc(sizeof(double)*gridInfoRead->nSpecies);
      for(i_us=0;i_us<gridInfoRead->nSpecies;i_us++)
        (*gp)[i_LL].nmol[i_us] = 0.0;
    }

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

    /* Read the TURBDPLR column:
    */
    fits_get_colnum(fptr, CASEINSEN, "TURBDPLR", &colNum, &status);
    processFitsError(status);

    dopb = malloc(sizeof(*dopb)*numGridCells);
    fits_read_col(fptr, TFLOAT, colNum, firstRow, firstElem, numGridCells, 0, dopb, &anynul, &status);
    processFitsError(status);

    for(i_LL=0;i_LL<numGridCells;i_LL++) {
      (*gp)[i_LL].dopb = (double)dopb[i_LL];
    }
    free(dopb);

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

    /* Read the TEMPKNTC column:
    */
    fits_get_colnum(fptr, CASEINSEN, "TEMPKNTC", &colNum, &status);
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
    processFitsError(status);

    fits_read_col(fptr, TFLOAT, colNum, firstRow, firstElem, numGridCells, 0, t, &anynul, &status);
    processFitsError(status);

    for(i_LL=0;i_LL<numGridCells;i_LL++) {
      (*gp)[i_LL].t[1] = (double)t[i_LL];
    }
    free(t);

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
  }

  /* Read kwds:
  */
  fits_read_key(fptr, TDOUBLE, "RADIUS  ", &modelRadius, NULL, &status);
  processFitsError(status);

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
void readLinksExtFromFits(fitsfile *fptr, struct gridInfoType *gridInfoRead\
  , struct grid *gp, struct linkType **links, const int dataStageI){
  /*
See the comment at the beginning of gridio.c for a description of how the LINKS extension relates to the grid struct.

The present function mallocs the pointer *links.
  */

  LONGLONG totalNumLinks, firstRow=1, firstElem=1, i_LL;
  int status=0, colNum, anynul=0;
  char colName[21];
  unsigned int *ids=NULL, totalNumGridPoints, ppi;
  double *aCoeffs=NULL;
  char message[80];
  unsigned short i_s;

  totalNumGridPoints = gridInfoRead->nInternalPoints + gridInfoRead->nSinkPoints; /* Just for a bit more brevity. */

  /* Go to the LINKS extension.
  */
  fits_movnam_hdu(fptr, BINARY_TBL, "LINKS", 0, &status);
  processFitsError(status);

  /* Find out how many rows there are, then malloc the array.
  */
  fits_get_num_rowsll(fptr, &totalNumLinks, &status);
  processFitsError(status);

  gridInfoRead->nLinks = (unsigned int)totalNumLinks;
  *links = malloc(sizeof(**links)*totalNumLinks);

  for(i_LL=0;i_LL<totalNumLinks;i_LL++) {
    (*links)[i_LL].id = (unsigned int)i_LL;
  }

  /* Read GRID_I_1 column.
  */
  fits_get_colnum(fptr, CASEINSEN, "GRID_I_1", &colNum, &status);
  processFitsError(status);

  ids = malloc(sizeof(*ids)*totalNumLinks);
  fits_read_col(fptr, TUINT, colNum, firstRow, firstElem, totalNumLinks, 0, ids, &anynul, &status);
  processFitsError(status);

  for(i_LL=0;i_LL<totalNumLinks;i_LL++) {
    ppi = ids[i_LL];
    if(ppi<0 || ppi>=totalNumGridPoints){
      if(!silent){
        sprintf(message, "GRID_I_1 %dth-row value %ud is outside range [0,%ud]", (int)i_LL, ppi, totalNumGridPoints);
        bail_out(message);
      }
      exit(1);
    }
    (*links)[i_LL].g[0] = &gp[ppi];
  }

  /* Read GRID_I_2 column.
  */
  fits_get_colnum(fptr, CASEINSEN, "GRID_I_2", &colNum, &status);
  processFitsError(status);

  fits_read_col(fptr, TUINT, colNum, firstRow, firstElem, totalNumLinks, 0, ids, &anynul, &status);
  processFitsError(status);

  for(i_LL=0;i_LL<totalNumLinks;i_LL++) {
    ppi = ids[i_LL];
    if(ppi<0 || ppi>=totalNumGridPoints){
      if(!silent){
        sprintf(message, "GRID_I_2 %dth-row value %ud is outside range [0,%ud]", (int)i_LL, ppi, totalNumGridPoints);
        bail_out(message);
      }
      exit(1);
    }
    (*links)[i_LL].g[1] = &gp[ppi];
  }
  free(ids);

  if(dataStageI>2){ /* Means we are in at least stage C. */
    /* Find out how many ACOEFF_* columns there are.
    */
    gridInfoRead->nACoeffs = (unsigned short)countColumns(fptr, "ACOEFF_");

    for(i_LL=0;i_LL<totalNumLinks;i_LL++)
      (*links)[i_LL].aCoeffs = malloc(sizeof(double)*gridInfoRead->nACoeffs);

    aCoeffs = malloc(sizeof(*aCoeffs)*totalNumLinks);
    for(i_s=0;i_s<gridInfoRead->nACoeffs;i_s++){
      /* Read the ACOEFF_n columns.
      */
      sprintf(colName, "ACOEFF_%d", (int)i_s+1);
      fits_get_colnum(fptr, CASEINSEN, colName, &colNum, &status);
      processFitsError(status);

      fits_read_col(fptr, TDOUBLE, colNum, firstRow, firstElem, totalNumLinks, 0, aCoeffs, &anynul, &status);
      processFitsError(status);

      for(i_LL=0;i_LL<totalNumLinks;i_LL++)
        (*links)[i_LL].aCoeffs[i_s] = aCoeffs[i_LL];
    }
    free(aCoeffs);

  }else{
    gridInfoRead->nACoeffs = 0;
    for(i_LL=0;i_LL<totalNumLinks;i_LL++)
      (*links)[i_LL].aCoeffs = NULL;
  } 
}

/*....................................................................*/
void readNnIndicesExtFromFits(fitsfile *fptr, struct linkType *links\
  , struct linkType ***nnLinks, struct gridInfoType *gridInfoRead){
  /*
See the comment at the beginning of gridio.c for a description of how the NN_INDICES extension relates to the grid struct.

The function mallocs the pointer *nnLinks.
  */

  LONGLONG totalNumNeigh, firstRow=1, firstElem=1, i_LL;
  int status=0, colNum, anynul=0;
  unsigned int *linkIs=NULL;

  /* Go to the NN_INDICES extension.
  */
  fits_movnam_hdu(fptr, BINARY_TBL, "NN_INDICES", 0, &status);
  processFitsError(status);

  /* Find out how many rows there are, then malloc the array.
  */
  fits_get_num_rowsll(fptr, &totalNumNeigh, &status);
  processFitsError(status);

  gridInfoRead->nNNIndices = (unsigned int)totalNumNeigh;
  *nnLinks = malloc(sizeof(**nnLinks)*totalNumNeigh);

  /* Read LINK_I column.
  */
  fits_get_colnum(fptr, CASEINSEN, "LINK_I", &colNum, &status);
  processFitsError(status);

  linkIs = malloc(sizeof(*linkIs)*totalNumNeigh);
  fits_read_col(fptr, TUINT, colNum, firstRow, firstElem, totalNumNeigh, 0, linkIs, &anynul, &status);
  processFitsError(status);

  for(i_LL=0;i_LL<totalNumNeigh;i_LL++)
    (*nnLinks)[i_LL] = &links[linkIs[i_LL]];

  free(linkIs);
}

/*....................................................................*/
_Bool checkPopsFitsExtExists(fitsfile *fptr, const unsigned short speciesI){
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
void readPopsExtFromFits(fitsfile *fptr, const unsigned short speciesI\
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
  unsigned int numGridPoints, i_u;

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
  for(i_u=0;i_u<numGridPoints;i_u++){
    fpixels[0]=1;
    fpixels[1]=(int)i_u+1;
    lpixels[0]=gridInfoRead->mols[speciesI].nLevels;
    lpixels[1]=(int)i_u+1;

    fits_read_subset(fptr, TFLOAT, fpixels, lpixels, inc, 0, row, &anynul, &status);
    processFitsError(status);

    gp[i_u].mol[speciesI].pops = malloc(sizeof(double)*gridInfoRead->mols[speciesI].nLevels);
    for(xi=0;xi<gridInfoRead->mols[speciesI].nLevels;xi++)
      gp[i_u].mol[speciesI].pops[xi] = (double)row[xi];
  }

  free(row);

  /* Read kwds:
  */
  fits_read_key(fptr, TSTRING, "MOL_NAME", molNameRead, NULL, &status);
  gridInfoRead->mols[speciesI].molName = molNameRead;//*****??
  processFitsError(status);
}


