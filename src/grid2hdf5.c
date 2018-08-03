/*
 *  grid2hdf5.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
*/

#include "lime.h"
#include "gridio.h"

/*
The present module contains routines for transferring the LIME grid point data to or from an HDF5 format file. The purpose of the present comment block is to describe the HDF5 file format. Note that the amount and type of information stored depends on the 'data stage' of the grid struct, as described in the header remarks to module gridio.c.

In the description below, the data stage bit associated with the presence of a particular extension, column or keyword is given on the leftmost place of each line. When the file is read, all the objects associated with a bit must be present for the bit to be set.

Note that all extensions are binary table except where indicated. The letter in the second row for column descriptions indicates the HDF5 data type according to the following key:

	A	H5T_NATIVE_CHAR
	I	H5T_NATIVE_SHORT
	U	H5T_NATIVE_USHORT
	V	H5T_NATIVE_UINT
	E	H5T_NATIVE_FLOAT
	D	H5T_NATIVE_DOUBLE

Where column names contain a lower-case letter, this is a placeholder for a digit as explained in the respective comment.

0    The file.
         Attributes:
0            RADIUS             D	# The model radius in metres.

0        Group 'GRID'
             Attributes:
0                COLLPARn       A	# 1 for each nth collision partner.

0            Group 'columns'
                 Number of elements in each of the datasets = number of grid points.

                 Datasets:
0                    ID         V
0                    Xj         D	# Cartesian components of the point location, 1 col per jth dimension.
0                    IS_SINK    I	# =True iff the point lies on the edge of the model.
1                    NUMNEIGH   U
1                    FIRST_NN   V	# See explanation in section 3 below.
2                    VELj       D	# 1 col per jth dimension.
3                    DENSITYn   D	# 1 per nth collision partner.
4                    DENSMOLm   E	# 1 per mth molecular species. (Note that we will allow reading of ABUNMOLm for backwards compatibility.)
5                    TURBDPLR   E	# Given Gaussian lineshape exp(-v^2/[B^2 + 2*k*T/m]), this is B.
6                    TEMPKNTC   E	# From t[0].
6                    TEMPDUST   E	# From t[1].

1        Group 'NN_INDICES' (see explanation in the header to gridio.c)
1            Group 'columns'
                 Number of elements in each of the datasets = number of grid points * average number of Delaunay links per point.

                 Datasets:
1                    LINK_I     V

1        Group 'LINKS' (see explanation in the header to gridio.c)
1            Group 'columns'
                 Number of elements in each of the datasets = number of Delaunay links.

                 Datasets:
1                    GRID_I_1   V
1                    GRID_I_2   V
7                    V_p_j      D	# 1 per pth velocity sample per jth dimension.

8        Group 'LEVEL_POPS_m' (1 per mth molecular species)
             Attributes:
8                MOL_NAME       A
8            Dataset 'array'    E
             dimensions = [(number of energy levels this species),(number of grid cells)].


Note that at present the data in the 'partner' element of grid.mol is *NOT* being stored.
*/

/*....................................................................*/
hid_t
openHDF5FileForWrite(char *outFileName){
  hid_t fptr;

  fptr = H5Fcreate(outFileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  return fptr;
}

/*....................................................................*/
void
closeHDF5File(hid_t file){
  herr_t status=0;

  (void)status; // just to stop compiler warnings because this return value is currently unused.

  status = H5Fclose(file);
//*** check status
}

/*....................................................................*/
void
_convertKwdDataTypesToHDF5(struct keywordType *kwds, const int numKeywords){
  int i;
  char message[80];

  for(i=0;i<numKeywords;i++){
    if(     kwds[i].datatype==lime_CHAR)
      kwds[i].datatype = H5T_NATIVE_CHAR;
    else if(kwds[i].datatype==lime_INT)
      kwds[i].datatype = H5T_NATIVE_INT;
    else if(kwds[i].datatype==lime_FLOAT)
      kwds[i].datatype = H5T_NATIVE_FLOAT;
    else if(kwds[i].datatype==lime_DOUBLE)
      kwds[i].datatype = H5T_NATIVE_DOUBLE;
    else{
      if(!silent){
        sprintf(message, "Keyword %d dataype %d is not currently accepted.", i, kwds[i].datatype);
        bail_out(message);
      }
      exit(1);
    }
//*** check status
  }

}

/*....................................................................*/
void
writeKeywordsToHDF5(hid_t parent, struct keywordType *kwds\
  , const int numKeywords){

  hid_t kwdSpace,datatype=0,kwdAttr;
  herr_t status=0;
  int i;
  char message[80];

  (void)status; // just to stop compiler warnings because this return value is currently unused.

  _convertKwdDataTypesToHDF5(kwds, numKeywords);

  kwdSpace = H5Screate(H5S_SCALAR);

  for(i=0;i<numKeywords;i++){
    if(kwds[i].datatype==H5T_NATIVE_CHAR){
      datatype = H5Tcopy(H5T_C_S1);
      status = H5Tset_size(datatype, 1+strlen(kwds[i].charValue));
//*** test status
      status = H5Tset_strpad(datatype, H5T_STR_NULLTERM);
//*** test status
      kwdAttr = H5Acreate(parent, kwds[i].keyname, datatype, kwdSpace, H5P_DEFAULT, H5P_DEFAULT);
    }else
      kwdAttr = H5Acreate(parent, kwds[i].keyname, kwds[i].datatype, kwdSpace, H5P_DEFAULT, H5P_DEFAULT);

    if(     kwds[i].datatype==H5T_NATIVE_CHAR)
      status = H5Awrite(kwdAttr, datatype, kwds[i].charValue); /* datatype is correct here, **NOT** kwds[i].datatype */
    else if(kwds[i].datatype==H5T_NATIVE_INT)
      status = H5Awrite(kwdAttr, kwds[i].datatype, &kwds[i].intValue); 
    else if(kwds[i].datatype==H5T_NATIVE_FLOAT)
      status = H5Awrite(kwdAttr, kwds[i].datatype, &kwds[i].floatValue); 
    else if(kwds[i].datatype==H5T_NATIVE_DOUBLE)
      status = H5Awrite(kwdAttr, kwds[i].datatype, &kwds[i].doubleValue); 
    else{
      if(!silent){
        sprintf(message, "Keyword %d dataype %d is not currently accepted.", i, kwds[i].datatype);
        bail_out(message);
      }
      exit(1);
    }
//*** test status

    if(kwds[i].datatype==H5T_NATIVE_CHAR){
      status = H5Tclose(datatype);
//*** test status
    }

    /* Close attribute. */
    status = H5Aclose(kwdAttr);
//*** test status
  }

  /* Close attribute dataspace. */
  status = H5Sclose(kwdSpace); 
//*** test status

}

/*....................................................................*/
void
_defineAndLoadColumns(struct gridInfoType gridInfo\
  , const int dataFlags, const unsigned short numColNameChars, char ***allColNames, int **allColNumbers\
  , int *maxNumCols, int *numValidCols, int **colDataTypes, char ***colUnits){
  /*
To understand what is going on in the present function one needs to be aware that there is a set of all possible columns, corresponding to the set of quasi-scalar elements of the grid struct; and also a subset of that (the 'columns to write'), which is the columns for which data currently exists in the grid struct. The latter set is defined by the bits set in dataFlags.

The function returns 4 vectors: allColNames, allColNumbers and colDataTypes, described as follows:

      allColNames: pretty self-explanatory.

      allColNumbers: this is a vector of size equal to the total number of possible columns, but its values are indices in the range {1,...,numColsToWrite}. If one of the possible columns is not scheduled to be written (because the appropriate bit of dataFlags was not set), then the matching value of allColNumbers will be 0. The number of columns scheduled to be written is therefore the number of non-zero elements of allColNumbers.

      colDataTypes: this contains data types, but NOT for all columns, just for those to be written.

      colUnits: this contains unit strings, but NOT for all columns, just for those to be written.

NOTES:
  - The calling routine needs to free allColNames, allColNumbers and colDataTypes after it is finshed with them.
  - Data types in all capitals are defined in fitsio.h.
  */

  const unsigned short maxNumDims=9, maxNumSpecies=9, maxNumDensities=9;
  unsigned short i_us;
  int colI,i,colToWriteI;
  char message[80];
  char **unitAllCols=NULL;
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
  unitAllCols     = malloc(sizeof(*unitAllCols)    *(*maxNumCols));
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
  unitAllCols[colI] = "\0";
  dataTypeAllCols[colI] = H5T_NATIVE_UINT;

  /* should rather have a vector column? */
  for(i_us=0;i_us<gridInfo.nDims;i_us++){
    colI++;
    sprintf((*allColNames)[colI], "X%d", (int)i_us+1);
    if(bitIsSet(dataFlags, DS_bit_x)){
      colToWriteI++;
      (*allColNumbers)[colI] = colToWriteI;
    }
    unitAllCols[colI] = "m";
    dataTypeAllCols[colI] = H5T_NATIVE_DOUBLE;
  }

  colI++;
  sprintf((*allColNames)[colI], "IS_SINK");
  if(bitIsSet(dataFlags, DS_bit_x)){
    colToWriteI++;
    (*allColNumbers)[colI] = colToWriteI;
  }
  unitAllCols[colI] = "\0";
  dataTypeAllCols[colI] = H5T_NATIVE_CHAR;//H5T_NATIVE_HBOOL?;

  colI++;
  sprintf((*allColNames)[colI], "NUMNEIGH");
  if(bitIsSet(dataFlags, DS_bit_neighbours)){
    colToWriteI++;
    (*allColNumbers)[colI] = colToWriteI;
  }
  unitAllCols[colI] = "\0";
  dataTypeAllCols[colI] = H5T_NATIVE_USHORT;

  colI++;
  sprintf((*allColNames)[colI], "FIRST_NN");
  if(bitIsSet(dataFlags, DS_bit_neighbours)){
    colToWriteI++;
    (*allColNumbers)[colI] = colToWriteI;
  }
  unitAllCols[colI] = "\0";
  dataTypeAllCols[colI] = H5T_NATIVE_UINT;

  /* should rather have a vector column? */
  for(i_us=0;i_us<gridInfo.nDims;i_us++){
    colI++;
    sprintf((*allColNames)[colI], "VEL%d", (int)i_us+1);
    if(bitIsSet(dataFlags, DS_bit_velocity)){
      colToWriteI++;
      (*allColNumbers)[colI] = colToWriteI;
    }
    unitAllCols[colI] = "m/s";
    dataTypeAllCols[colI] = H5T_NATIVE_DOUBLE;
  }

  /* should rather have a vector column? */
  for(i_us=0;i_us<gridInfo.nDensities;i_us++){
    colI++;
    sprintf((*allColNames)[colI], "DENSITY%d", (int)i_us+1);
    if(bitIsSet(dataFlags, DS_bit_density)){
      colToWriteI++;
      (*allColNumbers)[colI] = colToWriteI;
    }
    unitAllCols[colI] = "kg/m^3";
    dataTypeAllCols[colI] = H5T_NATIVE_DOUBLE;
  }

  /* should rather have a vector column? */
  for(i_us=0;i_us<gridInfo.nSpecies;i_us++){
    colI++;
    sprintf((*allColNames)[colI], "DENSMOL%d", (int)i_us+1);
    if(bitIsSet(dataFlags, DS_bit_abundance)){
      colToWriteI++;
      (*allColNumbers)[colI] = colToWriteI;
    }
    unitAllCols[colI] = "\0";
    dataTypeAllCols[colI] = H5T_NATIVE_FLOAT;
  }

  colI++;
  sprintf((*allColNames)[colI], "TURBDPLR");
  if(bitIsSet(dataFlags, DS_bit_turb_doppler)){
    colToWriteI++;
    (*allColNumbers)[colI] = colToWriteI;
  }
  unitAllCols[colI] = "m/s";
  dataTypeAllCols[colI] = H5T_NATIVE_FLOAT;

  colI++;
  sprintf((*allColNames)[colI], "TEMPKNTC");
  if(bitIsSet(dataFlags, DS_bit_temperatures)){
    colToWriteI++;
    (*allColNumbers)[colI] = colToWriteI;
  }
  unitAllCols[colI] = "K";
  dataTypeAllCols[colI] = H5T_NATIVE_FLOAT;

  colI++;
  sprintf((*allColNames)[colI], "TEMPDUST");
  if(bitIsSet(dataFlags, DS_bit_temperatures)){
    colToWriteI++;
    (*allColNumbers)[colI] = colToWriteI;
  }
  unitAllCols[colI] = "K";
  dataTypeAllCols[colI] = H5T_NATIVE_FLOAT;

  /* should rather have a vector column? */
  for(i=0;i<3;i++){/* **** should rather loop to gridInfo.nDims but only entre here if it ==3?? */
    colI++;
    sprintf((*allColNames)[colI], "B_FIELD%d", i+1);
    if(bitIsSet(dataFlags, DS_bit_magfield)){
      colToWriteI++;
      (*allColNumbers)[colI] = colToWriteI;
    }
    unitAllCols[colI] = "T";
    dataTypeAllCols[colI] = H5T_NATIVE_FLOAT;
  }

  /* Define the datatype and physical units for the columns which will actually be written.
  */
  *colDataTypes = malloc(sizeof(**colDataTypes)*colToWriteI);
  *colUnits     = malloc(sizeof(**colUnits    )*colToWriteI);
  for(i=0;i<colToWriteI;i++){
    (*colUnits)[i] = malloc(sizeof(char)*STRLEN_KCHAR);
  }

  for(colI=0;colI<(*maxNumCols);colI++){
    i = (*allColNumbers)[colI];
    if(i>0){
      (*colDataTypes)[i-1] = dataTypeAllCols[colI];
      strcpy((*colUnits)[i-1], unitAllCols[colI]);
    }
  }

  *numValidCols = colToWriteI;

  free(unitAllCols);
  free(dataTypeAllCols);
}

/*....................................................................*/
int
_getColIndex(char **allColNames, const int maxNumCols, char *colName){
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
_writeKwdsToHDF5Col(hid_t dset, const int colI, char *colName, char *colUnit){
  int numKwds,i;
  struct keywordType *kwds=NULL;

  numKwds = 4;
  kwds = malloc(sizeof(*kwds)*numKwds);

  i = 0;
  initializeKeyword(&kwds[i]);
  kwds[i].datatype = lime_CHAR;
  sprintf(kwds[i].keyname, "CLASS");
  sprintf(kwds[i].charValue, "COLUMN");

  i++;
  initializeKeyword(&kwds[i]);
  kwds[i].datatype = lime_CHAR;
  sprintf(kwds[i].keyname, "UNIT");
  strcpy(kwds[i].charValue, colUnit);

  i++;
  initializeKeyword(&kwds[i]);
  kwds[i].datatype = lime_INT;
  sprintf(kwds[i].keyname, "POSITION");
  kwds[i].intValue = colI;

  i++;
  initializeKeyword(&kwds[i]);
  kwds[i].datatype = lime_CHAR;
  sprintf(kwds[i].keyname, "COL_NAME");
  strcpy(kwds[i].charValue, colName);

  writeKeywordsToHDF5(dset, kwds, numKwds);
  freeKeywords(kwds, numKwds);
}

/*....................................................................*/
void
_setUpHDF5Column(hid_t dataGroup, char *colName, const unsigned int numEntries\
  , const int dataType, hid_t *space, hid_t *dset){

  hsize_t dims[1] = {(hsize_t)numEntries};

  (*space) = H5Screate_simple(1, dims, NULL);
  (*dset) = H5Dcreate(dataGroup, colName, dataType, (*space), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
}

/*....................................................................*/
herr_t
_closeHDF5Column(hid_t dset, hid_t space, const int colI, char *colName, char *colUnit){

  herr_t status = 0;

  _writeKwdsToHDF5Col(dset, colI, colName, colUnit);

  status = H5Dclose(dset);
//*** check status?
  status = H5Sclose(space);
//*** check status?

  return status;
}

/*....................................................................*/
herr_t
_writeColumnToHDF5_ui(hid_t dataGroup, const int colI, char *colName\
  , char *colUnit, const unsigned int numEntries, unsigned int *colValues){

  herr_t status = 0;
  hid_t space,dset;

  _setUpHDF5Column(dataGroup, colName, numEntries, H5T_NATIVE_UINT, &space, &dset);

  status = H5Dwrite(dset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, colValues);
//*** check status?

  status = _closeHDF5Column(dset, space, colI, colName, colUnit);
//*** check status?

  return status;
}

/*....................................................................*/
herr_t
_writeColumnToHDF5_double(hid_t dataGroup, const int colI, char *colName\
  , char *colUnit, const unsigned int numEntries, double *colValues){

  herr_t status = 0;
  hid_t space,dset;

  _setUpHDF5Column(dataGroup, colName, numEntries, H5T_NATIVE_DOUBLE, &space, &dset);

  status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, colValues);
//*** check status?

  status = _closeHDF5Column(dset, space, colI, colName, colUnit);
//*** check status?

  return status;
}

/*....................................................................*/
herr_t
_writeColumnToHDF5_bool(hid_t dataGroup, const int colI, char *colName\
  , char *colUnit, const unsigned int numEntries, _Bool *colValues){

  herr_t status = 0;
  hid_t space,dset;
  short *shortValues=NULL;
  unsigned int ui;

  shortValues = malloc(sizeof(*shortValues)*numEntries);
  for(ui=0;ui<numEntries;ui++){
    if(colValues[ui])
      shortValues[ui] = 1;
    else
      shortValues[ui] = 0;
  }

  _setUpHDF5Column(dataGroup, colName, numEntries, H5T_NATIVE_SHORT, &space, &dset);

  status = H5Dwrite(dset, H5T_NATIVE_SHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, shortValues);
//*** check status?

  status = _closeHDF5Column(dset, space, colI, colName, colUnit);
//*** check status?

  free(shortValues);

  return status;
}

/*....................................................................*/
herr_t
_writeColumnToHDF5_us(hid_t dataGroup, const int colI, char *colName\
  , char *colUnit, const unsigned int numEntries, unsigned short *colValues){

  herr_t status = 0;
  hid_t space,dset;

  _setUpHDF5Column(dataGroup, colName, numEntries, H5T_NATIVE_USHORT, &space, &dset);

  status = H5Dwrite(dset, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, colValues);
//*** check status?

  status = _closeHDF5Column(dset, space, colI, colName, colUnit);
//*** check status?

  return status;
}

/*....................................................................*/
herr_t
_writeColumnToHDF5_float(hid_t dataGroup, const int colI, char *colName\
  , char *colUnit, const unsigned int numEntries, float *colValues){

  herr_t status = 0;
  hid_t space,dset;

  _setUpHDF5Column(dataGroup, colName, numEntries, H5T_NATIVE_FLOAT, &space, &dset);

  status = H5Dwrite(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, colValues);
//*** check status?

  status = _closeHDF5Column(dset, space, colI, colName, colUnit);
//*** check status?

  return status;
}

/*....................................................................*/
void
writeGridExtToHDF5(hid_t file, struct gridInfoType gridInfo\
  , struct grid *gp, unsigned int *firstNearNeigh\
  , char **collPartNames, const int dataFlags){
  /*
This writes whatever information is in the grid struct (as specified by the dataFlags) and which also has a dimensionality which is some simple multiple of the number of grid points, to a single FITS binary table extension called GRID. The function tries to be fairly forgiving of screwy situations but it will exit if the minimum information is not present (defined as allBitsSet(dataFlags, DS_mask_x), which implies that elements .id, .x and .sink should all contain valid values).

Note that data types in all capitals are defined in fitsio.h.
  */

  const unsigned int totalNumGridPoints = gridInfo.nInternalPoints+gridInfo.nSinkPoints;
  const unsigned short maxNumCollPart = 9;
  char extname[] = "GRID";
  unsigned int *ids=NULL,i_ui;
  double *xj=NULL;
  _Bool *sink=NULL;/* Probably should be char* but this seems to work. */
  unsigned short *numNeigh=NULL,i_us,localNumCollPart;
  double *velj=NULL,*densn=NULL;
  float *dopb=NULL,*t=NULL,*densm=NULL,*bField=NULL;
  int colI=0,i,di,maxNumCols,numValidCols,numKwds;
  char message[80];
  char colName[STRLEN_KCHAR];
  char **allColNames=NULL,**colUnits=NULL;
  int *allColNumbers=NULL,*colDataTypes=NULL;
  hid_t hduGroup,dataGroup;
  herr_t status=0;
  struct keywordType *kwds=NULL;

  (void)status; /* just to stop compiler warnings because this return value is currently unused. */

  if(!allBitsSet(dataFlags, DS_mask_x)){
    if(!silent) bail_out("Data stage indicates no grid data!");
    exit(1);
  }

  /*
Ok we have a bit of a tricky situation here in that the number of columns we write is going to depend on the information available in gp, as encoded in the dataFlags. We need to work out which columns we are going to write ahead of time because we need the appropriate data on ALL the columns to set up the table size before we can start to write their individual data values. I also want to avoid checking dataFlags twice in two different contexts - that is how errors arise. So I've arranged that the following routine will do all the donkey work of setting up only those columns we can write and then using that information to define the table size. The routine also returns four more vectors:

	allColNames   - This is returned to remove what would otherwise be a hard-wired dependence that the column ordering was the same in the present routine as in defineAndLoadColumns(). With this vector, the present routine can search for a column name in it and then use the returned vector index to access (from the next vector) the number of the column in the (smaller) sequence of valid columns.

	allColNumbers - This contains the number of a column in the sequence (beginning at 1) of those for which gp has data. If gp contains no data for a given column name, its entry in allColNumbers will be 0.

	colDataTypes  - This just contains the data types for the valid columns.

	colUnits      - This just contains the unit strings for the valid columns.
  */
  _defineAndLoadColumns(gridInfo, dataFlags, STRLEN_KCHAR, &allColNames\
     , &allColNumbers, &maxNumCols, &numValidCols, &colDataTypes, &colUnits);

  hduGroup = H5Gcreate(file, extname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dataGroup = H5Gcreate(hduGroup, "columns", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  /* Write the columns:
  */
  colI = allColNumbers[_getColIndex(allColNames, maxNumCols, "ID")];
  if(colI<=0){
    if(!silent) bail_out("This should not occur, it is some sort of bug.");
    exit(1);
  }
  ids = malloc(sizeof(*ids)*totalNumGridPoints);
  for(i_ui=0;i_ui<totalNumGridPoints;i_ui++)
    ids[i_ui] = (unsigned int)gp[i_ui].id;
  status = _writeColumnToHDF5_ui(dataGroup, colI, "ID", colUnits[colI-1], totalNumGridPoints, ids);
  free(ids);

  xj = malloc(sizeof(*xj)*totalNumGridPoints);
  for(i_us=0;i_us<gridInfo.nDims;i_us++){
    sprintf(colName, "X%d", (int)i_us+1);
    colI = allColNumbers[_getColIndex(allColNames, maxNumCols, colName)];
    if(colI<=0){
      if(!silent) bail_out("This should not occur, it is some sort of bug.");
      exit(1);
    }

    for(i_ui=0;i_ui<totalNumGridPoints;i_ui++)
      xj[i_ui] = gp[i_ui].x[i_us];
    status = _writeColumnToHDF5_double(dataGroup, colI, colName, colUnits[colI-1], totalNumGridPoints, xj);
  }
  free(xj);

  colI = allColNumbers[_getColIndex(allColNames, maxNumCols, "IS_SINK")];
  if(colI<=0){
    if(!silent) bail_out("This should not occur, it is some sort of bug.");
    exit(1);
  }
  sink = malloc(sizeof(*sink)*totalNumGridPoints);
  for(i_ui=0;i_ui<totalNumGridPoints;i_ui++)
    sink[i_ui] = (_Bool)gp[i_ui].sink;
  status = _writeColumnToHDF5_bool(dataGroup, colI, "IS_SINK", colUnits[colI-1], totalNumGridPoints, sink);
  free(sink);

  colI = allColNumbers[_getColIndex(allColNames, maxNumCols, "NUMNEIGH")];
  if(colI>0){
    numNeigh = malloc(sizeof(*numNeigh)*totalNumGridPoints);
    for(i_ui=0;i_ui<totalNumGridPoints;i_ui++)
      numNeigh[i_ui] = (unsigned short)gp[i_ui].numNeigh;
    status = _writeColumnToHDF5_us(dataGroup, colI, "NUMNEIGH", colUnits[colI-1], totalNumGridPoints, numNeigh);
    free(numNeigh);

    colI = allColNumbers[_getColIndex(allColNames, maxNumCols, "FIRST_NN")];
    if(colI<=0){
      if(!silent) bail_out("This should not occur, it is some sort of bug.");
      exit(1);
    }
    status = _writeColumnToHDF5_ui(dataGroup, colI, "FIRST_NN", colUnits[colI-1], totalNumGridPoints, firstNearNeigh);
  }

  /* Check if first VEL column has info:
  */
  colI = allColNumbers[_getColIndex(allColNames, maxNumCols, "VEL1")];
  if(colI>0){
    velj = malloc(sizeof(*velj)*totalNumGridPoints);
    for(i_us=0;i_us<gridInfo.nDims;i_us++){
      sprintf(colName, "VEL%d", (int)i_us+1);
      colI = allColNumbers[_getColIndex(allColNames, maxNumCols, colName)];
      if(colI<=0){
        if(!silent) bail_out("This should not occur, it is some sort of bug.");
        exit(1);
      }

      for(i_ui=0;i_ui<totalNumGridPoints;i_ui++)
        velj[i_ui] = gp[i_ui].vel[i_us];
      status = _writeColumnToHDF5_double(dataGroup, colI, colName, colUnits[colI-1], totalNumGridPoints, velj);
    }
    free(velj);
  }

  /* Check if first DENSITY column has info:
  */
  colI = allColNumbers[_getColIndex(allColNames, maxNumCols, "DENSITY1")];
  if(colI>0){
    densn = malloc(sizeof(*densn)*totalNumGridPoints);
    for(i_us=0;i_us<gridInfo.nDensities;i_us++){
      sprintf(colName, "DENSITY%d", (int)i_us+1);
      colI = allColNumbers[_getColIndex(allColNames, maxNumCols, colName)];
      if(colI<=0){
        if(!silent) bail_out("This should not occur, it is some sort of bug.");
        exit(1);
      }

      for(i_ui=0;i_ui<totalNumGridPoints;i_ui++)
        densn[i_ui] = gp[i_ui].dens[i_us];
      status = _writeColumnToHDF5_double(dataGroup, colI, colName, colUnits[colI-1], totalNumGridPoints, densn);
    }
    free(densn);
  }

  /* Check if first DENSMOL column has info:
  */
  colI = allColNumbers[_getColIndex(allColNames, maxNumCols, "DENSMOL1")];
  if(colI>0){
    densm = malloc(sizeof(*densm)*totalNumGridPoints);
    for(i_us=0;i_us<gridInfo.nSpecies;i_us++){
      sprintf(colName, "DENSMOL%d", (int)i_us+1);
      colI = allColNumbers[_getColIndex(allColNames, maxNumCols, colName)];
      if(colI<=0){
        if(!silent) bail_out("This should not occur, it is some sort of bug.");
        exit(1);
      }

      for(i_ui=0;i_ui<totalNumGridPoints;i_ui++)
        densm[i_ui] = (float)gp[i_ui].mol[i_us].nmol;
      status = _writeColumnToHDF5_float(dataGroup, colI, colName, colUnits[colI-1], totalNumGridPoints, densm);
    }
    free(densm);
  }

  colI = allColNumbers[_getColIndex(allColNames, maxNumCols, "TURBDPLR")];
  if(colI>0){
    dopb = malloc(sizeof(*dopb)*totalNumGridPoints);
    for(i_ui=0;i_ui<totalNumGridPoints;i_ui++)
      dopb[i_ui] = (float)gp[i_ui].dopb_turb;
    status = _writeColumnToHDF5_float(dataGroup, colI, "TURBDPLR", colUnits[colI-1], totalNumGridPoints, dopb);
    free(dopb);
  }

  colI = allColNumbers[_getColIndex(allColNames, maxNumCols, "TEMPKNTC")];
  if(colI>0){
    t = malloc(sizeof(*t)*totalNumGridPoints);
    for(i_ui=0;i_ui<totalNumGridPoints;i_ui++)
      t[i_ui] = (float)gp[i_ui].t[0];
    status = _writeColumnToHDF5_float(dataGroup, colI, "TEMPKNTC", colUnits[colI-1], totalNumGridPoints, t);

    colI = allColNumbers[_getColIndex(allColNames, maxNumCols, "TEMPDUST")];
    if(colI<=0){
      if(!silent) bail_out("This should not occur, it is some sort of bug.");
      exit(1);
    }

    for(i_ui=0;i_ui<totalNumGridPoints;i_ui++)
      t[i_ui] = (float)gp[i_ui].t[1];
    status = _writeColumnToHDF5_float(dataGroup, colI, "TEMPDUST", colUnits[colI-1], totalNumGridPoints, t);
    free(t);
  }

  /* Check if first B_FIELD column has info:
  */
  colI = allColNumbers[_getColIndex(allColNames, maxNumCols, "B_FIELD1")];
  if(colI>0){
    bField = malloc(sizeof(*bField)*totalNumGridPoints);
    for(di=0;di<3;di++){//**** keep this hard-wired or rather test that gridInfo.nDims==3??
      sprintf(colName, "B_FIELD%d", di+1);
      colI = allColNumbers[_getColIndex(allColNames, maxNumCols, colName)];
      if(colI<=0){
        if(!silent) bail_out("This should not occur, it is some sort of bug.");
        exit(1);
      }

      for(i_ui=0;i_ui<totalNumGridPoints;i_ui++)
        bField[i_ui] = (float)gp[i_ui].B[di];
      status = _writeColumnToHDF5_float(dataGroup, colI, colName, colUnits[colI-1], totalNumGridPoints, bField);
    }
    free(bField);
  }

  /* Write keywords.
  */
  numKwds = 1;
  kwds = malloc(sizeof(*kwds)*numKwds);

  i = 0;
  initializeKeyword(&kwds[i]);
  kwds[i].datatype = lime_CHAR;
  sprintf(kwds[i].keyname, "CLASS");
  sprintf(kwds[i].charValue, "DATA_GROUP");

  writeKeywordsToHDF5(dataGroup, kwds, numKwds);
  freeKeywords(kwds, numKwds);

  status = H5Gclose(dataGroup);
//*** check status?

  localNumCollPart = 0; /* Default. */
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
  }

  numKwds = localNumCollPart + 3;
  kwds = malloc(sizeof(*kwds)*numKwds);
  for(i=0;i<numKwds;i++)
    initializeKeyword(&kwds[i]);

  i = 0;
  kwds[i].datatype = lime_CHAR;
  sprintf(kwds[i].keyname, "CLASS");
  sprintf(kwds[i].charValue, "HDU");

  i++;
  kwds[i].datatype = lime_CHAR;
  sprintf(kwds[i].keyname, "EXTNAME");
  strcpy(kwds[i].charValue, extname);

  i++;
  kwds[i].datatype = lime_INT;
  sprintf(kwds[i].keyname, "HDUNUM");
  kwds[i].intValue = 0;

  for(i_us=0;i_us<localNumCollPart;i_us++){
    i = i_us+3;
    kwds[i].datatype = lime_CHAR;
    sprintf(kwds[i].keyname, "COLLPAR%d", (int)i_us+1);
    strcpy(kwds[i].charValue, collPartNames[i_us]);
  }

  writeKeywordsToHDF5(hduGroup, kwds, numKwds);
  freeKeywords(kwds, numKwds);

  status = H5Gclose(hduGroup);
//*** check status?

  for(i=0;i<maxNumCols;i++) free(allColNames[i]);
  free(allColNames);
  for(i=0;i<numValidCols;i++) free(colUnits[i]);
  free(colUnits);
  free(allColNumbers);
  free(colDataTypes);
}

/*....................................................................*/
void
writeNnIndicesExtToHDF5(hid_t file, struct gridInfoType gridInfo\
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
  int i,numKwds;
  char extname[] = "NN_INDICES";
  hid_t hduGroup,dataGroup;
  herr_t status=0;
  struct keywordType *kwds=NULL;

  (void)status; /* just to stop compiler warnings because this return value is currently unused. */

  if (nnLinks==NULL){
    if(!silent) bail_out("No link or near-neighbour data!");
    exit(1);
  }

  /* Append a new empty binary table onto the FITS file.
  */
  hduGroup = H5Gcreate(file, extname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dataGroup = H5Gcreate(hduGroup, "columns", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  linkIs = malloc(sizeof(*linkIs)*gridInfo.nNNIndices);
  for(i_ui=0;i_ui<gridInfo.nNNIndices;i_ui++)
    linkIs[i_ui] = nnLinks[i_ui]->id;
  status = _writeColumnToHDF5_ui(dataGroup, 0, "LINK_I", "\0", gridInfo.nNNIndices, linkIs);
  free(linkIs);

  numKwds = 1;
  kwds = malloc(sizeof(*kwds)*numKwds);

  i = 0;
  initializeKeyword(&kwds[i]);
  kwds[i].datatype = lime_CHAR;
  sprintf(kwds[i].keyname, "CLASS");
  sprintf(kwds[i].charValue, "DATA_GROUP");

  writeKeywordsToHDF5(dataGroup, kwds, numKwds);
  freeKeywords(kwds, numKwds);

  status = H5Gclose(dataGroup);
//*** check status?

  numKwds = 3;
  kwds = malloc(sizeof(*kwds)*numKwds);

  i = 0;
  initializeKeyword(&kwds[i]);
  kwds[i].datatype = lime_CHAR;
  sprintf(kwds[i].keyname, "CLASS");
  sprintf(kwds[i].charValue, "HDU");

  i++;
  initializeKeyword(&kwds[i]);
  kwds[i].datatype = lime_CHAR;
  sprintf(kwds[i].keyname, "EXTNAME");
  strcpy(kwds[i].charValue, extname);

  i++;
  initializeKeyword(&kwds[i]);
  kwds[i].datatype = lime_INT;
  sprintf(kwds[i].keyname, "HDUNUM");
  kwds[i].intValue = 1;

  writeKeywordsToHDF5(hduGroup, kwds, numKwds);
  freeKeywords(kwds, numKwds);

  status = H5Gclose(hduGroup);
//*** check status?
}

/*....................................................................*/
void
writeLinksExtToHDF5(hid_t file, struct gridInfoType gridInfo\
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
  unsigned int *ids=NULL,i_ui;
  unsigned short i_us,j_us;
  double *vels=NULL;
  int colI=0,i,numKwds;
  int numCols;
  char extname[] = "LINKS";
  char **tempColNames=NULL;
  hid_t hduGroup,dataGroup;
  herr_t status=0;
  struct keywordType *kwds=NULL;

  (void)status; /* just to stop compiler warnings because this return value is currently unused. */

  if(links==NULL){
    if(!silent) bail_out("No link data!");
    exit(1);
  }

  if(links[0].vels==NULL)
    numCols = 2;
  else
    numCols = 2 + gridInfo.nDims*gridInfo.nLinkVels;

  tempColNames = malloc(sizeof(*tempColNames)*numCols);
  for(i=0;i<numCols;i++) tempColNames[i]=malloc(sizeof(char)*STRLEN_KCHAR);

  colI = 0;
  sprintf(tempColNames[colI], "GRID_I_1");

  colI++;
  sprintf(tempColNames[colI], "GRID_I_2");

  if(links[0].vels!=NULL){
    /* Should rather have vector columns? */
    for(i_us=0;i_us<gridInfo.nLinkVels;i_us++){
      for(j_us=0;j_us<gridInfo.nDims;j_us++){
        colI++;
        sprintf(tempColNames[colI], "V_%d_%d", (int)i_us+1, (int)j_us+1);
      }
    }
  }

  /* Append a new empty binary table onto the FITS file.
  */
  hduGroup = H5Gcreate(file, extname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dataGroup = H5Gcreate(hduGroup, "columns", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  /* Write columns.
  */
  colI = 0;
  ids = malloc(sizeof(*ids)*gridInfo.nLinks);
  for(i_ui=0;i_ui<gridInfo.nLinks;i_ui++)
    ids[i_ui] = links[i_ui].gis[0];
  status = _writeColumnToHDF5_ui(dataGroup, colI, tempColNames[colI], "\0", gridInfo.nLinks, ids);

  colI++;
  for(i_ui=0;i_ui<gridInfo.nLinks;i_ui++)
    ids[i_ui] = links[i_ui].gis[1];
  status = _writeColumnToHDF5_ui(dataGroup, colI, tempColNames[colI], "\0", gridInfo.nLinks, ids);
  free(ids);

  if(links[0].vels!=NULL){
    vels = malloc(sizeof(*vels)*gridInfo.nLinks);
    for(i_us=0;i_us<gridInfo.nLinkVels;i_us++){
      for(j_us=0;j_us<gridInfo.nDims;j_us++){
        colI++;
        for(i_ui=0;i_ui<gridInfo.nLinks;i_ui++)
          vels[i_ui] = links[i_ui].vels[gridInfo.nDims*i_us + j_us];
        status = _writeColumnToHDF5_double(dataGroup, colI, tempColNames[colI], "m/s", gridInfo.nLinks, vels);
      }
    }
    free(vels);
  }

  for(i=0;i<numCols;i++)
    free(tempColNames[i]);
  free(tempColNames);

  numKwds = 1;
  kwds = malloc(sizeof(*kwds)*numKwds);

  i = 0;
  initializeKeyword(&kwds[i]);
  kwds[i].datatype = lime_CHAR;
  sprintf(kwds[i].keyname, "CLASS");
  sprintf(kwds[i].charValue, "DATA_GROUP");

  writeKeywordsToHDF5(dataGroup, kwds, numKwds);
  freeKeywords(kwds, numKwds);

  status = H5Gclose(dataGroup);
//*** check status?

  numKwds = 3;
  kwds = malloc(sizeof(*kwds)*numKwds);

  i = 0;
  initializeKeyword(&kwds[i]);
  kwds[i].datatype = lime_CHAR;
  sprintf(kwds[i].keyname, "CLASS");
  sprintf(kwds[i].charValue, "HDU");

  i++;
  initializeKeyword(&kwds[i]);
  kwds[i].datatype = lime_CHAR;
  sprintf(kwds[i].keyname, "EXTNAME");
  strcpy(kwds[i].charValue, extname);

  i++;
  initializeKeyword(&kwds[i]);
  kwds[i].datatype = lime_INT;
  sprintf(kwds[i].keyname, "HDUNUM");
  kwds[i].intValue = 2;

  writeKeywordsToHDF5(hduGroup, kwds, numKwds);
  freeKeywords(kwds, numKwds);

  status = H5Gclose(hduGroup);
//*** check status?
}

/*....................................................................*/
void
writePopsGroupToHDF5(hid_t file, struct gridInfoType gridInfo\
  , const unsigned short speciesI, struct grid *gp){
  /*
	Extension name: LEVEL_POPS_m (1 per molecular species)
	This is an image extension of size (number of grid cells)*(number of energy levels this species)

	Keywords:
		MOL_NAME

Note that data types in all capitals are defined in fitsio.h.
  */

  const unsigned int totalNumGridPoints = gridInfo.nInternalPoints+gridInfo.nSinkPoints;
  const int naxis = 2;  /* i.e. 2-dimensional image */    
  char groupName[20];
  int numEnergyLevels = (int)gridInfo.mols[speciesI].nLevels;
  hid_t hduGroup,fileSpace,dset,memSpace;
  hsize_t naxes[2] = {(hsize_t)numEnergyLevels, (hsize_t)totalNumGridPoints};
  hsize_t memBufDim[] = {(hsize_t)totalNumGridPoints};
  hsize_t memStart[1],memStride[1],memCount[1],memBlock[1];
  hsize_t fileStart[2],fileStride[2],fileCount[2],fileBlock[2];
  herr_t status=0;
  float *row=NULL;
  int xi,numKwds,i;
  unsigned int i_ui;
  struct keywordType *kwds=NULL;

  (void)status; /* just to stop compiler warnings because this return value is currently unused. */

  sprintf(groupName, "LEVEL_POPS_%d", (int)speciesI+1);

  hduGroup = H5Gcreate(file, groupName, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  fileSpace = H5Screate_simple(naxis, naxes, NULL);
  dset = H5Dcreate(hduGroup, "array", H5T_NATIVE_FLOAT, fileSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  memSpace = H5Screate_simple(1, memBufDim, NULL);
  memStart[0]  = 0;
  memStride[0] = 1;
  memCount[0]  = (hsize_t)totalNumGridPoints;
  memBlock[0]  = 1;
  status = H5Sselect_hyperslab(memSpace, H5S_SELECT_SET, memStart, memStride, memCount, memBlock);
//*** check status?

  fileStart[1]  = 0;
  fileStride[0] = 1; fileStride[1] = 1;
  fileCount[0]  = 1; fileCount[1]  = (hsize_t)totalNumGridPoints;    
  fileBlock[0]  = 1; fileBlock[1]  = 1;

  row = malloc(sizeof(*row)*totalNumGridPoints);

  /* Write data.
  */
  for(xi=0;xi<numEnergyLevels;xi++){
    for(i_ui=0;i_ui<totalNumGridPoints;i_ui++)
      row[i_ui] = (float)gp[i_ui].mol[speciesI].pops[xi]; 

    fileStart[0] = (hsize_t)xi;
    status = H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, fileStart, fileStride, fileCount, fileBlock);
//*** check status?
    status = H5Dwrite(dset, H5T_NATIVE_FLOAT, memSpace, fileSpace, H5P_DEFAULT, row);
//*** check status?
  }

  free(row);

  status = H5Sclose(memSpace);
//*** check status?
  status = H5Dclose(dset);
//*** check status?
  status = H5Sclose(fileSpace);
//*** check status?

  numKwds = 4;
  kwds = malloc(sizeof(*kwds)*numKwds);

  i = 0;
  initializeKeyword(&kwds[i]);
  kwds[i].datatype = lime_CHAR;
  sprintf(kwds[i].keyname, "CLASS");
  sprintf(kwds[i].charValue, "HDU");

  i++;
  initializeKeyword(&kwds[i]);
  kwds[i].datatype = lime_CHAR;
  sprintf(kwds[i].keyname, "EXTNAME");
  sprintf(kwds[i].charValue, "LEVEL_POPS_%d", (int)speciesI+1);

  i++;
  initializeKeyword(&kwds[i]);
  kwds[i].datatype = lime_INT;
  sprintf(kwds[i].keyname, "HDUNUM");
  kwds[i].intValue = 3 + speciesI;

  i++;
  initializeKeyword(&kwds[i]);
  kwds[i].datatype = lime_CHAR;
  sprintf(kwds[i].keyname, "MOL_NAME");
  strcpy(kwds[i].charValue, gridInfo.mols[speciesI].molName);

  writeKeywordsToHDF5(hduGroup, kwds, numKwds);
  freeKeywords(kwds, numKwds);

  status = H5Gclose(hduGroup);
//*** check status?
}

/*....................................................................*/
hid_t
openHDF5FileForRead(char *inFileName){
  hid_t fptr;

  fptr = H5Fopen(inFileName, H5F_ACC_RDWR, H5P_DEFAULT);

  return fptr;
}

/*....................................................................*/
int
countDataSetNamePlusInt(hid_t parent, char *baseName){
  const int maxI = 999;
  char message[80],datasetName[20];
  int i;
  htri_t myi;

  i = 0;
  while(i<maxI){
    sprintf(datasetName, "%s%d", baseName, i+1);
    myi = H5Lexists(parent, datasetName, H5P_DEFAULT);
    if((int)myi<=0)
  break;

    i++;
  }

  if(i>=maxI){
    if(!silent){
      sprintf(message, "Max num datasets %d exceeded.", maxI);
      bail_out(message);
    }
    exit(1);
  }

  return i;
}

/*....................................................................*/
int
countAttributeNamePlusInt(hid_t parent, char *baseName){
  const int maxI = 999;
  char message[80],attrName[20];
  int i;
  htri_t myi;

  i = 0;
  while(i<maxI){
    sprintf(attrName, "%s%d", baseName, i+1);
    myi = H5Aexists(parent, attrName);
    if((int)myi<=0)
  break;

    i++;
  }

  if(i>=maxI){
    if(!silent){
      sprintf(message, "Max num attributes %d exceeded.", maxI);
      bail_out(message);
    }
    exit(1);
  }

  return i;
}

/*....................................................................*/
void
readKeywordsFromHDF5(hid_t parent, struct keywordType *kwds\
  , const int numKeywords){

  int i,status=0;
  hid_t kwdAttr,datatype;
  hsize_t lenStr;
  char message[80];

  (void)status; /* just to stop compiler warnings because this return value is currently unused. */

  for(i=0;i<numKeywords;i++){
    kwdAttr = H5Aopen(parent, kwds[i].keyname, H5P_DEFAULT);

    if(     kwds[i].datatype==lime_CHAR){
      datatype = H5Tcopy(H5T_C_S1);
      lenStr = H5Aget_storage_size(kwdAttr);
      status = H5Tset_size(datatype, lenStr);
//*** check status
      status = H5Aread(kwdAttr, datatype, kwds[i].charValue);
//*** check status
      status = H5Tclose(datatype);
//*** test status
    }else if(kwds[i].datatype==lime_INT){
      status = H5Aread(kwdAttr, H5T_NATIVE_INT, &kwds[i].intValue);
//*** check status
    }else if(kwds[i].datatype==lime_FLOAT){
      status = H5Aread(kwdAttr, H5T_NATIVE_FLOAT, &kwds[i].floatValue);
//*** check status
    }else if(kwds[i].datatype==lime_DOUBLE){
      status = H5Aread(kwdAttr, H5T_NATIVE_DOUBLE, &kwds[i].doubleValue);
//*** check status
    }else{
      if(!silent){
        sprintf(message, "Keyword %d datatype %d is not currently accepted.", i, kwds[i].datatype);
        bail_out(message);
      }
      exit(1);
    }

    kwds[i].comment = NULL;

    /* Close attribute. */
    status = H5Aclose(kwdAttr);
//*** test status
  }
}

/*....................................................................*/
void
readGridExtFromHDF5(hid_t file, struct gridInfoType *gridInfoRead\
  , struct grid **gp, unsigned int **firstNearNeigh, char ***collPartNames\
  , int *numCollPartRead, int *dataFlags, _Bool *densMolColsExists){
  /*
The present function mallocs 'gp' and sets defaults for all the simple or first-level struct elements.

If COLLPARn attributes are found in the GRID group then collPartNames is malloc'd to the number of these which have 'n' values forming a consecutive sequence increasing from 1.

Note that the calling routine needs to free gp, firstNearNeigh and collPartNames after use.
  */

  hid_t hduGroup,dataGroup,dset,space,attr,atype,atype_mem;
  int i,numSpaceDims,numBFieldCols;
  hsize_t spaceDims[1],maxSpaceDims[1];
  herr_t status = 0;
  unsigned int numGridCells,i_ui;
  unsigned int *ids=NULL;
  char message[80],genericKwd[10];
  double *xj=NULL;
  char datasetName[20];
  short *sink=NULL;
  unsigned short *numNeigh=NULL, i_us;
  double *velj=NULL,*densn=NULL;
  float *abunm=NULL,*densm=NULL,*dopb=NULL,*t=NULL,*bField=NULL;
  htri_t myi;

  (void)status; /* just to stop compiler warnings because this return value is currently unused. */

  hduGroup = H5Gopen(file, "GRID", H5P_DEFAULT);
  dataGroup = H5Gopen(hduGroup, "columns", H5P_DEFAULT);

  /* Find out how many rows there are, then malloc the array.
  */
  dset = H5Dopen(dataGroup, "ID", H5P_DEFAULT);
//*** test and error if return <0.
  space = H5Dget_space(dset);
  numSpaceDims = H5Sget_simple_extent_dims(space, spaceDims, maxSpaceDims);
  status = H5Sclose(space);
//*** check status?

  if(numSpaceDims<0){
    if(!silent) bail_out("ID column numSpaceDims<0");
    exit(1);
  }

  if(numSpaceDims!=1){
    if(!silent){
      sprintf(message, "ID column: expected numSpaceDims==1 but got %d", numSpaceDims);
      bail_out(message);
    }
    exit(1);
  }

  numGridCells = (unsigned int)spaceDims[0];
  if(numGridCells<=0){
    if(!silent) warning("No rows found in grid dataset.");

    status = H5Dclose(dset);
//*** check status?
    status = H5Gclose(dataGroup);
//*** check status?
    status = H5Gclose(hduGroup);
//*** check status?

    return; /* I.e. with dataFlags left unchanged. */
  }

  /* Find out if the user has supplied ABUNMOLn or DENSMOLn columns.
  */
  *densMolColsExists = FALSE;
  dset = H5Dopen(dataGroup, "DENSMOL1", H5P_DEFAULT);
  if(dset==0)
    *densMolColsExists = TRUE;
  status = H5Dclose(dset);

  /* Count the numbers of ABUNMOLn/DENSMOLn columns to get the number of species:
  */
  if(*densMolColsExists)
    gridInfoRead->nSpecies = (unsigned short)countDataSetNamePlusInt(dataGroup, "DENSMOL");
  else
    gridInfoRead->nSpecies = (unsigned short)countDataSetNamePlusInt(dataGroup, "ABUNMOL");
  mallocAndSetDefaultGrid(gp, (size_t)numGridCells, gridInfoRead->nSpecies);

  /* Read the columns.
  */
  ids = malloc(sizeof(*ids)*numGridCells);
  status = H5Dread(dset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ids);
//*** check status?
  status = H5Dclose(dset);
//*** check status?

  for(i_ui=0;i_ui<numGridCells;i_ui++)
    (*gp)[i_ui].id = (int)ids[i_ui];
  free(ids);

  gridInfoRead->nDims = (unsigned short)countDataSetNamePlusInt(dataGroup, "X");
  if(gridInfoRead->nDims<=0){
    if(!silent) warning("No X columns found in grid dataset.");
    status = H5Gclose(dataGroup);
//*** check status?
    status = H5Gclose(hduGroup);
//*** check status?
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
    sprintf(datasetName, "X%d", (int)i_us+1);
    dset = H5Dopen(dataGroup, datasetName, H5P_DEFAULT);
    status = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, xj);
//*** check status?
    status = H5Dclose(dset);
//*** check status?

    for(i_ui=0;i_ui<numGridCells;i_ui++) {
      (*gp)[i_ui].x[i_us] = xj[i_ui];
    }
  }
  free(xj);

  dset = H5Dopen(dataGroup, "IS_SINK", H5P_DEFAULT);
  if(dset<0){
    if(!silent) warning("No IS_SINK column found in grid dataset.");

    status = H5Dclose(dset);
//*** check status?
    status = H5Gclose(dataGroup);
//*** check status?
    status = H5Gclose(hduGroup);
//*** check status?

    return; /* I.e. with dataFlags left unchanged. */
  }

  sink = malloc(sizeof(*sink)*numGridCells);
  status = H5Dread(dset, H5T_NATIVE_SHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, sink);
//*** check status?
  status = H5Dclose(dset);
//*** check status?

  gridInfoRead->nSinkPoints = 0;
  for(i_ui=0;i_ui<numGridCells;i_ui++) {
    (*gp)[i_ui].sink = (int)sink[i_ui];
    if((*gp)[i_ui].sink)
      gridInfoRead->nSinkPoints++;
  }
  free(sink);

  gridInfoRead->nInternalPoints = numGridCells - gridInfoRead->nSinkPoints;

  /* If we have made it this far, we can set the first bit of dataFlags. Woot!
  */
  (*dataFlags) |= (1 << DS_bit_x);

  myi = H5Lexists(dataGroup, "NUMNEIGH", H5P_DEFAULT);
  if((int)myi>0){
    dset = H5Dopen(dataGroup, "NUMNEIGH", H5P_DEFAULT);
    numNeigh = malloc(sizeof(*numNeigh)*numGridCells);
    status = H5Dread(dset, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, numNeigh);
//*** check status?
    status = H5Dclose(dset);
//*** check status?

    for(i_ui=0;i_ui<numGridCells;i_ui++) {
      (*gp)[i_ui].numNeigh = (int)numNeigh[i_ui];
    }
    free(numNeigh);

    myi = H5Lexists(dataGroup, "FIRST_NN", H5P_DEFAULT);
    if((int)myi>0){
      dset = H5Dopen(dataGroup, "FIRST_NN", H5P_DEFAULT);

      *firstNearNeigh = malloc(sizeof(**firstNearNeigh)*numGridCells);
      status = H5Dread(dset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, *firstNearNeigh);
//*** check status?
      status = H5Dclose(dset);
//*** check status?

      /* If we made it to here, we can set the neighbour bit of dataFlags. Note however that this bit is only on trial til we check for the LINKS extension. So it had better behave itself.
      */
      (*dataFlags) |= (1 << DS_bit_neighbours);
    }
  }

  /* See if there are any VEL columns:
  */
  if(countDataSetNamePlusInt(dataGroup, "VEL")>0){
    /* Read the VEL columns:
    */
    velj = malloc(sizeof(*velj)*numGridCells);
    for(i_us=0;i_us<gridInfoRead->nDims;i_us++){
      sprintf(datasetName, "VEL%d", (int)i_us+1);
      dset = H5Dopen(dataGroup, datasetName, H5P_DEFAULT);

      status = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, velj);
//*** check status?
      status = H5Dclose(dset);
//*** check status?

      for(i_ui=0;i_ui<numGridCells;i_ui++) {
        (*gp)[i_ui].vel[i_us] = velj[i_ui];
      }
    }
    free(velj);

    (*dataFlags) |= (1 << DS_bit_velocity);
  }

  /* Count the numbers of DENSITYn columns:
  */
  gridInfoRead->nDensities = (unsigned short)countDataSetNamePlusInt(dataGroup, "DENSITY");
  if(gridInfoRead->nDensities > 0){
    for(i_ui=0;i_ui<numGridCells;i_ui++)
      (*gp)[i_ui].dens = malloc(sizeof(double)*gridInfoRead->nDensities);

    /* Read the DENSITY columns:
    */
    densn = malloc(sizeof(*densn)*numGridCells);
    for(i_us=0;i_us<gridInfoRead->nDensities;i_us++){
      sprintf(datasetName, "DENSITY%d", (int)i_us+1);
      dset = H5Dopen(dataGroup, datasetName, H5P_DEFAULT);

      status = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, densn);
//*** check status?
      status = H5Dclose(dset);
//*** check status?

      for(i_ui=0;i_ui<numGridCells;i_ui++) {
        (*gp)[i_ui].dens[i_us] = densn[i_ui];
      }
    }
    free(densn);

    (*dataFlags) |= (1 << DS_bit_density);
  }

  if(gridInfoRead->nSpecies > 0){
    if(*densMolColsExists){
      /* Read the DENSMOL columns:
      */
      densm = malloc(sizeof(*densm)*numGridCells);
      for(i_us=0;i_us<gridInfoRead->nSpecies;i_us++){
        sprintf(datasetName, "DENSMOL%d", (int)i_us+1);
        dset = H5Dopen(dataGroup, datasetName, H5P_DEFAULT);

        status = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, densm);
//*** check status?
        status = H5Dclose(dset);
//*** check status?

        for(i_ui=0;i_ui<numGridCells;i_ui++) {
          (*gp)[i_ui].mol[i_us].nmol = (double)densm[i_ui];
        }
      }
      free(densm);
    }else{
      /* Read the ABUNMOL columns:
      */
      abunm = malloc(sizeof(*abunm)*numGridCells);
      for(i_us=0;i_us<gridInfoRead->nSpecies;i_us++){
        sprintf(datasetName, "ABUNMOL%d", (int)i_us+1);
        dset = H5Dopen(dataGroup, datasetName, H5P_DEFAULT);

        status = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, abunm);
//*** check status?
        status = H5Dclose(dset);
//*** check status?

        for(i_ui=0;i_ui<numGridCells;i_ui++) {
          (*gp)[i_ui].mol[i_us].abun = (double)abunm[i_ui];
        }
      }
      free(abunm);
    }

//    par->useAbun = !densMolColsExists;
    (*dataFlags) |= (1 << DS_bit_abundance);
  }

  /* Read the TURBDPLR column:
  */
  myi = H5Lexists(dataGroup, "TURBDPLR", H5P_DEFAULT);
  if((int)myi>0){
    dset = H5Dopen(dataGroup, "TURBDPLR", H5P_DEFAULT);

    dopb = malloc(sizeof(*dopb)*numGridCells);
    status = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, dopb);
//*** check status?
    status = H5Dclose(dset);
//*** check status?

    for(i_ui=0;i_ui<numGridCells;i_ui++) {
      (*gp)[i_ui].dopb_turb = (double)dopb[i_ui];
    }
    free(dopb);

    (*dataFlags) |= (1 << DS_bit_turb_doppler);
  }

//*** there is probably a neater way to do the temperatures. Resetting the t0 defaults is a bit ugly.
  /* Read the TEMPKNTC column:
  */
  myi = H5Lexists(dataGroup, "TEMPKNTC", H5P_DEFAULT);
  if((int)myi>0){
    dset = H5Dopen(dataGroup, "TEMPKNTC", H5P_DEFAULT);

    t = malloc(sizeof(*t)*numGridCells);
    status = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, t);
//*** check status?
    status = H5Dclose(dset);
//*** check status?

    for(i_ui=0;i_ui<numGridCells;i_ui++) {
      (*gp)[i_ui].t[0] = (double)t[i_ui];
    }

    /* Read the TEMPDUST column:
    */
    myi = H5Lexists(dataGroup, "TEMPDUST", H5P_DEFAULT);
    if((int)myi<=0){
      for(i_ui=0;i_ui<numGridCells;i_ui++) {
        (*gp)[i_ui].t[0] = -1;
      }
    }else{
      dset = H5Dopen(dataGroup, "TEMPDUST", H5P_DEFAULT);

      status = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, t);
//*** check status?
      status = H5Dclose(dset);
//*** check status?

      for(i_ui=0;i_ui<numGridCells;i_ui++) {
        (*gp)[i_ui].t[1] = (double)t[i_ui];
      }

      (*dataFlags) |= (1 << DS_bit_temperatures);
    }
    free(t);
  }

  /* See if there are any B_FIELD columns:
  */
  numBFieldCols = (int)countDataSetNamePlusInt(dataGroup, "B_FIELD");
  if(numBFieldCols>0){
    if(numBFieldCols!=3){
      if(!silent){
        sprintf(message, "Wrong number (%d) of B_FIELD columns - there should be %d.", numBFieldCols, 3);
        bail_out(message);
      }
      exit(1);
    }

    /* Read the B_FIELD columns:
    */
    bField = malloc(sizeof(*bField)*numGridCells);
    for(i=0;i<3;i++){//**** keep this hard-wired or rather test that gridInfo.nDims==3??
      sprintf(datasetName, "B_FIELD%d", i+1);
      dset = H5Dopen(dataGroup, datasetName, H5P_DEFAULT);
      status = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, bField);
//*** check status?
      status = H5Dclose(dset);
//*** check status?

      for(i_ui=0;i_ui<numGridCells;i_ui++) {
        (*gp)[i_ui].B[i] = (double)bField[i_ui];
      }
    }
    free(bField);

    (*dataFlags) |= (1 << DS_bit_magfield);
  }

  status = H5Gclose(dataGroup);
//*** check status?

  /* Check if there are any COLLPAR keywords.
  */
  *numCollPartRead = countAttributeNamePlusInt(hduGroup, "COLLPAR");
  if(*numCollPartRead <= 0){
    *collPartNames = NULL;
  }else{
    *collPartNames = malloc(sizeof(**collPartNames)*(*numCollPartRead));
    for(i=0;i<(*numCollPartRead);i++){
      sprintf(genericKwd, "COLLPAR%d", i+1);
      (*collPartNames)[i] = malloc(sizeof(char)*100);
      attr = H5Aopen(hduGroup, genericKwd, H5P_DEFAULT);
      atype = H5Aget_type(attr);
      atype_mem = H5Tget_native_type(atype, H5T_DIR_ASCEND);
      status = H5Aread(attr, atype_mem, (*collPartNames)[i]);
//*** check status?
      status = H5Tclose(atype_mem);
//*** check status?
      status = H5Tclose(atype);
//*** check status?
      status = H5Aclose(attr);
//*** check status?
    }
  }

  status = H5Gclose(hduGroup);
//*** check status?

}

/*....................................................................*/
void
readLinksExtFromHDF5(hid_t file, struct gridInfoType *gridInfoRead\
  , struct grid *gp, struct linkType **links, int *dataFlags){
  /*
See the comment at the beginning of gridio.c for a description of how the LINKS extension relates to the grid struct.

The present function mallocs the pointer *links.
  */

  unsigned int totalNumGridPoints,i_ui,j_ui,*ids=NULL;
  hid_t hduGroup,dataGroup,dset,space;
  herr_t status = 0;
  int i,numSpaceDims;
  hsize_t spaceDims[1],maxSpaceDims[1];
  char message[80],datasetName[22];
  double *vels=NULL;
  unsigned short i_us,j_us;

  (void)status; // just to stop compiler warnings because this return value is currently unused.

  if(!bitIsSet(*dataFlags, DS_bit_neighbours))
    return;

  totalNumGridPoints = gridInfoRead->nInternalPoints + gridInfoRead->nSinkPoints; /* Just for a bit more brevity. */

  hduGroup = H5Gopen(file, "LINKS", H5P_DEFAULT);
  if(hduGroup<0){
    if(!silent) warning("No LINKS group found in grid dataset.");

    /* Unset the DS_mask_neighbours bit: */
    *dataFlags &= ~(1 << DS_bit_neighbours);

    status = H5Gclose(hduGroup);
//*** check status?

    return;
  }
  dataGroup = H5Gopen(hduGroup, "columns", H5P_DEFAULT);

  /* Find out how many rows there are, then malloc the array.
  */
  dset = H5Dopen(dataGroup, "GRID_I_1", H5P_DEFAULT);
//*** test and error if return <0. ********* other H5Dopen too.
  space = H5Dget_space(dset);
  numSpaceDims = H5Sget_simple_extent_dims(space, spaceDims, maxSpaceDims);
  status = H5Sclose(space);
//*** check status?

  if(numSpaceDims<0){
    if(!silent) bail_out("GRID_I_1 column numSpaceDims<0");
    exit(1);
  }

  if(numSpaceDims!=1){
    if(!silent){
      sprintf(message, "GRID_I_1 column: expected numSpaceDims==1 but got %d", numSpaceDims);
      bail_out(message);
    }
    exit(1);
  }

  gridInfoRead->nLinks = (unsigned int)spaceDims[0];
/*  free(spaceDims); */
  if(gridInfoRead->nLinks<=0){
    if(!silent) warning("No rows in LINKS group of grid dataset.");

    /* Unset the neighbours bit: */
    *dataFlags &= ~(1 << DS_bit_neighbours);

    status = H5Gclose(dset);
//*** check status?
    status = H5Gclose(dataGroup);
//*** check status?
    status = H5Gclose(hduGroup);
//*** check status?

    return;
  }
  *links = malloc(sizeof(**links)*gridInfoRead->nLinks);

  for(i_ui=0;i_ui<gridInfoRead->nLinks;i_ui++) {
    (*links)[i_ui].id = i_ui;
  }

  /* Read GRID_I_1 column.
  */
  ids = malloc(sizeof(*ids)*gridInfoRead->nLinks);
  status = H5Dread(dset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ids);
//*** check status?
  status = H5Dclose(dset);
//*** check status?

  for(i_ui=0;i_ui<gridInfoRead->nLinks;i_ui++) {
    j_ui = ids[i_ui];
    if(j_ui<0 || j_ui>=totalNumGridPoints){
      if(!silent){
        sprintf(message, "GRID_I_1 %udth-row value %ud is outside range [0,%ud]", i_ui, j_ui, totalNumGridPoints);
        bail_out(message);
      }
      exit(1);
    }
    (*links)[i_ui].gis[0] = gp[j_ui].id;
  }

  /* Read GRID_I_2 column.
  */
  dset = H5Dopen(dataGroup, "GRID_I_2", H5P_DEFAULT);
//*** test and error if return <0.
  status = H5Dread(dset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ids);
//*** check status?
  status = H5Dclose(dset);
//*** check status?

  for(i_ui=0;i_ui<gridInfoRead->nLinks;i_ui++) {
    j_ui = ids[i_ui];
    if(j_ui<0 || j_ui>=totalNumGridPoints){
      if(!silent){
        sprintf(message, "GRID_I_2 %udth-row value %ud is outside range [0,%ud]", i_ui, j_ui, totalNumGridPoints);
        bail_out(message);
      }
      exit(1);
    }
    (*links)[i_ui].gis[1] = gp[j_ui].id;
  }
  free(ids);

  /* Find out how many V_* columns there are.
  */
  i = 0;
  status = 0;
  while(!status){
    sprintf(datasetName, "V_%d_", i+1);
    if((unsigned short)countDataSetNamePlusInt(dataGroup, datasetName)!=gridInfoRead->nDims)
      status = 1;

    i++;
  }
  gridInfoRead->nLinkVels = (unsigned short)(i - 1);
  status = 0;

  if(gridInfoRead->nLinkVels<=0){
    for(i_ui=0;i_ui<gridInfoRead->nLinks;i_ui++)
      (*links)[i_ui].vels = NULL;

    status = H5Gclose(dataGroup);
//*** check status?
    status = H5Gclose(hduGroup);
//*** check status?

    return;
  }

  for(i_ui=0;i_ui<gridInfoRead->nLinks;i_ui++)
    (*links)[i_ui].vels = malloc(sizeof(double)*gridInfoRead->nLinkVels*gridInfoRead->nDims);

  vels = malloc(sizeof(*vels)*gridInfoRead->nLinks);
  for(i_us=0;i_us<gridInfoRead->nLinkVels;i_us++){
    for(j_us=0;j_us<gridInfoRead->nDims;j_us++){
      /* Read the V_n_d columns.
      */
      sprintf(datasetName, "V_%d_%d", (int)i_us+1, (int)j_us+1);
      dset = H5Dopen(dataGroup, datasetName, H5P_DEFAULT);
//*** test and error if return <0.
      status = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vels);
//*** check status?
      status = H5Dclose(dset);
//*** check status?

      for(i_ui=0;i_ui<gridInfoRead->nLinks;i_ui++)
        (*links)[i_ui].vels[gridInfoRead->nDims*i_us + j_us] = vels[i_ui];
    }
  }
  free(vels);

  (*dataFlags) |= (1 << DS_bit_ACOEFF);

  status = H5Gclose(dataGroup);
//*** check status?
  status = H5Gclose(hduGroup);
//*** check status?
}

/*....................................................................*/
void
readNnIndicesExtFromHDF5(hid_t file, struct linkType *links\
  , struct linkType ***nnLinks, struct gridInfoType *gridInfoRead, int *dataFlags){
  /*
See the comment at the beginning of gridio.c for a description of how the NN_INDICES extension relates to the grid struct.

The function mallocs the pointer *nnLinks.
  */

  hid_t hduGroup,dataGroup,dset,space;
  herr_t status = 0;
  int totalNumNeigh,i,numSpaceDims;
  hsize_t spaceDims[1],maxSpaceDims[1];
  unsigned int *linkIs=NULL;
  char message[80];

  (void)status; // just to stop compiler warnings because this return value is currently unused.

  if(!bitIsSet(*dataFlags, DS_bit_neighbours))
    return;

  /* Go to the NN_INDICES extension.
  */
  hduGroup = H5Gopen(file, "NN_INDICES", H5P_DEFAULT);
  if(hduGroup<0){
    if(!silent) warning("No NN_INDICES group found in grid dataset.");

    /* Unset the DS_mask_neighbours bit: */
    *dataFlags &= ~(1 << DS_bit_neighbours);

    status = H5Gclose(hduGroup);
//*** check status?

    return;
  }
  dataGroup = H5Gopen(hduGroup, "columns", H5P_DEFAULT);

  /* Find out how many rows there are, then malloc the array.
  */
  dset = H5Dopen(dataGroup, "LINK_I", H5P_DEFAULT);
//*** test and error if return <0.
  space = H5Dget_space(dset);
  numSpaceDims = H5Sget_simple_extent_dims(space, spaceDims, maxSpaceDims);
  status = H5Sclose(space);
//*** check status?

  if(numSpaceDims<0){
    if(!silent) bail_out("LINK_I column numSpaceDims<0");
    exit(1);
  }

  if(numSpaceDims!=1){
    if(!silent){
      sprintf(message, "LINK_I column: expected numSpaceDims==1 but got %d", numSpaceDims);
      bail_out(message);
    }
    exit(1);
  }

  totalNumNeigh = spaceDims[0];
  gridInfoRead->nNNIndices = (unsigned int)totalNumNeigh;
  if(gridInfoRead->nNNIndices<=0){
    if(!silent) warning("No rows in NN_INDICES group of grid dataset.");

    /* Unset the neighbours bit: */
    *dataFlags &= ~(1 << DS_bit_neighbours);

    status = H5Gclose(dset);
//*** check status?
    status = H5Gclose(dataGroup);
//*** check status?
    status = H5Gclose(hduGroup);
//*** check status?

    return;
  }
  *nnLinks = malloc(sizeof(**nnLinks)*totalNumNeigh);

  /* Read LINK_I column.
  */
  linkIs = malloc(sizeof(*linkIs)*totalNumNeigh);
  status = H5Dread(dset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, linkIs);
//*** check status?
  status = H5Dclose(dset);
//*** check status?

  for(i=0;i<totalNumNeigh;i++)
    (*nnLinks)[i] = &links[linkIs[i]];

  free(linkIs);

  status = H5Gclose(dataGroup);
//*** check status?
  status = H5Gclose(hduGroup);
//*** check status?
}

/*....................................................................*/
_Bool
checkPopsHDF5GroupExists(hid_t file, const unsigned short speciesI){
  const unsigned short maxNumSpecies = 9;
  char message[80];
  char groupName[14];
  int returnStatus=0;
  htri_t myi;

  if(speciesI+1>maxNumSpecies){
    if(!silent){
      sprintf(message, "Species block %d is greater than the limit %d", (int)speciesI+1, (int)maxNumSpecies);
      bail_out(message);
    }
    exit(1);
  }

  sprintf(groupName, "LEVEL_POPS_%d", (int)speciesI+1);

  /* Try to move to extension LEVEL_POPS_<speciesI>:
  */
  myi = H5Lexists(file, groupName, H5P_DEFAULT);
  if((int)myi<=0)
    returnStatus = 0;
  else
    returnStatus = 1;

  return returnStatus;
}

/*....................................................................*/
void
readPopsGroupFromHDF5(hid_t file, const unsigned short speciesI\
  , struct grid *gp, struct gridInfoType *gridInfoRead){
  /*
See the comment at the beginning of the present module for a description of how the LEVEL_POPS_m extensions relate to the grid struct.

The function mallocs gp[yi].mol[speciesI].pops for each grid point yi and species.
  */

  const int maxLenMolName = 8;
  const unsigned short maxNumSpecies = 9;
  char message[80],groupName[14];
  char molNameRead[maxLenMolName+1];
  hid_t hduGroup,dset,memSpace,fileSpace,attr,atype,atype_mem;
  int numEnergyLevels;
  hsize_t memBufDim[1];
  hsize_t memStart[1],memStride[1],memCount[1],memBlock[1];
  hsize_t fileStart[2],fileStride[2],fileCount[2],fileBlock[2];
  herr_t status = 0;
  int xi,numSpaceDims;
  hsize_t spaceDims[2],maxSpaceDims[2];
  unsigned int totalNumGridPoints,i_ui;
  float *row=NULL;

  (void)status; // just to stop compiler warnings because this return value is currently unused.

  if(speciesI+1>maxNumSpecies){
    if(!silent){
      sprintf(message, "Species block %d is greater than the limit %d", (int)speciesI+1, (int)maxNumSpecies);
      bail_out(message);
    }
    exit(1);
  }
  sprintf(groupName, "LEVEL_POPS_%d", (int)speciesI+1);

  /* Try to move to extension LEVEL_POPS_<speciesI>:
  */
  hduGroup = H5Gopen(file, groupName, H5P_DEFAULT);
  dset = H5Dopen(hduGroup, "array", H5P_DEFAULT);

  fileSpace = H5Dget_space(dset);
  numSpaceDims = H5Sget_simple_extent_dims(fileSpace, spaceDims, maxSpaceDims);

  if(numSpaceDims<0){
    if(!silent){
      sprintf(message, "%s/array numSpaceDims<0", groupName);
      bail_out(message);
    }
    exit(1);
  }

  if(numSpaceDims!=2){
    if(!silent){
      sprintf(message, "%s/array: expected numSpaceDims==2 but got %d", groupName, numSpaceDims);
      bail_out(message);
    }
    exit(1);
  }

  totalNumGridPoints = gridInfoRead->nInternalPoints + gridInfoRead->nSinkPoints;
  if((hsize_t)totalNumGridPoints != spaceDims[1]){
    if(!silent){
      sprintf(message, "Expected %ld grid points but extension %s has %d"\
        , (long)totalNumGridPoints, groupName, (int)spaceDims[1]);
      bail_out(message);
    }
    exit(1);
  }
  gridInfoRead->mols[speciesI].nLevels = spaceDims[0];
  gridInfoRead->mols[speciesI].nLines = -1;
  numEnergyLevels = gridInfoRead->mols[speciesI].nLevels;

  memBufDim[0] = (hsize_t)totalNumGridPoints;

  memSpace = H5Screate_simple(1, memBufDim, NULL);
  memStart[0]  = 0;
  memStride[0] = 1;
  memCount[0]  = (hsize_t)totalNumGridPoints;
  memBlock[0]  = 1;
  status = H5Sselect_hyperslab(memSpace, H5S_SELECT_SET, memStart, memStride, memCount, memBlock);
//*** check status?

  fileStart[1]  = 0;
  fileStride[0] = 1; fileStride[1] = 1;
  fileCount[0]  = 1; fileCount[1]  = (hsize_t)totalNumGridPoints;    
  fileBlock[0]  = 1; fileBlock[1]  = 1;

  row = malloc(sizeof(*row)*totalNumGridPoints);

  for(i_ui=0;i_ui<totalNumGridPoints;i_ui++)
    gp[i_ui].mol[speciesI].pops = malloc(sizeof(double)*gridInfoRead->mols[speciesI].nLevels);

  /* Read data.
  */
  for(xi=0;xi<numEnergyLevels;xi++){
    fileStart[0] = (hsize_t)xi;
    status = H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, fileStart, fileStride, fileCount, fileBlock);
//*** check status?
    status = H5Dread(dset, H5T_NATIVE_FLOAT, memSpace, fileSpace, H5P_DEFAULT, row);
//*** check status?

    for(i_ui=0;i_ui<totalNumGridPoints;i_ui++)
      gp[i_ui].mol[speciesI].pops[xi] = (double)row[i_ui]; 
  }

  free(row);

  status = H5Sclose(memSpace);
//*** check status?
  status = H5Sclose(fileSpace);
//*** check status?
  status = H5Dclose(dset);
//*** check status?

  /* Read kwds:
  */
  attr = H5Aopen(hduGroup, "MOL_NAME", H5P_DEFAULT);
  atype = H5Aget_type(attr);
  atype_mem = H5Tget_native_type(atype, H5T_DIR_ASCEND);
  status = H5Aread(attr, atype_mem, molNameRead);

  gridInfoRead->mols[speciesI].molName = malloc(sizeof(char)*(STR_LEN_0+1));
  sprintf(gridInfoRead->mols[speciesI].molName, "%s", molNameRead);

//*** check status?
  status = H5Tclose(atype_mem);
//*** check status?
  status = H5Tclose(atype);
//*** check status?
  status = H5Aclose(attr);
//*** check status?

  status = H5Gclose(hduGroup);
//*** check status?
}


/*....................................................................*/
int
countDensityColsHDF5(char *inFileName){
  int numDensities;
  hid_t file,hduGroup,dataGroup;
  int status=0;

  (void)status; // just to stop compiler warnings because this return value is currently unused.

  file = openHDF5FileForRead(inFileName);
  hduGroup = H5Gopen(file, "GRID", H5P_DEFAULT);
  dataGroup = H5Gopen(hduGroup, "columns", H5P_DEFAULT);
  numDensities = countDataSetNamePlusInt(dataGroup, "DENSITY");
  status = H5Gclose(dataGroup);
//*** check status?
  status = H5Gclose(hduGroup);
//*** check status?
  closeHDF5File(file);

  return numDensities;
}

