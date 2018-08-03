/*
 *  grid2hdf5.h
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#ifndef GRID2HDF5_H
#define GRID2HDF5_H

/* The following are required to be defined for any file format: */
#define STRLEN_KNAME	20
#define STRLEN_KCOMM	200
#define STRLEN_KCHAR	200
#define lime_fptr	hid_t
#define lime_init	0
#define _FAILED_TO_OPEN	<0

void
_defineAndLoadColumns(struct gridInfoType gridInfo\
  , const int dataFlags, const unsigned short numColNameChars, char ***allColNames, int **allColNumbers\
  , int *maxNumCols, int *numValidCols, int **colDataTypes, char ***colUnits);
herr_t
_writeColumnToHDF5_ui(hid_t dataGroup, const int colI, char *colName\
  , char *colUnit, const unsigned int numEntries, unsigned int *colValues);
int
_getColIndex(char **allColNames, const int maxNumCols, char *colName);


_Bool	checkPopsHDF5GroupExists(hid_t file, const unsigned short);
int	countDensityColsHDF5(char *inFileName);
void	closeHDF5File(hid_t file);
void	readGridExtFromHDF5(hid_t file, struct gridInfoType*, struct grid**, unsigned int**, char***, int*, int*, _Bool *densMolColsExists);
void	readKeywordsFromHDF5(hid_t parent, struct keywordType *kwds, const int numKeywords);
void	readLinksExtFromHDF5(hid_t file, struct gridInfoType*, struct grid*, struct linkType**, int*);
void	readNnIndicesExtFromHDF5(hid_t file, struct linkType*, struct linkType***, struct gridInfoType*, int*);
void	readPopsGroupFromHDF5(hid_t file, const unsigned short, struct grid*, struct gridInfoType*);
hid_t	openHDF5FileForRead(char*);
hid_t	openHDF5FileForWrite(char *outFileName);
void	writeKeywordsToHDF5(hid_t parent, struct keywordType *kwds, const int numKeywords);
void	writeGridExtToHDF5(hid_t file, struct gridInfoType, struct grid*, unsigned int*, char**, const int);
void	writeLinksExtToHDF5(hid_t file, struct gridInfoType, struct linkType*);
void	writeNnIndicesExtToHDF5(hid_t file, struct gridInfoType, struct linkType**);
void	writePopsGroupToHDF5(hid_t file, struct gridInfoType, const unsigned short, struct grid*);

#endif /* GRID2HDF5_H */

