/*
 *  grid2fits.h
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *

NOTE! This file is not stand-alone, in needs to be included in an environment which defines (on top of the usual types) the following data types:
  fitsfile
  struct gridInfoType
  struct grid
  struct keywordType
  struct linkType

 */

#ifndef GRID2FITS_H
#define GRID2FITS_H

/* The following are required to be defined for any file format: */
#define STRLEN_KNAME	9
#define STRLEN_KCOMM	48
#define STRLEN_KCHAR	21
#define lime_fptr	fitsfile*
#define lime_init	NULL
#define _FAILED_TO_OPEN	==NULL

void	processFitsError(int);
_Bool	checkPopsFITSExtExists(fitsfile*, const unsigned short);
void	closeFITSFile(fitsfile*);
void	readGridExtFromFITS(fitsfile*, struct gridInfoType*, struct grid**, unsigned int**, char***, int*, int*, _Bool *densMolColsExists);
void	readKeywordsFromFITS(fitsfile*, struct keywordType*, const int);
void	readLinksExtFromFITS(fitsfile*, struct gridInfoType*, struct grid*, struct linkType**, int*);
void	readNnIndicesExtFromFITS(fitsfile*, struct linkType*, struct linkType***, struct gridInfoType*, int*);
void	readPopsExtFromFITS(fitsfile*, const unsigned short, struct grid*, struct gridInfoType*);
fitsfile *openFITSFileForRead(char*);
fitsfile *openFITSFileForWrite(char*);
void	writeKeywordsToFITS(fitsfile*, struct keywordType*, const int);
void	writeGridExtToFITS(fitsfile*, struct gridInfoType, struct grid*, unsigned int*, char**, const int);
void	writeLinksExtToFITS(fitsfile*, struct gridInfoType, struct linkType*);
void	writeNnIndicesExtToFITS(fitsfile*, struct gridInfoType, struct linkType**);//, struct linkType*);
void	writePopsExtToFITS(fitsfile*, struct gridInfoType, const unsigned short, struct grid*);
int	countDensityColsFITS(char *inFileName);

#endif /* GRID2FITS_H */

