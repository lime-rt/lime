/*
 *  grid2fits.h
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2017 The LIME development team
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

_Bool	checkPopsFITSExtExists(fitsfile*, const unsigned short);
void	closeFITSFile(fitsfile*);
void	readGridExtFromFITS(fitsfile*, const int, struct gridInfoType*, struct grid**, unsigned int**, char***, int*, int*);
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

