/*
 *  grid2fits.h
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2016 The LIME development team
 *
 */

#ifndef GRID2FITS_H
#define GRID2FITS_H

#define lime_fptr	fitsfile

struct keywordType{
  int datatype; /* Accepted: TSTRING, TINT, TFLOAT, TDOUBLE */
  char *keyname, *comment;
  int intValue;
  float floatValue;
  double doubleValue;
  char *charValue;
};


_Bool	checkPopsFitsExtExists(fitsfile*, const unsigned short);
void	closeFITSFile(fitsfile*);
void	initializeKeyword(struct keywordType*);
void	readGridExtFromFits(fitsfile*, struct gridInfoType*, struct grid**, unsigned int**, char***, int*, int*);
void	readKeywordsFromFits(lime_fptr*, struct keywordType*, const int);
void	readLinksExtFromFits(fitsfile*, struct gridInfoType*, struct grid*, struct linkType**, int*);
void	readNnIndicesExtFromFits(fitsfile*, struct linkType*, struct linkType***, struct gridInfoType*, int*);
void	readPopsExtFromFits(fitsfile*, const unsigned short, struct grid*, struct gridInfoType*);
fitsfile *openFITSFileForRead(char*);
fitsfile *openFITSFileForWrite(char*);
void	writeKeywordsToFits(lime_fptr*, struct keywordType*, const int);
void	writeGridExtToFits(fitsfile*, struct gridInfoType, struct grid*, unsigned int*, char**, const int);
void	writeLinksExtToFits(fitsfile*, struct gridInfoType, struct linkType*);
void	writeNnIndicesExtToFits(fitsfile*, struct gridInfoType, struct linkType**);//, struct linkType*);
void	writePopsExtToFits(fitsfile*, struct gridInfoType, const unsigned short, struct grid*);

#endif /* GRID2FITS_H */

