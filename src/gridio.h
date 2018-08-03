/*
 *  gridio.h
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#ifndef GRIDIO_H
#define GRIDIO_H

#define lime_FITS	1
#define lime_HDF5	2

#ifdef USEHDF5
#include <hdf5.h>
#define lime_IO		lime_HDF5
#else
#define lime_IO		lime_FITS
#endif

#define lime_CHAR	0
#define lime_INT	1
#define lime_FLOAT	2
#define lime_DOUBLE	3
#define lime_BOOL	4

struct linkType {
  unsigned int id, gis[2];
  double *vels;
};

struct molInfoType{
  char *molName;
  int nLevels, nLines;
};

struct gridInfoType{
  unsigned int nInternalPoints, nSinkPoints, nLinks, nNNIndices;
  unsigned short nDims, nSpecies, nDensities, nLinkVels;
  struct molInfoType *mols;
};

struct keywordType{
  int datatype; /* Codes given above. */
  char *keyname,*comment,*charValue;
  int intValue;
  float floatValue;
  double doubleValue;
  _Bool boolValue;
};

#if defined(lime_IO) && lime_IO==lime_HDF5
  #include "grid2hdf5.h"
#else
  #include "grid2fits.h"
#endif

void	initializeKeyword(struct keywordType*);
void	freeGridInfo(struct gridInfoType*);
void	freeKeywords(struct keywordType*, const int);
int	readGrid(char*, struct gridInfoType*, struct keywordType*, const int, struct grid**, char***, int*, int*, _Bool *densMolColsExists);
int	writeGrid(char*, struct gridInfoType, struct keywordType*, const int, struct grid*, char**, const int);
int	countDensityCols(char*, int*);

#endif /* GRIDIO_H */

