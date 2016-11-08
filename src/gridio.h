/*
 *  gridio.h
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2016 The LIME development team
 *
 */

#ifndef GRIDIO_H
#define GRIDIO_H

#define lime_FITS	1

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

#include "grid2fits.h"

int	readGrid(char*, const int, struct gridInfoType*, struct keywordType*, const int, struct grid**, char***, int*, int*);
int	writeGrid(char*, const int, struct gridInfoType, struct keywordType*, const int, struct grid*, char**, const int);

#endif /* GRIDIO_H */

