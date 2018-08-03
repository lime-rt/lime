/*
 *  aux.h
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#ifndef AUX_H
#define AUX_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

void	reportInfAtOrigin(const double, const char*);
void	reportInfsAtOrigin(const int, const double*, const char*);
void	sigintHandler(int sigI);
void	checkFgets(char *fgetsResult, char *message);
void	checkFscanf(const int fscanfResult, const int expectedNum, char *message);
void	checkFread(const size_t freadResult, const size_t expectedNum, char *message);
void	checkFwrite(const size_t fwriteResult, const size_t expectedNum, char *message);
double	gaussline(const double, const double);
double	dotProduct3D(const double*, const double*);
void	copyInparStr(const char*, char**);
_Bool	charPtrIsNullOrEmpty(const char *inStr);
_Bool	allBitsSet(const int flags, const int mask);
_Bool	anyBitSet(const int flags, const int mask);
_Bool	bitIsSet(const int flags, const int bitI);
_Bool	onlyBitsSet(const int flags, const int mask);

#endif /* AUX_H */

