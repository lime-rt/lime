/*
 *  pytypes.h
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2016 The LIME development team
 *
 */

#ifndef PYTYPES_H
#define PYTYPES_H

#include "lime.h"

typedef struct {
  char name[STR_LEN_0+1], type[STR_LEN_0+1];
  _Bool mandatory,isList;
} parTemplateType;

struct tempType{
  int intValue;
  double doubleValue;
  char strValue[STR_LEN_0+1];
  _Bool boolValue;
};

/* Defaults for user-specifiable functions */
void density_default(double,double,double,double *);
void temperature_default(double,double,double,double *);
void abundance_default(double,double,double,double *);
void doppler_default(double,double,double, double *);
void velocity_default(double,double,double,double *);
void magfield_default(double,double,double,double *);
void gasIIdust_default(double,double,double,double *);
double gridDensity_default(configInfo*, double*);

int	checkAttributes(PyObject *pParInstance, parTemplateType *parTemplate\
  , const int nPars);
_Bool	checkAttrType(PyObject *pObj, const char *parTemplateType);
void	extractValue(PyObject *pPars, const char *attrName, const char *attrType\
  , const int maxStrLen, struct tempType *tempValue);
int	extractValues(PyObject *pPars, const char *attrName, const char *attrType\
  , const int maxStrLen, struct tempType **tempValues);
int	getParTemplates(const char *headerModuleName, parTemplateType **parTemplate, int *nPars\
  , parTemplateType **imgParTemplate, int *nImgPars);
void	getPythonFunc(PyObject *pModule, const char *funcName, const int verbosity\
  , PyObject **pFunc);
int	initParImg(PyObject *pModule, PyObject *pMacros, parTemplateType *parTemplate\
  , const int nPars, parTemplateType *imgParTemplate, const int nImgPars\
  , inputPars *par, image **img, int *nImages);
int	myStrCpy(const char source[], char destination[], const int strlenDest);
void	printOrClearPyError(void);
void	pyerror(char *message);
void	pyToC(PyObject *pAttr, const char *attrType, const int maxStrLen\
  , struct tempType *tempValue);

#endif /* PYTYPES_H */


