/*
 *  pyutils.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2016 The LIME development team
 *
TODO:
 */

#include <Python.h>
#include "lime.h"
#include "pytypes.h"

/*....................................................................*/
int
myStrCpy(const char source[], char destination[], const int strlenDest){
  int strlenSrc,i,status=0;

  strlenSrc = strlen(source);
  if(strlenSrc>strlenDest)
    return 1;

  for(i=0;i<strlenSrc;i++)
    destination[i] = source[i];
  destination[strlenSrc] = '\0';

  return status;
}

/*....................................................................*/
void
printOrClearPyError(){
  /* This can cause a problem if called when !PyErr_Occurred(). */
  if(silent)
    PyErr_Clear();
  else
    PyErr_Print(); /* Also clears. */
}

/*....................................................................*/
void
pyerror(char *message){
  if(!silent) bail_out(message);
  Py_Exit(1);
}

/*....................................................................*/
int
getParTemplates(const char *headerModuleName, parTemplateType **parTemplate, int *nPars\
  , parTemplateType **imgParTemplate, int *nImgPars){
  /*
This reads templates for the 'ordinary' parameter list and the 'image' parameter lists. The templates, which are read from a python module called here 'headerModuleName' (note this should be given minus the '.py'), define the name, type and default value of each parameter, as well as whether they are a list parameter and whether the parameter is mandatory.
  */

  PyObject *pName,*pModule,*pParClass,*pResult,*pListItem,*pTupleItem;
  int i,status=0;
  const char *tempStr;

  pName = PyString_FromString(headerModuleName);
  if(pName==NULL){
    printOrClearPyError();
    return 1;
  }

  pModule = PyImport_Import(pName);
  Py_DECREF(pName);

  if(pModule==NULL){
    printOrClearPyError();
    return 2;
  }

  /* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .*/
  /* Get the parameter template:
  */
  pParClass = PyObject_GetAttrString(pModule, "ModelParameters");
  if(pParClass==NULL){
    printOrClearPyError();
    Py_DECREF(pModule);
    return 3;
  }

  if (!PyClass_Check(pParClass)){
    if(PyErr_Occurred())
      printOrClearPyError();
    Py_DECREF(pParClass);
    Py_DECREF(pModule);
    return 4;
  }

  pResult = PyObject_GetAttrString(pParClass, "_listOfAttrs");
  if(pResult==NULL){
    printOrClearPyError();
    Py_DECREF(pParClass);
    Py_DECREF(pModule);
    return 5;
  }

  if(!PyList_Check(pResult)){
    if(PyErr_Occurred())
      printOrClearPyError();
    Py_DECREF(pResult);
    Py_DECREF(pParClass);
    Py_DECREF(pModule);
    return 6;
  }

  *nPars = (int)PyList_Size(pResult);
  if(PyErr_Occurred()){
    printOrClearPyError();
    Py_DECREF(pResult);
    Py_DECREF(pParClass);
    Py_DECREF(pModule);
    return 7;
  }

  *parTemplate = malloc(sizeof(**parTemplate)*(*nPars));

  for(i=0;i<(*nPars);i++){
    pListItem = PyList_GetItem(pResult, (Py_ssize_t)i);
    /* Not going to check for errors, since I think I have covered all the possibilities already...? */

    if(!PyTuple_CheckExact(pListItem)){
      if(PyErr_Occurred())
        printOrClearPyError();
      Py_DECREF(pResult);
      Py_DECREF(pParClass);
      Py_DECREF(pModule);
      return 8;
    }

    pTupleItem = PyTuple_GetItem(pListItem, (Py_ssize_t)0);
    if(pTupleItem==NULL){
      printOrClearPyError();
      Py_DECREF(pResult);
      Py_DECREF(pParClass);
      Py_DECREF(pModule);
      return 9;
    }

    tempStr = PyString_AsString(pTupleItem);
    if(tempStr==NULL){
      printOrClearPyError();
      Py_DECREF(pResult);
      Py_DECREF(pParClass);
      Py_DECREF(pModule);
      return 10;
    }

    if(myStrCpy(tempStr, (*parTemplate)[i].name, STR_LEN_0)){
      Py_DECREF(pResult);
      Py_DECREF(pParClass);
      Py_DECREF(pModule);
      return 11;
    }

    pTupleItem = PyTuple_GetItem(pListItem, (Py_ssize_t)1);
    if(pTupleItem==NULL){
      printOrClearPyError();
      Py_DECREF(pResult);
      Py_DECREF(pParClass);
      Py_DECREF(pModule);
      return 12;
    }

    tempStr = PyString_AsString(pTupleItem);
    if(tempStr==NULL){
      printOrClearPyError();
      Py_DECREF(pResult);
      Py_DECREF(pParClass);
      Py_DECREF(pModule);
      return 13;
    }

    if(myStrCpy(tempStr, (*parTemplate)[i].type, STR_LEN_0)){
      Py_DECREF(pResult);
      Py_DECREF(pParClass);
      Py_DECREF(pModule);
      return 14;
    }

    pTupleItem = PyTuple_GetItem(pListItem, (Py_ssize_t)2);
    if(pTupleItem==NULL){
      printOrClearPyError();
      Py_DECREF(pResult);
      Py_DECREF(pParClass);
      Py_DECREF(pModule);
      return 15;
    }

    (*parTemplate)[i].isList = (_Bool)PyInt_AsLong(pTupleItem);
    if(PyErr_Occurred()){
      printOrClearPyError();
      Py_DECREF(pResult);
      Py_DECREF(pParClass);
      Py_DECREF(pModule);
      return 16;
    }

    pTupleItem = PyTuple_GetItem(pListItem, (Py_ssize_t)3);
    if(pTupleItem==NULL){
      printOrClearPyError();
      Py_DECREF(pResult);
      Py_DECREF(pParClass);
      Py_DECREF(pModule);
      return 17;
    }

    (*parTemplate)[i].mandatory = (_Bool)PyInt_AsLong(pTupleItem);
    if(PyErr_Occurred()){
      printOrClearPyError();
      Py_DECREF(pResult);
      Py_DECREF(pParClass);
      Py_DECREF(pModule);
      return 18;
    }

    /* Don't need to call Py_DECREF for pListItem or pTupleItem because they are 'borrowed' references. */
  }

  Py_DECREF(pResult);
  Py_DECREF(pParClass);

  /* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .*/
  /* Now get the image-parameter template:
  */
  pParClass = PyObject_GetAttrString(pModule, "ImageParameters");
  if(pParClass==NULL){
    printOrClearPyError();
    Py_DECREF(pModule);
    return 19;
  }

  if (!PyClass_Check(pParClass)){
    if(PyErr_Occurred())
      printOrClearPyError();
    Py_DECREF(pParClass);
    Py_DECREF(pModule);
    return 20;
  }

  pResult = PyObject_GetAttrString(pParClass, "_listOfAttrs");
  if(pResult==NULL){
    printOrClearPyError();
    Py_DECREF(pParClass);
    Py_DECREF(pModule);
    return 21;
  }

  if(!PyList_Check(pResult)){
    if(PyErr_Occurred())
      printOrClearPyError();
    Py_DECREF(pResult);
    Py_DECREF(pParClass);
    Py_DECREF(pModule);
    return 22;
  }

  *nImgPars = (int)PyList_Size(pResult);
  if(PyErr_Occurred()){
    printOrClearPyError();
    Py_DECREF(pResult);
    Py_DECREF(pParClass);
    Py_DECREF(pModule);
    return 23;
  }

  *imgParTemplate = malloc(sizeof(**imgParTemplate)*(*nImgPars));

  for(i=0;i<(*nImgPars);i++){
    pListItem = PyList_GetItem(pResult, (Py_ssize_t)i);
    /* Not going to check for errors, since I think I have covered all the possibilities already...? */

    if(!PyTuple_CheckExact(pListItem)){
      if(PyErr_Occurred())
        printOrClearPyError();
      Py_DECREF(pResult);
      Py_DECREF(pParClass);
      Py_DECREF(pModule);
      return 24;
    }

    pTupleItem = PyTuple_GetItem(pListItem, (Py_ssize_t)0);
    if(pTupleItem==NULL){
      printOrClearPyError();
      Py_DECREF(pResult);
      Py_DECREF(pParClass);
      Py_DECREF(pModule);
      return 25;
    }

    tempStr = PyString_AsString(pTupleItem);
    if(tempStr==NULL){
      printOrClearPyError();
      Py_DECREF(pResult);
      Py_DECREF(pParClass);
      Py_DECREF(pModule);
      return 26;
    }

    if(myStrCpy(tempStr, (*imgParTemplate)[i].name, STR_LEN_0)){
      Py_DECREF(pResult);
      Py_DECREF(pParClass);
      Py_DECREF(pModule);
      return 27;
    }

    pTupleItem = PyTuple_GetItem(pListItem, (Py_ssize_t)1);
    if(pTupleItem==NULL){
      printOrClearPyError();
      Py_DECREF(pResult);
      Py_DECREF(pParClass);
      Py_DECREF(pModule);
      return 28;
    }

    tempStr = PyString_AsString(pTupleItem);
    if(tempStr==NULL){
      printOrClearPyError();
      Py_DECREF(pResult);
      Py_DECREF(pParClass);
      Py_DECREF(pModule);
      return 29;
    }

    if(myStrCpy(tempStr, (*imgParTemplate)[i].type, STR_LEN_0)){
      Py_DECREF(pResult);
      Py_DECREF(pParClass);
      Py_DECREF(pModule);
      return 30;
    }

    pTupleItem = PyTuple_GetItem(pListItem, (Py_ssize_t)2);
    if(pTupleItem==NULL){
      printOrClearPyError();
      Py_DECREF(pResult);
      Py_DECREF(pParClass);
      Py_DECREF(pModule);
      return 31;
    }

    (*imgParTemplate)[i].isList = (_Bool)PyInt_AsLong(pTupleItem);
    if(PyErr_Occurred()){
      printOrClearPyError();
      Py_DECREF(pResult);
      Py_DECREF(pParClass);
      Py_DECREF(pModule);
      return 32;
    }

    pTupleItem = PyTuple_GetItem(pListItem, (Py_ssize_t)3);
    if(pTupleItem==NULL){
      printOrClearPyError();
      Py_DECREF(pResult);
      Py_DECREF(pParClass);
      Py_DECREF(pModule);
      return 33;
    }

    (*imgParTemplate)[i].mandatory = (_Bool)PyInt_AsLong(pTupleItem);
    if(PyErr_Occurred()){
      printOrClearPyError();
      Py_DECREF(pResult);
      Py_DECREF(pParClass);
      Py_DECREF(pModule);
      return 34;
    }

    /* Don't need to call Py_DECREF for pListItem or pTupleItem because they are 'borrowed' references. */
  }

  Py_DECREF(pResult);
  Py_DECREF(pParClass);
  Py_DECREF(pModule);

  return status;
}

/*....................................................................*/
_Bool
checkAttrType(PyObject *pObj, const char *parTemplateType){
  /* This compares the type of the parameter attribute read from the user's file to the type specified in the parameter template. */

  _Bool typesMatch=0;

  if(      strcmp(parTemplateType,"int"  )==0){
    typesMatch = (_Bool)PyInt_CheckExact(pObj);
  }else if(strcmp(parTemplateType,"float")==0){
    typesMatch = (_Bool)(PyFloat_CheckExact(pObj) || PyInt_CheckExact(pObj)); /* We also allow ints for this */
  }else if(strcmp(parTemplateType,"bool" )==0){
    typesMatch = (_Bool)(PyBool_Check(pObj) || PyInt_CheckExact(pObj)); /* We also allow ints for this */;
  }else if(strcmp(parTemplateType,"str"  )==0){
    typesMatch = (_Bool)PyString_CheckExact(pObj);
  }else if(strcmp(parTemplateType,"obj" )==0){
    typesMatch = (_Bool)1;
  }

  return typesMatch;
}

/*....................................................................*/
int
checkAttributes(PyObject *pParInstance, parTemplateType *parTemplate\
  , const int nPars){
  /*
The user should set and return both 'ordinary' and 'image' parameters as attributes of instances of par_classes.ModelParameters and lime_h.ImageParameters classes respectively. These classes define all the possible/permittable user-settable parameters, and are also used to set up the templates. The templates are used here basically to check that the user has not gone off and returned some object without some of the necessary parameters, or parameters of the wrong type.

Note that this routine can be, and is, used both for 'ordinary' and 'image' parameters.
  */

  int i,nItems,status=0;
  PyObject *pAttr,*pFirstItem;
  _Bool typesMatch=1;

  for(i=0;i<nPars;i++){
    pAttr = PyObject_GetAttrString(pParInstance, parTemplate[i].name);
    if(pAttr==NULL){
      printOrClearPyError();
      return 1;
    }

    if(PyList_CheckExact(pAttr)){ /* Is it a list? */
      nItems = (int)PyList_Size(pAttr);
      if(PyErr_Occurred()){
        printOrClearPyError();
       Py_DECREF(pAttr);
       return 2;
      }
      if(nItems>0){
        pFirstItem = PyList_GetItem(pAttr, (Py_ssize_t)0); /* Don't have to DECREF pFirstItem, it is a borrowed reference. */
        if(pFirstItem==NULL){
          printOrClearPyError();
          Py_DECREF(pAttr);
          return 3;
        }
        typesMatch = checkAttrType(pFirstItem, parTemplate[i].type);
      }
    }else{ /* ...no, just a scalar. */
      typesMatch = checkAttrType(pAttr, parTemplate[i].type);
    }

    Py_DECREF(pAttr);

    if(!typesMatch)
      return 4;
  }

  return status;
}

/*....................................................................*/
void
pyToC(PyObject *pAttr, const char *attrType, const int maxStrLen\
  , struct tempType *tempValue){
  /*
This does the heavy lifting of converting the parameter attribute read from the user's supplied file from a python-format object to a C one.
 
Note: type 'obj' is possible, but results in no action here. Such attributes need to be dealt with elsewhere.
  */

  (*tempValue).intValue = 0;
  (*tempValue).doubleValue = 0.0;
  (*tempValue).boolValue = 0;
  (*tempValue).strValue[0] = '\0';

  if(      strcmp(attrType,"int"  )==0){
    (*tempValue).intValue = (int)PyInt_AS_LONG(pAttr);
  }else if(strcmp(attrType,"float")==0){
    if(PyInt_CheckExact(pAttr))
      (*tempValue).doubleValue = (double)PyInt_AS_LONG(pAttr);
    else
      (*tempValue).doubleValue = PyFloat_AS_DOUBLE(pAttr);
  }else if(strcmp(attrType,"bool" )==0){
    (*tempValue).boolValue = (_Bool)PyInt_AS_LONG(pAttr);
  }else if(strcmp(attrType,"str"  )==0){
    (void)myStrCpy(PyString_AsString(pAttr), (*tempValue).strValue, maxStrLen);
  }
}

/*....................................................................*/
void
extractValue(PyObject *pPars, const char *attrName, const char *attrType\
  , const int maxStrLen, struct tempType *tempValue){
  /*
A slight wrapper for pyToC() in the case of scalar parameters, justified because extractValues() is a much more extensive wrapper to pyToC() for list parameters.

Note: NO ERROR CHECKING is performed, since (we hope) all errors possible here have already been caught.
  */
  PyObject *pAttr = PyObject_GetAttrString(pPars, attrName);
  pyToC(pAttr, attrType, maxStrLen, tempValue);
  Py_DECREF(pAttr);
}

/*....................................................................*/
int
extractValues(PyObject *pPars, const char *attrName, const char *attrType\
  , const int maxStrLen, struct tempType **tempValues){
  /*
A wrapper for pyToC() to deal with list parameters.

Notes:
  - pAttr should be a list.
  - NO ERROR CHECKING is performed, since (we hope) all errors possible here have already been caught.
  - tempValues should be either malloc'd or NULL before the call.
  */
  int i,nValues=0;
  PyObject *pAttr,*pItem;

  free(*tempValues);
  *tempValues = NULL;

  pAttr = PyObject_GetAttrString(pPars, attrName);

  nValues = (int)PyList_Size(pAttr);
  if(nValues>0){
    *tempValues = malloc(sizeof(**tempValues)*nValues);
    for(i=0;i<nValues;i++){
      pItem = PyList_GetItem(pAttr, (Py_ssize_t)i); /* Don't have to DECREF pItem, it is a borrowed reference. */
      pyToC(pItem, attrType, maxStrLen, &(*tempValues)[i]);
    }
  }

  Py_DECREF(pAttr);

  return nValues;
}

/*....................................................................*/
void
getPythonFunc(PyObject *pModule, const char *funcName, const int verbosity\
  , PyObject **pFunc){
  /* Wraps some error checking around code to return a function (i.e., a callable) attribute from a python object. */

  *pFunc = PyObject_GetAttrString(pModule, funcName);
  if(*pFunc==NULL){
    if(verbosity>0)
      PyErr_Print();
    else
      PyErr_Clear();
    return;
  }

  if (!PyCallable_Check(*pFunc)){
    if(verbosity>0 && PyErr_Occurred())
      PyErr_Print();
    else
      PyErr_Clear();
    Py_DECREF(*pFunc);
    *pFunc=NULL;
    return;
  }
}

/*....................................................................*/
int
initParImg(PyObject *pModule, PyObject *pMacros, parTemplateType *parTemplate\
  , const int nPars, parTemplateType *imgParTemplate, const int nImgPars\
  , inputPars *par, image **img, int *nImages){
  /*
Here we extract par and img from the user's supplied 'model' module by calling their python function 'input(macros)'.
  */

  int status=0,i,j,nValues;
  PyObject *pFunc,*pArgs,*pPars,*pImgList,*pImgPars;
  struct tempType tempValue,*tempValues=NULL;

  getPythonFunc(pModule, "input", 1, &pFunc);
  if(pFunc==NULL){
    printOrClearPyError();
    return 1;
  }

  pArgs = Py_BuildValue("(O)", pMacros);
  if(pArgs==NULL){
    printOrClearPyError();
    Py_DECREF(pFunc);
    return 3;
  }

  pPars = PyObject_CallObject(pFunc, pArgs); /* Should return an instance of a par_classes.ModelParameters object. */
  Py_DECREF(pArgs);
  Py_DECREF(pFunc);
  if(pPars==NULL){
    printOrClearPyError();
    return 4;
  }
  if(!PyInstance_Check(pPars)){
    Py_DECREF(pPars);
    return 5;
  }

  /* Check that all the required attributes are there and have the correct types.
  */
  status = checkAttributes(pPars, parTemplate, nPars);
  if(status!=0){
    if(status>0)
      status += 5;
    Py_DECREF(pPars);
    return status;
  }

  /* Get the number of images and check that the image attributes, for all images, have the correct types.
  */
  pImgList = PyObject_GetAttrString(pPars, "img");
  if(pImgList==NULL){
    printOrClearPyError();
    Py_DECREF(pPars);
    return 10;
  }

  *nImages = (int)PyList_Size(pImgList);
  if(PyErr_Occurred()){
    printOrClearPyError();
    Py_DECREF(pImgList);
    Py_DECREF(pPars);
    return 11;
  }

  for(j=0;j<(*nImages);j++){
    pImgPars = PyList_GetItem(pImgList, (Py_ssize_t)j); /* Don't have to DECREF pImgPars, it is a borrowed reference. */
    if(pImgPars==NULL){
      printOrClearPyError();
      Py_DECREF(pImgList);
      Py_DECREF(pPars);
      return 12;
    }

    status = checkAttributes(pImgPars, imgParTemplate, nImgPars);
    if(status!=0){
      if(status>0)
        status += 12;
      Py_DECREF(pImgList);
      Py_DECREF(pPars);
      return status;
    }
  }

  Py_DECREF(pImgList);

  /* Finally, copy the parameters to the par struct.
  */
  i = 0;
  extractValue(pPars, "radius",            parTemplate[i++].type, STR_LEN_0, &tempValue);
  par->radius            = tempValue.doubleValue;
  extractValue(pPars, "minScale",          parTemplate[i++].type, STR_LEN_0, &tempValue);
  par->minScale          = tempValue.doubleValue;
  extractValue(pPars, "tcmb",              parTemplate[i++].type, STR_LEN_0, &tempValue);
  par->tcmb              = tempValue.doubleValue;
  extractValue(pPars, "sinkPoints",        parTemplate[i++].type, STR_LEN_0, &tempValue);
  par->sinkPoints        = tempValue.intValue;
  extractValue(pPars, "pIntensity",        parTemplate[i++].type, STR_LEN_0, &tempValue);
  par->pIntensity        = tempValue.intValue;
  extractValue(pPars, "blend",             parTemplate[i++].type, STR_LEN_0, &tempValue);
  par->blend             = tempValue.boolValue;
  extractValue(pPars, "traceRayAlgorithm", parTemplate[i++].type, STR_LEN_0, &tempValue);
  par->traceRayAlgorithm = tempValue.intValue;
  extractValue(pPars, "sampling",          parTemplate[i++].type, STR_LEN_0, &tempValue);
  par->sampling          = tempValue.intValue;
  extractValue(pPars, "lte_only",          parTemplate[i++].type, STR_LEN_0, &tempValue);
  par->lte_only          = tempValue.boolValue;
  extractValue(pPars, "init_lte",          parTemplate[i++].type, STR_LEN_0, &tempValue);
  par->init_lte          = tempValue.boolValue;
  extractValue(pPars, "antialias",         parTemplate[i++].type, STR_LEN_0, &tempValue);
  par->antialias         = tempValue.intValue;
  extractValue(pPars, "polarization",      parTemplate[i++].type, STR_LEN_0, &tempValue);
  par->polarization      = tempValue.boolValue;
  extractValue(pPars, "nThreads",          parTemplate[i++].type, STR_LEN_0, &tempValue);
  par->nThreads          = tempValue.intValue;

  extractValue(pPars, "outputfile",        parTemplate[i++].type, STR_LEN_0, &tempValue);
  if(strlen(tempValue.strValue)>0) strcpy(par->outputfile,    tempValue.strValue);
  extractValue(pPars, "binoutputfile",     parTemplate[i++].type, STR_LEN_0, &tempValue);
  if(strlen(tempValue.strValue)>0) strcpy(par->binoutputfile, tempValue.strValue);
  extractValue(pPars, "gridfile",          parTemplate[i++].type, STR_LEN_0, &tempValue);
  if(strlen(tempValue.strValue)>0) strcpy(par->gridfile,      tempValue.strValue);
  extractValue(pPars, "pregrid",           parTemplate[i++].type, STR_LEN_0, &tempValue);
  if(strlen(tempValue.strValue)>0) strcpy(par->pregrid,       tempValue.strValue);
  extractValue(pPars, "restart",           parTemplate[i++].type, STR_LEN_0, &tempValue);
  if(strlen(tempValue.strValue)>0) strcpy(par->restart,       tempValue.strValue);
  extractValue(pPars, "dust",              parTemplate[i++].type, STR_LEN_0, &tempValue);
  if(strlen(tempValue.strValue)>0) strcpy(par->dust,          tempValue.strValue);

  nValues = extractValues(pPars, "nMolWeights", parTemplate[i++].type, STR_LEN_0, &tempValues);
  if(nValues>0){
    if(nValues>MAX_N_COLL_PART){
      nValues = MAX_N_COLL_PART;
      if(!silent) warning("Too many nMolWeights supplied!");
    }
    for(j=0;j<nValues;j++)
      par->nMolWeights[j] = tempValues[j].doubleValue;
  }
  nValues = extractValues(pPars, "dustWeights", parTemplate[i++].type, STR_LEN_0, &tempValues);
  if(nValues>0){
    if(nValues>MAX_N_COLL_PART){
      nValues = MAX_N_COLL_PART;
      if(!silent) warning("Too many dustWeights supplied!");
    }
    for(j=0;j<nValues;j++)
      par->dustWeights[j] = tempValues[j].doubleValue;
  }
  nValues = extractValues(pPars, "collPartIds", parTemplate[i++].type, STR_LEN_0, &tempValues);
  if(nValues>0){
    if(nValues>MAX_N_COLL_PART){
      nValues = MAX_N_COLL_PART;
      if(!silent) warning("Too many collPartIds supplied!");
    }
    for(j=0;j<nValues;j++)
      par->collPartIds[j] = tempValues[j].intValue;
  }
  nValues = extractValues(pPars, "moldatfile",  parTemplate[i++].type, STR_LEN_0, &tempValues);
  if(nValues>0){
    if(nValues>MAX_NSPECIES){
      nValues = MAX_NSPECIES;
      if(!silent) warning("Too many moldatfile supplied!");
    }
    for(j=0;j<nValues;j++)
      strcpy(par->moldatfile[j], tempValues[j].strValue);
  }
  free(tempValues);

  /* Now we process the images.
  */
  if(*nImages>0){
    if(*nImages>MAX_NIMAGES){
      *nImages = MAX_NIMAGES;
      if(!silent) warning("Too many images supplied!");
    }
    (*img) = malloc(sizeof(**img)*(*nImages));
    for(j=0;j<(*nImages);j++){
      (*img)[j].filename = malloc(sizeof(char)*(STR_LEN_0+1)); /* Have to do this here because we want to strcpy() rather than copy the pointer to the start of a read-only string, as in traditional Lime. */
      (*img)[j].filename[0] = '\0';
    }

    /* No error checking, because (in theory anyway) we have already checked all possible errors. */
    pImgList = PyObject_GetAttrString(pPars, "img");
    for(j=0;j<(*nImages);j++){
      pImgPars = PyList_GET_ITEM(pImgList, (Py_ssize_t)j); /* Don't have to DECREF pImgPars, it is a borrowed reference. */

      i = 0;
      extractValue(pImgPars, "nchan",      imgParTemplate[i++].type, STR_LEN_0, &tempValue);
      (*img)[j].nchan = tempValue.intValue;
      extractValue(pImgPars, "trans",      imgParTemplate[i++].type, STR_LEN_0, &tempValue);
      (*img)[j].trans = tempValue.intValue;
      extractValue(pImgPars, "molI",       imgParTemplate[i++].type, STR_LEN_0, &tempValue);
      (*img)[j].molI = tempValue.intValue;
      extractValue(pImgPars, "velres",     imgParTemplate[i++].type, STR_LEN_0, &tempValue);
      (*img)[j].velres = tempValue.doubleValue;
      extractValue(pImgPars, "imgres",     imgParTemplate[i++].type, STR_LEN_0, &tempValue);
      (*img)[j].imgres = tempValue.doubleValue;
      extractValue(pImgPars, "pxls",       imgParTemplate[i++].type, STR_LEN_0, &tempValue);
      (*img)[j].pxls = tempValue.intValue;
      extractValue(pImgPars, "unit",       imgParTemplate[i++].type, STR_LEN_0, &tempValue);
      (*img)[j].unit = tempValue.intValue;
      extractValue(pImgPars, "freq",       imgParTemplate[i++].type, STR_LEN_0, &tempValue);
      (*img)[j].freq = tempValue.doubleValue;
      extractValue(pImgPars, "bandwidth",  imgParTemplate[i++].type, STR_LEN_0, &tempValue);
      (*img)[j].bandwidth = tempValue.doubleValue;
      extractValue(pImgPars, "source_vel", imgParTemplate[i++].type, STR_LEN_0, &tempValue);
      (*img)[j].source_vel = tempValue.doubleValue;
      extractValue(pImgPars, "theta",      imgParTemplate[i++].type, STR_LEN_0, &tempValue);
      (*img)[j].theta = tempValue.doubleValue;
      extractValue(pImgPars, "phi",        imgParTemplate[i++].type, STR_LEN_0, &tempValue);
      (*img)[j].phi = tempValue.doubleValue;
      extractValue(pImgPars, "incl",       imgParTemplate[i++].type, STR_LEN_0, &tempValue);
      (*img)[j].incl = tempValue.doubleValue;
      extractValue(pImgPars, "posang",     imgParTemplate[i++].type, STR_LEN_0, &tempValue);
      (*img)[j].posang = tempValue.doubleValue;
      extractValue(pImgPars, "azimuth",    imgParTemplate[i++].type, STR_LEN_0, &tempValue);
      (*img)[j].azimuth = tempValue.doubleValue;
      extractValue(pImgPars, "distance",   imgParTemplate[i++].type, STR_LEN_0, &tempValue);
      (*img)[j].distance = tempValue.doubleValue;

      extractValue(pImgPars, "filename",   imgParTemplate[i++].type, STR_LEN_0, &tempValue);
      if(strlen(tempValue.strValue)>0)
        strcpy((*img)[j].filename, tempValue.strValue);
    }
    Py_DECREF(pImgList);
  }

  Py_DECREF(pPars);

  return status;
}

