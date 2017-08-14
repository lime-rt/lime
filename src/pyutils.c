/*
 *  pyutils.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2017 The LIME development team
 *
TODO:
 */

#include "lime.h"
#include "pytypes.h"

_Bool _doTest=FALSE;

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
The user should set and return both 'ordinary' and 'image' parameters as attributes of instances of par_classes.ModelParameters and par_classes.ImageParameters classes respectively. These classes define all the possible/permittable user-settable parameters, and are also used to set up the templates. The templates are used here basically to check that the user has not gone off and returned some object without some of the necessary parameters, or parameters of the wrong type.

Note that this routine can be, and is, used both for 'ordinary' and 'image' parameters.
  */

  int i,nItems,nItems1,status=0;
  PyObject *pAttr,*pFirstItem,*pFirstItem1;
  _Bool typesMatch=TRUE;

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

        if(PyList_CheckExact(pFirstItem)){ /* Is it a list? */
          nItems1 = (int)PyList_Size(pFirstItem);
          if(PyErr_Occurred()){
            printOrClearPyError();
            Py_DECREF(pAttr);
            return 4;
          }
          if(nItems1>0){
            pFirstItem1 = PyList_GetItem(pFirstItem, (Py_ssize_t)0); /* Don't have to DECREF pFirstItem1, it is a borrowed reference. */
            if(pFirstItem1==NULL){
              printOrClearPyError();
              Py_DECREF(pAttr);
              return 5;
            }
            typesMatch = checkAttrType(pFirstItem1, parTemplate[i].type);
          }
        }else{ /* ...no, just a scalar. */
          typesMatch = checkAttrType(pFirstItem, parTemplate[i].type);
        }
      }
    }else{ /* ...no, just a scalar. */
      typesMatch = checkAttrType(pAttr, parTemplate[i].type);
    }

    Py_DECREF(pAttr);

if(_doTest) printf("%d  %s  %s\n", i, parTemplate[i].name, parTemplate[i].type);

    if(!typesMatch)
      return 6;
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
extractScalarValue(PyObject *pPars, const char *attrName, const char *attrType\
  , const int maxStrLen, struct tempType *tempValue){
  /*
A slight wrapper for pyToC() in the case of scalar parameters, justified because extractListValues() is a much more extensive wrapper to pyToC() for list parameters.

Note: NO ERROR CHECKING is performed, since (we hope) all errors possible here have already been caught.
  */
  PyObject *pAttr = PyObject_GetAttrString(pPars, attrName);
  pyToC(pAttr, attrType, maxStrLen, tempValue);
  Py_DECREF(pAttr);
}

/*....................................................................*/
int
extractListValues(PyObject *pPars, const char *attrName, const char *attrType\
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
extractListListValues(PyObject *pPars, const char *attrName, const char *attrType\
  , const int maxStrLen, struct tempType ***tempValues, int (*nValues)[2]){
  /*
A wrapper for pyToC() to deal with list of lists parameters (actually pythonized 2D arrays).

Notes:
  - pAttr0 and pAttr1 should be lists.
  - The number of elements in pAttr1 is assumed to be constant.
  - NO ERROR CHECKING is performed, since (we hope) all errors possible here have already been caught.
  - tempValues should be NULL before the call.
  */
  int i,j,nValues0=0,nValues1=0;
  PyObject *pAttr0,*pAttr1,*pItem;

  pAttr0 = PyObject_GetAttrString(pPars, attrName);
  nValues0 = (int)PyList_Size(pAttr0);

  if(nValues0>0){
    *tempValues = malloc(sizeof(**tempValues)*nValues0);

    for(i=0;i<nValues0;i++){
      pAttr1 = PyList_GetItem(pAttr0, (Py_ssize_t)i); /* Don't have to DECREF pAttr1, it is a borrowed reference. */
      nValues1 = (int)PyList_Size(pAttr1);

      if(nValues1>0){
        (*tempValues)[i] = malloc(sizeof(*(*tempValues)[i])*nValues1);

        for(j=0;j<nValues1;j++){
          pItem = PyList_GetItem(pAttr1, (Py_ssize_t)j); /* Don't have to DECREF pItem, it is a borrowed reference. */
          pyToC(pItem, attrType, maxStrLen, &((*tempValues)[i])[j]);
        }
      }else
        (*tempValues)[i] = NULL;
    }
  }

  Py_DECREF(pAttr0);

  (*nValues)[0] = nValues0;
  (*nValues)[1] = nValues1;
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

***NOTE*** that the number, order and type of the parameters must be the same as given in ../python/par_classes.py.
  */

  int status=0,i,j,k,nValues,dims[2];
  PyObject *pFunc,*pArgs,*pPars,*pImgList,*pImgPars;
  struct tempType tempValue,*tempValues=NULL,**tempValues2=NULL;

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
  extractScalarValue(pPars, "radius",            parTemplate[i++].type, STR_LEN_0, &tempValue);
  par->radius            = tempValue.doubleValue;
  extractScalarValue(pPars, "minScale",          parTemplate[i++].type, STR_LEN_0, &tempValue);
  par->minScale          = tempValue.doubleValue;
  extractScalarValue(pPars, "pIntensity",        parTemplate[i++].type, STR_LEN_0, &tempValue);
  par->pIntensity        = tempValue.intValue;
  extractScalarValue(pPars, "sinkPoints",        parTemplate[i++].type, STR_LEN_0, &tempValue);
  par->sinkPoints        = tempValue.intValue;

  extractScalarValue(pPars, "dust",              parTemplate[i++].type, STR_LEN_0, &tempValue);
  if(strlen(tempValue.strValue)>0) strcpy(par->dust,          tempValue.strValue);
  extractScalarValue(pPars, "outputfile",        parTemplate[i++].type, STR_LEN_0, &tempValue);
  if(strlen(tempValue.strValue)>0) strcpy(par->outputfile,    tempValue.strValue);
  extractScalarValue(pPars, "binoutputfile",     parTemplate[i++].type, STR_LEN_0, &tempValue);
  if(strlen(tempValue.strValue)>0) strcpy(par->binoutputfile, tempValue.strValue);
  extractScalarValue(pPars, "gridfile",          parTemplate[i++].type, STR_LEN_0, &tempValue);
  if(strlen(tempValue.strValue)>0) strcpy(par->gridfile,      tempValue.strValue);
  extractScalarValue(pPars, "pregrid",           parTemplate[i++].type, STR_LEN_0, &tempValue);
  if(strlen(tempValue.strValue)>0) strcpy(par->pregrid,       tempValue.strValue);
  extractScalarValue(pPars, "restart",           parTemplate[i++].type, STR_LEN_0, &tempValue);
  if(strlen(tempValue.strValue)>0) strcpy(par->restart,       tempValue.strValue);
  extractScalarValue(pPars, "gridInFile",        parTemplate[i++].type, STR_LEN_0, &tempValue);
  if(strlen(tempValue.strValue)>0) strcpy(par->gridInFile,    tempValue.strValue);

  nValues = extractListValues(pPars, "collPartIds", parTemplate[i++].type, STR_LEN_0, &tempValues);
  if(nValues>0){
    if(nValues>MAX_N_COLL_PART){
      nValues = MAX_N_COLL_PART;
      if(!silent) warning("Too many collPartIds supplied!");
    }
    for(j=0;j<nValues;j++)
      par->collPartIds[j] = tempValues[j].intValue;
  }
  nValues = extractListValues(pPars, "nMolWeights", parTemplate[i++].type, STR_LEN_0, &tempValues);
  if(nValues>0){
    if(nValues>MAX_N_COLL_PART){
      nValues = MAX_N_COLL_PART;
      if(!silent) warning("Too many nMolWeights supplied!");
    }
    for(j=0;j<nValues;j++)
      par->nMolWeights[j] = tempValues[j].doubleValue;
  }
  nValues = extractListValues(pPars, "dustWeights", parTemplate[i++].type, STR_LEN_0, &tempValues);
  if(nValues>0){
    if(nValues>MAX_N_COLL_PART){
      nValues = MAX_N_COLL_PART;
      if(!silent) warning("Too many dustWeights supplied!");
    }
    for(j=0;j<nValues;j++)
      par->dustWeights[j] = tempValues[j].doubleValue;
  }
  nValues = extractListValues(pPars, "collPartMolWeights", parTemplate[i++].type, STR_LEN_0, &tempValues);
  if(nValues>0){
    if(nValues>MAX_N_COLL_PART){
      nValues = MAX_N_COLL_PART;
      if(!silent) warning("Too many collPartMolWeights supplied!");
    }
    for(j=0;j<nValues;j++)
      par->collPartMolWeights[j] = tempValues[j].doubleValue;
  }

  nValues = extractListValues(pPars, "gridDensMaxValues", parTemplate[i++].type, STR_LEN_0, &tempValues);
  if(nValues>0){
    if(nValues>MAX_N_HIGH){
      nValues = MAX_N_HIGH;
      if(!silent) warning("Too many gridDensMaxValues supplied!");
    }
    for(j=0;j<nValues;j++)
      par->gridDensMaxValues[j] = tempValues[j].doubleValue;
  }

  extractListListValues(pPars, "gridDensMaxLoc", parTemplate[i++].type, STR_LEN_0, &tempValues2, &dims);
  if(dims[0]>0){
    if(dims[0]>MAX_N_HIGH){
      dims[0] = MAX_N_HIGH;
      if(!silent) warning("Too many gridDensMaxLoc supplied!");
    }
    if(dims[1]>DIM){
      dims[1] = DIM;
      if(!silent) warning("Too many spatial dimensions in gridDensMaxLoc!");
    }
    for(j=0;j<dims[0];j++){
      for(k=0;k<dims[1];k++)
        par->gridDensMaxLoc[j][k] = tempValues2[j][k].doubleValue;
      free(tempValues2[j]);
    }
  }
  free(tempValues2);

  extractScalarValue(pPars, "tcmb",              parTemplate[i++].type, STR_LEN_0, &tempValue);
  par->tcmb              = tempValue.doubleValue;
  extractScalarValue(pPars, "lte_only",          parTemplate[i++].type, STR_LEN_0, &tempValue);
  par->lte_only          = tempValue.boolValue;
  extractScalarValue(pPars, "init_lte",          parTemplate[i++].type, STR_LEN_0, &tempValue);
  par->init_lte          = tempValue.boolValue;
  extractScalarValue(pPars, "samplingAlgorithm", parTemplate[i++].type, STR_LEN_0, &tempValue);
  par->samplingAlgorithm = tempValue.intValue;
  extractScalarValue(pPars, "sampling",          parTemplate[i++].type, STR_LEN_0, &tempValue);
  par->sampling          = tempValue.intValue;
  extractScalarValue(pPars, "blend",             parTemplate[i++].type, STR_LEN_0, &tempValue);
  par->blend             = tempValue.boolValue;
  extractScalarValue(pPars, "antialias",         parTemplate[i++].type, STR_LEN_0, &tempValue);
  par->antialias         = tempValue.intValue;
  extractScalarValue(pPars, "polarization",      parTemplate[i++].type, STR_LEN_0, &tempValue);
  par->polarization      = tempValue.boolValue;
  extractScalarValue(pPars, "nThreads",          parTemplate[i++].type, STR_LEN_0, &tempValue);
  par->nThreads          = tempValue.intValue;
  extractScalarValue(pPars, "nSolveIters",       parTemplate[i++].type, STR_LEN_0, &tempValue);
  par->nSolveIters       = tempValue.intValue;
  extractScalarValue(pPars, "traceRayAlgorithm", parTemplate[i++].type, STR_LEN_0, &tempValue);
  par->traceRayAlgorithm = tempValue.intValue;
  extractScalarValue(pPars, "resetRNG",          parTemplate[i++].type, STR_LEN_0, &tempValue);
  par->resetRNG          = tempValue.boolValue;
  extractScalarValue(pPars, "doSolveRTE",        parTemplate[i++].type, STR_LEN_0, &tempValue);
  par->doSolveRTE        = tempValue.boolValue;

  nValues = extractListValues(pPars, "gridOutFiles",  parTemplate[i++].type, STR_LEN_0, &tempValues);
  if(nValues>0){
    if(nValues>NUM_GRID_STAGES){
      nValues = NUM_GRID_STAGES;
      if(!silent) warning("Too many gridOutFiles supplied!");
    }
    for(j=0;j<nValues;j++)
      strcpy(par->gridOutFiles[j], tempValues[j].strValue);
  }

  nValues = extractListValues(pPars, "moldatfile",  parTemplate[i++].type, STR_LEN_0, &tempValues);
  if(nValues>0){
    if(nValues>MAX_NSPECIES){
      nValues = MAX_NSPECIES;
      if(!silent) warning("Too many moldatfile supplied!");
    }
    for(j=0;j<nValues;j++)
      strcpy(par->moldatfile[j], tempValues[j].strValue);
  }

  nValues = extractListValues(pPars, "girdatfile",  parTemplate[i++].type, STR_LEN_0, &tempValues);
  if(nValues>0){
    if(nValues>MAX_NSPECIES){
      nValues = MAX_NSPECIES;
      if(!silent) warning("Too many girdatfile supplied!");
    }
    for(j=0;j<nValues;j++)
      strcpy(par->girdatfile[j], tempValues[j].strValue);
  }

  nValues = extractListValues(pPars, "collPartNames",  parTemplate[i++].type, STR_LEN_0, &tempValues);
  if(nValues>0){
    if(nValues>MAX_N_COLL_PART){
      nValues = MAX_N_COLL_PART;
      if(!silent) warning("Too many collPartNames supplied!");
    }
    for(j=0;j<nValues;j++)
      strcpy(par->collPartNames[j], tempValues[j].strValue);
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
      (*img)[j].units = malloc(sizeof(char)*(STR_LEN_0+1)); /* Have to do this here because we want to strcpy() rather than copy the pointer to the start of a read-only string, as in traditional Lime. */
      (*img)[j].units[0] = '\0';
    }

    /* No error checking, because (in theory anyway) we have already checked all possible errors. */
    pImgList = PyObject_GetAttrString(pPars, "img");
    for(j=0;j<(*nImages);j++){
      pImgPars = PyList_GET_ITEM(pImgList, (Py_ssize_t)j); /* Don't have to DECREF pImgPars, it is a borrowed reference. */

      i = 0;
      extractScalarValue(pImgPars, "nchan",      imgParTemplate[i++].type, STR_LEN_0, &tempValue);
      (*img)[j].nchan = tempValue.intValue;
      extractScalarValue(pImgPars, "trans",      imgParTemplate[i++].type, STR_LEN_0, &tempValue);
      (*img)[j].trans = tempValue.intValue;
      extractScalarValue(pImgPars, "molI",       imgParTemplate[i++].type, STR_LEN_0, &tempValue);
      (*img)[j].molI = tempValue.intValue;
      extractScalarValue(pImgPars, "velres",     imgParTemplate[i++].type, STR_LEN_0, &tempValue);
      (*img)[j].velres = tempValue.doubleValue;
      extractScalarValue(pImgPars, "imgres",     imgParTemplate[i++].type, STR_LEN_0, &tempValue);
      (*img)[j].imgres = tempValue.doubleValue;
      extractScalarValue(pImgPars, "pxls",       imgParTemplate[i++].type, STR_LEN_0, &tempValue);
      (*img)[j].pxls = tempValue.intValue;
      extractScalarValue(pImgPars, "unit",       imgParTemplate[i++].type, STR_LEN_0, &tempValue);
      (*img)[j].unit = tempValue.intValue;
      extractScalarValue(pImgPars, "freq",       imgParTemplate[i++].type, STR_LEN_0, &tempValue);
      (*img)[j].freq = tempValue.doubleValue;
      extractScalarValue(pImgPars, "bandwidth",  imgParTemplate[i++].type, STR_LEN_0, &tempValue);
      (*img)[j].bandwidth = tempValue.doubleValue;
      extractScalarValue(pImgPars, "source_vel", imgParTemplate[i++].type, STR_LEN_0, &tempValue);
      (*img)[j].source_vel = tempValue.doubleValue;
      extractScalarValue(pImgPars, "theta",      imgParTemplate[i++].type, STR_LEN_0, &tempValue);
      (*img)[j].theta = tempValue.doubleValue;
      extractScalarValue(pImgPars, "phi",        imgParTemplate[i++].type, STR_LEN_0, &tempValue);
      (*img)[j].phi = tempValue.doubleValue;
      extractScalarValue(pImgPars, "incl",       imgParTemplate[i++].type, STR_LEN_0, &tempValue);
      (*img)[j].incl = tempValue.doubleValue;
      extractScalarValue(pImgPars, "posang",     imgParTemplate[i++].type, STR_LEN_0, &tempValue);
      (*img)[j].posang = tempValue.doubleValue;
      extractScalarValue(pImgPars, "azimuth",    imgParTemplate[i++].type, STR_LEN_0, &tempValue);
      (*img)[j].azimuth = tempValue.doubleValue;
      extractScalarValue(pImgPars, "distance",   imgParTemplate[i++].type, STR_LEN_0, &tempValue);
      (*img)[j].distance = tempValue.doubleValue;
      extractScalarValue(pImgPars, "doInterpolateVels",imgParTemplate[i++].type, STR_LEN_0, &tempValue);
      (*img)[j].doInterpolateVels = tempValue.boolValue;

      extractScalarValue(pImgPars, "filename",   imgParTemplate[i++].type, STR_LEN_0, &tempValue);
      if(strlen(tempValue.strValue)>0)
        strcpy((*img)[j].filename, tempValue.strValue);
      extractScalarValue(pImgPars, "units",      imgParTemplate[i++].type, STR_LEN_0, &tempValue);
      if(strlen(tempValue.strValue)>0)
        strcpy((*img)[j].units, tempValue.strValue);
    }
    Py_DECREF(pImgList);
  }

  Py_DECREF(pPars);

  return status;
}

