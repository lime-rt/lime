/*
 *  py_utils.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2017 The LIME development team
 *
TODO:
 */

#include "lime_config.h"
#include "inpars.h" /* for inputPars */
#include "error_codes.h"
#include "py_utils.h"
#include "constants.h"

PyObject *pMacros_global = NULL;

/* Globals for user-supplied result functions:
*/
PyObject *pDensity       = NULL,\
         *pTemperature   = NULL,\
         *pAbundance     = NULL,\
         *pMolNumDensity = NULL,\
         *pDoppler       = NULL,\
         *pVelocity      = NULL,\
         *pMagfield      = NULL,\
         *pGasIIdust     = NULL,\
         *pGridDensity   = NULL;

_Bool userFuncsInitialized=FALSE;

/*....................................................................*/
parTemplateType*
setTemplateDefaults(void){
//**** unused at present
  parTemplateType *template=NULL;

  template->name[0] = '\0';
  template->type[0] = '\0';
  template->mandatory = FALSE;
  template->isList    = FALSE;

  return template;
}

/*....................................................................*/
void
myStrCpy(const char *source, char *destination, const int maxStrlenDest){
  /*
**NOTE** that:
  - the source pointer must either be NULL or point to a 0-terminated series of characters in memory;
  - the destination pointer must either be allocated or a character array, with size maxStrlenDest+1.
  */

  int strlenSrc;

  if(source==NULL)
    destination[0] = '\0';
  else{
    strlenSrc = strlen(source);
    if(strlenSrc>maxStrlenDest)
      strncpy(destination, source, (size_t)maxStrlenDest); /* I expect this to write '\0' to destination[maxStrlenDest]. */
    else
      strcpy(destination, source);
  }
}

/*....................................................................*/
int
getModuleFromName(char *moduleNameNoSuffix, PyObject **pModule){
  /* Calling routine is expected to decref pModule */

  int status=0;
  PyObject *pName;

  pName = PyString_FromString(moduleNameNoSuffix);
  if(pName==NULL){
    *pModule = NULL;
return PY_STRING_READ_FAIL;
  }

  *pModule = PyImport_Import(pName);
  Py_DECREF(pName);
  if(*pModule==NULL){
return PY_IMPORT_FAIL;
  }

return status;
}

/*....................................................................*/
int
getParTemplates(PyObject *pParsClassOrInstance, parTemplateType **parTemplates\
  , int *nPars){
  /*
This converts templates for either the 'ordinary' parameter list or the 'image' parameter list. The templates, which are contained in the _listOfAttrs attribute of the relevant class or instance, define the name, type and default value of each parameter, as well as whether they are a list parameter and whether the parameter is mandatory.
  */

  PyObject *pListOfPars,*pListItem,*pTupleItem;
  int i,status=0;
  const char *tempStr;

  if (!PyClass_Check(pParsClassOrInstance) && !PyInstance_Check(pParsClassOrInstance)){
return GPT_NOT_CLASS;
  }

  pListOfPars = PyObject_GetAttrString(pParsClassOrInstance, "_listOfAttrs");
  if(pListOfPars==NULL){
return GPT_LIST_ATTR_NOT_FND;
  }

  if(!PyList_Check(pListOfPars)){
    Py_DECREF(pListOfPars);
return GPT_BAD_LIST;
  }

  *nPars = -1;
  *nPars = (int)PyList_Size(pListOfPars);
  if(*nPars<0){
    Py_DECREF(pListOfPars);
return GPT_BAD_N_LIST_ELEM;
  }

  *parTemplates = malloc(sizeof(**parTemplates)*(*nPars));
  if(*parTemplates==NULL){ /* malloc failed */
    Py_DECREF(pListOfPars);
    PyErr_NoMemory();
return GPT_MALLOC_FAIL;
  }

  for(i=0;i<(*nPars);i++){
    pListItem = PyList_GetItem(pListOfPars, (Py_ssize_t)i);
    /* Not going to check for errors, since I think I have covered all the possibilities already...? */

    if(!PyTuple_CheckExact(pListItem)){
      status = GPT_NOT_TUPLE;
  break;
    }

    pTupleItem = PyTuple_GetItem(pListItem, (Py_ssize_t)0);
    if(pTupleItem==NULL){
      status = GPT_BAD_TUPLE_ITEM;
  break;
    }

    tempStr = PyString_AsString(pTupleItem);
    if(tempStr==NULL){
      status = GPT_ITEM_NOT_STR;
  break;
    }

    myStrCpy(tempStr, (*parTemplates)[i].name, PY_MAX_LEN_PAR_NAME);

    pTupleItem = PyTuple_GetItem(pListItem, (Py_ssize_t)1);
    if(pTupleItem==NULL){
      status = GPT_BAD_TUPLE_ITEM;
  break;
    }

    tempStr = PyString_AsString(pTupleItem);
    if(tempStr==NULL){
      status = GPT_ITEM_NOT_STR;
  break;
    }

    myStrCpy(tempStr, (*parTemplates)[i].type, PY_MAX_LEN_PAR_TYPE);

    pTupleItem = PyTuple_GetItem(pListItem, (Py_ssize_t)2);
    if(pTupleItem==NULL){
      status = GPT_BAD_TUPLE_ITEM;
  break;
    }

    if(PyObject_IsTrue(pTupleItem)){
     (*parTemplates)[i].isList = TRUE;
    }else{
     (*parTemplates)[i].isList = FALSE;
    }

    pTupleItem = PyTuple_GetItem(pListItem, (Py_ssize_t)3);
    if(pTupleItem==NULL){
      status = GPT_BAD_TUPLE_ITEM;
  break;
    }

    if(PyObject_IsTrue(pTupleItem)){
     (*parTemplates)[i].mandatory = TRUE;
    }else{
     (*parTemplates)[i].mandatory = FALSE;
    }

    /* Don't need to call Py_DECREF for pListItem or pTupleItem because they are 'borrowed' references. */
  }

  Py_DECREF(pListOfPars);

return status;
}

/*....................................................................*/
int
mallocInputParStrs(inputPars *inpar){
  int i,j,status=0;

  #define RETURN_NO_MEM(status){\
            PyErr_NoMemory();\
            return status;\
           }

  inpar->collPartIds  = malloc(sizeof(int)*MAX_N_COLL_PART);
  if(inpar->collPartIds==NULL) RETURN_NO_MEM(1)
  for(i=0;i<MAX_N_COLL_PART;i++) inpar->collPartIds[i] = 0; /* Possible values start at 1. */

  inpar->nMolWeights  = malloc(sizeof(double)*MAX_N_COLL_PART);
  if(inpar->nMolWeights==NULL) RETURN_NO_MEM(2)
  for(i=0;i<MAX_N_COLL_PART;i++) inpar->nMolWeights[i] = -1.0;

  inpar->dustWeights  = malloc(sizeof(double)*MAX_N_COLL_PART); /* This param no longer has any effect. */
  if(inpar->dustWeights==NULL) RETURN_NO_MEM(3)
  for(i=0;i<MAX_N_COLL_PART;i++) inpar->dustWeights[i] = -1.0;

  inpar->collPartMolWeights = malloc(sizeof(double)*MAX_N_COLL_PART);
  if(inpar->collPartMolWeights==NULL) RETURN_NO_MEM(4)
  for(i=0;i<MAX_N_COLL_PART;i++) inpar->collPartMolWeights[i] = -1.0;

  inpar->gridDensMaxValues = malloc(sizeof(*(inpar->gridDensMaxValues))*MAX_N_HIGH);
  if(inpar->gridDensMaxValues==NULL) RETURN_NO_MEM(5)

  inpar->gridDensMaxLoc    = malloc(sizeof(*(inpar->gridDensMaxLoc))*MAX_N_HIGH);
  if(inpar->gridDensMaxLoc==NULL) RETURN_NO_MEM(6)
  for(i=0;i<MAX_N_HIGH;i++){
    inpar->gridDensMaxValues[i] = -1.0; /* Impossible default value. */
    for(j=0;j<DIM;j++) inpar->gridDensMaxLoc[i][j] = 0.0;
  }

  /* We have to malloc the strings in 'par' here (even though it means we will have to free them again afterward) because we want to strcpy() values into them rather than have them point to read-only strings.
  */
  inpar->outputfile    = malloc(sizeof(char)*(PY_STR_LEN_0+1));
  if(inpar->outputfile   ==NULL) RETURN_NO_MEM(7)
  inpar->outputfile[0]    = '\0';

  inpar->binoutputfile = malloc(sizeof(char)*(PY_STR_LEN_0+1));
  if(inpar->binoutputfile==NULL) RETURN_NO_MEM(8)
  inpar->binoutputfile[0] = '\0';

  inpar->gridfile      = malloc(sizeof(char)*(PY_STR_LEN_0+1));
  if(inpar->gridfile     ==NULL) RETURN_NO_MEM(9)
  inpar->gridfile[0]      = '\0';

  inpar->pregrid       = malloc(sizeof(char)*(PY_STR_LEN_0+1));
  if(inpar->pregrid      ==NULL) RETURN_NO_MEM(10)
  inpar->pregrid[0]       = '\0';

  inpar->restart       = malloc(sizeof(char)*(PY_STR_LEN_0+1));
  if(inpar->restart      ==NULL) RETURN_NO_MEM(11)
  inpar->restart[0]       = '\0';

  inpar->dust          = malloc(sizeof(char)*(PY_STR_LEN_0+1));
  if(inpar->dust         ==NULL) RETURN_NO_MEM(12)
  inpar->dust[0]          = '\0';

  inpar->gridInFile    = malloc(sizeof(char)*(PY_STR_LEN_0+1));
  if(inpar->gridInFile   ==NULL) RETURN_NO_MEM(13)
  inpar->gridInFile[0]    = '\0';

  inpar->moldatfile = malloc(sizeof(char *)*MAX_NSPECIES);
  if(inpar->moldatfile==NULL) RETURN_NO_MEM(14)

  inpar->girdatfile = malloc(sizeof(char *)*MAX_NSPECIES);
  if(inpar->girdatfile==NULL) RETURN_NO_MEM(15)

  for(i=0;i<MAX_NSPECIES;i++){
    inpar->moldatfile[i] = malloc(sizeof(char)*(PY_STR_LEN_0+1));
    if(inpar->moldatfile[i]==NULL) RETURN_NO_MEM(16)
    inpar->moldatfile[i][0] = '\0';

    inpar->girdatfile[i] = malloc(sizeof(char)*(PY_STR_LEN_0+1));
    if(inpar->girdatfile[i]==NULL) RETURN_NO_MEM(17)
    inpar->girdatfile[i][0] = '\0';
  }

  inpar->gridOutFiles = malloc(sizeof(char *)*NUM_GRID_STAGES);
  if(inpar->gridOutFiles==NULL) RETURN_NO_MEM(18)
  for(i=0;i<NUM_GRID_STAGES;i++){
    inpar->gridOutFiles[i] = malloc(sizeof(char)*(PY_STR_LEN_0+1));
    if(inpar->gridOutFiles[i]==NULL) RETURN_NO_MEM(19)
    inpar->gridOutFiles[i][0] = '\0';
  }

  /* Allocate initial space for (non-LAMDA) collision partner names */
  inpar->collPartNames = malloc(sizeof(char *)*MAX_N_COLL_PART);
  if(inpar->collPartNames==NULL) RETURN_NO_MEM(20)
  for(i=0;i<MAX_N_COLL_PART;i++){
    inpar->collPartNames[i] = malloc(sizeof(char)*(PY_STR_LEN_0+1));
    if(inpar->collPartNames[i]==NULL) RETURN_NO_MEM(21)
    inpar->collPartNames[i][0] = '\0';
  }

return status;
}

/*....................................................................*/
_Bool
_checkAttrType(PyObject *pObj, const char *parTemplateType){
  /* This compares the type of the parameter attribute read from the user's file to the type specified in the parameter template. */

  _Bool typesMatch=0;

  if(      strcmp(parTemplateType,"int"  )==0){
    typesMatch = (_Bool)PyInt_CheckExact(pObj);
  }else if(strcmp(parTemplateType,"float")==0){
    typesMatch = (_Bool)(PyFloat_CheckExact(pObj) || PyInt_CheckExact(pObj)); /* We also allow ints for this */
  }else if(strcmp(parTemplateType,"bool" )==0){
    typesMatch = (_Bool)(PyBool_Check(pObj) || PyInt_CheckExact(pObj)); /* We also allow ints for this */;
  }else if(strcmp(parTemplateType,"str"  )==0){
    typesMatch = ((_Bool)PyString_CheckExact(pObj) || pObj==Py_None); /* None is allowed for strings; should result in a value of NULL being given to the respective C character pointer. */
  }else if(strcmp(parTemplateType,"obj" )==0){
    typesMatch = (_Bool)1;
  }

  return typesMatch;
}

/*....................................................................*/
int
_checkAttributes(PyObject *pParInstance, parTemplateType *parTemplates\
  , const int nPars){
  /*
The user should set and return both 'ordinary' and 'image' parameters as attributes of instances of par_classes.ModelParameters and par_classes.ImageParameters classes respectively. These classes define all the possible/permittable user-settable parameters, and are also used to set up the templates. The templates are used here basically to check that the user has not gone off and returned some object without some of the necessary parameters, or parameters of the wrong type.

Note that this routine can be, and is, used both for 'ordinary' and 'image' parameters.
  */

  int i,nItems,nItems1,status=0;
  PyObject *pAttr,*pFirstItem,*pFirstItem1;
  _Bool typesMatch=TRUE;

  for(i=0;i<nPars;i++){
    pAttr = PyObject_GetAttrString(pParInstance, parTemplates[i].name);
    if(pAttr==NULL){
return CA_ATTR_READ_FAIL;
    }

    if(PyList_CheckExact(pAttr)){ /* Is it a list? */
      nItems = (int)PyList_Size(pAttr);
      if(nItems>0){
        pFirstItem = PyList_GetItem(pAttr, (Py_ssize_t)0); /* Don't have to DECREF pFirstItem, it is a borrowed reference. */
        if(pFirstItem==NULL){
          Py_DECREF(pAttr);
return CA_LIST_ITEM_READ_FAIL;
        }

        if(PyList_CheckExact(pFirstItem)){ /* Is it a list? */
          nItems1 = (int)PyList_Size(pFirstItem);
          if(nItems1>0){
            pFirstItem1 = PyList_GetItem(pFirstItem, (Py_ssize_t)0); /* Don't have to DECREF pFirstItem1, it is a borrowed reference. */
            if(pFirstItem1==NULL){
              Py_DECREF(pAttr);
return CA_LISTLIST_ITEM_READ_FAIL;
            }
            typesMatch = _checkAttrType(pFirstItem1, parTemplates[i].type);
          }
        }else{ /* ...no, just a scalar. */
          typesMatch = _checkAttrType(pFirstItem, parTemplates[i].type);
        }
      }
    }else{ /* ...no, just a scalar. */
      typesMatch = _checkAttrType(pAttr, parTemplates[i].type);
    }

    Py_DECREF(pAttr);

    if(!typesMatch)
return CA_TYPES_DONT_MATCH;
  }

return status;
}

/*....................................................................*/
void
_pyToC(PyObject *pAttr, const char *attrType, struct tempType *tempValue){
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
    if(pAttr==Py_None)
      (*tempValue).isNone = TRUE;
    else{
      (*tempValue).isNone = FALSE;
      myStrCpy(PyString_AsString(pAttr), (*tempValue).strValue, PY_STR_LEN_0);
    }
  }
}

/*....................................................................*/
void
_extractScalarValue(PyObject *pPars, const char *attrName, const char *attrType\
  , struct tempType *tempValue){
  /*
A slight wrapper for _pyToC() in the case of scalar parameters, justified because _extractListValues() is a much more extensive wrapper to _pyToC() for list parameters.

Note: NO ERROR CHECKING is performed, since (we hope) all errors possible here have already been caught.
  */
  PyObject *pAttr = PyObject_GetAttrString(pPars, attrName);
  _pyToC(pAttr, attrType, tempValue);
  Py_DECREF(pAttr);
}

/*....................................................................*/
int
_extractListValues(PyObject *pPars, const char *attrName, const char *attrType\
  , struct tempType **tempValues){
  /*
A wrapper for _pyToC() to deal with list parameters.

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
      _pyToC(pItem, attrType, &(*tempValues)[i]);
    }
  }

  Py_DECREF(pAttr);

  return nValues;
}

/*....................................................................*/
void
_extractListListValues(PyObject *pPars, const char *attrName, const char *attrType\
  , struct tempType ***tempValues, int (*nValues)[2]){
  /*
A wrapper for _pyToC() to deal with list of lists parameters (actually pythonized 2D arrays).

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
          _pyToC(pItem, attrType, &((*tempValues)[i])[j]);
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
pyFreeInputImgPars(image *inimg, int nImages){
  int i;

  if(inimg!=NULL){
    for(i=0;i<nImages;i++){
      free(inimg[i].filename);
      free(inimg[i].units);
    }
    free(inimg);
  }
}

/*....................................................................*/
int
readParImg(PyObject *pPars, parTemplateType *parTemplates\
  , const int nPars, parTemplateType *imgParTemplates, const int nImgPars\
  , inputPars *inpar, image **inimg, int *nImages, void (*warning)(char *message)){
  /*
Here we convert the 'ordinary' parameters from the supplied python object to its C struct form.

***NOTE TO DEVELOPERS: the number, order and type of the parameters listed in the present function must be the same as given in ../python/par_classes.py.
  */

  int status=0,i,j,k,nValues,dims[2];
  PyObject *pImgList,*pImgPars;
  struct tempType tempValue,*tempValues=NULL,**tempValues2=NULL;

  /* Check that all the required attributes are there and have the correct types.
  */
  status = _checkAttributes(pPars, parTemplates, nPars);
  if(status)
return status;

  /* Get the number of images and check that the image attributes, for all images, have the correct types.
  */
  pImgList = PyObject_GetAttrString(pPars, "img");
  if(pImgList==NULL)
return RPI_ATTR_READ_FAIL;

  *nImages = (int)PyList_Size(pImgList);

  for(j=0;j<(*nImages);j++){
    pImgPars = PyList_GetItem(pImgList, (Py_ssize_t)j); /* Don't have to DECREF pImgPars, it is a borrowed reference. */
    if(pImgPars==NULL){
      Py_DECREF(pImgList);
return RPI_LIST_ITEM_READ_FAIL;
    }

    status = _checkAttributes(pImgPars, imgParTemplates, nImgPars);
    if(status){
      if(status>0)
        status += 9;
      Py_DECREF(pImgList);
return status;
    }
  }

  Py_DECREF(pImgList);

  /* Finally, copy the parameters to the par struct.
  */
  i = 0;
  _extractScalarValue(pPars, "radius",            parTemplates[i++].type, &tempValue);
  inpar->radius            = tempValue.doubleValue;
  _extractScalarValue(pPars, "minScale",          parTemplates[i++].type, &tempValue);
  inpar->minScale          = tempValue.doubleValue;
  _extractScalarValue(pPars, "pIntensity",        parTemplates[i++].type, &tempValue);
  inpar->pIntensity        = tempValue.intValue;
  _extractScalarValue(pPars, "sinkPoints",        parTemplates[i++].type, &tempValue);
  inpar->sinkPoints        = tempValue.intValue;

  /* The following can have None as their python value:
  */
  _extractScalarValue(pPars, "dust",              parTemplates[i++].type, &tempValue);
  if(tempValue.isNone){
    free(inpar->dust);
    inpar->dust = NULL;
  }else if(strlen(tempValue.strValue)>0) /* otherwise leave the destination string as initialized to point to '\0' */
    strcpy(inpar->dust,          tempValue.strValue);

  _extractScalarValue(pPars, "outputfile",        parTemplates[i++].type, &tempValue);
  if(tempValue.isNone){
    free(inpar->outputfile);
    inpar->outputfile = NULL;
  }else if(strlen(tempValue.strValue)>0) /* otherwise leave the destination string as initialized to point to '\0' */
    strcpy(inpar->outputfile,    tempValue.strValue);

  _extractScalarValue(pPars, "binoutputfile",     parTemplates[i++].type, &tempValue);
  if(tempValue.isNone){
    free(inpar->binoutputfile);
    inpar->binoutputfile = NULL;
  }else if(strlen(tempValue.strValue)>0) /* otherwise leave the destination string as initialized to point to '\0' */
    strcpy(inpar->binoutputfile, tempValue.strValue);

  _extractScalarValue(pPars, "gridfile",          parTemplates[i++].type, &tempValue);
  if(tempValue.isNone){
    free(inpar->gridfile);
    inpar->gridfile = NULL;
  }else if(strlen(tempValue.strValue)>0) /* otherwise leave the destination string as initialized to point to '\0' */
    strcpy(inpar->gridfile,      tempValue.strValue);

  _extractScalarValue(pPars, "pregrid",           parTemplates[i++].type, &tempValue);
  if(tempValue.isNone){
    free(inpar->pregrid);
    inpar->pregrid = NULL;
  }else if(strlen(tempValue.strValue)>0) /* otherwise leave the destination string as initialized to point to '\0' */
    strcpy(inpar->pregrid,       tempValue.strValue);

  _extractScalarValue(pPars, "restart",           parTemplates[i++].type, &tempValue);
  if(tempValue.isNone){
    free(inpar->restart);
    inpar->restart = NULL;
  }else if(strlen(tempValue.strValue)>0) /* otherwise leave the destination string as initialized to point to '\0' */
    strcpy(inpar->restart,       tempValue.strValue);

  _extractScalarValue(pPars, "gridInFile",        parTemplates[i++].type, &tempValue);
  if(tempValue.isNone){
    free(inpar->gridInFile);
    inpar->gridInFile = NULL;
  }else if(strlen(tempValue.strValue)>0) /* otherwise leave the destination string as initialized to point to '\0' */
    strcpy(inpar->gridInFile,    tempValue.strValue);


  nValues = _extractListValues(pPars, "collPartIds", parTemplates[i++].type, &tempValues);
  if(nValues>0){
    if(nValues>MAX_N_COLL_PART){
      nValues = MAX_N_COLL_PART;
      if(!silent) warning("Too many collPartIds supplied!");//****** where do warnings fit in the picture if this is being called ultimately from python?
    }
    for(j=0;j<nValues;j++)
      inpar->collPartIds[j] = tempValues[j].intValue;
  }
  nValues = _extractListValues(pPars, "nMolWeights", parTemplates[i++].type, &tempValues);
  if(nValues>0){
    if(nValues>MAX_N_COLL_PART){
      nValues = MAX_N_COLL_PART;
      if(!silent) warning("Too many nMolWeights supplied!");//******
    }
    for(j=0;j<nValues;j++)
      inpar->nMolWeights[j] = tempValues[j].doubleValue;
  }
  nValues = _extractListValues(pPars, "dustWeights", parTemplates[i++].type, &tempValues);
  if(nValues>0){
    if(nValues>MAX_N_COLL_PART){
      nValues = MAX_N_COLL_PART;
      if(!silent) warning("Too many dustWeights supplied!");//******
    }
    for(j=0;j<nValues;j++)
      inpar->dustWeights[j] = tempValues[j].doubleValue;
  }
  nValues = _extractListValues(pPars, "collPartMolWeights", parTemplates[i++].type, &tempValues);
  if(nValues>0){
    if(nValues>MAX_N_COLL_PART){
      nValues = MAX_N_COLL_PART;
      if(!silent) warning("Too many collPartMolWeights supplied!");//******
    }
    for(j=0;j<nValues;j++)
      inpar->collPartMolWeights[j] = tempValues[j].doubleValue;
  }

  nValues = _extractListValues(pPars, "gridDensMaxValues", parTemplates[i++].type, &tempValues);
  if(nValues>0){
    if(nValues>MAX_N_HIGH){
      nValues = MAX_N_HIGH;
      if(!silent) warning("Too many gridDensMaxValues supplied!");//******
    }
    for(j=0;j<nValues;j++)
      inpar->gridDensMaxValues[j] = tempValues[j].doubleValue;
  }

  _extractListListValues(pPars, "gridDensMaxLoc", parTemplates[i++].type, &tempValues2, &dims);
  if(dims[0]>0){
    if(dims[0]>MAX_N_HIGH){
      dims[0] = MAX_N_HIGH;
      if(!silent) warning("Too many gridDensMaxLoc supplied!");//******
    }
    if(dims[1]>DIM){
      dims[1] = DIM;
      if(!silent) warning("Too many spatial dimensions in gridDensMaxLoc!");//******
    }
    for(j=0;j<dims[0];j++){
      for(k=0;k<dims[1];k++)
        inpar->gridDensMaxLoc[j][k] = tempValues2[j][k].doubleValue;
      free(tempValues2[j]);
    }
  }
  free(tempValues2);

  _extractScalarValue(pPars, "tcmb",              parTemplates[i++].type, &tempValue);
  inpar->tcmb              = tempValue.doubleValue;
  _extractScalarValue(pPars, "lte_only",          parTemplates[i++].type, &tempValue);
  inpar->lte_only          = tempValue.boolValue;
  _extractScalarValue(pPars, "init_lte",          parTemplates[i++].type, &tempValue);
  inpar->init_lte          = tempValue.boolValue;
  _extractScalarValue(pPars, "samplingAlgorithm", parTemplates[i++].type, &tempValue);
  inpar->samplingAlgorithm = tempValue.intValue;
  _extractScalarValue(pPars, "sampling",          parTemplates[i++].type, &tempValue);
  inpar->sampling          = tempValue.intValue;
  _extractScalarValue(pPars, "blend",             parTemplates[i++].type, &tempValue);
  inpar->blend             = tempValue.boolValue;
  _extractScalarValue(pPars, "antialias",         parTemplates[i++].type, &tempValue);
  inpar->antialias         = tempValue.intValue;
  _extractScalarValue(pPars, "polarization",      parTemplates[i++].type, &tempValue);
  inpar->polarization      = tempValue.boolValue;
  _extractScalarValue(pPars, "nThreads",          parTemplates[i++].type, &tempValue);
  inpar->nThreads          = tempValue.intValue;
  _extractScalarValue(pPars, "nSolveIters",       parTemplates[i++].type, &tempValue);
  inpar->nSolveIters       = tempValue.intValue;
  _extractScalarValue(pPars, "traceRayAlgorithm", parTemplates[i++].type, &tempValue);
  inpar->traceRayAlgorithm = tempValue.intValue;
  _extractScalarValue(pPars, "resetRNG",          parTemplates[i++].type, &tempValue);
  inpar->resetRNG          = tempValue.boolValue;
  _extractScalarValue(pPars, "doSolveRTE",        parTemplates[i++].type, &tempValue);
  inpar->doSolveRTE        = tempValue.boolValue;

  nValues = _extractListValues(pPars, "gridOutFiles",  parTemplates[i++].type, &tempValues);
  if(nValues>0){
    if(nValues>NUM_GRID_STAGES){
      nValues = NUM_GRID_STAGES;
      if(!silent) warning("Too many gridOutFiles supplied!");//******
    }
    for(j=0;j<nValues;j++)
      strcpy(inpar->gridOutFiles[j], tempValues[j].strValue);
  }

  nValues = _extractListValues(pPars, "moldatfile",  parTemplates[i++].type, &tempValues);
  if(nValues>0){
    if(nValues>MAX_NSPECIES){
      nValues = MAX_NSPECIES;
      if(!silent) warning("Too many moldatfile supplied!");//******
    }
    for(j=0;j<nValues;j++)
      strcpy(inpar->moldatfile[j], tempValues[j].strValue);
  }

  nValues = _extractListValues(pPars, "girdatfile",  parTemplates[i++].type, &tempValues);
  if(nValues>0){
    if(nValues>MAX_NSPECIES){
      nValues = MAX_NSPECIES;
      if(!silent) warning("Too many girdatfile supplied!");//******
    }
    for(j=0;j<nValues;j++)
      strcpy(inpar->girdatfile[j], tempValues[j].strValue);
  }

  nValues = _extractListValues(pPars, "collPartNames",  parTemplates[i++].type, &tempValues);
  if(nValues>0){
    if(nValues>MAX_N_COLL_PART){
      nValues = MAX_N_COLL_PART;
      if(!silent) warning("Too many collPartNames supplied!");//******
    }
    for(j=0;j<nValues;j++)
      strcpy(inpar->collPartNames[j], tempValues[j].strValue);
  }
  free(tempValues);

  /* Now we process the images.
  */
  if(*nImages>0){
    if(*nImages>MAX_NIMAGES){
      *nImages = MAX_NIMAGES;
      if(!silent) warning("Too many images supplied!");//******
    }
    *inimg = malloc(sizeof(**inimg)*(*nImages));
    if(*inimg==NULL){ /* malloc failed. */
      PyErr_NoMemory();
return RPI_MALLOC_FAIL;
    }

    /* Malloc the character pointers:
    */
    for(j=0;j<(*nImages);j++){
      (*inimg)[j].filename = malloc(sizeof(char)*(PY_STR_LEN_0+1)); /* Have to do this here because we want to strcpy() rather than copy the pointer to the start of a read-only string, as in traditional Lime. */
      if((*inimg)[j].filename==NULL){ /* malloc failed. */
        PyErr_NoMemory();
return RPI_MALLOC_FAIL;
      }
      (*inimg)[j].filename[0] = '\0';

      (*inimg)[j].units = malloc(sizeof(char)*(PY_STR_LEN_0+1)); /* Have to do this here because we want to strcpy() rather than copy the pointer to the start of a read-only string, as in traditional Lime. */
      if((*inimg)[j].units==NULL){ /* malloc failed. */
        PyErr_NoMemory();
return RPI_MALLOC_FAIL;
      }
      (*inimg)[j].units[0] = '\0';
    }

    /* No error checking, because (in theory anyway) we have already checked all possible errors. */
    pImgList = PyObject_GetAttrString(pPars, "img");
    for(j=0;j<(*nImages);j++){
      pImgPars = PyList_GET_ITEM(pImgList, (Py_ssize_t)j); /* Don't have to DECREF pImgPars, it is a borrowed reference. */

      i = 0;
      _extractScalarValue(pImgPars, "nchan",      imgParTemplates[i++].type, &tempValue);
      (*inimg)[j].nchan = tempValue.intValue;
      _extractScalarValue(pImgPars, "trans",      imgParTemplates[i++].type, &tempValue);
      (*inimg)[j].trans = tempValue.intValue;
      _extractScalarValue(pImgPars, "molI",       imgParTemplates[i++].type, &tempValue);
      (*inimg)[j].molI = tempValue.intValue;
      _extractScalarValue(pImgPars, "velres",     imgParTemplates[i++].type, &tempValue);
      (*inimg)[j].velres = tempValue.doubleValue;
      _extractScalarValue(pImgPars, "imgres",     imgParTemplates[i++].type, &tempValue);
      (*inimg)[j].imgres = tempValue.doubleValue;
      _extractScalarValue(pImgPars, "pxls",       imgParTemplates[i++].type, &tempValue);
      (*inimg)[j].pxls = tempValue.intValue;
      _extractScalarValue(pImgPars, "unit",       imgParTemplates[i++].type, &tempValue);
      (*inimg)[j].unit = tempValue.intValue;
      _extractScalarValue(pImgPars, "freq",       imgParTemplates[i++].type, &tempValue);
      (*inimg)[j].freq = tempValue.doubleValue;
      _extractScalarValue(pImgPars, "bandwidth",  imgParTemplates[i++].type, &tempValue);
      (*inimg)[j].bandwidth = tempValue.doubleValue;
      _extractScalarValue(pImgPars, "source_vel", imgParTemplates[i++].type, &tempValue);
      (*inimg)[j].source_vel = tempValue.doubleValue;
      _extractScalarValue(pImgPars, "theta",      imgParTemplates[i++].type, &tempValue);
      (*inimg)[j].theta = tempValue.doubleValue;
      _extractScalarValue(pImgPars, "phi",        imgParTemplates[i++].type, &tempValue);
      (*inimg)[j].phi = tempValue.doubleValue;
      _extractScalarValue(pImgPars, "incl",       imgParTemplates[i++].type, &tempValue);
      (*inimg)[j].incl = tempValue.doubleValue;
      _extractScalarValue(pImgPars, "posang",     imgParTemplates[i++].type, &tempValue);
      (*inimg)[j].posang = tempValue.doubleValue;
      _extractScalarValue(pImgPars, "azimuth",    imgParTemplates[i++].type, &tempValue);
      (*inimg)[j].azimuth = tempValue.doubleValue;
      _extractScalarValue(pImgPars, "distance",   imgParTemplates[i++].type, &tempValue);
      (*inimg)[j].distance = tempValue.doubleValue;
      _extractScalarValue(pImgPars, "doInterpolateVels",imgParTemplates[i++].type, &tempValue);
      (*inimg)[j].doInterpolateVels = tempValue.boolValue;

      /* The following can have None as their python value:
      */
      _extractScalarValue(pImgPars, "filename",   imgParTemplates[i++].type, &tempValue);
      if(tempValue.isNone){
        free((*inimg)[j].filename);
        (*inimg)[j].filename = NULL;
      }else if(strlen(tempValue.strValue)>0) /* otherwise leave the destination string as initialized to point to '\0' */
        strcpy((*inimg)[j].filename, tempValue.strValue);

      _extractScalarValue(pImgPars, "units",      imgParTemplates[i++].type, &tempValue);
      if(tempValue.isNone){
        free((*inimg)[j].units);
        (*inimg)[j].units = NULL;
      }else if(strlen(tempValue.strValue)>0) /* otherwise leave the destination string as initialized to point to '\0' */
        strcpy((*inimg)[j].units, tempValue.strValue);
    }
    Py_DECREF(pImgList);
  }

  return status;
}

/*....................................................................*/
int
setMacros(void){
  const int nDblMacros=14,nIntMacros=7;
  struct {char *name; double value;} macrosDbl[nDblMacros];
  struct {char *name;    int value;} macrosInt[nIntMacros];
  int status=0,i;
  PyObject *pValue;

  i = 0;
  macrosDbl[i++].name = "AMU";
  macrosDbl[i++].name = "CLIGHT";
  macrosDbl[i++].name = "HPLANCK";
  macrosDbl[i++].name = "KBOLTZ";
  macrosDbl[i++].name = "YJULIAN";
  macrosDbl[i++].name = "STEFANB";
  macrosDbl[i++].name = "GRAV";
  macrosDbl[i++].name = "AU";
  macrosDbl[i++].name = "LOCAL_CMB_TEMP";
  macrosDbl[i++].name = "PC";
  macrosDbl[i++].name = "MSUN";
  macrosDbl[i++].name = "RSUN";
  macrosDbl[i++].name = "PI";
  macrosDbl[i++].name = "SQRT_PI";
  i = 0;
  macrosDbl[i++].value = AMU;
  macrosDbl[i++].value = CLIGHT;
  macrosDbl[i++].value = HPLANCK;
  macrosDbl[i++].value = KBOLTZ;
  macrosDbl[i++].value = YJULIAN;
  macrosDbl[i++].value = STEFANB;
  macrosDbl[i++].value = GRAV;
  macrosDbl[i++].value = AU;
  macrosDbl[i++].value = LOCAL_CMB_TEMP;
  macrosDbl[i++].value = PC;
  macrosDbl[i++].value = MSUN;
  macrosDbl[i++].value = RSUN;
  macrosDbl[i++].value = M_PI;
  macrosDbl[i++].value = sqrt(M_PI);

  i = 0;
  macrosInt[i++].name = "CP_H2";
  macrosInt[i++].name = "CP_p_H2";
  macrosInt[i++].name = "CP_o_H2";
  macrosInt[i++].name = "CP_e";
  macrosInt[i++].name = "CP_H";
  macrosInt[i++].name = "CP_He";
  macrosInt[i++].name = "CP_Hplus";

  for(i=0;i<nIntMacros;i++)
    macrosInt[i].value = i+1;

  /* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .*/
  /* Construct the 'macro' list argument:
  */
  pMacros_global = PyDict_New();
  if(pMacros_global==NULL){
return PY_MACROS_FAIL;
  }

  for(i=0;i<nDblMacros;i++){
    pValue = PyFloat_FromDouble(macrosDbl[i].value);
    if(pValue==NULL){
      Py_DECREF(pMacros_global);
return PY_CONVERT_FAIL;
    }

    if(PyDict_SetItemString(pMacros_global, macrosDbl[i].name, pValue)){
      Py_DECREF(pValue); /* This is correct because PyDict_SetItemString() does NOT steal the reference (unlike with lists or tuples). */
      Py_DECREF(pMacros_global);
return PY_DICT_SET_FAIL;
    }

    Py_DECREF(pValue);
  }

  for(i=0;i<nIntMacros;i++){
    pValue = Py_BuildValue("i", macrosInt[i].value);
    if(pValue==NULL){
      Py_DECREF(pMacros_global);
return PY_CONVERT_FAIL;
    }

    if(PyDict_SetItemString(pMacros_global, macrosInt[i].name, pValue)){
      Py_DECREF(pValue); /* This is correct because PyDict_SetItemString() does NOT steal the reference (unlike with lists or tuples). */
      Py_DECREF(pMacros_global);
return PY_DICT_SET_FAIL;
    }

    Py_DECREF(pValue);
  }

return status;
}

/*....................................................................*/
void
unsetMacros(void){
  Py_XDECREF(pMacros_global); //************ replace the one with the other - a 1-line function is a bit silly.
}

/*....................................................................*/
void
getPythonFunc(PyObject *pModule, const char *funcName, PyObject **pFunc){
  /*
Wraps some error checking around code to return a function (i.e., a callable) attribute from a python object. Note I have changed this so it no longer returns a status: instead failure is indicated simply by a NULL value of pFunc. This is because it is usual for a less than complete set of these functions to be defined in an external python module - usually some are left to go to default.
  */

  *pFunc = PyObject_GetAttrString(pModule, funcName);
  if(*pFunc==NULL){
    PyErr_Clear(); //**** modellib seems to need this. :-/
return;
  }

  if (!PyCallable_Check(*pFunc)){
    Py_DECREF(*pFunc);
    *pFunc=NULL;
return;
  }
}

/*....................................................................*/
void
setUpUserPythonFuncs(PyObject *pModule){

  getPythonFunc(pModule, "density",            &pDensity);
  getPythonFunc(pModule, "temperature",    &pTemperature);
  getPythonFunc(pModule, "abundance",        &pAbundance);
  getPythonFunc(pModule, "molNumDensity",&pMolNumDensity);
  getPythonFunc(pModule, "doppler",            &pDoppler);
  getPythonFunc(pModule, "velocity",          &pVelocity);
  getPythonFunc(pModule, "magfield",          &pMagfield);
  getPythonFunc(pModule, "gasIIdust",        &pGasIIdust);
  getPythonFunc(pModule, "gridDensity",    &pGridDensity);

  userFuncsInitialized = TRUE;
}

/*....................................................................*/
void
decrefAllUserFuncs(void){
  Py_XDECREF(pDensity);
  Py_XDECREF(pTemperature);
  Py_XDECREF(pAbundance);
  Py_XDECREF(pMolNumDensity);
  Py_XDECREF(pDoppler);
  Py_XDECREF(pVelocity);
  Py_XDECREF(pMagfield);
  Py_XDECREF(pGasIIdust);
  Py_XDECREF(pGridDensity);
}

/*....................................................................*/
int
readScalarDouble(PyObject *pValue, double *value, _Bool *isNull){ 
  int status=0;

  *isNull = FALSE; /* default */

  if(PyFloat_CheckExact(pValue))
    *value = PyFloat_AS_DOUBLE(pValue);
  else if(PyInt_CheckExact(pValue))
    *value = (double)PyInt_AS_LONG(pValue);
  else if(pValue==Py_None)
    *isNull = TRUE;
  else{ /* raise exception */
return 1;
  }

return status;
}

/*....................................................................*/
int
userFuncWrapper(PyObject *pFunc, const char *funcName\
  , double x, double y, double z, double *resultBuffer, int *numElemInUserFuncReturn){ 
  /*
This is a generic wrapper for LIME functions (density, temperature etc) supplied by the user in a python module. 

The status return codes are in error_codes.h.
  */
  int status=0,i;
  PyObject *pArgs,*pResult,*pListItem;

  pArgs = Py_BuildValue("(Offf)", pMacros_global, x, y, z);
  if(pArgs==NULL)
return PY_BUILDVALUE_NULL;
  else{

    pResult = PyObject_CallObject(pFunc, pArgs);
    Py_DECREF(pArgs);
    if(pResult==NULL)
return PY_CALLOBJECT_NULL;
  }

  /* The returned value can be either a scalar or a list, but the datatypes should be floats or ints.
  */
  if(PyList_Check(pResult)){
    *numElemInUserFuncReturn = (int)PyList_Size(pResult);
    if(PyErr_Occurred() || (*numElemInUserFuncReturn)<=0){
      Py_DECREF(pResult);
      if(PyErr_Occurred())
return PY_LISTCHECK_ERR;
      else{
return PY_EMPTY_LIST;
      }
    }

    if(*numElemInUserFuncReturn>UFUNC_BUFFER_SIZE){
      Py_DECREF(pResult);
return PY_LIST_OVERFLOW;
    }

    for(i=0;i<(*numElemInUserFuncReturn);i++){
      pListItem = PyList_GetItem(pResult, (Py_ssize_t)i); /* pListItem is a borrowed reference, don't DECREF it. */
      /* Not going to check for errors, since I think I have covered all the possibilities already...? */

      if(PyFloat_CheckExact(pListItem))
        resultBuffer[i] = PyFloat_AS_DOUBLE(pListItem);
      else if(PyInt_CheckExact(pListItem))
        resultBuffer[i] = (double)PyInt_AS_LONG(pListItem);
      else{ /* raise exception */
        Py_DECREF(pResult);
return PY_NON_NUMERIC;
      }
    }

  }else{ /* Scalar returned. */
    *numElemInUserFuncReturn = 1;

    /* The returned value should be a float or an int.
    */
    i = 0;

    if(PyFloat_CheckExact(pResult))
      resultBuffer[i] = PyFloat_AS_DOUBLE(pResult);
    else if(PyInt_CheckExact(pResult))
      resultBuffer[i] = (double)PyInt_AS_LONG(pResult);
    else{
      Py_DECREF(pResult);
return PY_NON_NUMERIC;
    }
  }

  Py_DECREF(pResult);

return status;
}

