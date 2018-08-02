/*
 *  py_utils.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#include "lime_config.h"
#include "inpars.h" /* for inputPars */
#include "error_codes.h"
#include "local_err.h"
#include "py_lime.h"
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
void myStrNCpy(char *destination, const char *source, const size_t maxStrlenDest){
  /*
This is designed to behave like strncpy, i.e. to copy the first N characters contained in 'source' into 'destination', but to be a bit more forgiving. In particular it always returns 'destination' with '\0' somewhere in its N+1 allowed characters.

**NOTE** that:
  - the source pointer must either be NULL or point to a 0-terminated series of characters in memory;
  - the destination pointer must either be allocated or a character array, with size maxStrlenDest+1.
  */

  size_t strlenSrc;

  if(source==NULL)
    destination[0] = '\0';
  else{
    strlenSrc = strlen(source);
    if(strlenSrc>maxStrlenDest){
      strncpy(destination, source, maxStrlenDest);
      destination[maxStrlenDest] = '\0';
    }else
      strcpy(destination, source);
  }
}

/*....................................................................*/
errType
getModuleFromName(const char *moduleNameNoSuffix, PyObject **pModule){
  /* Calling routine is expected to decref pModule if and only if the return status is == 0. */

  errType err=init_local_err();
  PyObject *pName;
  char message[ERR_STR_LEN];

  pName = PyString_FromString(moduleNameNoSuffix);
  if(pName==NULL){
    *pModule = NULL;
    snprintf(message, ERR_STR_LEN-1, "Could not extract module object for %s", moduleNameNoSuffix);
return write_local_err(PY_STRING_READ_FAIL, message);
  }

  *pModule = PyImport_Import(pName);
  Py_DECREF(pName);
  if(*pModule==NULL){
    snprintf(message, ERR_STR_LEN-1, "Could not import %s", moduleNameNoSuffix);
return write_local_err(PY_IMPORT_FAIL, message);
  }

return err;
}

/*....................................................................*/
errType
getParTemplates(PyObject *pParsClassOrInstance, parTemplateType **parTemplates\
  , int *nPars){
  /*
This converts templates for either the 'ordinary' parameter list or the 'image' parameter list. The templates, which are contained in the _listOfAttrs attribute of the relevant class or instance, define the name, type and default value of each parameter, as well as whether they are a list parameter and whether the parameter is mandatory.
  */

  errType err=init_local_err();
  PyObject *pListOfPars,*pListItem,*pTupleItem;
  int i;
  const char *tempStr;
  char message[ERR_STR_LEN];

  if (!PyClass_Check(pParsClassOrInstance) && !PyInstance_Check(pParsClassOrInstance)){
return write_local_err(GPT_NOT_CLASS, "Object is neither a class nor an instance.");
  }

  pListOfPars = PyObject_GetAttrString(pParsClassOrInstance, "_listOfAttrs");
  if(pListOfPars==NULL){
return write_local_err(GPT_LIST_ATTR_NOT_FND, "Attribute '_listOfAttrs' not found in the instance.");
  }

  if(!PyList_Check(pListOfPars)){
    Py_DECREF(pListOfPars);
return write_local_err(GPT_BAD_LIST, "Attribute '_listOfAttrs' does not actually seem to be a list.");
  }

  *nPars = -1;
  *nPars = (int)PyList_Size(pListOfPars);
  if(*nPars<0){
    Py_DECREF(pListOfPars);
return write_local_err(GPT_BAD_N_LIST_ELEM, "Zero elements counted in attribute '_listOfAttrs'.");
  }

  *parTemplates = malloc(sizeof(**parTemplates)*(*nPars));
  if(*parTemplates==NULL){ /* malloc failed */
    Py_DECREF(pListOfPars);
    PyErr_NoMemory();
return write_local_err(GPT_MALLOC_FAIL, "Malloc failed for '*parTemplates'.");
  }

  for(i=0;i<(*nPars);i++){
    pListItem = PyList_GetItem(pListOfPars, (Py_ssize_t)i);
    /* Not going to check for errors, since I think I have covered all the possibilities already...? */

    if(!PyTuple_CheckExact(pListItem)){
      snprintf(message, ERR_STR_LEN-1, "List item %d is not a tuple.", i);
      err = write_local_err(GPT_NOT_TUPLE, message);
  break;
    }

    pTupleItem = PyTuple_GetItem(pListItem, (Py_ssize_t)0);
    if(pTupleItem==NULL){
      snprintf(message, ERR_STR_LEN-1, "Tuple %d element 0 is bad.", i);
      err = write_local_err(GPT_BAD_TUPLE_ITEM, message);
  break;
    }

    tempStr = PyString_AsString(pTupleItem);
    if(tempStr==NULL){
      snprintf(message, ERR_STR_LEN-1, "Tuple %d element 0 is not a string.", i);
      err = write_local_err(GPT_ITEM_NOT_STR, message);
  break;
    }

    myStrNCpy((*parTemplates)[i].name, tempStr, PY_MAX_LEN_PAR_NAME);

    pTupleItem = PyTuple_GetItem(pListItem, (Py_ssize_t)1);
    if(pTupleItem==NULL){
      snprintf(message, ERR_STR_LEN-1, "Tuple %d element 1 is bad.", i);
      err = write_local_err(GPT_BAD_TUPLE_ITEM, message);
  break;
    }

    tempStr = PyString_AsString(pTupleItem);
    if(tempStr==NULL){
      snprintf(message, ERR_STR_LEN-1, "Tuple %d element 1 is not a string.", i);
      err = write_local_err(GPT_ITEM_NOT_STR, message);
  break;
    }

    myStrNCpy((*parTemplates)[i].type, tempStr, PY_MAX_LEN_PAR_TYPE);

    pTupleItem = PyTuple_GetItem(pListItem, (Py_ssize_t)2);
    if(pTupleItem==NULL){
      snprintf(message, ERR_STR_LEN-1, "Tuple %d element 2 is bad.", i);
      err = write_local_err(GPT_BAD_TUPLE_ITEM, message);
  break;
    }

    if(PyObject_IsTrue(pTupleItem)){
     (*parTemplates)[i].isList = TRUE;
    }else{
     (*parTemplates)[i].isList = FALSE;
    }

    pTupleItem = PyTuple_GetItem(pListItem, (Py_ssize_t)3);
    if(pTupleItem==NULL){
      snprintf(message, ERR_STR_LEN-1, "Tuple %d element 3 is bad.", i);
      err = write_local_err(GPT_BAD_TUPLE_ITEM, message);
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

return err;
}

/*....................................................................*/
errType
getParTemplatesWrapper(const char *headerModuleName, parTemplateType **parTemplates\
  , int *nPars, parTemplateType **imgParTemplates, int *nImgPars){

  errType err=init_local_err();
  PyObject *pModule,*pParClass,*pImgParClass;
  char message[ERR_STR_LEN];

  err = getModuleFromName(headerModuleName, &pModule);
  if(err.status){
    printOrClearPyError();
return err;
  }

  /* Get the parameter templates:
  */
  pParClass = PyObject_GetAttrString(pModule, "ModelParameters");
  if(pParClass==NULL){
    Py_DECREF(pModule);
    printOrClearPyError();
    snprintf(message, ERR_STR_LEN-1, "Attribute 'ModelParameters' not found in module.");
return write_local_err(3, message);
  }

  err = getParTemplates(pParClass, parTemplates, nPars);
  if(err.status){
    Py_DECREF(pParClass);
    Py_DECREF(pModule);
return err;
  }

  Py_DECREF(pParClass);

  /* Get the image parameter templates:
  */
  pImgParClass = PyObject_GetAttrString(pModule, "ImageParameters");
  if(pImgParClass==NULL){
    Py_DECREF(pModule);
    printOrClearPyError();
    snprintf(message, ERR_STR_LEN-1, "Attribute 'ImageParameters' not found in module.");
return write_local_err(5, message);
  }

  Py_DECREF(pModule);

  err = getParTemplates(pImgParClass, imgParTemplates, nImgPars);
  if(err.status){
    Py_DECREF(pImgParClass);
return err;
  }

  Py_DECREF(pImgParClass);

return err;
}

/*....................................................................*/
errType
mallocInputParStrs(inputPars *inpar){
  int i,j;
  errType err=init_local_err();

  #define RETURN_NO_MEM(status, message){\
            PyErr_NoMemory();\
            return write_local_err(status, message);\
           }

  inpar->collPartIds  = malloc(sizeof(int)*MAX_N_COLL_PART);
  if(inpar->collPartIds==NULL) RETURN_NO_MEM(1, "Malloc failed for inpar->collPartIds.")
  for(i=0;i<MAX_N_COLL_PART;i++) inpar->collPartIds[i] = 0; /* Possible values start at 1. */

  inpar->nMolWeights  = malloc(sizeof(double)*MAX_N_COLL_PART);
  if(inpar->nMolWeights==NULL) RETURN_NO_MEM(2, "Malloc failed for inpar->nMolWeights.")
  for(i=0;i<MAX_N_COLL_PART;i++) inpar->nMolWeights[i] = -1.0;

  inpar->dustWeights  = malloc(sizeof(double)*MAX_N_COLL_PART); /* This param no longer has any effect. */
  if(inpar->dustWeights==NULL) RETURN_NO_MEM(3, "Malloc failed for inpar->dustWeights.")
  for(i=0;i<MAX_N_COLL_PART;i++) inpar->dustWeights[i] = -1.0;

  inpar->collPartMolWeights = malloc(sizeof(double)*MAX_N_COLL_PART);
  if(inpar->collPartMolWeights==NULL) RETURN_NO_MEM(4, "Malloc failed for inpar->collPartMolWeights.")
  for(i=0;i<MAX_N_COLL_PART;i++) inpar->collPartMolWeights[i] = -1.0;

  inpar->gridDensMaxValues = malloc(sizeof(*(inpar->gridDensMaxValues))*MAX_N_HIGH);
  if(inpar->gridDensMaxValues==NULL) RETURN_NO_MEM(5, "Malloc failed for inpar->gridDensMaxValues.")

  inpar->gridDensMaxLoc    = malloc(sizeof(*(inpar->gridDensMaxLoc))*MAX_N_HIGH);
  if(inpar->gridDensMaxLoc==NULL) RETURN_NO_MEM(6, "Malloc failed for inpar->gridDensMaxLoc.")
  for(i=0;i<MAX_N_HIGH;i++){
    inpar->gridDensMaxValues[i] = -1.0; /* Impossible default value. */
    for(j=0;j<DIM;j++) inpar->gridDensMaxLoc[i][j] = 0.0;
  }

  /* We have to malloc the strings in 'par' here (even though it means we will have to free them again afterward) because we want to strcpy() values into them rather than have them point to read-only strings.
  */
  inpar->outputfile    = malloc(sizeof(char)*(PY_STR_LEN_0+1));
  if(inpar->outputfile   ==NULL) RETURN_NO_MEM(7, "Malloc failed for inpar->outputfile.")
  inpar->outputfile[0]    = '\0';

  inpar->binoutputfile = malloc(sizeof(char)*(PY_STR_LEN_0+1));
  if(inpar->binoutputfile==NULL) RETURN_NO_MEM(8, "Malloc failed for inpar->binoutputfile.")
  inpar->binoutputfile[0] = '\0';

  inpar->gridfile      = malloc(sizeof(char)*(PY_STR_LEN_0+1));
  if(inpar->gridfile     ==NULL) RETURN_NO_MEM(9, "Malloc failed for inpar->gridfile.")
  inpar->gridfile[0]      = '\0';

  inpar->pregrid       = malloc(sizeof(char)*(PY_STR_LEN_0+1));
  if(inpar->pregrid      ==NULL) RETURN_NO_MEM(10, "Malloc failed for inpar->pregrid.")
  inpar->pregrid[0]       = '\0';

  inpar->restart       = malloc(sizeof(char)*(PY_STR_LEN_0+1));
  if(inpar->restart      ==NULL) RETURN_NO_MEM(11, "Malloc failed for inpar->restart.")
  inpar->restart[0]       = '\0';

  inpar->dust          = malloc(sizeof(char)*(PY_STR_LEN_0+1));
  if(inpar->dust         ==NULL) RETURN_NO_MEM(12, "Malloc failed for inpar->dust.")
  inpar->dust[0]          = '\0';

  inpar->gridInFile    = malloc(sizeof(char)*(PY_STR_LEN_0+1));
  if(inpar->gridInFile   ==NULL) RETURN_NO_MEM(13, "Malloc failed for inpar->gridInFile.")
  inpar->gridInFile[0]    = '\0';

  inpar->moldatfile = malloc(sizeof(char *)*MAX_NSPECIES);
  if(inpar->moldatfile==NULL) RETURN_NO_MEM(14, "Malloc failed for inpar->moldatfile.")

  inpar->girdatfile = malloc(sizeof(char *)*MAX_NSPECIES);
  if(inpar->girdatfile==NULL) RETURN_NO_MEM(15, "Malloc failed for inpar->girdatfile.")

  for(i=0;i<MAX_NSPECIES;i++){
    inpar->moldatfile[i] = malloc(sizeof(char)*(PY_STR_LEN_0+1));
    if(inpar->moldatfile[i]==NULL) RETURN_NO_MEM(16, "Malloc failed for inpar->moldatfile.")
    inpar->moldatfile[i][0] = '\0';

    inpar->girdatfile[i] = malloc(sizeof(char)*(PY_STR_LEN_0+1));
    if(inpar->girdatfile[i]==NULL) RETURN_NO_MEM(17, "Malloc failed for inpar->girdatfile.")
    inpar->girdatfile[i][0] = '\0';
  }

  inpar->gridOutFiles = malloc(sizeof(char *)*NUM_GRID_STAGES);
  if(inpar->gridOutFiles==NULL) RETURN_NO_MEM(18, "Malloc failed for inpar->gridOutFiles.")
  for(i=0;i<NUM_GRID_STAGES;i++){
    inpar->gridOutFiles[i] = malloc(sizeof(char)*(PY_STR_LEN_0+1));
    if(inpar->gridOutFiles[i]==NULL) RETURN_NO_MEM(19, "Malloc failed for inpar->gridOutFiles.")
    inpar->gridOutFiles[i][0] = '\0';
  }

  /* Allocate initial space for (non-LAMDA) collision partner names */
  inpar->collPartNames = malloc(sizeof(char *)*MAX_N_COLL_PART);
  if(inpar->collPartNames==NULL) RETURN_NO_MEM(20, "Malloc failed for inpar->collPartNames.")
  for(i=0;i<MAX_N_COLL_PART;i++){
    inpar->collPartNames[i] = malloc(sizeof(char)*(PY_STR_LEN_0+1));
    if(inpar->collPartNames[i]==NULL) RETURN_NO_MEM(21, "Malloc failed for inpar->collPartNames.")
    inpar->collPartNames[i][0] = '\0';
  }

return err;
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
errType
_checkAttributes(PyObject *pParInstance, parTemplateType *parTemplates\
  , const int nPars){
  /*
The user should set and return both 'ordinary' and 'image' parameters as attributes of instances of limepar_classes.ModelParameters and limepar_classes.ImageParameters classes respectively. These classes define all the possible/permittable user-settable parameters, and are also used to set up the templates. The templates are used here basically to check that the user has not gone off and returned some object without some of the necessary parameters, or parameters of the wrong type.

Note that this routine can be, and is, used both for 'ordinary' and 'image' parameters.
  */

  errType err=init_local_err();
  int i,j,nItems,nItems1;
  PyObject *pAttr,*pFirstItem,*pFirstItem1,*pListOfAttrs,*pString;
  _Bool typesMatch=TRUE,parFound;
  char *tempStr;
  char message[ERR_STR_LEN];

  for(i=0;i<nPars;i++){
    pAttr = PyObject_GetAttrString(pParInstance, parTemplates[i].name);
    if(pAttr==NULL){
      snprintf(message, ERR_STR_LEN-1, "Could not read attribute %s.", parTemplates[i].name);
return write_local_err(CA_ATTR_READ_FAIL, message);
    }

    if(PyList_CheckExact(pAttr)){ /* Is it a list? */
      nItems = (int)PyList_Size(pAttr);
      if(nItems>0){
        pFirstItem = PyList_GetItem(pAttr, (Py_ssize_t)0); /* Don't have to DECREF pFirstItem, it is a borrowed reference. */
        if(pFirstItem==NULL){
          Py_DECREF(pAttr);
          snprintf(message, ERR_STR_LEN-1, "Could not read first item in list attribute %s.", parTemplates[i].name);
return write_local_err(CA_LIST_ITEM_READ_FAIL, message);
        }

        if(PyList_CheckExact(pFirstItem)){ /* Is it a list? */
          nItems1 = (int)PyList_Size(pFirstItem);
          if(nItems1>0){
            pFirstItem1 = PyList_GetItem(pFirstItem, (Py_ssize_t)0); /* Don't have to DECREF pFirstItem1, it is a borrowed reference. */
            if(pFirstItem1==NULL){
              Py_DECREF(pAttr);
              snprintf(message, ERR_STR_LEN-1, "Could not read item [0][0] in attribute %s.", parTemplates[i].name);
return write_local_err(CA_LISTLIST_ITEM_READ_FAIL, message);
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

    if(!typesMatch){
      snprintf(message, ERR_STR_LEN-1, "Attribute does not match the template %s.", parTemplates[i].name);
return write_local_err(CA_TYPES_DONT_MATCH, message);
    }
  }

  /* We want now to check that the user has not specified any extra, non-recognized attributes. These will have no effect on the running of the program but the user should be warned that such things are simply thrown away.
  */
  pListOfAttrs = PyObject_Dir(pParInstance);
  nItems = (int)PyList_Size(pListOfAttrs);
  for(i=0;i<nItems;i++){
    pString = PyList_GetItem(pListOfAttrs, (Py_ssize_t)i); /* Don't have to DECREF pItem, it is a borrowed reference. */
    tempStr = PyString_AsString(pString);
    if(tempStr[0]!='_'){
      parFound = FALSE; /* default */

      for(j=0;j<nPars;j++){
        if(strcmp(tempStr, parTemplates[j].name)==0){
          parFound = TRUE;
      break;
        }
      }

      if(!parFound){
        snprintf(message, ERR_STR_LEN-1, "User has added non-canonical attribute %s. This *won't* be read!", tempStr);
        err=write_local_err(CA_ILLEGAL_PAR, message);
      }
    }
  }

return err;
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
      myStrNCpy((*tempValue).strValue, PyString_AsString(pAttr), PY_STR_LEN_0);
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
errType
readParImg(PyObject *pPars, parTemplateType *parTemplates\
  , const int nPars, parTemplateType *imgParTemplates, const int nImgPars\
  , inputPars *inpar, image **inimg, int *nImages, void (*warning)(char *message)){
  /*
Here we convert the 'ordinary' parameters from the supplied python object to its C struct form.

***NOTE TO DEVELOPERS: the number, order and type of the parameters listed in the present function must be the same as given in ../python/limepar_classes.py.
  */

  errType err=init_local_err();
  int i,j,k,nValues,dims[2];
  PyObject *pImgList,*pImgPars;
  struct tempType tempValue,*tempValues=NULL,**tempValues2=NULL;
  char message[ERR_STR_LEN];

  /* Check that all the required attributes are there and have the correct types.
  */
  err = _checkAttributes(pPars, parTemplates, nPars);
  if(err.status)
return err;

  /* Get the number of images and check that the image attributes, for all images, have the correct types.
  */
  pImgList = PyObject_GetAttrString(pPars, "img");
  if(pImgList==NULL)
return write_local_err(RPI_ATTR_READ_FAIL, "Could not read attribute 'img'.");

  *nImages = (int)PyList_Size(pImgList);

  for(j=0;j<(*nImages);j++){
    pImgPars = PyList_GetItem(pImgList, (Py_ssize_t)j); /* Don't have to DECREF pImgPars, it is a borrowed reference. */
    if(pImgPars==NULL){
      Py_DECREF(pImgList);
      snprintf(message, ERR_STR_LEN-1, "Could not read 'img' item %d.", j);
return write_local_err(RPI_LIST_ITEM_READ_FAIL, message);
    }

    err = _checkAttributes(pImgPars, imgParTemplates, nImgPars);
    if(err.status){
      Py_DECREF(pImgList);
return err;
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
  }else if(strlen(tempValue.strValue)>0){ /* otherwise leave the destination string as initialized to point to '\0' */
    strcpy(inpar->dust,          tempValue.strValue);
  }

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
return write_local_err(RPI_MALLOC_FAIL, "Malloc of '*inimg' failed.");
    }

    /* Malloc the character pointers:
    */
    for(j=0;j<(*nImages);j++){
      (*inimg)[j].filename = malloc(sizeof(char)*(PY_STR_LEN_0+1)); /* Have to do this here because we want to strcpy() rather than copy the pointer to the start of a read-only string, as in traditional Lime. */
      if((*inimg)[j].filename==NULL){ /* malloc failed. */
        PyErr_NoMemory();
        snprintf(message, ERR_STR_LEN-1, "Malloc of '(*inimg)[%d].filename' failed.", j);
return write_local_err(RPI_MALLOC_FAIL, message);
      }
      (*inimg)[j].filename[0] = '\0';

      (*inimg)[j].units = malloc(sizeof(char)*(PY_STR_LEN_0+1)); /* Have to do this here because we want to strcpy() rather than copy the pointer to the start of a read-only string, as in traditional Lime. */
      if((*inimg)[j].units==NULL){ /* malloc failed. */
        PyErr_NoMemory();
        snprintf(message, ERR_STR_LEN-1, "Malloc of '(*inimg)[%d].units' failed.", j);
return write_local_err(RPI_MALLOC_FAIL, message);
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

return err;
}

/*....................................................................*/
errType
setMacros(void){
  const int nDblMacros=14,nIntMacros=7;
  struct {char *name; double value;} macrosDbl[nDblMacros];
  struct {char *name;    int value;} macrosInt[nIntMacros];
  int i;
  PyObject *pValue;
  errType err=init_local_err();
  char message[ERR_STR_LEN];

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
return write_local_err(PY_MACROS_FAIL, "Failed to create python dict for macros.");
  }

  for(i=0;i<nDblMacros;i++){
    pValue = PyFloat_FromDouble(macrosDbl[i].value);
    if(pValue==NULL){
      Py_DECREF(pMacros_global);
      snprintf(message, ERR_STR_LEN-1, "Failed to build python object for double macro %s.", macrosDbl[i].name);
return write_local_err(PY_CONVERT_FAIL, message);
    }

    if(PyDict_SetItemString(pMacros_global, macrosDbl[i].name, pValue)){
      Py_DECREF(pValue); /* This is correct because PyDict_SetItemString() does NOT steal the reference (unlike with lists or tuples). */
      Py_DECREF(pMacros_global);
      snprintf(message, ERR_STR_LEN-1, "Failed to set double macro %s.", macrosDbl[i].name);
return write_local_err(PY_DICT_SET_FAIL, message);
    }

    Py_DECREF(pValue);
  }

  for(i=0;i<nIntMacros;i++){
    pValue = Py_BuildValue("i", macrosInt[i].value);
    if(pValue==NULL){
      Py_DECREF(pMacros_global);
      snprintf(message, ERR_STR_LEN-1, "Failed to build python object for integer macro %s.", macrosInt[i].name);
return write_local_err(PY_CONVERT_FAIL, message);
    }

    if(PyDict_SetItemString(pMacros_global, macrosInt[i].name, pValue)){
      Py_DECREF(pValue); /* This is correct because PyDict_SetItemString() does NOT steal the reference (unlike with lists or tuples). */
      Py_DECREF(pMacros_global);
      snprintf(message, ERR_STR_LEN-1, "Failed to set integer macro %s.", macrosInt[i].name);
return write_local_err(PY_DICT_SET_FAIL, message);
    }

    Py_DECREF(pValue);
  }

return err;
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

//*** does not seem to be used.

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
errType
_userFuncWrapper(PyObject *pFunc, const char *funcName\
  , double x, double y, double z, double *resultBuffer, int *numElemInUserFuncReturn){ 
  /*
This is a generic wrapper for LIME functions (density, temperature etc) supplied by the user in a python module. 

The status return codes are in error_codes.h.
  */

//***** can this be moved to a separate module? It seems only to be used by ml_models.c and py_models.c, neither of which use anything else from py_utils.c.

  int i;
  PyObject *pArgs,*pResult,*pListItem;
  errType err=init_local_err();
  char message[ERR_STR_LEN];

  pArgs = Py_BuildValue("(Offf)", pMacros_global, x, y, z);
  if(pArgs==NULL){
    snprintf(message, ERR_STR_LEN-1, "Could not build arguments object for user function %s", funcName);
return write_local_err(PY_BUILDVALUE_NULL, message);
  }else{

    pResult = PyObject_CallObject(pFunc, pArgs);
    Py_DECREF(pArgs);
    if(pResult==NULL){
      snprintf(message, ERR_STR_LEN-1, "User function %s returned a null result.", funcName);
return write_local_err(PY_CALLOBJECT_NULL, message);
    }
  }

  /* The returned value can be either a scalar or a list, but the datatypes should be floats or ints.
  */
  if(PyList_Check(pResult)){
    *numElemInUserFuncReturn = (int)PyList_Size(pResult);
    if(PyErr_Occurred() || (*numElemInUserFuncReturn)<=0){
      Py_DECREF(pResult);
      if(PyErr_Occurred()){
        snprintf(message, ERR_STR_LEN-1, "Error in calling user function %s", funcName);
return write_local_err(PY_LISTCHECK_ERR, message);
      }else{
        snprintf(message, ERR_STR_LEN-1, "User function %s returned an empty list", funcName);
return write_local_err(PY_EMPTY_LIST, message);
      }
    }

    if(*numElemInUserFuncReturn>UFUNC_BUFFER_SIZE){
      Py_DECREF(pResult);
      snprintf(message, ERR_STR_LEN-1, "User function %s returned %d items, more than the maximum allowed of %d."\
        , funcName, *numElemInUserFuncReturn, UFUNC_BUFFER_SIZE);
return write_local_err(PY_LIST_OVERFLOW, message);
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
        snprintf(message, ERR_STR_LEN-1, "Can't interpret item %d in the list returned by user function %s.", i, funcName);
return write_local_err(PY_NON_NUMERIC, message);
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
      snprintf(message, ERR_STR_LEN-1, "Can't interpret the value returned by user function %s.", funcName);
return write_local_err(PY_NON_NUMERIC, message);
    }
  }

  Py_DECREF(pResult);

return err;
}

/*....................................................................*/
errType
userFuncWrapper(PyObject *pFunc, const char *funcName, const int requiredNumElements\
  , double x, double y, double z, double *resultBuffer, int *numElemInUserFuncReturn){

  errType err=init_local_err();
  char message[ERR_STR_LEN];

  err = _userFuncWrapper(pFunc, funcName, x, y, z, resultBuffer, numElemInUserFuncReturn);

  if(err.status!=0)
return err;

  /* requiredNumElements set to <0 flags that we don't want to do this check. */
  if(requiredNumElements>=0 && *numElemInUserFuncReturn!=requiredNumElements){
    if(requiredNumElements==1)
      snprintf(message, ERR_STR_LEN-1, "User %s() function should return a scalar result.", funcName);
    else
      snprintf(message, ERR_STR_LEN-1, "User %s() function should return %d results.", funcName, requiredNumElements);
return write_local_err(1, message);
  }

return err;
}

/*....................................................................*/
void
pyFreeInputImgPars(image *inimg, const int nImages){
  /*
In 'standard' LIME the user copies the location of a (read-only) string to img[i].filename. Trying to free img[i].filename then results in an error. However the projected python version will malloc img[i].filename and copy string characters into that memory space. It is for this purpose that we preserve the present function.
  */
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
void
pyFreeInputPars(inputPars *par){
  /*
In 'standard' LIME the user copies the location of a (read-only) string to each of the char* elements of inputPars. Trying to free that element then results in an error. However the projected python version will malloc these elements and copy string characters into that memory space. It is for this purpose that we preserve the present function.
  */
  int i;

  free(par->collPartIds);
  free(par->nMolWeights);
  free(par->dustWeights);
  free(par->collPartMolWeights);

  free(par->gridDensMaxValues);
  free(par->gridDensMaxLoc);

  free(par->outputfile);
  free(par->binoutputfile);
  free(par->gridfile);
  free(par->pregrid);
  free(par->restart);
  free(par->dust);
  free(par->gridInFile);

  if(par->moldatfile!= NULL){
    for(i=0;i<MAX_NSPECIES;i++)
      free(par->moldatfile[i]);
    free(par->moldatfile);
  }
  if(par->girdatfile!= NULL){
    for(i=0;i<MAX_NSPECIES;i++)
      free(par->girdatfile[i]);
    free(par->girdatfile);
  }

  if(par->gridOutFiles!= NULL){
    for(i=0;i<NUM_GRID_STAGES;i++)
      free(par->gridOutFiles[i]);
    free(par->gridOutFiles);
  }

  if(par->collPartNames!= NULL){
    for(i=0;i<MAX_N_COLL_PART;i++)
      free(par->collPartNames[i]);
    free(par->collPartNames);
  }
}
