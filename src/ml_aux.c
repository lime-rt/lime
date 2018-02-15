/*
 *  ml_aux.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2017 The LIME development team
 *
TODO:
 */

#include "constants.h"
#include "py_utils.h"
#include "ml_types.h"
#include "ml_funcs.h"
#include "ml_models.h"

_Bool mla_doTest = FALSE;

/*....................................................................*/
int
extractParams(PyObject *pModelObj, char *errStr){
  /*
The input argument pModelObj should be an instance of the python class 'modellib_classes._Model'. This is expected to have an attribute 'paramList' which is a simple python list of objects, each of which is expected to be an instance of the class 'modellib_classes._ModelParam'. The purpose of the present function is to unpack these parameters and store them in lists modelIntPars, modelDblPars and modelParams, all of which are declared in ml_models.h.

These are parameters of the model-library recipes. Since we should not assume that they have the same order in the python code which sets them and supplies them via the Python object pModelObj to the present function, a key connecting the name of each parameter with its index in modelDblPars etc is stored in the array modelParams. The key is read in function ml_models.getParamI(), which is called from the *_onFinalizeConfiguration() function of the appropriate model recipe in ml_recipes/.

Note that errStr should be malloc'd (or declared as an array) by the calling routine with STR_LEN_0 bytes.
  */
  char *attrName="paramList",dType[MAX_LEN_PAR_TYPE+1],idStr[MAX_LEN_PAR_NAME+1],*attrNameValue="_value",*tempStr=NULL;
  PyObject *pAttr,*pParamObj,*pParamValue;
  int numElements=-1,numDoubles=0,numInts=0,numStrings=0,i,maxStrLen,status=0;

  /* Extract the parameters:
  */
  pAttr = PyObject_GetAttrString(pModelObj, attrName);
  if(pAttr==NULL){
    snprintf(errStr, STR_LEN_0, "Attribute %s not found in model instance.", attrName);
return 1;
  }

  if(!PyList_CheckExact(pAttr)){
    Py_DECREF(pAttr);
    snprintf(errStr, STR_LEN_0, "Attribute %s in model instance is not a list.", attrName);
return 2;
  }

  /* Find out first how many int, double and string parameters there are:
  */
  numElements = (int)PyList_Size(pAttr);
  if(numElements<0){
    Py_DECREF(pAttr);
    snprintf(errStr, STR_LEN_0, "Bad size of list %s in model instance.", attrName);
return 3;
  }

if(mla_doTest) printf("  +++ In ml_aux.extractParams(). Num elements in %s is %d\n", attrName, numElements);

  for(i=0;i<numElements;i++){
    pParamObj = PyList_GetItem(pAttr, (Py_ssize_t)i); /* Don't have to DECREF pParamObj, it is a borrowed reference. */

    if(!PyInstance_Check(pParamObj)){
      Py_DECREF(pAttr);
      snprintf(errStr, STR_LEN_0, "Element %d of %s is not an instance.", i, attrName);
return 4;
    }

    status = getStringAttribute(pParamObj, "dType", MAX_LEN_PAR_TYPE, dType);
    if(status!=0){
      Py_DECREF(pAttr);
      if(status==1){
        snprintf(errStr, STR_LEN_0, "Model parameter attribute 'dType' not read.");
        status = 5;
      }else if(status==2){
        snprintf(errStr, STR_LEN_0, "Model parameter attribute 'dType' is not a string.");
        status = 6;
      }else{
        snprintf(errStr, STR_LEN_0, "Call to read model parameter attribute 'dType' returned status %d", status);
        status = 7;
      }
return status;
    }

    if(strcmp(dType, "double")==0)
      numDoubles++;
    else if(strcmp(dType, "int")==0)
      numInts++;
    else if(strcmp(dType, "string")==0)
      numStrings++;
    else{
      Py_DECREF(pAttr);
      snprintf(errStr, STR_LEN_0, "Datatype %s not recognized.", dType);
return 8;
    }
  } /* end loop over paramList. */

  if(numStrings>1){
    Py_DECREF(pAttr);
    snprintf(errStr, STR_LEN_0, "Can only handle %d string parameter, not the %d supplied.", 1, numStrings);
return 9;
  }

  numModelParams = numElements; /* numModelParams is defined extern in the header of ml_models.h */

  modelParams = malloc(sizeof(*modelParams)*numModelParams); /* modelParams is defined extern in ml_models.h */
  if(modelParams==NULL){ /* malloc failed */
    Py_DECREF(pAttr);
    snprintf(errStr, STR_LEN_0, "Malloc of modelParams failed.");
return 10;
  }

  modelDblPars = malloc(sizeof(*modelDblPars)*numDoubles); // defined extern in ml_models
  if(modelDblPars==NULL){ /* malloc failed */
    Py_DECREF(pAttr);
    snprintf(errStr, STR_LEN_0, "Malloc of modelDblPars failed.");
return 11;
  }

  modelIntPars = malloc(sizeof(*modelIntPars)*numDoubles); // defined extern in ml_models
  if(modelIntPars==NULL){ /* malloc failed */
    Py_DECREF(pAttr);
    snprintf(errStr, STR_LEN_0, "Malloc of modelIntPars failed.");
return 12;
  }

  /* Now we know how many parameters of what type we have, we can copy them over. We don't have to duplicate error checks we've already done.
  */
  numDoubles = 0;
  numInts = 0;
  for(i=0;i<numElements;i++){
    pParamObj = PyList_GetItem(pAttr, (Py_ssize_t)i); /* Don't have to DECREF pParamObj, it is a borrowed reference. */

    status = getStringAttribute(pParamObj, "dType", MAX_LEN_PAR_TYPE, dType);
    if(status!=0){
      Py_DECREF(pAttr);
      if(status==1){
        snprintf(errStr, STR_LEN_0, "Model parameter attribute 'dType' not read.");
        status = 13;
      }else if(status==2){
        snprintf(errStr, STR_LEN_0, "Model parameter attribute 'dType' is not a string.");
        status = 14;
      }else{
        snprintf(errStr, STR_LEN_0, "Call to read model parameter attribute 'dType' returned status %d", status);
        status = 15;
      }
return status;
    }

    status = getStringAttribute(pParamObj, "idStr", MAX_LEN_PAR_NAME, idStr);
    if(status!=0){
      Py_DECREF(pAttr);
      if(status==1){
        snprintf(errStr, STR_LEN_0, "Model parameter attribute 'idStr' not read.");
        status = 16;
      }else if(status==2){
        snprintf(errStr, STR_LEN_0, "Model parameter attribute 'idStr' is not a string.");
        status = 17;
      }else{
        snprintf(errStr, STR_LEN_0, "Call to read model parameter attribute 'idStr' returned status %d", status);
        status = 18;
      }
return status;
    }

    pParamValue = PyObject_GetAttrString(pParamObj, attrNameValue);
    if(pParamValue==NULL){
      Py_DECREF(pAttr);
      snprintf(errStr, STR_LEN_0, "Attribute %s not found in model parameter instance.", attrNameValue);
return 19;
    }

    myStrNCpy(modelParams[i].name,  idStr, MAX_LEN_PAR_NAME); /* Note that the .name and .dtype attributes */
    myStrNCpy(modelParams[i].dtype, dType, MAX_LEN_PAR_TYPE); /* are dimensioned to MAX_LEN_PAR_TYPE+1. */

    if(strcmp(dType, "double")==0){
      if(!(PyFloat_CheckExact(pParamValue) || PyInt_CheckExact(pParamValue))){ /* We also allow ints for this */
        Py_DECREF(pParamValue);
        Py_DECREF(pAttr);
        snprintf(errStr, STR_LEN_0, "Parameter %s flagged as double is not convertable to double.", idStr);
return 20;
      }

      if(PyInt_CheckExact(pParamValue))
        modelDblPars[numDoubles] = (double)PyInt_AS_LONG(pParamValue);
      else
        modelDblPars[numDoubles] = PyFloat_AS_DOUBLE(pParamValue);

      modelParams[i].index = numDoubles;
      numDoubles++;

    }else if(strcmp(dType, "int")==0){
      if(!PyInt_CheckExact(pParamValue)){
        Py_DECREF(pParamValue);
        Py_DECREF(pAttr);
        snprintf(errStr, STR_LEN_0, "Parameter %s flagged as int is not convertable to int.", idStr);
return 21;
      }

      if(PyInt_CheckExact(pParamValue))
        modelIntPars[numInts] = (int)PyInt_AS_LONG(pParamValue);

      modelParams[i].index = numInts;
      numInts++;

    }else{ /* assume dType is "string". */
      if(!PyString_CheckExact(pParamValue)){
        Py_DECREF(pParamValue);
        Py_DECREF(pAttr);
        snprintf(errStr, STR_LEN_0, "Parameter %s flagged as string is not in fact a string.", idStr);
return 22;
      }
      tempStr = PyString_AsString(pParamValue); /* This should NOT be deallocated or written to. */
      if(tempStr==NULL){
        Py_DECREF(pParamValue);
        Py_DECREF(pAttr);
        snprintf(errStr, STR_LEN_0, "String parameter returned NULL.");
return 23;
      }

      maxStrLen = strlen(tempStr);
      modelStrPar = malloc(sizeof(*modelStrPar)*(maxStrLen+1)); /* modelStrPar defined extern in ml_models. */
      if(modelStrPar==NULL){ /* malloc failed */
        Py_DECREF(pParamValue);
        Py_DECREF(pAttr);
        snprintf(errStr, STR_LEN_0, "Malloc of modelStrPar failed.");
return 24;
      }

      myStrNCpy(modelStrPar, tempStr, maxStrLen); /* Last use of tempStr. */

      modelParams[i].index = 0;
    }

    Py_DECREF(pParamValue); /* .paramList[i]._value */
  } /* end loop over paramList. */

  Py_DECREF(pAttr); /* .paramList */

return 0;
}

/*....................................................................*/
int
_extractFuncParams(PyObject *pFuncObj, const int resultI, char *funcName\
  , const int funcI, const int funcTypeI, char *errStr){
  /*
These are parameters of the functions listed in ml_funcs.c as _scalar*() or _vector(). Since we should not assume that they have the same order in the python code which sets them and supplies them via the Python object pFuncObj to the present function, a key connecting the name of each parameter with its index in modelDblPars etc is available via the function getFuncParIndexFromName().

Note that errStr should be malloc'd (or declared as an array) by the calling routine with STR_LEN_0 bytes.
  */
  char *attrName="_argsDict",*attrNameValue="_value",funcParName[MAX_LEN_PAR_NAME+1];
  PyObject *pAttr,*pFuncParamValue,*pParamName,*pParamObj;
  int numPars,expectedNumPars,i,parIndex,status=0;
  Py_ssize_t pos=0;

  pAttr = PyObject_GetAttrString(pFuncObj, attrName);
  if(pAttr==NULL){
    snprintf(errStr, STR_LEN_0, "Attribute %s not found in function instance.", attrName);
return 1;
  }

  if(!PyDict_CheckExact(pAttr)){
    Py_DECREF(pAttr);
    snprintf(errStr, STR_LEN_0, "Attribute %s in function instance is not a dictionary.", attrName);
return 2;
  }

  /* Check that the dictionary has the right number of key:value pairs:
  */
  numPars = (int)PyDict_Size(pAttr);

  if(funcTypeI==FUNC_scalar)
    expectedNumPars = numScalarFuncPars[funcI]; /* defined extern in ml_funcs.c */
  else if(funcTypeI==FUNC_vector)
    expectedNumPars = numVectorFuncPars[funcI]; /* defined extern in ml_funcs.c */
  else{
    Py_DECREF(pAttr);
    snprintf(errStr, STR_LEN_0, "Unrecognized function type integer %d", funcTypeI);
return 3;
  }

  if(numPars!=expectedNumPars){
    Py_DECREF(pAttr);
    snprintf(errStr, STR_LEN_0, "Expected %d parameters for function %s but found %d.", expectedNumPars, funcName, numPars);
return 4;
  }

  /* Global variable defined in ml_funcs.h.
  */
  funcPars[currentModelI][resultI] = malloc(sizeof(***funcPars)*numPars);
  if(funcPars[currentModelI][resultI]==NULL){ /* malloc failed */
    Py_DECREF(pAttr);
    snprintf(errStr, STR_LEN_0, "Malloc of funcPars[%d][%d] failed.", currentModelI, resultI);
return 5;
  }

  pos = 0;
  while(PyDict_Next(pAttr, &pos, &pParamName, &pParamObj)){ /* pParamName, pParamObj are borrowed references - don't need to be DECREF'd. */
    /* Get the key string - should be the name of a result.
    */
    if(!PyString_CheckExact(pParamName)){
      Py_DECREF(pAttr);
      snprintf(errStr, STR_LEN_0, "Non-string key found for %s.", attrName);
return 6;
    }

    if(!PyInstance_Check(pParamObj)){
      Py_DECREF(pAttr);
      snprintf(errStr, STR_LEN_0, "Element %d of %s is not an instance.", i, attrName);
return 7;
    }

    status = getStringAttribute(pParamObj, "name", MAX_LEN_PAR_NAME, funcParName);
    if(status!=0){
      Py_DECREF(pAttr);
      if(status==1){
        snprintf(errStr, STR_LEN_0, "Model parameter attribute 'name' not read.");
        status = 8;
      }else if(status==2){
        snprintf(errStr, STR_LEN_0, "Model parameter attribute 'name' is not a string.");
        status = 9;
      }else{
        snprintf(errStr, STR_LEN_0, "Call to read model parameter attribute 'name' returned status %d", status);
        status = 10;
      }
return status;
    }

    pFuncParamValue = PyObject_GetAttrString(pParamObj, attrNameValue);
    if(pFuncParamValue==NULL){
      Py_DECREF(pAttr);
      snprintf(errStr, STR_LEN_0, "Attribute %s not found in function parameter instance.", attrNameValue);
return 11;
    }

    if(!(PyFloat_CheckExact(pFuncParamValue) || PyInt_CheckExact(pFuncParamValue))){ /* We also allow ints for this */
      Py_DECREF(pFuncParamValue);
      Py_DECREF(pAttr);
      snprintf(errStr, STR_LEN_0, "Function parameter type is not convertable to double.");
return 12;
    }

    status = getFuncParIndexFromName(funcI, funcTypeI, funcParName, &parIndex);
    if(status){
      Py_DECREF(pFuncParamValue);
      Py_DECREF(pAttr);
      snprintf(errStr, STR_LEN_0, "getFuncParIndexFromName() returned status %d", status);
return 13;
    }

    if(PyInt_CheckExact(pFuncParamValue))
      funcPars[currentModelI][resultI][parIndex] = (double)PyInt_AS_LONG(pFuncParamValue);
    else
      funcPars[currentModelI][resultI][parIndex] = PyFloat_AS_DOUBLE(pFuncParamValue);

    Py_DECREF(pFuncParamValue); /* .funcDict[resultName]._argList[i]._value */
  } /* end loop over .funcDict[resultName]._argList. */

  Py_DECREF(pAttr); /* .funcDict[resultName]._argList */

return 0;
}

/*....................................................................*/
int
extractFuncs(PyObject *pModelObj, char *errStr){
  char *attrName="funcDict",resultName[MAX_LEN_RESULT_NAME+1],funcTypeName[MAX_LEN_PAR_TYPE+1],funcName[MAX_LEN_PAR_NAME+1];
  PyObject *pAttr,*pResultName,*pFuncObj;
  Py_ssize_t pos=0;
  int resultI,funcTypeI,funcI=0,status=0;
  char message[STR_LEN_0];

if(mla_doTest) printf("  >>> Entering ml_aux.extractFuncs()\n");

  /* Time now to set up any functions.
  */
  pAttr = PyObject_GetAttrString(pModelObj, attrName);
  if(pAttr==NULL){
    snprintf(errStr, STR_LEN_0, "Attribute %s not found in model instance.", attrName);
return 1;
  }

  if(!PyDict_CheckExact(pAttr)){
    Py_DECREF(pAttr);
    snprintf(errStr, STR_LEN_0, "Attribute %s in model instance is not a dictionary.", attrName);
return 2;
  }

  pos = 0;
  while(PyDict_Next(pAttr, &pos, &pResultName, &pFuncObj)){ /* pResultName, pFuncObj are borrowed references - don't need to be DECREF'd. */
    /* Get the key string - should be the name of a result.
    */
    if(!PyString_CheckExact(pResultName)){
      Py_DECREF(pAttr);
      snprintf(errStr, STR_LEN_0, "Non-string key found for %s.", attrName);
return 3;
    }

    myStrNCpy(resultName, PyString_AsString(pResultName), MAX_LEN_RESULT_NAME);
    resultI = getResultIFromName(resultName);
    if(resultI<0){
      Py_DECREF(pAttr);
      snprintf(errStr, STR_LEN_0, "Result name %s is not recognized.", resultName);
return 4;
    }

    /* Now get the value - should be a _Function object.
    */
    if(!PyInstance_Check(pFuncObj)){
      Py_DECREF(pAttr);
      snprintf(errStr, STR_LEN_0, "Value for key %s is not an instance.", resultName);
return 5;
    }

    /* Get the type ('scalar' or 'vector') of the function.
    */
    status = getStringAttribute(pFuncObj, "typeStr", MAX_LEN_PAR_TYPE, funcTypeName);
    if(status!=0){
      Py_DECREF(pAttr);
      if(status==1){
        snprintf(errStr, STR_LEN_0, "Function object attribute 'typeStr' not read.");
        status = 6;
      }else if(status==2){
        snprintf(errStr, STR_LEN_0, "Function object attribute 'typeStr' is not a string.");
        status = 7;
      }else{
        snprintf(errStr, STR_LEN_0, "Call to read function object attribute 'typeStr' returned status %d", status);
        status = 8;
      }
return status;
    }

    funcTypeI = getFuncTypeIFromName(funcTypeName);
    if(funcTypeI<0){
      Py_DECREF(pAttr);
      snprintf(errStr, STR_LEN_0, "Unrecognized function type %s", funcTypeName);
return 9;
    }

    /* Get the name of the function.
    */
    status = getStringAttribute(pFuncObj, "idStr", MAX_LEN_PAR_NAME, funcName);
    if(status!=0){
      Py_DECREF(pAttr);
      if(status==1){
        snprintf(errStr, STR_LEN_0, "Function object attribute 'idStr' not read.");
        status = 10;
      }else if(status==2){
        snprintf(errStr, STR_LEN_0, "Function object attribute 'idStr' is not a string.");
        status = 11;
      }else{
        snprintf(errStr, STR_LEN_0, "Call to read function object attribute 'idStr' returned status %d", status);
        status = 12;
      }
return status;
    }

    /* Convert the function name to an ID integer.
    */
    if(funcTypeI==FUNC_scalar)
      funcI = getScalarFuncIFromName(funcName);
    else if(funcTypeI==FUNC_vector)
      funcI = getVectorFuncIFromName(funcName);
    else{
      Py_DECREF(pAttr);
      snprintf(errStr, STR_LEN_0, "Unrecognized function type integer %d", funcTypeI);
return 13;
    }

    if(funcI<0){
      Py_DECREF(pAttr);
      snprintf(errStr, STR_LEN_0, "Unrecognized function %s", funcName);
return 14;
    }

    /* Global variables defined in ml_funcs.h.
    */
    funcIs[    currentModelI][resultI] = funcI;
    funcTypeIs[currentModelI][resultI] = funcTypeI;

    /* Now we can copy over the function parameters. These are in the _argsDict attribute of the function object.
    */
    status = _extractFuncParams(pFuncObj, resultI, funcName, funcI, funcTypeI, message);
    if(status!=0){
      snprintf(errStr, STR_LEN_0, "_extractFuncParams() returned status %d", status);
return 15;
    }
  } /* end loop over funcDict keys. */

  Py_DECREF(pAttr); /* .funcDict */

if(mla_doTest) printf("  <<< Leaving ml_aux.extractFuncs()\n");
return 0;
}

/*....................................................................*/
int
getModelI(PyObject *pModelObj, int *modelI, char *errMsg){
  /* Note that errMsg is expected to have been malloc'd (or to be an array) to length STR_LEN_0. */
  char modelName[MAX_LEN_MODEL_NAME+1];
  int status=0;

  *modelI = -1;

  status = getStringAttribute(pModelObj, "idStr", MAX_LEN_MODEL_NAME, modelName); /* in py_utils.c */
  if(status!=0){
    if(status==1){
      snprintf(errMsg, STR_LEN_0, "Model object attribute 'idStr' not read.");
    }else if(status==2){
      snprintf(errMsg, STR_LEN_0, "Model object attribute 'idStr' is not a string.");
    }else{
      snprintf(errMsg, STR_LEN_0, "Call to read model object attribute 'idStr' returned status %d", status);
    }
return status;
  }

  *modelI = getModelIFromName(modelName); /* in ml_models.c */
  if(*modelI<0){
    snprintf(errMsg, STR_LEN_0, "Model name %s was not recognized.", modelName);
return 3;
  }

return 0;
}

