/*
 *  ml_aux.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#include "constants.h"
#include "py_utils.h" /* For myStrNCpy */
#include "ml_types.h"
#include "ml_funcs.h"
#include "ml_models.h"
#include "local_err.h"

/*....................................................................*/
errType
_getStringAttribute(PyObject *pInstance, char *attrName, const int maxStrLenAttr\
  , char *attrStrValue){

  PyObject *pAttr;
  errType err=init_local_err();
  char message[ERR_STR_LEN];

  attrStrValue[0] = '\0';

  pAttr = PyObject_GetAttrString(pInstance, attrName);
  if(pAttr==NULL){
    snprintf(message, ERR_STR_LEN-1, "Model parameter attribute '%s' not read.", attrName);
return write_local_err(1, message);
  }

  if(!PyString_CheckExact(pAttr)){
    Py_DECREF(pAttr);
    snprintf(message, ERR_STR_LEN-1, "Model parameter attribute '%s' is not a string.", attrName);
return write_local_err(2, message);
  }

  myStrNCpy(attrStrValue, PyString_AsString(pAttr), maxStrLenAttr);
  Py_DECREF(pAttr);

return err;
}

/*....................................................................*/
errType
extractParams(PyObject *pModelObj){
  /*
The input argument pModelObj should be an instance of the python class 'modellib_classes._Model'. This is expected to have an attribute 'paramList' which is a simple python list of objects, each of which is expected to be an instance of the class 'modellib_classes._ModelParam'. The purpose of the present function is to unpack these parameters and store them in lists modelIntPars, modelDblPars and modelParams, all of which are declared in ml_models.h.

These are parameters of the model-library recipes. Since we should not assume that they have the same order in the python code which sets them and supplies them via the Python object pModelObj to the present function, a key connecting the name of each parameter with its index in modelDblPars etc is stored in the array modelParams. The key is read in function ml_models.getParamI(), which is called from the *_onFinalizeConfiguration() function of the appropriate model recipe in ml_recipes/.

Note that errStr should be malloc'd (or declared as an array) by the calling routine with STR_LEN_0 bytes.
  */
  char *attrName="paramList",dType[MAX_LEN_PAR_TYPE+1],idStr[MAX_LEN_PAR_NAME+1],*attrNameValue="_value",*tempStr=NULL;
  PyObject *pAttr,*pParamObj,*pParamValue;
  int numElements=-1,numDoubles=0,numInts=0,numStrings=0,i,maxStrLen,status=0;
  errType err=init_local_err();
  char message[ERR_STR_LEN];

  /* Extract the parameters:
  */
  pAttr = PyObject_GetAttrString(pModelObj, attrName);
  if(pAttr==NULL){
    snprintf(message, ERR_STR_LEN-1, "Attribute %s not found in model instance.", attrName);
return write_local_err(1, message);
  }

  if(!PyList_CheckExact(pAttr)){
    Py_DECREF(pAttr);
    snprintf(message, ERR_STR_LEN-1, "Attribute %s in model instance is not a list.", attrName);
return write_local_err(2, message);
  }

  /* Find out first how many int, double and string parameters there are:
  */
  numElements = (int)PyList_Size(pAttr);
  if(numElements<0){
    Py_DECREF(pAttr);
    snprintf(message, ERR_STR_LEN-1, "Bad size of list %s in model instance.", attrName);
return write_local_err(3, message);
  }

  for(i=0;i<numElements;i++){
    pParamObj = PyList_GetItem(pAttr, (Py_ssize_t)i); /* Don't have to DECREF pParamObj, it is a borrowed reference. */

    if(!PyInstance_Check(pParamObj)){
      Py_DECREF(pAttr);
      snprintf(message, ERR_STR_LEN-1, "Element %d of %s is not an instance.", i, attrName);
return write_local_err(4, message);
    }

    err = _getStringAttribute(pParamObj, "dType", MAX_LEN_PAR_TYPE, dType);
    if(err.status!=0){
      Py_DECREF(pAttr);
return err;
    }

    if(strcmp(dType, "double")==0)
      numDoubles++;
    else if(strcmp(dType, "int")==0)
      numInts++;
    else if(strcmp(dType, "string")==0)
      numStrings++;
    else{
      Py_DECREF(pAttr);
      snprintf(message, ERR_STR_LEN-1, "Datatype %s not recognized.", dType);
return write_local_err(8, message);
    }
  } /* end loop over paramList. */

  if(numStrings>1){
    Py_DECREF(pAttr);
    snprintf(message, ERR_STR_LEN-1, "Can only handle %d string parameter, not the %d supplied.", 1, numStrings);
return write_local_err(9, message);
  }

  numModelParams = numElements; /* numModelParams is defined extern in the header of ml_models.h */

  modelParams = malloc(sizeof(*modelParams)*numModelParams); /* modelParams is defined extern in ml_models.h */
  if(modelParams==NULL){ /* malloc failed */
    Py_DECREF(pAttr);
    snprintf(message, ERR_STR_LEN-1, "Malloc of modelParams failed.");
return write_local_err(10, message);
  }

  modelDblPars = malloc(sizeof(*modelDblPars)*numDoubles); // defined extern in ml_models
  if(modelDblPars==NULL){ /* malloc failed */
    Py_DECREF(pAttr);
    snprintf(message, ERR_STR_LEN-1, "Malloc of modelDblPars failed.");
return write_local_err(11, message);
  }

  modelIntPars = malloc(sizeof(*modelIntPars)*numDoubles); // defined extern in ml_models
  if(modelIntPars==NULL){ /* malloc failed */
    Py_DECREF(pAttr);
    snprintf(message, ERR_STR_LEN-1, "Malloc of modelIntPars failed.");
return write_local_err(12, message);
  }

  /* Now we know how many parameters of what type we have, we can copy them over. We don't have to duplicate error checks we've already done.
  */
  numDoubles = 0;
  numInts = 0;
  for(i=0;i<numElements;i++){
    pParamObj = PyList_GetItem(pAttr, (Py_ssize_t)i); /* Don't have to DECREF pParamObj, it is a borrowed reference. */

    err = _getStringAttribute(pParamObj, "dType", MAX_LEN_PAR_TYPE, dType);
    if(err.status!=0){
      Py_DECREF(pAttr);
return err;
    }

    err = _getStringAttribute(pParamObj, "idStr", MAX_LEN_PAR_NAME, idStr);
    if(status!=0){
      Py_DECREF(pAttr);
return err;
    }

    pParamValue = PyObject_GetAttrString(pParamObj, attrNameValue);
    if(pParamValue==NULL){
      Py_DECREF(pAttr);
      snprintf(message, ERR_STR_LEN-1, "Attribute %s not found in model parameter instance.", attrNameValue);
return write_local_err(19, message);
    }

    myStrNCpy(modelParams[i].name,  idStr, MAX_LEN_PAR_NAME); /* Note that the .name and .dtype attributes */
    myStrNCpy(modelParams[i].dtype, dType, MAX_LEN_PAR_TYPE); /* are dimensioned to MAX_LEN_PAR_TYPE+1. */

    if(strcmp(dType, "double")==0){
      if(!(PyFloat_CheckExact(pParamValue) || PyInt_CheckExact(pParamValue))){ /* We also allow ints for this */
        Py_DECREF(pParamValue);
        Py_DECREF(pAttr);
        snprintf(message, ERR_STR_LEN-1, "Parameter %s flagged as double is not convertable to double.", idStr);
return write_local_err(20, message);
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
        snprintf(message, ERR_STR_LEN-1, "Parameter %s flagged as int is not convertable to int.", idStr);
return write_local_err(21, message);
      }

      if(PyInt_CheckExact(pParamValue))
        modelIntPars[numInts] = (int)PyInt_AS_LONG(pParamValue);

      modelParams[i].index = numInts;
      numInts++;

    }else{ /* assume dType is "string". */
      if(!PyString_CheckExact(pParamValue)){
        Py_DECREF(pParamValue);
        Py_DECREF(pAttr);
        snprintf(message, ERR_STR_LEN-1, "Parameter %s flagged as string is not in fact a string.", idStr);
return write_local_err(22, message);
      }
      tempStr = PyString_AsString(pParamValue); /* This should NOT be deallocated or written to. */
      if(tempStr==NULL){
        Py_DECREF(pParamValue);
        Py_DECREF(pAttr);
        snprintf(message, ERR_STR_LEN-1, "String parameter returned NULL.");
return write_local_err(23, message);
      }

      maxStrLen = strlen(tempStr);
      modelStrPar = malloc(sizeof(*modelStrPar)*(maxStrLen+1)); /* modelStrPar defined extern in ml_models. */
      if(modelStrPar==NULL){ /* malloc failed */
        Py_DECREF(pParamValue);
        Py_DECREF(pAttr);
        snprintf(message, ERR_STR_LEN-1, "Malloc of modelStrPar failed.");
return write_local_err(24, message);
      }

      myStrNCpy(modelStrPar, tempStr, maxStrLen); /* Last use of tempStr. */

      modelParams[i].index = 0;
    }

    Py_DECREF(pParamValue); /* .paramList[i]._value */
  } /* end loop over paramList. */

  Py_DECREF(pAttr); /* .paramList */

return err;
}

/*....................................................................*/
errType
_extractFuncParams(PyObject *pFuncObj, const int resultI, char *funcName\
  , const int funcI, const int funcTypeI){
  /*
These are parameters of the functions listed in ml_funcs.c as _scalar*() or _vector(). Since we should not assume that they have the same order in the python code which sets them and supplies them via the Python object pFuncObj to the present function, a key connecting the name of each parameter with its index in modelDblPars etc is available via the function getFuncParIndexFromName().

Note that errStr should be malloc'd (or declared as an array) by the calling routine with STR_LEN_0 bytes.
  */
  char *attrName="_argsDict",*attrNameValue="_value",funcParName[MAX_LEN_PAR_NAME+1];
  PyObject *pAttr,*pFuncParamValue,*pParamName,*pParamObj;
  int numPars,expectedNumPars,i,parIndex,status=0;
  Py_ssize_t pos=0;
  errType err=init_local_err();
  char message[ERR_STR_LEN];

  pAttr = PyObject_GetAttrString(pFuncObj, attrName);
  if(pAttr==NULL){
    snprintf(message, ERR_STR_LEN-1, "Attribute %s not found in function instance.", attrName);
return write_local_err(1, message);
  }

  if(!PyDict_CheckExact(pAttr)){
    Py_DECREF(pAttr);
    snprintf(message, ERR_STR_LEN-1, "Attribute %s in function instance is not a dictionary.", attrName);
return write_local_err(2, message);
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
    snprintf(message, ERR_STR_LEN-1, "Unrecognized function type integer %d", funcTypeI);
return write_local_err(3, message);
  }

  if(numPars!=expectedNumPars){
    Py_DECREF(pAttr);
    snprintf(message, ERR_STR_LEN-1, "Expected %d parameters for function %s but found %d.", expectedNumPars, funcName, numPars);
return write_local_err(4, message);
  }

  /* Global variable defined in ml_funcs.h.
  */
  funcPars[currentModelI][resultI] = malloc(sizeof(***funcPars)*numPars);
  if(funcPars[currentModelI][resultI]==NULL){ /* malloc failed */
    Py_DECREF(pAttr);
    snprintf(message, ERR_STR_LEN-1, "Malloc of funcPars[%d][%d] failed.", currentModelI, resultI);
return write_local_err(5, message);
  }

  pos = 0;
  while(PyDict_Next(pAttr, &pos, &pParamName, &pParamObj)){ /* pParamName, pParamObj are borrowed references - don't need to be DECREF'd. */
    /* Get the key string - should be the name of a result.
    */
    if(!PyString_CheckExact(pParamName)){
      Py_DECREF(pAttr);
      snprintf(message, ERR_STR_LEN-1, "Non-string key found for %s.", attrName);
return write_local_err(6, message);
    }

    if(!PyInstance_Check(pParamObj)){
      Py_DECREF(pAttr);
      snprintf(message, ERR_STR_LEN-1, "Element %d of %s is not an instance.", i, attrName);
return write_local_err(7, message);
    }

    err = _getStringAttribute(pParamObj, "name", MAX_LEN_PAR_NAME, funcParName);
    if(err.status!=0){
      Py_DECREF(pAttr);
return err;
    }

    pFuncParamValue = PyObject_GetAttrString(pParamObj, attrNameValue);
    if(pFuncParamValue==NULL){
      Py_DECREF(pAttr);
      snprintf(message, ERR_STR_LEN-1, "Attribute %s not found in function parameter instance.", attrNameValue);
return write_local_err(11, message);
    }

    if(!(PyFloat_CheckExact(pFuncParamValue) || PyInt_CheckExact(pFuncParamValue))){ /* We also allow ints for this */
      Py_DECREF(pFuncParamValue);
      Py_DECREF(pAttr);
      snprintf(message, ERR_STR_LEN-1, "Function parameter type is not convertable to double.");
return write_local_err(12, message);
    }

    status = getFuncParIndexFromName(funcI, funcTypeI, funcParName, &parIndex);
    if(status){
      Py_DECREF(pFuncParamValue);
      Py_DECREF(pAttr);
      snprintf(message, ERR_STR_LEN-1, "getFuncParIndexFromName() returned status %d", status);
return write_local_err(13, message);
    }

    if(PyInt_CheckExact(pFuncParamValue))
      funcPars[currentModelI][resultI][parIndex] = (double)PyInt_AS_LONG(pFuncParamValue);
    else
      funcPars[currentModelI][resultI][parIndex] = PyFloat_AS_DOUBLE(pFuncParamValue);

    Py_DECREF(pFuncParamValue); /* .funcDict[resultName]._argList[i]._value */
  } /* end loop over .funcDict[resultName]._argList. */

  Py_DECREF(pAttr); /* .funcDict[resultName]._argList */

return err;
}

/*....................................................................*/
errType
extractFuncs(PyObject *pModelObj){
  char *attrName="funcDict",resultName[MAX_LEN_RESULT_NAME+1],funcTypeName[MAX_LEN_PAR_TYPE+1],funcName[MAX_LEN_PAR_NAME+1];
  PyObject *pAttr,*pResultName,*pFuncObj;
  Py_ssize_t pos=0;
  int resultI,funcTypeI,funcI=0,status=0;
  char message[ERR_STR_LEN];
  errType err=init_local_err();

  /* Time now to set up any functions.
  */
  pAttr = PyObject_GetAttrString(pModelObj, attrName);
  if(pAttr==NULL){
    snprintf(message, ERR_STR_LEN-1, "Attribute %s not found in model instance.", attrName);
return write_local_err(1, message);
  }

  if(!PyDict_CheckExact(pAttr)){
    Py_DECREF(pAttr);
    snprintf(message, ERR_STR_LEN-1, "Attribute %s in model instance is not a dictionary.", attrName);
return write_local_err(2, message);
  }

  pos = 0;
  while(PyDict_Next(pAttr, &pos, &pResultName, &pFuncObj)){ /* pResultName, pFuncObj are borrowed references - don't need to be DECREF'd. */
    /* Get the key string - should be the name of a result.
    */
    if(!PyString_CheckExact(pResultName)){
      Py_DECREF(pAttr);
      snprintf(message, ERR_STR_LEN-1, "Non-string key found for %s.", attrName);
return write_local_err(3, message);
    }

    myStrNCpy(resultName, PyString_AsString(pResultName), MAX_LEN_RESULT_NAME);
    resultI = getResultIFromName(resultName);
    if(resultI<0){
      Py_DECREF(pAttr);
      snprintf(message, ERR_STR_LEN-1, "Result name %s is not recognized.", resultName);
return write_local_err(4, message);
    }

    /* Now get the value - should be a _Function object.
    */
    if(!PyInstance_Check(pFuncObj)){
      Py_DECREF(pAttr);
      snprintf(message, ERR_STR_LEN-1, "Value for key %s is not an instance.", resultName);
return write_local_err(5, message);
    }

    /* Get the type ('scalar' or 'vector') of the function.
    */
    err = _getStringAttribute(pFuncObj, "typeStr", MAX_LEN_PAR_TYPE, funcTypeName);
    if(err.status!=0){
      Py_DECREF(pAttr);
return err;
    }

    funcTypeI = getFuncTypeIFromName(funcTypeName);
    if(funcTypeI<0){
      Py_DECREF(pAttr);
      snprintf(message, ERR_STR_LEN-1, "Unrecognized function type %s", funcTypeName);
return write_local_err(9, message);
    }

    /* Get the name of the function.
    */
    err = _getStringAttribute(pFuncObj, "idStr", MAX_LEN_PAR_NAME, funcName);
    if(status!=0){
      Py_DECREF(pAttr);
return err;
    }

    /* Convert the function name to an ID integer.
    */
    if(funcTypeI==FUNC_scalar)
      funcI = getScalarFuncIFromName(funcName);
    else if(funcTypeI==FUNC_vector)
      funcI = getVectorFuncIFromName(funcName);
    else{
      Py_DECREF(pAttr);
      snprintf(message, ERR_STR_LEN-1, "Unrecognized function type integer %d", funcTypeI);
return write_local_err(13, message);
    }

    if(funcI<0){
      Py_DECREF(pAttr);
      snprintf(message, ERR_STR_LEN-1, "Unrecognized function %s", funcName);
return write_local_err(14, message);
    }

    /* Global variables defined in ml_funcs.h.
    */
    funcIs[    currentModelI][resultI] = funcI;
    funcTypeIs[currentModelI][resultI] = funcTypeI;

    /* Now we can copy over the function parameters. These are in the _argsDict attribute of the function object.
    */
    err = _extractFuncParams(pFuncObj, resultI, funcName, funcI, funcTypeI);//, message);
    if(err.status!=0){
return write_local_err(15, message);
    }
  } /* end loop over funcDict keys. */

  Py_DECREF(pAttr); /* .funcDict */

return err;
}

/*....................................................................*/
errType
getModelI(PyObject *pModelObj, int *modelI){
  /* Note that errMsg is expected to have been malloc'd (or to be an array) to length STR_LEN_0. */
  char modelName[MAX_LEN_MODEL_NAME+1];
  errType err=init_local_err();
  char message[ERR_STR_LEN];

  *modelI = -1;

  err = _getStringAttribute(pModelObj, "idStr", MAX_LEN_MODEL_NAME, modelName);
  if(err.status!=0){
return err;
  }

  *modelI = getModelIFromName(modelName); /* in ml_models.c */
  if(*modelI<0){
    snprintf(message, ERR_STR_LEN-1, "Model name %s was not recognized.", modelName);
return write_local_err(3, message);
  }

return err;
}

