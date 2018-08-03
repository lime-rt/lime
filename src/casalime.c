/*
 *  casalime.c
 *  This file is part of LIME, the versatile line modeling engine.
 *
 *  See ../COPYRIGHT
 *
 */

#include "lime.h"
#include "error_codes.h"
#include "local_err.h"
#include "py_lime.h"
#include "py_utils.h"
#include "ml_types.h"
#include "ml_funcs.h"
#include "ml_models.h"

int copyTemp = 0; /* <0 means tgas should be set to equal tdust, >0 means the other way, ==0 means no action. */
int silent = 0;
int defaultFuncFlags = 0;

/*....................................................................*/
void
_readParsObjWrapper(const char *pklFileName, PyObject **pCurrentModel\
  , char *userModuleNameNoSuffix, const int maxLenName, int *copyTemperatureI, PyObject **pLimePars){
  /*
The user needs to pass four pieces of information to the present program:
  - The modellib model (an instance of the class modellib_classes.py:_Model).
  - The name (if any) of the python module which contains additional 'result' function specifications.
  - Which way to copy dust/gas temperatures, if at all.
  - The LIME parameters (an instance of the class limepar_classes.py:ModelParameters).

These are constructed as python objects by the top-layer python script and written by this script to a temporary pickle file. The script starts the present program (casalime) in a separate process and passes it the name of the pickle file. Casalime reads this file in the present routine via a function 'casalime_read.py:readPars()'. The three pieces of information required are returned as a tuple. They are then extracted and returned to the calling routine.

NOTE that I have used sequence-access routines here rather than tuple-access ones. This is because sequence access returns new references for the contents whereas tuple access only returns borrowed references. In practice this seems to mean that when I decref the tuple, the references to its contents disappear as well. Of course this means that the three objects extracted (*pCurrentModel, pUserModuleName, *pLimePars) must now be decrefed, 2 of them by the calling routine.
  */

  errType err=init_local_err();
  const int expectedNTupleMembers=4;
  const char *moduleNameNoSuffix="casalime_read",*funcName="readPars";
  int nTupleMembers=0;
  PyObject *pModule,*pFunc,*pArgs,*pResult,*pCopyTemperatureI,*pUserModuleName;
  char message[STR_LEN_1],*tempStr;

  err = getModuleFromName(moduleNameNoSuffix, &pModule); /* In py_utils.c */
  if(err.status!=0){ /* Don't need to decref pModule. */
pyerror(err.message);
  }

  getPythonFunc(pModule, funcName, &pFunc); /* In py_utils.c */
  if(pFunc==NULL){
    Py_DECREF(pModule);
    printOrClearPyError();
    snprintf(message, STR_LEN_1, "Error obtaining function %s from module %s.", funcName, moduleNameNoSuffix);
pyerror(message);
  }

  pArgs = Py_BuildValue("(s)", pklFileName);
  if(pArgs==NULL){
    Py_DECREF(pFunc);
    Py_DECREF(pModule);
    printOrClearPyError();
    snprintf(message, STR_LEN_1, "Error building arguments for %s function.", funcName);
pyerror(message);
  }

  pResult = PyObject_CallObject(pFunc, pArgs); /* Should return a tuple with three elements. */
  Py_DECREF(pArgs);
  Py_DECREF(pFunc);
  Py_DECREF(pModule);

  if(pResult==NULL){
    printOrClearPyError();
    snprintf(message, STR_LEN_1, "The %s function returned a NULL object rather than the expected tuple.", funcName);
pyerror(message);
  }

  if(!PySequence_Check(pResult)){
    Py_DECREF(pResult);
    snprintf(message, STR_LEN_1, "The %s function did not return a valid tuple as expected.", funcName);
pyerror(message);
  }

  nTupleMembers = (int)PySequence_Size(pResult);
  if(nTupleMembers!=expectedNTupleMembers){
    Py_DECREF(pResult);
    snprintf(message, STR_LEN_1, "The %s function returned a tuple with %d members (3 required).", funcName, nTupleMembers);
pyerror(message);
  }

  *pCurrentModel = PySequence_GetItem(pResult, (Py_ssize_t)0);
  if(*pCurrentModel==NULL){
    Py_DECREF(pResult);
    snprintf(message, STR_LEN_1, "The %s function returned a NULL object in the first tuple member.", funcName);
pyerror(message);
  }

  if(!PyInstance_Check(*pCurrentModel) && *pCurrentModel!=Py_None){
    /*
***NOTE*** that PyInstance_Check() only works with old-style classes. To check a new-style instance, it seems one would have to pass in the name of the module, read the class from that using PyObject_GetAttrString(), then compare the instance and the class via PyObject_IsInstance(<instance>, <class>).
    */
    Py_DECREF(*pCurrentModel);
    Py_DECREF(pResult);
    snprintf(message, STR_LEN_1, "The %s function did not return a valid modellib_classes._Model instance.", funcName);
pyerror(message);
  }


  pUserModuleName = PySequence_GetItem(pResult, (Py_ssize_t)1);
  if(pUserModuleName==NULL){
    Py_DECREF(pResult);
    snprintf(message, STR_LEN_1, "The %s function returned a NULL object in the second tuple member.", funcName);
pyerror(message);
  }

  if(pUserModuleName==Py_None){
    userModuleNameNoSuffix[0] = '\0';
  }else{
    tempStr = PyString_AsString(pUserModuleName);
    if(tempStr==NULL){
      Py_DECREF(pUserModuleName);
      Py_DECREF(pResult);
      snprintf(message, STR_LEN_1, "The %s function returned a NULL string in the second tuple member.", funcName);
pyerror(message);
    }
    myStrNCpy(userModuleNameNoSuffix, tempStr, maxLenName);
  }
  Py_DECREF(pUserModuleName);


  pCopyTemperatureI = PySequence_GetItem(pResult, (Py_ssize_t)2);
  if(pCopyTemperatureI==NULL){
    Py_DECREF(pResult);
    snprintf(message, STR_LEN_1, "The %s function returned a NULL object in the third tuple member.", funcName);
pyerror(message);
  }

  if(!PyInt_CheckExact(pCopyTemperatureI)){
    Py_DECREF(pCopyTemperatureI);
    Py_DECREF(pResult);
    snprintf(message, STR_LEN_1, "The %s function returned a non-integer in the third tuple member.", funcName);
pyerror(message);
  }

  *copyTemperatureI = (int)PyInt_AS_LONG(pCopyTemperatureI);

  Py_DECREF(pCopyTemperatureI);


  *pLimePars = PySequence_GetItem(pResult, (Py_ssize_t)3);
  if(*pLimePars==NULL){
    Py_DECREF(pResult);
    snprintf(message, STR_LEN_1, "The %s function returned a NULL object in the fourth tuple member.", funcName);
pyerror(message);
  }

  if(!PyInstance_Check(*pLimePars)){
    /*
***NOTE*** that PyInstance_Check() only works with old-style classes. To check a new-style instance, it seems one would have to pass in the name of the module, read the class from that using PyObject_GetAttrString(), then compare the instance and the class via PyObject_IsInstance(<instance>, <class>).
    */
    Py_DECREF(*pLimePars);
    Py_DECREF(pResult);
    snprintf(message, STR_LEN_1, "The %s function did not return a valid limepar_classes.ModelParameters instance.", funcName);
pyerror(message);
  }

  Py_DECREF(pResult);
}

/*....................................................................*/
int
main(int argc, char *argv[]){
  errType err=init_local_err();
  PyObject *pCurrentModel,*pLimePars,*pModule;
  const int maxLenName=100;
  char userModuleNameNoSuffix[maxLenName+1];
  int modelI=-1,nPars,nImgPars,nImages,status=0;
  char message[STR_LEN_1];
  const char *headerModuleName="limepar_classes";
  inputPars par;
  image *img = NULL;
  parTemplateType *parTemplates=NULL,*imgParTemplates=NULL;

  if (argc < 2){
    printf("Usage: casalime <name of file with pickled pars object>\n");
exit(1);
  }

  /* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .*/
  Py_Initialize();

  /* pCurrentModel should be a modellib_classes._Model instance, pLimePars should be a limepar_classes.ModelParameters instance. */
  _readParsObjWrapper(argv[1], &pCurrentModel, userModuleNameNoSuffix, maxLenName, &copyTemp, &pLimePars);

  /* Do some initialization */
  setDefaultFuncStuffs(); /* in ml_funcs */
  silent = 0;//********** pass it as argument?
  defaultFuncFlags = 0;

  if(pCurrentModel==Py_None){
    currentModelI = MODEL_None;
  }else{
    err = getModelI(pCurrentModel, &modelI); /* in ml_aux.c */
    if(err.status!=0){
      Py_DECREF(pCurrentModel);
      Py_DECREF(pLimePars);
pyerror(err.message);
    }

    currentModelI = modelI; /* global var. */

    /* Set some global arrays defined in the header of ml_models.c */
    err = extractParams(pCurrentModel); /* in ml_aux.c */
    if(err.status!=0){
      Py_DECREF(pCurrentModel);
      Py_DECREF(pLimePars);
pyerror(err.message);
    }

    /* Set some global arrays defined in the header of ml_funcs.c */
    err = extractFuncs(pCurrentModel); /* in ml_aux.c */
    if(err.status!=0){
      Py_DECREF(pCurrentModel);
      Py_DECREF(pLimePars);
pyerror(err.message);
    }
  }

  Py_DECREF(pCurrentModel);

  status = finalizeModelConfig(currentModelI); /* in ml_models.c */
  if(status!=0){
    Py_DECREF(pLimePars);
    sprintf(message, "finalizeModelConfig() returned status value %d", status);
pyerror(message);
  }

  /* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .*/
  /* Construct the 'macro' list argument:
  */
  err = setMacros(); /* in py_utils.c */
  if(err.status){
    Py_DECREF(pLimePars);
    unsetMacros(); /* in py_utils.c */
pyerror(err.message);
  }

  /* Set up any user-supplied result functions: */
  if(strlen(userModuleNameNoSuffix)>0){
    err = getModuleFromName(userModuleNameNoSuffix, &pModule); /* in py_utils.c */
    if(err.status!=0){ /* Don't need to decref pModule. */
      Py_DECREF(pLimePars);
pyerror(err.message);
    }

    /* Sets up global objects defined in the header of py_utils.c */
    setUpUserPythonFuncs(pModule); /* in py_utils.c */

    Py_DECREF(pModule);
  }

  /* Now get the lists of attribute names from the 2 classes in limepar_classes.py:
  */
  err = getParTemplatesWrapper(headerModuleName, &parTemplates, &nPars\
    , &imgParTemplates, &nImgPars); /* in py_utils.c */
  if(err.status!=0){
    Py_DECREF(pLimePars);
pyerror(err.message);
  }

  err = mallocInputParStrs(&par);
  if(err.status){
    unsetMacros();
    Py_DECREF(pLimePars);
pyerror(err.message);
  }

  /* Finally, unpack the LIME parameter values, following the templates: */
  err = readParImg(pLimePars, parTemplates, nPars, imgParTemplates\
    , nImgPars, &par, &img, &nImages, pywarning); /* in py_utils.c */

  if(err.status){
    Py_DECREF(pLimePars);
pyerror(err.message);
  }

  Py_DECREF(pLimePars);
  free(imgParTemplates);
  free(parTemplates);

  /* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .*/
  /* Now call the main bit of LIME:
  */
  status = run(par, img, nImages);

  /* Python-object clean up before status check and possible exit.
  */
  decrefAllUserFuncs(); /* in py_utils.c */
  free(modelDblPars); /* global in header of ml_models.c */
  free(modelIntPars); /* global in header of ml_models.c */
//************* why not the str pars??
  freeFuncsPars(); /* in ml_funcs.c */

  if(status){
    sprintf(message, "Function run() returned with status %d", status);
pyerror(message);
  }

  Py_Finalize();

  return 0;
}

