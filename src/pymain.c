/*
 *  pymain.c
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

#ifdef FASTEXP
double EXP_TABLE_2D[128][10];
double EXP_TABLE_3D[256][2][10];
/* I've hard-wired the dimensions of these arrays, but it would be better perhaps to declare them as pointers, and calculate the dimensions with the help of the function call:
  calcFastExpRange(FAST_EXP_MAX_TAYLOR, FAST_EXP_NUM_BITS, &numMantissaFields, &lowestExponent, &numExponentsUsed)
*/
#else
double EXP_TABLE_2D[1][1]; // nominal definitions so the fastexp.c module will compile.
double EXP_TABLE_3D[1][1][1];
#endif

PyObject *pModule_global = NULL,\
         *pMacros_global = NULL,\
         *pDensity     = NULL,\
         *pTemperature = NULL,\
         *pAbundance   = NULL,\
         *pDoppler     = NULL,\
         *pVelocity    = NULL,\
         *pMagfield    = NULL,\
         *pGasIIdust   = NULL,\
         *pGridDensity = NULL;

/*....................................................................*/
void
decrefAllGlobals(){
  Py_XDECREF(pDensity);
  Py_XDECREF(pTemperature);
  Py_XDECREF(pAbundance);
  Py_XDECREF(pDoppler);
  Py_XDECREF(pVelocity);
  Py_XDECREF(pMagfield);
  Py_XDECREF(pGasIIdust);
  Py_XDECREF(pGridDensity);
  Py_DECREF(pMacros_global);
  Py_DECREF(pModule_global);
}

/*....................................................................*/
void
userFuncWrapper(PyObject *pFunc, const char *funcName, PyObject *pMacros\
  , double x, double y, double z, PyObject **pResult){

  int status=0;
  PyObject *pArgs;

  pArgs = Py_BuildValue("(Offf)", pMacros, x, y, z);
  if(pArgs==NULL){
    status = 1;
  }else{
    *pResult = PyObject_CallObject(pFunc, pArgs);
    Py_DECREF(pArgs);
    if(*pResult==NULL){
      status = 2;
    }
  }

  if(status){
    char message[STR_LEN_0];
    if(!silent)
      PyErr_Print();
    decrefAllGlobals();
    sprintf(message, "User function %s() failed with status %d.", funcName, status);
    pyerror(message);
  }
}

/*....................................................................*/
void
density(double x, double y, double z, double *densValues){
  if(pDensity==NULL) /* User did not supply this function. */
    density_default(x, y, z, densValues);
  else{
    int nItems=0,i;
    PyObject *pResult=NULL,*pListItem;
    userFuncWrapper(pDensity, "density", pMacros_global, x, y, z, &pResult); /* Sets up and calls the function. pResult guaranteed non-NULL. */

    /* The returned value should be a list of floats (or ints).
    */
    if(!PyList_Check(pResult)){
      if(PyErr_Occurred())
        printOrClearPyError();
      Py_DECREF(pResult);
      decrefAllGlobals();
      pyerror("User function density() should return a list.");
    }

    nItems = (int)PyList_Size(pResult);
    if(PyErr_Occurred() || nItems<=0){
      if(PyErr_Occurred())
        printOrClearPyError();
      Py_DECREF(pResult);
      decrefAllGlobals();
      pyerror("User function density() returned an empty or defective list.");
    }

    for(i=0;i<nItems;i++){
      pListItem = PyList_GetItem(pResult, (Py_ssize_t)i); /* Borrowed reference, don't DECREF it. */
      /* Not going to check for errors, since I think I have covered all the possibilities already...? */

      if(PyFloat_CheckExact(pListItem))
        densValues[i] = PyFloat_AS_DOUBLE(pListItem);
      else if(PyInt_CheckExact(pListItem))
        densValues[i] = (double)PyInt_AS_LONG(pListItem);
      else{
        char message[STR_LEN_0];
        Py_DECREF(pResult);
        decrefAllGlobals();
        sprintf(message, "Item %d in return from user function density() is not numeric.", i);
        pyerror(message);
      }
    }
    Py_DECREF(pResult);
  }
}

/*....................................................................*/
void
temperature(double x, double y, double z, double *tValues){
  if(pTemperature==NULL) /* User did not supply this function. */
    temperature_default(x, y, z, tValues);
  else{
    int nItems=0,i;
    PyObject *pResult=NULL,*pTupleItem;
    userFuncWrapper(pTemperature, "temperature", pMacros_global, x, y, z, &pResult); /* Sets up and calls the function. pResult guaranteed non-NULL. */

    /* The returned value should be a tuple of 2 floats or ints (the second of these may also ==None however).
    */
    if(!PyTuple_Check(pResult)){
      if(PyErr_Occurred())
        printOrClearPyError();
      Py_DECREF(pResult);
      decrefAllGlobals();
      pyerror("User function temperature() should return a tuple.");
    }

    nItems = (int)PyTuple_Size(pResult);
    if(PyErr_Occurred() || nItems!=2){
      if(PyErr_Occurred())
        printOrClearPyError();
      Py_DECREF(pResult);
      decrefAllGlobals();
      pyerror("User function temperature() returned a defective tuple.");
    }

    for(i=0;i<nItems;i++){
      pTupleItem = PyTuple_GetItem(pResult, (Py_ssize_t)i); /* Borrowed reference, don't DECREF it. */
      /* Not going to check for errors, since I think I have covered all the possibilities already...? */

      if(PyFloat_CheckExact(pTupleItem))
        tValues[i] = PyFloat_AS_DOUBLE(pTupleItem);
      else if(PyInt_CheckExact(pTupleItem))
        tValues[i] = (double)PyInt_AS_LONG(pTupleItem);
      else if(i==1 && pTupleItem==Py_None)
        tValues[i] = -1.0; /* This is used in Lime as the default value for t[2]. ~:-/ */
      else{
        char message[STR_LEN_0];
        Py_DECREF(pResult);
        decrefAllGlobals();
        sprintf(message, "Item %d in return from user function temperature() is not numeric.", i);
        pyerror(message);
      }
    }
    Py_DECREF(pResult);
  }
}

/*....................................................................*/
void
abundance(double x, double y, double z, double *abunValues){
  if(pAbundance==NULL) /* User did not supply this function. */
    abundance_default(x, y, z, abunValues);
  else{
    int nItems=0,i;
    PyObject *pResult=NULL,*pListItem;
    userFuncWrapper(pAbundance, "abundance", pMacros_global, x, y, z, &pResult); /* Sets up and calls the function. pResult guaranteed non-NULL. */

    /* The returned value should be a list of floats (or ints).
    */
    if(!PyList_Check(pResult)){
      if(PyErr_Occurred())
        printOrClearPyError();
      Py_DECREF(pResult);
      decrefAllGlobals();
      pyerror("User function abundance() should return a list.");
    }

    nItems = (int)PyList_Size(pResult);
    if(PyErr_Occurred() || nItems<=0){
      if(PyErr_Occurred())
        printOrClearPyError();
      Py_DECREF(pResult);
      decrefAllGlobals();
      pyerror("User function abundance() returned an empty or defective list.");
    }

    for(i=0;i<nItems;i++){
      pListItem = PyList_GetItem(pResult, (Py_ssize_t)i); /* Borrowed reference, don't DECREF it. */
      /* Not going to check for errors, since I think I have covered all the possibilities already...? */

      if(PyFloat_CheckExact(pListItem))
        abunValues[i] = PyFloat_AS_DOUBLE(pListItem);
      else if(PyInt_CheckExact(pListItem))
        abunValues[i] = (double)PyInt_AS_LONG(pListItem);
      else{
        char message[STR_LEN_0];
        Py_DECREF(pResult);
        decrefAllGlobals();
        sprintf(message, "Item %d in return from user function abundance() is not numeric.", i);
        pyerror(message);
      }
    }
    Py_DECREF(pResult);
  }
}

/*....................................................................*/
void
doppler(double x, double y, double z, double *doppValue){
  if(pDoppler==NULL) /* User did not supply this function. */
    doppler_default(x, y, z, doppValue);
  else{
    PyObject *pResult=NULL;
    userFuncWrapper(pDoppler, "doppler", pMacros_global, x, y, z, &pResult); /* Sets up and calls the function. pResult guaranteed non-NULL. */

    /* The returned value should be a float or an int.
    */
    if(PyFloat_CheckExact(pResult))
      *doppValue = PyFloat_AS_DOUBLE(pResult);
    else if(PyInt_CheckExact(pResult))
      *doppValue = (double)PyInt_AS_LONG(pResult);
    else{
      char message[STR_LEN_0];
      Py_DECREF(pResult);
      decrefAllGlobals();
      sprintf(message, "Return from user function doppler() is not numeric.");
      pyerror(message);
    }

    Py_DECREF(pResult);
  }
}

/*....................................................................*/
void
velocity(double x, double y, double z, double *veloValues){
  if(pVelocity==NULL) /* User did not supply this function. */
    velocity_default(x, y, z, veloValues);
  else{
    int nItems=0,i;
    PyObject *pResult=NULL,*pListItem;
    userFuncWrapper(pVelocity, "velocity", pMacros_global, x, y, z, &pResult); /* Sets up and calls the function. pResult guaranteed non-NULL. */

    /* The returned value should be a list of 3 floats (or ints).
    */
    if(!PyList_Check(pResult)){
      if(PyErr_Occurred())
        printOrClearPyError();
      Py_DECREF(pResult);
      decrefAllGlobals();
      pyerror("User function velocity() should return a list.");
    }

    nItems = (int)PyList_Size(pResult);
    if(PyErr_Occurred() || nItems!=3){
      if(PyErr_Occurred())
        printOrClearPyError();
      Py_DECREF(pResult);
      decrefAllGlobals();
      pyerror("User function velocity() returned a defective list.");
    }

    for(i=0;i<nItems;i++){
      pListItem = PyList_GetItem(pResult, (Py_ssize_t)i); /* Borrowed reference, don't DECREF it. */
      /* Not going to check for errors, since I think I have covered all the possibilities already...? */

      if(PyFloat_CheckExact(pListItem))
        veloValues[i] = PyFloat_AS_DOUBLE(pListItem);
      else if(PyInt_CheckExact(pListItem))
        veloValues[i] = (double)PyInt_AS_LONG(pListItem);
      else{
        char message[STR_LEN_0];
        Py_DECREF(pResult);
        decrefAllGlobals();
        sprintf(message, "Item %d in return from user function velocity() is not numeric.", i);
        pyerror(message);
      }
    }
    Py_DECREF(pResult);
  }
}

/*....................................................................*/
void
magfield(double x, double y, double z, double *magfValues){
  if(pMagfield==NULL) /* User did not supply this function. */
    magfield_default(x, y, z, magfValues);
  else{
    int nItems=0,i;
    PyObject *pResult=NULL,*pListItem;
    userFuncWrapper(pMagfield, "magfield", pMacros_global, x, y, z, &pResult); /* Sets up and calls the function. pResult guaranteed non-NULL. */

    /* The returned value should be a list of 3 floats (or ints).
    */
    if(!PyList_Check(pResult)){
      if(PyErr_Occurred())
        printOrClearPyError();
      Py_DECREF(pResult);
      decrefAllGlobals();
      pyerror("User function magfield() should return a list.");
    }

    nItems = (int)PyList_Size(pResult);
    if(PyErr_Occurred() || nItems!=3){
      if(PyErr_Occurred())
        printOrClearPyError();
      Py_DECREF(pResult);
      decrefAllGlobals();
      pyerror("User function magfield() returned a defective list.");
    }

    for(i=0;i<nItems;i++){
      pListItem = PyList_GetItem(pResult, (Py_ssize_t)i); /* Borrowed reference, don't DECREF it. */
      /* Not going to check for errors, since I think I have covered all the possibilities already...? */

      if(PyFloat_CheckExact(pListItem))
        magfValues[i] = PyFloat_AS_DOUBLE(pListItem);
      else if(PyInt_CheckExact(pListItem))
        magfValues[i] = (double)PyInt_AS_LONG(pListItem);
      else{
        char message[STR_LEN_0];
        Py_DECREF(pResult);
        decrefAllGlobals();
        sprintf(message, "Item %d in return from user function magfield() is not numeric.", i);
        pyerror(message);
      }
    }
    Py_DECREF(pResult);
  }
}

/*....................................................................*/
void
gasIIdust(double x, double y, double z, double *gas2dustValue){
  if(pGasIIdust==NULL) /* User did not supply this function. */
    gasIIdust_default(x, y, z, gas2dustValue);
  else{
    PyObject *pResult=NULL;
    userFuncWrapper(pGasIIdust, "gasIIdust", pMacros_global, x, y, z, &pResult); /* Sets up and calls the function. pResult guaranteed non-NULL. */

    /* The returned value should be a float or an int.
    */
    if(PyFloat_CheckExact(pResult))
      *gas2dustValue = PyFloat_AS_DOUBLE(pResult);
    else if(PyInt_CheckExact(pResult))
      *gas2dustValue = (double)PyInt_AS_LONG(pResult);
    else{
      char message[STR_LEN_0];
      Py_DECREF(pResult);
      decrefAllGlobals();
      sprintf(message, "Return from user function gasIIdust() is not numeric.");
      pyerror(message);
    }

    Py_DECREF(pResult);
  }
}

/*....................................................................*/
double
gridDensity(configInfo *par, double *r){
  double fracDensity=0.0;

  if(pGridDensity==NULL) /* User did not supply this function. */
    fracDensity = gridDensity_default(par, r);
  else{
    /*  ***** NOTE ***** that we are throwing away the config info! */

    PyObject *pResult=NULL;
    userFuncWrapper(pGridDensity, "gridDensity", pMacros_global, r[0], r[1], r[2], &pResult); /* Sets up and calls the function. pResult guaranteed non-NULL. */

    /* The returned value should be a float or an int.
    */
    if(PyFloat_CheckExact(pResult))
      fracDensity = PyFloat_AS_DOUBLE(pResult);
    else if(PyInt_CheckExact(pResult))
      fracDensity = (double)PyInt_AS_LONG(pResult);
    else{
      char message[STR_LEN_0];
      Py_DECREF(pResult);
      decrefAllGlobals();
      sprintf(message, "Return from user function gridDensity() is not numeric.");
      pyerror(message);
    }

    Py_DECREF(pResult);
  }

  return fracDensity;
}

/*....................................................................*/
void
mallocInputParStrs(inputPars *par){
  int i;

  /* We have to malloc the strings in 'par' here (even though it means we will have to free them again afterward) because we want to strcpy() values into them rather than have them point to read-only strings.
  */
  par->outputfile    = malloc(sizeof(char)*(STR_LEN_0+1));
  par->binoutputfile = malloc(sizeof(char)*(STR_LEN_0+1));
  par->gridfile      = malloc(sizeof(char)*(STR_LEN_0+1));
  par->pregrid       = malloc(sizeof(char)*(STR_LEN_0+1));
  par->restart       = malloc(sizeof(char)*(STR_LEN_0+1));
  par->dust          = malloc(sizeof(char)*(STR_LEN_0+1));
  par->outputfile[0]    = '\0';
  par->binoutputfile[0] = '\0';
  par->gridfile[0]      = '\0';
  par->pregrid[0]       = '\0';
  par->restart[0]       = '\0';
  par->dust[0]          = '\0';

  for(i=0;i<MAX_NSPECIES;i++){
    par->moldatfile[i] = malloc(sizeof(char)*(STR_LEN_0+1));
    par->moldatfile[i][0] = '\0';
  }
}

/*....................................................................*/
int
main(int argc, char *argv[]){
  /*
Main program for stand-alone LIME with a python model file.

For pretty detailed documentation on embedding python in C, see

  https://docs.python.org/2/c-api/index.html
  */

  const int lenSuffix=3,maxLenNoSuffix=STR_LEN_0,nDblMacros=9,nIntMacros=7;
  const char *nameOfExecutable="lime", *headerModuleName="par_classes";
  const char *oldModulePath;
  char *modelName,modelNameNoSuffix[maxLenNoSuffix+1],message[STR_LEN_0],*newModulePath;
  char suffix[lenSuffix+1];
  int lenNoSuffix,nPars,nImgPars,nImages,i,status,strlenOMPath;
  PyObject *pName,*pValue;
  inputPars par;
  image *img = NULL;
  parTemplateType *parTemplate=NULL,*imgParTemplate=NULL;
  struct {char *name; double value;} macrosDbl[nDblMacros];
  struct {char *name;    int value;} macrosInt[nIntMacros];

  i = 0;
  macrosDbl[i++].name = "AMU";
  macrosDbl[i++].name = "CLIGHT";
  macrosDbl[i++].name = "HPLANCK";
  macrosDbl[i++].name = "KBOLTZ";
  macrosDbl[i++].name = "GRAV";
  macrosDbl[i++].name = "AU";
  macrosDbl[i++].name = "PC";
  macrosDbl[i++].name = "PI";
  macrosDbl[i++].name = "SPI";
  i = 0;
  macrosDbl[i++].value = AMU;
  macrosDbl[i++].value = CLIGHT;
  macrosDbl[i++].value = HPLANCK;
  macrosDbl[i++].value = KBOLTZ;
  macrosDbl[i++].value = GRAV;
  macrosDbl[i++].value = AU;
  macrosDbl[i++].value = PC;
  macrosDbl[i++].value = PI;
  macrosDbl[i++].value = SPI;

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

  silent = 0;

  mallocInputPars(&par);
  mallocInputParStrs(&par);

  /* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .*/
  /* Get the model name, then strip off the '.py' from it:
  */
  if(argc<2){
    if(!silent){
      sprintf(message, "Useage: %s <python model file>", nameOfExecutable);
      warning(message);
    }
    exit(1);
  }
  modelName = argv[1];

  lenNoSuffix = strlen(modelName) - lenSuffix;
  if(lenNoSuffix<1){
    sprintf(message, "Model file name must be more than %d characters long!", lenSuffix);
    error(message);
  }
  if(lenNoSuffix>maxLenNoSuffix){
    sprintf(message, "Model file name is longer than the permitted %d characters.\n", maxLenNoSuffix);
    error(message);
  }

  strncpy(suffix, modelName + lenNoSuffix, lenSuffix);
  suffix[lenSuffix] = '\0';

  if(strcmp(suffix,".py")!=0){
    sprintf(message, "Python files must end in '.py'");
    error(message);
  }

  strncpy(modelNameNoSuffix, modelName, strlen(modelName)-lenSuffix);
  modelNameNoSuffix[lenNoSuffix] = '\0';

  /* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .*/
  Py_Initialize();

  /* The first thing to do is add the PWD to sys.path, which doesn't happen by default when embedding.
  */
  oldModulePath = Py_GetPath();
  strlenOMPath = strlen(oldModulePath);
  newModulePath = malloc(sizeof(char)*strlenOMPath+3);
  if(myStrCpy(oldModulePath, newModulePath, strlenOMPath+2))
    pyerror("Could not copy existing sys.path to a new string.");

  newModulePath[strlenOMPath] = ':';
  newModulePath[strlenOMPath+1] = '.';
  newModulePath[strlenOMPath+2] = '\0';

  PySys_SetPath(newModulePath);
  free(newModulePath);
  if(PyErr_Occurred()){
    if(!silent)
      PyErr_Print();
    pyerror("Could not append PWD to sys.path");
  }

  /* Now get the lists of attribute names from the 2 classes in par_classes.py:
  */
  status = getParTemplates(headerModuleName, &parTemplate, &nPars\
    , &imgParTemplate, &nImgPars);
  if(status){
    sprintf(message, "Function getParTemplates() returned with status %d", status);
    pyerror(message);
  }

  /* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .*/
  /* Construct the 'macro' list argument:
  */
  pMacros_global = PyDict_New();
  if(pMacros_global==NULL){
    if(!silent)
      PyErr_Print();
    sprintf(message, "Failed to initialize macro dictionary.");
    pyerror(message);
  }

  for(i=0;i<nDblMacros;i++){
    pValue = PyFloat_FromDouble(macrosDbl[i].value);
    if(pValue==NULL){
      if(!silent)
        PyErr_Print();
      Py_DECREF(pMacros_global);
      sprintf(message, "Failed to convert type double macro %d", i);
      pyerror(message);
    }

    if(PyDict_SetItemString(pMacros_global, macrosDbl[i].name, pValue)){
      if(!silent)
        PyErr_Print();
      Py_DECREF(pValue);
      Py_DECREF(pMacros_global);
      sprintf(message, "Failed to set dictionary item for type double macro %d", i);
      pyerror(message);
    }

    Py_DECREF(pValue);
  }

  for(i=0;i<nIntMacros;i++){
    pValue = Py_BuildValue("i", macrosInt[i].value);

    if(pValue==NULL){
      if(!silent)
        PyErr_Print();
      Py_DECREF(pMacros_global);
      sprintf(message, "Failed to convert type int macro %d", i);
      pyerror(message);
    }

    if(PyDict_SetItemString(pMacros_global, macrosInt[i].name, pValue)){
      if(!silent)
        PyErr_Print();
      Py_DECREF(pValue);
      Py_DECREF(pMacros_global);
      sprintf(message, "Failed to set dictionary item for type int macro %d", i);
      pyerror(message);
    }

    Py_DECREF(pValue);
  }

  /* Now we open the user's 'model' module:
  */
  pName = PyString_FromString(modelNameNoSuffix);
  if(pName==NULL){
    if(!silent)
      PyErr_Print();
    Py_DECREF(pMacros_global);
    sprintf(message, "Could not convert module name to python string.");
    pyerror(message);
  }

  pModule_global = PyImport_Import(pName);
  Py_DECREF(pName);
  if(pModule_global==NULL){
    if(!silent)
      PyErr_Print();
    Py_DECREF(pMacros_global);
    sprintf(message, "Failed to load %s", modelName);
    pyerror(message);
  }

  /* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .*/
  /* Read user-supplied parameters from the 'model' module they supply:
  */
  status = initParImg(pModule_global, pMacros_global, parTemplate, nPars, imgParTemplate, nImgPars, &par, &img, &nImages);
  if(status){
    Py_DECREF(pMacros_global);
    Py_DECREF(pModule_global);
    sprintf(message, "Function initParImg() returned with status %d", status);
    pyerror(message);
  }

  free(imgParTemplate);
  free(parTemplate);

  /* Set up the 'user-supplied' functions (they are left at NULL if the user has not supplied them)
  */
  getPythonFunc(pModule_global, "density",     1,     &pDensity);
  getPythonFunc(pModule_global, "temperature", 0, &pTemperature);
  getPythonFunc(pModule_global, "abundance",   0,   &pAbundance);
  getPythonFunc(pModule_global, "doppler",     0,     &pDoppler);
  getPythonFunc(pModule_global, "velocity",    0,    &pVelocity);
  getPythonFunc(pModule_global, "magfield",    0,    &pMagfield);
  getPythonFunc(pModule_global, "gasIIdust",   0,   &pGasIIdust);
  getPythonFunc(pModule_global, "gridDensity", 0, &pGridDensity);

  /* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .*/
  /* Now call the main bit of LIME:
  */
  status = run(par, img, nImages);

  /* Python-object clean up before status check and possible exit.
  */
  decrefAllGlobals();

  if(status){
    sprintf(message, "Function run() returned with status %d", status);
    pyerror(message);
  }

  Py_Finalize();

//  free(par.outputfile);
//  free(par.binoutputfile);
//  free(par.gridfile);
//  free(par.pregrid);
//  free(par.restart);
//  free(par.dust);

  freeInputImg(nImages, img);
  freeInputPars(MAX_NSPECIES, &par);

  return 0;
}

