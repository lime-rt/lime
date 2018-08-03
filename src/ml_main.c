/*
 *  ml_main.c
 *  This file is part of LIME, the versatile line modeling engine.
 *
 *  See ../COPYRIGHT
 *
 */

#include "constants.h" /* for TRUE and FALSE */
#include "defaults.h"
#include "local_err.h"
#include "ml_types.h"
#include "ml_funcs.h"
#include "ml_models.h"
#include "py_utils.h"
#include "ufunc_types.h" /* for the 'result' function definitions. */


int copyTemp = 0; /* <0 means tgas should be set to equal tdust, >0 means the other way, ==0 means no action. */
_Bool isInitialized=FALSE,modelIsFinalized=FALSE;
double defaultDensyPower=DENSITY_POWER;

int silent = 0;
int defaultFuncFlags = 0;

int _verbosity=0;

/*....................................................................*/
static PyObject* py_initialize(PyObject* self, PyObject* args){
  errType err=init_local_err();

  if(!isInitialized){
    setDefaultFuncStuffs(); /* In ml_funcs.c */ 
    err = setMacros(); /* In py_utils.c */
    if(err.status){
      PyErr_SetString(PyExc_AttributeError, err.message);
return NULL;
    }

    silent = 0;//********** pass it as argument?
    defaultFuncFlags = 0;

    isInitialized = TRUE;
  }

Py_RETURN_NONE;
}

/*....................................................................*/
static PyObject* py_clean_up(PyObject* self, PyObject* args){

  if(isInitialized){
    unsetMacros(); /* In py_utils.c */
    decrefAllUserFuncs(); /* In py_utils.c */
    free(modelDblPars);
    free(modelIntPars);
    freeFuncsPars(); /* In ml_funcs.c */ 
    isInitialized = FALSE;
  }

Py_RETURN_NONE;
}

/*....................................................................*/
static PyObject* py_finalize_config(PyObject* self, PyObject* args){
  /* The calling routine is expected to pass in 'args' a tuple whose single element is model object which conforms to the pattern of class _Model in ../python/modellib.py.
  */
  errType err=init_local_err();
  char *functionName="finalize_config",message[STR_LEN_0+1],*userModuleNameNoSuffix;
  int modelI=-1,status=0;
  PyObject *pModelObj,*pModule;

if(_verbosity>0) printf(">>> Entering py_finalize_config()\n");

  if (!PyArg_ParseTuple(args, "(Osi)", &pModelObj, &userModuleNameNoSuffix, &copyTemp)) /* According to the docs, I don't have to decref pModelObj (I think). */
return NULL;

if(_verbosity>0) printf("ooo In py_finalize_config(). Arg tuple parsed\n");
if(_verbosity>0) printf("ooo In py_finalize_config(). User model is %s\n", userModuleNameNoSuffix);

  if(pModelObj==Py_None){ /* Means all the necessary 'result' functions should be supplied via userModuleNameNoSuffix. */
    currentModelI = MODEL_None;
  }else{
    /*
***NOTE*** that PyInstance_Check() only works with old-style classes. To check a new-style instance, it seems one would have to pass in the name of the module, read the class from that using PyObject_GetAttrString(), then compare the instance and the class via PyObject_IsInstance(<instance>, <class>).
    */
    if(!PyInstance_Check(pModelObj)){
      sprintf(message, "Object passed to %s is not an instance.", functionName);
      PyErr_SetString(PyExc_ValueError, message);
return NULL;
    }

    err = getModelI(pModelObj, &modelI); /* in ml_aux.c */
    if(err.status!=0){
      PyErr_SetString(PyExc_ValueError, err.message);
return NULL;
    }

    currentModelI = modelI;
if(_verbosity>0) printf("ooo In py_finalize_config(). Model I is %d\n", modelI);

    /* Clean up if these things were malloc'd:
    */
    free(modelDblPars);
    freeFuncsPars();

    err = extractParams(pModelObj); /* in ml_aux.c */
    if(err.status!=0){
      PyErr_SetString(PyExc_AttributeError, err.message);
return NULL;
    }

if(_verbosity>0) printf("ooo In py_finalize_config(). Extracted model parameters.\n");

    err = extractFuncs(pModelObj); /* in ml_aux.c */
    if(err.status!=0){
      PyErr_SetString(PyExc_AttributeError, err.message);
return NULL;
    }
if(_verbosity>0) printf("ooo In py_finalize_config(). Extracted model functions.\n");
  }
if(_verbosity>0) printf("ooo In py_finalize_config(). currentModelI = %d\n", currentModelI);

  status = finalizeModelConfig(currentModelI);
  if(status){
    sprintf(message, "finalizeModelConfig() returned status value %d", status);
    PyErr_SetString(PyExc_ValueError, message);
return NULL;
  }
if(_verbosity>0) printf("ooo In py_finalize_config(). Finalized model config.\n");

  /* Set up any user-supplied result functions: */
  if(strlen(userModuleNameNoSuffix)>0){
    err = getModuleFromName(userModuleNameNoSuffix, &pModule); /* in py_utils.c */
    if(err.status){
      PyErr_SetString(PyExc_ValueError, err.message);
return NULL;
    }

    /* Sets up global objects defined in the header of py_utils.c */
    setUpUserPythonFuncs(pModule); /* in py_utils.c */

    Py_DECREF(pModule);

  }else if(pModelObj==Py_None){
    sprintf(message, "No library model chosen, nor module with user functions provided.");
    PyErr_SetString(PyExc_ValueError, message);
return NULL;
  }

if(_verbosity>0) printf("ooo In py_finalize_config(). Set up any python functions.\n");

  modelIsFinalized = TRUE;
if(_verbosity>0) printf("<<< Leaving py_finalize_config()\n");

Py_RETURN_NONE;
}

/*....................................................................*/
static PyObject* py_is_config_finalized(PyObject* self, PyObject* args){
  if(modelIsFinalized)
    return Py_BuildValue("O", Py_True);
  else
    return Py_BuildValue("O", Py_False);
}

/*....................................................................*/
static PyObject* py_density(PyObject* self, PyObject* args){
  double x,y,z,values[0];

  if (!PyArg_ParseTuple(args, "ddd", &x, &y, &z))
return NULL;

  density(x, y, z, values);

return Py_BuildValue("d", values[0]);
}

/*....................................................................*/
static PyObject* py_doppler(PyObject* self, PyObject* args){
  double x,y,z,value=11.9;

if(_verbosity>0) printf("Entering py_doppler()\n");

  if (!PyArg_ParseTuple(args, "ddd", &x, &y, &z))
return NULL;

if(_verbosity>0) printf("About to call doppler()\n");

  doppler(x, y, z, &value);

if(_verbosity>0) printf("Leaving py_doppler()\n");
return Py_BuildValue("d", value);
}


/*....................................................................*/
/*
 * Bind Python function names to our C functions
 */
static PyMethodDef myModule_methods[] = {
  {"initialize",          py_initialize,          METH_VARARGS, "Not intended for direct user import"},
  {"clean_up",            py_clean_up,            METH_VARARGS, "Not intended for direct user import"},
  {"finalize_config",     py_finalize_config,     METH_VARARGS, "Not intended for direct user import"},
  {"is_config_finalized", py_is_config_finalized, METH_VARARGS, "Not intended for direct user import"},
  {"density",             py_density,             METH_VARARGS, "Not intended for direct user import"},
  {"doppler",             py_doppler,             METH_VARARGS, "Not intended for direct user import"},
  {NULL, NULL, 0, NULL}
};

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
initlibmodellib(void)
{
  (void) Py_InitModule("libmodellib", myModule_methods);
}

