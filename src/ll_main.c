/*
 *  ll_main.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2017 The LIME development team
 *
TODO:
 */

#include "lime.h"
#include "pyshared_io.h"
#include "py_utils.h"


int silent = 0;
statusType statusObj;
_Bool statusObjInitialized=FALSE;

int _cl_verbosity=0;

/*....................................................................*/
void initializeStatusObj(void){
  statusObj.progressGridBuilding      = 0.0;
  statusObj.progressGridSmoothing     = 0.0;
  statusObj.statusGrid                = 0;
  statusObj.progressConvergence       = 0.0;
  statusObj.numberIterations          = 0;
  statusObj.minsnr                    = 0.0;
  statusObj.median                    = 0.0;
  statusObj.progressPhotonPropagation = 0.0;
  statusObj.progressRayTracing        = 0.0;
  statusObj.statusRayTracing          = 0;
  statusObj.statusGlobal              = 0;
  statusObj.error                     = 0;
  statusObj.message[0]                = '\0';
}

/*....................................................................*/
static PyObject* py_set_silent(PyObject* self, PyObject* args){
  if (!PyArg_ParseTuple(args, "i", &silent))
    return NULL;
  return Py_BuildValue("i", 0);
}

/*....................................................................*/
static PyObject* py_get_silent(PyObject* self, PyObject* args){
  /* args should be empty */
  return Py_BuildValue("i", silent);
}

/*....................................................................*/
static PyObject* py_run_wrapper(PyObject* self, PyObject* args){
  /*
This is the function which is invoked to actually run LIME.

Note that 'args' should be an instance of type limepar_classes.py:ModelParameters.
  */

  int status=0,nPars,nImgPars=0,nImages=0;
  parTemplateType *parTemplates=NULL,*imgParTemplates=NULL;
  PyObject *pParClass,*pImgList,*pImgPars;
  inputPars inpars;
  image *inimg = NULL;
  char message[STR_LEN_0];

if(_cl_verbosity>0) printf("Entering py_run_wrapper()\n");

  if(!PyArg_ParseTuple(args, "O", &pParClass))
return NULL;

if(_cl_verbosity>0) printf("In py_run_wrapper(): extracted par class\n");

  status = getParTemplates(pParClass, &parTemplates, &nPars);
  if(status){
    snprintf(message, STR_LEN_0, "getParTemplates() returned status value %d", status);
if(_cl_verbosity>0) printf("getParTemplates() returned status value %d", status);
    PyErr_SetString(PyExc_ValueError, message);
    free(parTemplates);
return NULL;
  }

if(_cl_verbosity>0) printf("In py_run_wrapper(): got parameter templates\n");

  /* Get the number of images and check that the image attributes, for all images, have the correct types.
  */
  pImgList = PyObject_GetAttrString(pParClass, "img");
  if(pImgList==NULL){
    free(parTemplates);
return NULL;
  }

if(_cl_verbosity>0) printf("In py_run_wrapper(): got image list\n");

  nImages = (int)PyList_Size(pImgList);

if(_cl_verbosity>0) printf("In py_run_wrapper(): got number of images\n");

  if(nImages>0){
    pImgPars = PyList_GetItem(pImgList, (Py_ssize_t)0); /* Don't have to DECREF pImgPars, it is a borrowed reference. */

    status = getParTemplates(pImgPars, &imgParTemplates, &nImgPars);
    if(status){
      free(imgParTemplates);
      free(parTemplates);
      Py_DECREF(pImgList);
return NULL;
    }
  }

  Py_DECREF(pImgList);

if(_cl_verbosity>0) printf("In py_run_wrapper(): got image parameter templates\n");

  status = mallocInputParStrs(&inpars); /* Note that the inimg equivalents are malloc'd in readParImg(). This separation seems a little messy. */
  if(status){
    PyErr_NoMemory();
return NULL;
  }

if(_cl_verbosity>0) printf("In py_run_wrapper(): malloc'd input pars\n");

  status = readParImg(pParClass, parTemplates, nPars, imgParTemplates, nImgPars, &inpars, &inimg, &nImages, warning);

  free(imgParTemplates);
  free(parTemplates);

  if(status){
    snprintf(message, STR_LEN_0, "Function readParImg() returned with status %d.", status);
    PyErr_SetString(PyExc_AttributeError, message);
return NULL;
  }

if(_cl_verbosity>0) printf("In py_run_wrapper(): setup complete - ready to go\n");

  /* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .*/
  /* Now call the main bit of LIME:
  */
  status = run(inpars, inimg, nImages);
  if(status){
    snprintf(message, STR_LEN_0, "Function run() returned with status %d.", status);
    PyErr_SetString(PyExc_AttributeError, message);
return NULL;
  }

  pyFreeInputImgPars(inimg, nImages);
  freeInputPars(&inpars);

if(_cl_verbosity>0) printf("Leaving py_run_wrapper()\n");

Py_RETURN_NONE;
}

/*....................................................................*/
static PyObject* py_read_status(PyObject* self, PyObject* args){
  /* args should be empty */
  PyObject *pStatusTuple;

  if(!statusObjInitialized){
    initializeStatusObj();
    statusObjInitialized=TRUE;
  }

  pStatusTuple = Py_BuildValue("ffifiddffiiis"\
    , statusObj.progressGridBuilding\
    , statusObj.progressGridSmoothing\
    , statusObj.statusGrid\
    , statusObj.progressConvergence\
    , statusObj.numberIterations\
    , statusObj.minsnr\
    , statusObj.median\
    , statusObj.progressPhotonPropagation\
    , statusObj.progressRayTracing\
    , statusObj.statusRayTracing\
    , statusObj.statusGlobal\
    , statusObj.error\
    , statusObj.message\
  );

  return pStatusTuple;
}

/*....................................................................*/
/*
 * Bind Python function names to our C functions
 */
static PyMethodDef myModule_methods[] = {
  {"set_silent",   py_set_silent,  METH_VARARGS, "Not intended for direct user import"},
  {"get_silent",   py_get_silent,  METH_VARARGS, "Not intended for direct user import"},
  {"read_status",  py_read_status, METH_VARARGS, "Not intended for direct user import"},
  {"run_wrapper",  py_run_wrapper, METH_VARARGS, "Not intended for direct user import"},
  {NULL, NULL, 0, NULL}
};

/*
 * Python calls this to let us initialize our module
 */
#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
initliblime(void)
{
  (void) Py_InitModule("liblime", myModule_methods);
}

/*....................................................................*/
int
main(int argc, char *argv[]){
  // Test main function for testing liblime as an executable.

  printf("Testing liblime.\n");

  int testI=0;

  switch(testI){
    case 0: /* Test a few of the 'result' functions. */
//****
  break;
    default:
      printf("Unrecognized test code %d\n", testI);
  }

return 0;
}







