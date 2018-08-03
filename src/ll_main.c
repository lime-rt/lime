/*
 *  ll_main.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#include "lime.h"
#include "py_utils.h"
#include "local_err.h"


int silent = 0;
int _cl_verbosity=0;

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

  errType err=init_local_err();
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

  err = getParTemplates(pParClass, &parTemplates, &nPars);
  if(err.status){
    PyErr_SetString(PyExc_ValueError, err.message);
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

    err = getParTemplates(pImgPars, &imgParTemplates, &nImgPars);
    if(err.status){
      PyErr_SetString(PyExc_ValueError, err.message);
      free(imgParTemplates);
      free(parTemplates);
      Py_DECREF(pImgList);
return NULL;
    }
  }

  Py_DECREF(pImgList);

if(_cl_verbosity>0) printf("In py_run_wrapper(): got image parameter templates\n");

  err = mallocInputParStrs(&inpars); /* Note that the inimg equivalents are malloc'd in readParImg(). This separation seems a little messy. */
  if(err.status){
    PyErr_NoMemory();
return NULL;
  }

if(_cl_verbosity>0) printf("In py_run_wrapper(): malloc'd input pars\n");

  err = readParImg(pParClass, parTemplates, nPars, imgParTemplates, nImgPars, &inpars, &inimg, &nImages, warning);

  free(imgParTemplates);
  free(parTemplates);

  if(err.status){
    PyErr_SetString(PyExc_AttributeError, err.message);
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
/*
 * Bind Python function names to our C functions
 */
static PyMethodDef myModule_methods[] = {
  {"set_silent",   py_set_silent,  METH_VARARGS, "Not intended for direct user import"},
  {"get_silent",   py_get_silent,  METH_VARARGS, "Not intended for direct user import"},
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







