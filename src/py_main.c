/*
 *  py_main.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2017 The LIME development team
 *
TODO:
 */

#include "lime.h"
#include "error_codes.h"
#include "py_lime.h"
#include <argp.h>
#include "py_utils.h"

#ifdef NOVERBOSE
int silent = 1;
#else
int silent = 0;
#endif

int defaultFuncFlags = 0;

const char *argp_program_version = VERSION;
const char *argp_program_bug_address = "https://github.com/lime-rt/lime";
/* Program documentation. */
static char doc[] = "pylime - the version of LIME which accepts a model file written in python.";

/* A description of the arguments we accept. */
static char args_doc[] = "modelfilename";

/* The options we understand. */
static struct argp_option options[] = {
  {"silent",   's',     0, 0, "Suppress output messages." },
  {"testmode", 't',     0, 0, "Use fixed RNG seeds." },
  {"nthreads", 'p', "int", 0, "Run in parallel with NTHREADS threads (default: 1)" },
  { 0 }
};

/* Used by main to communicate with parse_opt. */
struct arguments
{
  _Bool doSilent,fixSeeds;
  char *modelfilename;
  int numThreads;
};

/*....................................................................*/
/* Parse a single option. */
static error_t
parse_opt(int key, char *arg, struct argp_state *state){
  /* Get the input argument from argp_parse, which we
     know is a pointer to our arguments structure. */
  struct arguments *arguments = state->input;

  switch (key)
    {
    case 's':
      arguments->doSilent = TRUE;
      break;
    case 't':
      arguments->fixSeeds = TRUE;
      break;
    case 'p':
      arguments->numThreads = atoi(arg);
      break;

    case ARGP_KEY_ARG:
      arguments->modelfilename = arg;
      break;

    case ARGP_KEY_END:
      break;

    default:
      return ARGP_ERR_UNKNOWN;
    }
  return 0;
}

/* Our argp parser. */
static struct argp argp = { options, parse_opt, args_doc, doc };

/*....................................................................*/
void
readParImgWrapper(PyObject *pModule, parTemplateType *parTemplates\
  , const int nPars, parTemplateType *imgParTemplates, const int nImgPars\
  , inputPars *par, image **img, int *nImages){

  int status=0;
  PyObject *pUserInputFunc=NULL,*pMacroArgs,*pPars;
  char message[STR_LEN_0];

  getPythonFunc(pModule, "input", &pUserInputFunc);
  if(pUserInputFunc==NULL){
    printOrClearPyError();
    sprintf(message, "getPythonFunc() didn't find the users input() function.");
pyerror(message);
  }

  pMacroArgs = Py_BuildValue("(O)", pMacros_global);
  if(pMacroArgs==NULL){
    unsetMacros();
    Py_DECREF(pUserInputFunc);
    printOrClearPyError();
    sprintf(message, "blah");
pyerror(message);
  }

  pPars = PyObject_CallObject(pUserInputFunc, pMacroArgs); /* Should return an instance of a limepar_classes.ModelParameters object. */
  Py_DECREF(pMacroArgs);
  Py_DECREF(pUserInputFunc);
  if(pPars==NULL){
    unsetMacros();
    printOrClearPyError();
    sprintf(message, "User's input() function returned a NULL ModelParameters object.");
pyerror(message);
  }
  if(!PyInstance_Check(pPars)){
    unsetMacros();
    Py_DECREF(pPars);
    sprintf(message, "User's input() function did not return a valid ModelParameters object.");
pyerror(message);
  }

  status = readParImg(pPars, parTemplates, nPars, imgParTemplates, nImgPars, par, img, nImages, pywarning);
  if(status){
    Py_DECREF(pPars);
    sprintf(message, "readParImg() returned status value %d", status);
pyerror(message);
  }

  Py_DECREF(pPars);
}

/*....................................................................*/
int
main(int argc, char *argv[]){
  /*
Main program for stand-alone LIME with a python model file.

For pretty detailed documentation on embedding python in C, see

  https://docs.python.org/2/c-api/index.html
  */

  const int lenSuffix=3,maxLenNoSuffix=STR_LEN_0;
  const char *nameOfExecutable="pylime",*headerModuleName="limepar_classes";
  const char *oldModulePath;
  char *modelName,modelNameNoSuffix[maxLenNoSuffix+1],message[STR_LEN_0],*newModulePath;
  char suffix[lenSuffix+1];
  int lenNoSuffix,nPars,nImgPars,nImages,status=0,strlenOMPath;
  PyObject *pModule = NULL;
  inputPars par;
  image *img = NULL;
  parTemplateType *parTemplates=NULL,*imgParTemplates=NULL;
  struct arguments arguments;

  /* Set defaults for argument returns:
  */
  arguments.doSilent      = FALSE;
  arguments.fixSeeds      = FALSE;
  arguments.numThreads    = -1;
  arguments.modelfilename = NULL;

  /* Parse our arguments; every option seen by parse_opt will be reflected in arguments.
  */
  argp_parse(&argp, argc, argv, 0, 0, &arguments);

  /* Interpret the options (numThreads is left until later):
  */
  if(arguments.doSilent)
    silent = 1;

  if(arguments.fixSeeds)
    fixRandomSeeds = TRUE;

  /* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .*/
  /* Get the model name, then strip off the '.py' from it:
  */
  if(arguments.modelfilename==NULL){
    sprintf(message, "Usage: %s [<arguments>] <python model file>", nameOfExecutable);
    error(message);
  }

  lenNoSuffix = strlen(arguments.modelfilename) - lenSuffix;
  if(lenNoSuffix<1){
    sprintf(message, "Model file name must be more than %d characters long!", lenSuffix);
    error(message);
  }
  if(lenNoSuffix>maxLenNoSuffix){
    sprintf(message, "Model file name is longer than the permitted %d characters.\n", maxLenNoSuffix);
    error(message);
  }

  strncpy(suffix, arguments.modelfilename + lenNoSuffix, lenSuffix);
  suffix[lenSuffix] = '\0';

  if(strcmp(suffix,".py")!=0){
    sprintf(message, "Python files must end in '.py'");
    error(message);
  }

  strncpy(modelNameNoSuffix, arguments.modelfilename, strlen(arguments.modelfilename)-lenSuffix);
  modelNameNoSuffix[lenNoSuffix] = '\0';

  /* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .*/
  Py_Initialize();

  /* The first thing to do is add the PWD to sys.path, which doesn't happen by default when embedding.
  */
  oldModulePath = Py_GetPath();
  strlenOMPath = strlen(oldModulePath);
  newModulePath = malloc(sizeof(char)*strlenOMPath+3);
  myStrNCpy(newModulePath, oldModulePath, strlenOMPath+2);

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

  /* Now get the lists of attribute names from the 2 classes in limepar_classes.py:
  */
  status = getParTemplatesWrapper(headerModuleName, &parTemplates, &nPars\
    , &imgParTemplates, &nImgPars, message); /* in py_utils.c */
  if(status!=0){
pyerror(message);
  }

  /* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .*/
  /* Construct the 'macro' list argument:
  */
  status = setMacros(); /* in py_utils.c */
  if(status){
    unsetMacros(); /* in py_utils.c */
    sprintf(message, "Function setMacros() returned with status %d", status);
pyerror(message);
  }

  /* Now we open the user's 'model' module:
  */
  status = getModuleFromName(modelNameNoSuffix, &pModule);
  if(status){
    if(status==PY_STRING_READ_FAIL){
      sprintf(message, "Could not convert module name to python string.");
    }else if(status==PY_IMPORT_FAIL){
      sprintf(message, "Failed to load %s", modelName);
    }else{
      sprintf(message, "getModuleFromName() returned with unknown status %d", status);
    }
    if(!silent)
      PyErr_Print();
    unsetMacros();
pyerror(message);
  }

  /* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .*/
  /* Read user-supplied parameters from the 'model' module they supply:
  */
  status = mallocInputParStrs(&par);
  if(status){
    unsetMacros();
    Py_DECREF(pModule);
    sprintf(message, "Function mallocInputParStrs() returned with status %d", status);
pyerror(message);
  }

  readParImgWrapper(pModule, parTemplates, nPars, imgParTemplates, nImgPars, &par, &img, &nImages);

  free(imgParTemplates);
  free(parTemplates);

  if(arguments.numThreads>0) /* Indicates that the user has set it to something sensible on the command line. */
    par.nThreads = arguments.numThreads;

  setUpUserPythonFuncs(pModule);

  Py_DECREF(pModule);

  /* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .*/
  /* Now call the main bit of LIME:
  */
  status = run(par, img, nImages);

  /* Python-object clean up before status check and possible exit.
  */
  decrefAllUserFuncs(); /* in py_utils.c */

  if(status){
    sprintf(message, "Function run() returned with status %d", status);
pyerror(message);
  }

  Py_Finalize();

  pyFreeInputImgPars(img, nImages);
  pyFreeInputPars(&par);

  return 0;
}

