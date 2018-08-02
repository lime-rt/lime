/*
 *  py_main.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#include "lime.h"
#include "error_codes.h"
#include "local_err.h"
#include "py_lime.h"
#include <getopt.h>
#include "py_utils.h"

#ifdef NOVERBOSE
int silent = 1;
#else
int silent = 0;
#endif

int defaultFuncFlags = 0;

typedef struct {
  _Bool doSilent,fixSeeds;
  char *modelfilename;
  int numThreads;
} infoType;

/*....................................................................*/
void
printUsageAndQuit(void){
  printf("\nUseage:\n  pylime [-s] [-t] [-p <num threads>] <model file>\n");

exit(1);
}

/*....................................................................*/
infoType
_initInfo(void){
  infoType info;

  info.doSilent = FALSE;
  info.fixSeeds = FALSE;
  info.numThreads = 0;
  info.modelfilename = NULL;

return info;
}

/*....................................................................*/
errType
strToInt(char *inStr, int *result){
  errType err=init_local_err();
  char *nonIntPtr;
  long trialLong;
  char message[ERR_STR_LEN]; /* Allow an extra character for the \0 in snprintf. */

  trialLong = strtol(inStr, &nonIntPtr, 0);
  if (*nonIntPtr != '\0'){
    snprintf(message, ERR_STR_LEN-1, "Supposed integer string field %s contains bad characters.", inStr);
return write_local_err(1, message);
  }

  /* The compiler directive is needed only when `int` and `long` have different ranges. */
#if LONG_MIN < INT_MIN || LONG_MAX > INT_MAX
  if (trialLong < INT_MIN || trialLong > INT_MAX){
    snprintf(message, ERR_STR_LEN-1, "Integer string field %s is out of range.", inStr);
return write_local_err(2, message);
  }
#endif

  *result = (int)trialLong;

return err;
}

/*....................................................................*/
infoType
_parseArgs(int argc, char *argv[]){
  /*
This extracts information from the command-line arguments, storing it in the returned struct 'info'.

  */
  infoType info;
  int optI;
  /* getopt_long stores the option index here. */
  int option_index=0,numNonOptionArgs;
  char message[STR_LEN_0+1]; /* Allow an extra character for the \0 in snprintf. */
  const int specifiedNumNonOptionArgs=1;
  errType err=init_local_err();

  info = _initInfo();

  while(TRUE){
    static struct option long_options[] = {
          {"version", no_argument,       0, 'v'},
          {"help",    no_argument,       0, 'h'},
          {"silent",  no_argument,       0, 's'},
          {"test",    no_argument,       0, 't'},
          {"threads", required_argument, 0, 'p'},
          {0, 0, 0, 0}
    };

    optI = getopt_long(argc, argv, "stp:",
                       long_options, &option_index);

    /* Detect the end of the options. */
    if (optI == -1)
  break;

    switch (optI){
        case 'v':
          printf("LIME (pylime) %s\n", VERSION);
exit(0);

        case 'h':
printUsageAndQuit();

        case 's':
          info.doSilent = TRUE;
    break;

        case 't':
          info.fixSeeds = TRUE;
    break;

        case 'p':
          err = strToInt(optarg, &info.numThreads);
          if(err.status){
error(1,err.message);
          }
    break;

        case '?':
          /* getopt_long already printed an error message. */
    break;

        default:
printUsageAndQuit();
    }
  }

  /* Parse the non-option arguments:
  */
  numNonOptionArgs = argc - optind;
  if(numNonOptionArgs<specifiedNumNonOptionArgs)
error(1,"Too few non-option args");

  if(numNonOptionArgs>specifiedNumNonOptionArgs){
    snprintf(message, STR_LEN_0, "You have %d non-option arguments but you only should have %d.", numNonOptionArgs, specifiedNumNonOptionArgs);
    warning(message);
  }

  info.modelfilename = argv[optind++];

return info;
}

/*....................................................................*/
void
_readParImgWrapper(PyObject *pModule, parTemplateType *parTemplates\
  , const int nPars, parTemplateType *imgParTemplates, const int nImgPars\
  , inputPars *par, image **img, int *nImages){

  errType err=init_local_err();
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

  /*
Unpack all the model and image parameters from the pPars object, returning them as (i) a pointer 'par' to an inputPars struct, and (ii) a pointer 'img' to a series of image structs.
  */
  err = readParImg(pPars, parTemplates, nPars, imgParTemplates, nImgPars, par, img, nImages, pywarning); /* In py_utils.c */
  if(err.status){
    Py_DECREF(pPars);
pyerror(err.message);
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

  infoType info;
  errType err=init_local_err();
  const int lenSuffix=3,maxLenNoSuffix=STR_LEN_0;
  const char *nameOfExecutable="pylime",*headerModuleName="limepar_classes";
  const char *oldModulePath;
  char modelNameNoSuffix[maxLenNoSuffix+1],message[STR_LEN_0],*newModulePath;
  char suffix[lenSuffix+1];
  int lenNoSuffix,nPars,nImgPars,nImages,status=0,strlenOMPath;
  PyObject *pModule = NULL;
  inputPars par;
  image *img = NULL;
  parTemplateType *parTemplates=NULL,*imgParTemplates=NULL;

  /* Parse our arguments; every option seen by parse_opt will be reflected in arguments.
  */
  info = _parseArgs(argc, argv);

  /* Interpret the options (numThreads is left until later):
  */
  if(info.doSilent)
    silent = 1;

  if(info.fixSeeds)
    fixRandomSeeds = TRUE;

  /* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .*/
  /* Get the model name, then strip off the '.py' from it:
  */
  if(info.modelfilename==NULL){
    sprintf(message, "Usage: %s [<arguments>] <python model file>", nameOfExecutable);
    error(1, message);
  }

  lenNoSuffix = strlen(info.modelfilename) - lenSuffix;
  if(lenNoSuffix<1){
    sprintf(message, "Model file name must be more than %d characters long!", lenSuffix);
    error(1, message);
  }
  if(lenNoSuffix>maxLenNoSuffix){
    sprintf(message, "Model file name is longer than the permitted %d characters.\n", maxLenNoSuffix);
    error(1, message);
  }

  strncpy(suffix, info.modelfilename + lenNoSuffix, lenSuffix);
  suffix[lenSuffix] = '\0';

  if(strcmp(suffix,".py")!=0){
    sprintf(message, "Python files must end in '.py'");
    error(1, message);
  }

  strncpy(modelNameNoSuffix, info.modelfilename, strlen(info.modelfilename)-lenSuffix);
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
  err = getParTemplatesWrapper(headerModuleName, &parTemplates, &nPars\
    , &imgParTemplates, &nImgPars); /* in py_utils.c */
  if(err.status!=0){
pyerror(err.message);
  }

  /* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .*/
  /* Construct the 'macro' list argument:
  */
  err = setMacros(); /* in py_utils.c */
  if(err.status){
    unsetMacros(); /* in py_utils.c */
pyerror(err.message);
  }

  /* Now we open the user's 'model' module:
  */
  err = getModuleFromName(modelNameNoSuffix, &pModule);
  if(err.status){
    if(!silent)
      PyErr_Print();
    unsetMacros();
pyerror(err.message);
  }

  /* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .*/
  /* Read user-supplied parameters from the 'model' module they supply:
  */
  err = mallocInputParStrs(&par);
  if(err.status){
    unsetMacros();
    Py_DECREF(pModule);
pyerror(err.message);
  }

  _readParImgWrapper(pModule, parTemplates, nPars, imgParTemplates, nImgPars, &par, &img, &nImages);

  free(imgParTemplates);
  free(parTemplates);

  if(info.numThreads>0) /* Indicates that the user has set it to something sensible on the command line. */
    par.nThreads = info.numThreads;

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

