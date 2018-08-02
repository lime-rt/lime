/*
 *  py_models.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#include "lime.h" /* for configInfo, and function definitions. */
#include "py_lime.h" /* for pyerror definitions. */
#include "defaults.h"
#include "py_utils.h"
#include "local_err.h"

#define PY_NUM_DIMS    3

/*....................................................................*/
void
density(double x, double y, double z, double *density){
  char *resultName="density";
  int numElemInUserFuncReturn,i;
  double resultsBuffer[UFUNC_BUFFER_SIZE];
  errType err=init_local_err();

  if(pDensity!=NULL){ /* User supplied this function. */
    err = userFuncWrapper(pDensity, resultName, -1, x, y, z, resultsBuffer, &numElemInUserFuncReturn);

    if(err.status!=0){ /* All we can do is quit. There's no 'nice' way, e.g. passing on the status value, since we're constrained by the function interface. */
      unsetMacros();
      decrefAllUserFuncs();
      pyerror(err.message);
    }

    for(i=0;i<numElemInUserFuncReturn;i++)
      density[i] = resultsBuffer[i];
  }else
    default_density(x, y, z, density);
}

/*....................................................................*/
void
temperature(double x, double y, double z, double *temperature){
  char *resultName="temperature";
  int numElemInUserFuncReturn,i;
  double resultsBuffer[UFUNC_BUFFER_SIZE];
  errType err=init_local_err();

  if(pTemperature!=NULL){ /* User supplied this function. */
    err = userFuncWrapper(pTemperature, resultName, 2, x, y, z, resultsBuffer, &numElemInUserFuncReturn);

    if(err.status!=0){ /* All we can do is quit. There's no 'nice' way, e.g. passing on the status value, since we're constrained by the function interface. */
      unsetMacros();
      decrefAllUserFuncs();
      pyerror(err.message);
    }

    for(i=0;i<numElemInUserFuncReturn;i++)
      temperature[i] = resultsBuffer[i];
  }else /* There's no library function specified. */
    default_temperature(x, y, z, temperature);
}

/*....................................................................*/
void
abundance(double x, double y, double z, double *abundance){
  char *resultName="abundance";
  int numElemInUserFuncReturn,i;
  double resultsBuffer[UFUNC_BUFFER_SIZE];
  errType err=init_local_err();

  if(pAbundance!=NULL){ /* User supplied this function. */
    err = userFuncWrapper(pAbundance, resultName, -1, x, y, z, resultsBuffer, &numElemInUserFuncReturn);

    if(err.status!=0){ /* All we can do is quit. There's no 'nice' way, e.g. passing on the status value, since we're constrained by the function interface. */
      unsetMacros();
      decrefAllUserFuncs();
      pyerror(err.message);
    }

    for(i=0;i<numElemInUserFuncReturn;i++)
      abundance[i] = resultsBuffer[i];
  }else /* There's no library function specified. */
    default_abundance(x, y, z, abundance);
}

/*....................................................................*/
void
molNumDensity(double x, double y, double z, double *molNumDens){
  char *resultName="molNumDensity";
  int numElemInUserFuncReturn,i;
  double resultsBuffer[UFUNC_BUFFER_SIZE];
  errType err=init_local_err();

  if(pMolNumDensity!=NULL){ /* User supplied this function. */
    err = userFuncWrapper(pMolNumDensity, resultName, -1, x, y, z, resultsBuffer, &numElemInUserFuncReturn);

    if(err.status!=0){ /* All we can do is quit. There's no 'nice' way, e.g. passing on the status value, since we're constrained by the function interface. */
      unsetMacros();
      decrefAllUserFuncs();
      pyerror(err.message);
    }

    for(i=0;i<numElemInUserFuncReturn;i++)
      molNumDens[i] = resultsBuffer[i];
  }else /* There's no library function specified. */
    default_molNumDensity(x, y, z, molNumDens);
}

/*....................................................................*/
void
doppler(double x, double y, double z, double *doppler){
  char *resultName="doppler";
  int numElemInUserFuncReturn,i;
  double resultsBuffer[UFUNC_BUFFER_SIZE];
  errType err=init_local_err();

  if(pDoppler!=NULL){ /* User supplied this function. */
    err = userFuncWrapper(pDoppler, resultName, 1, x, y, z, resultsBuffer, &numElemInUserFuncReturn);

    if(err.status!=0){ /* All we can do is quit. There's no 'nice' way, e.g. passing on the status value, since we're constrained by the function interface. */
      unsetMacros();
      decrefAllUserFuncs();
      pyerror(err.message);
    }

    i=0;
    *doppler = resultsBuffer[i];
  }else /* There's no library function specified. */
    default_doppler(x, y, z, doppler);
}

/*....................................................................*/
void
velocity(double x, double y, double z, double *vel){
  char *resultName="velocity";
  int numElemInUserFuncReturn,i;
  double resultsBuffer[UFUNC_BUFFER_SIZE];
  errType err=init_local_err();

  if(pVelocity!=NULL){ /* User supplied this function. */
    err = userFuncWrapper(pVelocity, resultName, PY_NUM_DIMS, x, y, z, resultsBuffer, &numElemInUserFuncReturn);

    if(err.status!=0){ /* All we can do is quit. There's no 'nice' way, e.g. passing on the status value, since we're constrained by the function interface. */
      unsetMacros();
      decrefAllUserFuncs();
      pyerror(err.message);
    }

    for(i=0;i<numElemInUserFuncReturn;i++)
      vel[i] = resultsBuffer[i];
  }else /* There's no library function specified. */
    default_velocity(x, y, z, vel);
}

/*....................................................................*/
void
magfield(double x, double y, double z, double *B){
  char *resultName="magfield";
  int numElemInUserFuncReturn,i;
  double resultsBuffer[UFUNC_BUFFER_SIZE];
  errType err=init_local_err();

  if(pMagfield!=NULL){ /* User supplied this function. */
    err = userFuncWrapper(pMagfield, resultName, PY_NUM_DIMS, x, y, z, resultsBuffer, &numElemInUserFuncReturn);

    if(err.status!=0){ /* All we can do is quit. There's no 'nice' way, e.g. passing on the status value, since we're constrained by the function interface. */
      unsetMacros();
      decrefAllUserFuncs();
      pyerror(err.message);
    }

    for(i=0;i<numElemInUserFuncReturn;i++)
      B[i] = resultsBuffer[i];
  }else /* There's no library function specified. */
    default_magfield(x, y, z, B);
}

/*....................................................................*/
void
gasIIdust(double x, double y, double z, double *gas2dust){
  char *resultName="gasIIdust";
  int numElemInUserFuncReturn,i;
  double resultsBuffer[UFUNC_BUFFER_SIZE];
  errType err=init_local_err();

  if(pGasIIdust!=NULL){ /* User supplied this function. */
    err = userFuncWrapper(pGasIIdust, resultName, 1, x, y, z, resultsBuffer, &numElemInUserFuncReturn);

    if(err.status!=0){ /* All we can do is quit. There's no 'nice' way, e.g. passing on the status value, since we're constrained by the function interface. */
      unsetMacros();
      decrefAllUserFuncs();
      pyerror(err.message);
    }

    i=0;
    *gas2dust = resultsBuffer[i];
  }else /* There's no library function specified. */
    default_gasIIdust(x, y, z, gas2dust);
}

/*....................................................................*/
double
gridDensity(configInfo *par, double *rVec){
  char *resultName="gridDensity";
  int numElemInUserFuncReturn,i;
  double resultsBuffer[UFUNC_BUFFER_SIZE];
  double value;
  errType err=init_local_err();

  if(pGridDensity!=NULL){ /* User supplied this function. */
    /*  ***** NOTE ***** that we are throwing away the config info! */
    err = userFuncWrapper(pGridDensity, resultName, 1, rVec[0], rVec[1], rVec[2], resultsBuffer, &numElemInUserFuncReturn);

    if(err.status!=0){ /* All we can do is quit. There's no 'nice' way, e.g. passing on the status value, since we're constrained by the function interface. */
      unsetMacros();
      decrefAllUserFuncs();
      pyerror(err.message);
    }

    i=0;
    value = resultsBuffer[i];
  }else /* There's no library function specified. */
    value = default_gridDensity(par, rVec, density);

  return value;
}


