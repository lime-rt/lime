/*
 *  py_models.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2017 The LIME development team
 *
 */

#include "lime.h" /* for configInfo, and function definitions. */
#include "py_lime.h" /* for pyerror definitions. */
#include "defaults.h"
#include "py_utils.h"

#define PY_NUM_DIMS    3

/*....................................................................*/
void
density(double x, double y, double z, double *density){
  char *resultName="density",message[STR_LEN_0+1];
  int status=0,numElemInUserFuncReturn,i;
  double resultsBuffer[UFUNC_BUFFER_SIZE];

  if(pDensity!=NULL){ /* User supplied this function. */
    status = userFuncWrapper(pDensity, resultName, x, y, z, resultsBuffer, &numElemInUserFuncReturn);

    if(status!=0){ /* All we can do is quit. There's no 'nice' way, e.g. passing on the status value, since we're constrained by the function interface. */
      unsetMacros();
      decrefAllUserFuncs();
      sprintf(message, "User %s() function generated error status %d", resultName, status);
      pyerror(message);
    }

    for(i=0;i<numElemInUserFuncReturn;i++)
      density[i] = resultsBuffer[i];
  }else
    default_density(x, y, z, density);
}

/*....................................................................*/
void
temperature(double x, double y, double z, double *temperature){
  char *resultName="temperature",message[STR_LEN_0+1];
  int status=0,numElemInUserFuncReturn,i;
  double resultsBuffer[UFUNC_BUFFER_SIZE];

  if(pTemperature!=NULL){ /* User supplied this function. */
    status = userFuncWrapper(pTemperature, resultName, x, y, z, resultsBuffer, &numElemInUserFuncReturn);

    if(status!=0){ /* All we can do is quit. There's no 'nice' way, e.g. passing on the status value, since we're constrained by the function interface. */
      unsetMacros();
      decrefAllUserFuncs();
      sprintf(message, "User %s() function generated error status %d", resultName, status);
      pyerror(message);
    }

    if(numElemInUserFuncReturn!=2){
      unsetMacros();
      decrefAllUserFuncs();
      sprintf(message, "User %s() function should return %d results.", resultName, 2);
      pyerror(message);
    }

    for(i=0;i<numElemInUserFuncReturn;i++)
      temperature[i] = resultsBuffer[i];
  }else /* There's no library function specified. */
    default_temperature(x, y, z, temperature);
}

/*....................................................................*/
void
abundance(double x, double y, double z, double *abundance){
  char *resultName="abundance",message[STR_LEN_0+1];
  int status=0,numElemInUserFuncReturn,i;
  double resultsBuffer[UFUNC_BUFFER_SIZE];

  if(pAbundance!=NULL){ /* User supplied this function. */
    status = userFuncWrapper(pAbundance, resultName, x, y, z, resultsBuffer, &numElemInUserFuncReturn);

    if(status!=0){ /* All we can do is quit. There's no 'nice' way, e.g. passing on the status value, since we're constrained by the function interface. */
      unsetMacros();
      decrefAllUserFuncs();
      sprintf(message, "User %s() function generated error status %d", resultName, status);
      pyerror(message);
    }

    for(i=0;i<numElemInUserFuncReturn;i++)
      abundance[i] = resultsBuffer[i];
  }else /* There's no library function specified. */
    default_abundance(x, y, z, abundance);
}

/*....................................................................*/
void
molNumDensity(double x, double y, double z, double *molNumDens){
  char *resultName="molNumDensity",message[STR_LEN_0+1];
  int status=0,numElemInUserFuncReturn,i;
  double resultsBuffer[UFUNC_BUFFER_SIZE];

  if(pMolNumDensity!=NULL){ /* User supplied this function. */
    status = userFuncWrapper(pMolNumDensity, resultName, x, y, z, resultsBuffer, &numElemInUserFuncReturn);

    if(status!=0){ /* All we can do is quit. There's no 'nice' way, e.g. passing on the status value, since we're constrained by the function interface. */
      unsetMacros();
      decrefAllUserFuncs();
      sprintf(message, "User %s() function generated error status %d", resultName, status);
      pyerror(message);
    }

    for(i=0;i<numElemInUserFuncReturn;i++)
      molNumDens[i] = resultsBuffer[i];
  }else /* There's no library function specified. */
    default_molNumDensity(x, y, z, molNumDens);
}

/*....................................................................*/
void
doppler(double x, double y, double z, double *doppler){
  char *resultName="doppler",message[STR_LEN_0+1];
  int status=0,numElemInUserFuncReturn,i;
  double resultsBuffer[UFUNC_BUFFER_SIZE];

  if(pDoppler!=NULL){ /* User supplied this function. */
    status = userFuncWrapper(pDoppler, resultName, x, y, z, resultsBuffer, &numElemInUserFuncReturn);

    if(status!=0){ /* All we can do is quit. There's no 'nice' way, e.g. passing on the status value, since we're constrained by the function interface. */
      unsetMacros();
      decrefAllUserFuncs();
      sprintf(message, "User %s() function generated error status %d", resultName, status);
      pyerror(message);
    }

    if(numElemInUserFuncReturn>1){
      unsetMacros();
      decrefAllUserFuncs();
      sprintf(message, "User %s() function should return a scalar result.", resultName);
      pyerror(message);
    }

    i=0;
    *doppler = resultsBuffer[i];
  }else /* There's no library function specified. */
    default_doppler(x, y, z, doppler);
}

/*....................................................................*/
void
velocity(double x, double y, double z, double *vel){
  char *resultName="velocity",message[STR_LEN_0+1];
  int status=0,numElemInUserFuncReturn,i;
  double resultsBuffer[UFUNC_BUFFER_SIZE];

  if(pVelocity!=NULL){ /* User supplied this function. */
    status = userFuncWrapper(pVelocity, resultName, x, y, z, resultsBuffer, &numElemInUserFuncReturn);

    if(status!=0){ /* All we can do is quit. There's no 'nice' way, e.g. passing on the status value, since we're constrained by the function interface. */
      unsetMacros();
      decrefAllUserFuncs();
      sprintf(message, "User %s() function generated error status %d", resultName, status);
      pyerror(message);
    }

    if(numElemInUserFuncReturn!=PY_NUM_DIMS){
      unsetMacros();
      decrefAllUserFuncs();
      sprintf(message, "User %s() function should return %d results.", resultName, PY_NUM_DIMS);
      pyerror(message);
    }

    for(i=0;i<numElemInUserFuncReturn;i++)
      vel[i] = resultsBuffer[i];
  }else /* There's no library function specified. */
    default_velocity(x, y, z, vel);
}

/*....................................................................*/
void
magfield(double x, double y, double z, double *B){
  char *resultName="magfield",message[STR_LEN_0+1];
  int status=0,numElemInUserFuncReturn,i;
  double resultsBuffer[UFUNC_BUFFER_SIZE];

  if(pMagfield!=NULL){ /* User supplied this function. */
    status = userFuncWrapper(pMagfield, resultName, x, y, z, resultsBuffer, &numElemInUserFuncReturn);

    if(status!=0){ /* All we can do is quit. There's no 'nice' way, e.g. passing on the status value, since we're constrained by the function interface. */
      unsetMacros();
      decrefAllUserFuncs();
      sprintf(message, "User %s() function generated error status %d", resultName, status);
      pyerror(message);
    }

    if(numElemInUserFuncReturn!=PY_NUM_DIMS){
      unsetMacros();
      decrefAllUserFuncs();
      sprintf(message, "User %s() function should return %d results.", resultName, PY_NUM_DIMS);
      pyerror(message);
    }

    for(i=0;i<numElemInUserFuncReturn;i++)
      B[i] = resultsBuffer[i];
  }else /* There's no library function specified. */
    default_magfield(x, y, z, B);
}

/*....................................................................*/
void
gasIIdust(double x, double y, double z, double *gas2dust){
  char *resultName="gasIIdust",message[STR_LEN_0+1];
  int status=0,numElemInUserFuncReturn,i;
  double resultsBuffer[UFUNC_BUFFER_SIZE];

  if(pGasIIdust!=NULL){ /* User supplied this function. */
    status = userFuncWrapper(pGasIIdust, resultName, x, y, z, resultsBuffer, &numElemInUserFuncReturn);

    if(status!=0){ /* All we can do is quit. There's no 'nice' way, e.g. passing on the status value, since we're constrained by the function interface. */
      unsetMacros();
      decrefAllUserFuncs();
      sprintf(message, "User %s() function generated error status %d", resultName, status);
      pyerror(message);
    }

    if(numElemInUserFuncReturn>1){
      unsetMacros();
      decrefAllUserFuncs();
      sprintf(message, "User %s() function should return a scalar result.", resultName);
      pyerror(message);
    }

    i=0;
    *gas2dust = resultsBuffer[i];
  }else /* There's no library function specified. */
    default_gasIIdust(x, y, z, gas2dust);
}

/*....................................................................*/
double
gridDensity(configInfo *par, double *rVec){
  char *resultName="gridDensity",message[STR_LEN_0+1];
  int status=0,numElemInUserFuncReturn,i;
  double resultsBuffer[UFUNC_BUFFER_SIZE];
  double value;

  if(pGridDensity!=NULL){ /* User supplied this function. */
    /*  ***** NOTE ***** that we are throwing away the config info! */
    status = userFuncWrapper(pGridDensity, resultName, rVec[0], rVec[1], rVec[2], resultsBuffer, &numElemInUserFuncReturn);

    if(status!=0){ /* All we can do is quit. There's no 'nice' way, e.g. passing on the status value, since we're constrained by the function interface. */
      unsetMacros();
      decrefAllUserFuncs();
      sprintf(message, "User %s() function generated error status %d", resultName, status);
      pyerror(message);
    }

    if(numElemInUserFuncReturn>1){
      unsetMacros();
      decrefAllUserFuncs();
      sprintf(message, "User %s() function should return a scalar result.", resultName);
      pyerror(message);
    }

    i=0;
    value = resultsBuffer[i];
  }else /* There's no library function specified. */
    value = default_gridDensity(par, rVec, density);

  return value;
}


