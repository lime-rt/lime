/*
 *  ml_models.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
************ change sprintfs to snprintfs
 */

#include "defaults.h" /* includes lime_config.h for configInfo */
#include "ml_recipes/ml_recipes.h"
#include "ml_models.h" /* for modelParamType */
#include "py_utils.h" /* for getPythonFunc() etc */
#include "ufunc_types.h" /* for all user-function templates */
#include "messages.h" /* for all user-function templates */
#include "local_err.h"

int currentModelI = -1; /* <0 signifies it has not yet been set. */
int *modelIntPars=NULL,numModelParams=0;
double *modelDblPars=NULL;
char *modelStrPar=NULL;

modelParamType *modelParams=NULL;

/*....................................................................*/
int
getParamI(char *idStr, int *index){
  int status=1,i;

  *index = -1; /* Nonsensical value */
  for(i=0;i<numModelParams;i++){
    if(strcmp(idStr, modelParams[i].name)==0){
      status = 0;
      *index = modelParams[i].index;
      break;
    }
  }

  return status;
}

/*....................................................................*/
int getModelIFromName(char *modelName){
  int modelI=-1;

  if(strcmp(modelName,"BonnorEbert56")==0)
    modelI = MODEL_BoEb56;
  else if(strcmp(modelName,"CG97")==0)
    modelI = MODEL_CG97;
  else if(strcmp(modelName,"DDN01")==0)
    modelI = MODEL_DDN01;
  else if(strcmp(modelName,"LiShu96")==0)
    modelI = MODEL_LiSh96;
  else if(strcmp(modelName,"Mamon88")==0)
    modelI = MODEL_Ma88;
  else if(strcmp(modelName,"Mendoza09")==0)
    modelI = MODEL_Me09;
  else if(strcmp(modelName,"Shu77")==0)
    modelI = MODEL_Shu77;
  else if(strcmp(modelName,"Ulrich76")==0)
    modelI = MODEL_Ul76;
  else if(strcmp(modelName,"allen03a")==0)
    modelI = MODEL_Al03;
  else
    modelI = -1; /* signals that model was not recognized. */

  return modelI;
}

/*....................................................................*/
int getModelNameFromI(const int modelI, char *modelName){
  int status=0;

  switch(modelI){
    case MODEL_BoEb56:
      modelName = "BonnorEbert56";
      break;
    case MODEL_CG97:
      modelName = "CG97";
      break;
    case MODEL_DDN01:
      modelName = "DDN01";
      break;
    case MODEL_LiSh96:
      modelName = "LiShu96";
      break;
    case MODEL_Ma88:
      modelName = "Mamon88";
      break;
    case MODEL_Me09:
      modelName = "Mendoza09";
      break;
    case MODEL_Shu77:
      modelName = "Shu77";
      break;
    case MODEL_Ul76:
      modelName = "Ulrich76";
      break;
    case MODEL_Al03:
      modelName = "allen03a";
      break;
    default:
      modelName = '\0'; /* signals that model was not recognized. */
  }

  return status;
}

/*....................................................................*/
int finalizeModelConfig(const int modelI){
  int status=0;

  switch(modelI){
    case MODEL_Al03:
      status = Al03_onFinalizeConfiguration();
      break;
    case MODEL_BoEb56:
      status = BoEb56_onFinalizeConfiguration();
      break;
    case MODEL_CG97:
      status = CG97_onFinalizeConfiguration();
      break;
    case MODEL_DDN01:
      status = DDN01_onFinalizeConfiguration();
      break;
    case MODEL_LiSh96:
      status = LiSh96_onFinalizeConfiguration();
      break;
    case MODEL_Ma88:
      status = Ma88_onFinalizeConfiguration();
      break;
    case MODEL_Me09:
      status = Me09_onFinalizeConfiguration();
      break;
    case MODEL_Shu77:
      status = Shu77_onFinalizeConfiguration();
      break;
    case MODEL_Ul76:
      status = Ul76_onFinalizeConfiguration();
      break;
    case MODEL_None:
      /* no action needed. */
      break;
    default:
      status = ML_UNRECOG_MODEL;
  }

  return status;
}

/*....................................................................*/
void
density(double x, double y, double z, double *density){
  char *ufuncName="density";
  int numElemInUserFuncReturn,i;
  double resultsBuffer[UFUNC_BUFFER_SIZE];
  errType err=init_local_err();

  switch(currentModelI){
    case MODEL_Al03:
      density[0] = Al03_density(x,y,z);
      break;
    case MODEL_BoEb56:
      density[0] = BoEb56_density(x,y,z);
      break;
    case MODEL_CG97:
      density[0] = CG97_density(x,y,z);
      break;
    case MODEL_DDN01:
      density[0] = DDN01_density(x,y,z);
      break;
    case MODEL_LiSh96:
      density[0] = LiSh96_density(x,y,z);
      break;
    case MODEL_Ma88:
      density[0] = Ma88_density(x,y,z);
      break;
    case MODEL_Me09:
      density[0] = Me09_density(x,y,z);
      break;
    case MODEL_Shu77:
      density[0] = Shu77_density(x,y,z);
      break;
    case MODEL_Ul76:
      density[0] = Ul76_density(x,y,z);
      break;
    default:
      if(pDensity!=NULL){ /* User supplied this function. */
        err = userFuncWrapper(pDensity, ufuncName, -1, x, y, z, resultsBuffer, &numElemInUserFuncReturn);

        if(err.status!=0){ /* All we can do is quit. There's no 'nice' way, e.g. passing on the status value, since we're constrained by the function interface. */
          unsetMacros();
          decrefAllUserFuncs();
          if(!silent){
            bail_out(err.message);
          }
exit(1);
        }

        for(i=0;i<numElemInUserFuncReturn;i++)
          density[i] = resultsBuffer[i];
      }else{
        default_density(x, y, z, density);
      }
  } 

}

/*....................................................................*/
void
temperature(double x, double y, double z, double *temperature){
  char *ufuncName="temperature";
  double rVec[ML_NUM_DIMS],*funcParA=NULL,*funcParB=NULL;
  int funcIa=-1,funcIb=-1;
  double resultsBuffer[UFUNC_BUFFER_SIZE];
  errType err=init_local_err();
  int numElemInUserFuncReturn,i;

  if(currentModelI>=0){//********************* shouldn't it dump out if no model has been set?
    funcIa   = funcIs[  currentModelI][RESULT_temperature];
    funcParA = funcPars[currentModelI][RESULT_temperature];
    funcIb   = funcIs[  currentModelI][RESULT_tdust];
    funcParB = funcPars[currentModelI][RESULT_tdust];
  }

  rVec[0] = x;
  rVec[1] = y;
  rVec[2] = z;

  switch(currentModelI){
    case MODEL_Al03:
      temperature[0] = Al03_temperature(x,y,z);
      if(copyTemp>0) temperature[1] = temperature[0];
      break;
    case MODEL_BoEb56:
      temperature[0] = BoEb56_temperature(x,y,z);
      if(copyTemp>0) temperature[1] = temperature[0];
      break;
    case MODEL_CG97:
      temperature[1] = CG97_t_dust(x,y,z);
      if(copyTemp<0) temperature[0] = temperature[1];
      break;
    case MODEL_DDN01:
      temperature[1] = DDN01_t_dust(x,y,z);
      if(copyTemp<0) temperature[0] = temperature[1];
      break;
    case MODEL_LiSh96:
      temperature[0] = LiSh96_temperature(x,y,z);
      if(copyTemp>0) temperature[1] = temperature[0];
      break;
    case MODEL_Ma88:
      temperature[0] = Ma88_temperature(x,y,z);
      if(copyTemp>0) temperature[1] = temperature[0];
      break;
//    case MODEL_Me09:
//      break;
    case MODEL_Shu77:
      temperature[0] = Shu77_temperature(x,y,z);
      temperature[1] = Shu77_t_dust(x,y,z);
      break;
//    case MODEL_Ul76:
//      break;
    default: /* Catches MODEL_Me09 and MODEL_Ul76. */
      if(pTemperature!=NULL){ /* User supplied this function. */
        err = userFuncWrapper(pTemperature, ufuncName, 2, x, y, z, resultsBuffer, &numElemInUserFuncReturn);

        if(err.status!=0){ /* All we can do is quit. There's no 'nice' way, e.g. passing on the status value, since we're constrained by the function interface. */
          unsetMacros();
          decrefAllUserFuncs();
          if(!silent){
            bail_out(err.message);
          }
exit(1);
        }

        for(i=0;i<numElemInUserFuncReturn;i++)
          temperature[i] = resultsBuffer[i];

      }else if(funcIa<0 && funcIb<0){ /* There's no library function specified. */
        default_temperature(x, y, z, temperature);

      }else if(funcIb<0){ /* then funcIa must be >= 0 */
        temperature[0] = scalarFunctionSwitch(funcIa, funcParA, rVec);
        if(copyTemp>0) temperature[1] = temperature[0];
      }else if(funcIa<0){ /* then funcIb must be >= 0 */
        temperature[1] = scalarFunctionSwitch(funcIb, funcParB, rVec);
        if(copyTemp<0) temperature[0] = temperature[1];
      }else{ /* then both funcIa and funcIb must be >= 0 */
        temperature[0] = scalarFunctionSwitch(funcIa, funcParA, rVec);
        temperature[1] = scalarFunctionSwitch(funcIb, funcParB, rVec);
      }
  }
}

/*....................................................................*/
void
abundance(double x, double y, double z, double *abundance){
  //**** I'm going to assume there is only 1 radiating species (thus only 1 abundance element) for the time being. Note that abundance() should never be called in the case that there are zero species.
  char *ufuncName="abundance";
  double rVec[ML_NUM_DIMS],*funcPar=NULL;
  const int resultI=RESULT_abundance;
  int funcI=-1;
  double resultsBuffer[UFUNC_BUFFER_SIZE];
  errType err=init_local_err();
  int numElemInUserFuncReturn,i;

  if(currentModelI>=0){
    funcI   = funcIs[  currentModelI][resultI];
    funcPar = funcPars[currentModelI][resultI];
  }

  rVec[0] = x;
  rVec[1] = y;
  rVec[2] = z;

  switch(currentModelI){
//    case MODEL_BoEb56:
//      break;
//    case MODEL_CG97:
//      break;
//    case MODEL_DDN01:
//      break;
//    case MODEL_LiSh96:
//      break;
    case MODEL_Ma88:
      abundance[0] = Ma88_abundance(x,y,z);
      break;
//    case MODEL_Me09:
//      break;
//    case MODEL_Shu77:
//      break;
//    case MODEL_Ul76:
//      break;
//    case MODEL_Al03:
//      break;
    default:
      if(pAbundance!=NULL){ /* User supplied this function. */
        err = userFuncWrapper(pAbundance, ufuncName, -1, x, y, z, resultsBuffer, &numElemInUserFuncReturn);

        if(err.status!=0){ /* All we can do is quit. There's no 'nice' way, e.g. passing on the status value, since we're constrained by the function interface. */
          unsetMacros();
          decrefAllUserFuncs();
          if(!silent){
            bail_out(err.message);
          }
exit(1);
        }

        for(i=0;i<numElemInUserFuncReturn;i++)
          abundance[i] = resultsBuffer[i];
      }else if(funcI<0){ /* There's no library function specified. */
        default_abundance(x, y, z, abundance);
      }else{
        abundance[0] = scalarFunctionSwitch(funcI, funcPar, rVec);
      }
  }
}

/*....................................................................*/
void
molNumDensity(double x, double y, double z, double *molNumDens){
  char *ufuncName="molNumDensity";
  double resultsBuffer[UFUNC_BUFFER_SIZE];
  errType err=init_local_err();
  int numElemInUserFuncReturn,i;

  if(pMolNumDensity!=NULL){ /* User supplied this function. */
    err = userFuncWrapper(pMolNumDensity, ufuncName, -1, x, y, z, resultsBuffer, &numElemInUserFuncReturn);

    if(err.status!=0){ /* All we can do is quit. There's no 'nice' way, e.g. passing on the status value, since we're constrained by the function interface. */
      unsetMacros();
      decrefAllUserFuncs();
      if(!silent){
        bail_out(err.message);
      }
exit(1);
    }

    for(i=0;i<numElemInUserFuncReturn;i++)
      molNumDens[i] = resultsBuffer[i];
  }else /* There's no library function specified. */
    default_molNumDensity(x, y, z, molNumDens);
}

/*....................................................................*/
void
doppler(double x, double y, double z, double *doppler){
  char *ufuncName="doppler";
  double rVec[ML_NUM_DIMS],*funcPar=NULL;
  const int resultI=RESULT_doppler;
  int funcI=-1;
  double resultsBuffer[UFUNC_BUFFER_SIZE];
  errType err=init_local_err();
  int numElemInUserFuncReturn,i;

  if(currentModelI>=0){
    funcI   = funcIs[  currentModelI][resultI];
    funcPar = funcPars[currentModelI][resultI];
  }

  rVec[0] = x;
  rVec[1] = y;
  rVec[2] = z;

  if(pDoppler!=NULL){ /* User supplied this function. */
    err = userFuncWrapper(pDoppler, ufuncName, 1, x, y, z, resultsBuffer, &numElemInUserFuncReturn);

    if(err.status!=0){ /* All we can do is quit. There's no 'nice' way, e.g. passing on the status value, since we're constrained by the function interface. */
      unsetMacros();
      decrefAllUserFuncs();
      if(!silent){
        bail_out(err.message);
      }
exit(1);
    }

    i=0;
    *doppler = resultsBuffer[i];
  }else if(funcI<0){ /* There's no library function specified. */
    default_doppler(x, y, z, doppler);
  }else{
    *doppler = scalarFunctionSwitch(funcI, funcPar, rVec);
  }
}

/*....................................................................*/
void
velocity(double x, double y, double z, double *vel){
  char *ufuncName="velocity";
  double rVec[ML_NUM_DIMS],*funcPar=NULL;
  const int resultI=RESULT_velocity;
  int funcI=-1;
  double resultsBuffer[UFUNC_BUFFER_SIZE];
  int numElemInUserFuncReturn,i;
  errType err=init_local_err();

  if(currentModelI>=0){
    funcI   = funcIs[  currentModelI][resultI];
    funcPar = funcPars[currentModelI][resultI];
  }

  rVec[0] = x;
  rVec[1] = y;
  rVec[2] = z;

  switch(currentModelI){
    case MODEL_Al03:
      Al03_velocity(x,y,z,vel);
      break;
    case MODEL_BoEb56:
      BoEb56_velocity(x,y,z,vel);
      break;
    case MODEL_CG97:
      CG97_velocity(x,y,z,vel);
      break;
    case MODEL_DDN01:
      DDN01_velocity(x,y,z,vel);
      break;
    case MODEL_LiSh96:
      LiSh96_velocity(x,y,z,vel);
      break;
    case MODEL_Ma88:
      Ma88_velocity(x,y,z,vel);
      break;
    case MODEL_Me09:
      Me09_velocity(x,y,z,vel);
      break;
    case MODEL_Shu77:
      Shu77_velocity(x,y,z,vel);
      break;
    case MODEL_Ul76:
      Ul76_velocity(x,y,z,vel);
      break;
    default:
      if(pVelocity!=NULL){ /* User supplied this function. */
        err = userFuncWrapper(pVelocity, ufuncName, ML_NUM_DIMS, x, y, z, resultsBuffer, &numElemInUserFuncReturn);

        if(err.status!=0){ /* All we can do is quit. There's no 'nice' way, e.g. passing on the status value, since we're constrained by the function interface. */
          unsetMacros();
          decrefAllUserFuncs();
          if(!silent){
            bail_out(err.message);
          }
exit(1);
        }

        for(i=0;i<numElemInUserFuncReturn;i++)
          vel[i] = resultsBuffer[i];
      }else if(funcI<0){ /* There's no library function specified. */
        default_velocity(x, y, z, vel);
      }else{
        vectorFunctionSwitch(funcI, funcPar, rVec, vel);
      }
  }
}

/*....................................................................*/
void
magfield(double x, double y, double z, double *B){
  char *ufuncName="magfield";
  double rVec[ML_NUM_DIMS],*funcPar=NULL;
  const int resultI=RESULT_bmag;
  int funcI=-1;
  double resultsBuffer[UFUNC_BUFFER_SIZE];
  errType err=init_local_err();
  int numElemInUserFuncReturn,i;

  if(currentModelI>=0){
    funcI   = funcIs[  currentModelI][resultI];
    funcPar = funcPars[currentModelI][resultI];
  }

  rVec[0] = x;
  rVec[1] = y;
  rVec[2] = z;

  switch(currentModelI){
    case MODEL_Al03:
      Al03_bmag(x, y, z, B);
      break;
//    case MODEL_BoEb56:
//      break;
//    case MODEL_CG97:
//      break;
//    case MODEL_DDN01:
//      break;
    case MODEL_LiSh96:
      LiSh96_bmag(x, y, z, B);
      break;
//    case MODEL_Ma88:
//      break;
//    case MODEL_Me09:
//      break;
//    case MODEL_Shu77:
//      break;
//    case MODEL_Ul76:
//      break;
    default:
      if(pMagfield!=NULL){ /* User supplied this function. */
        err = userFuncWrapper(pMagfield, ufuncName, ML_NUM_DIMS, x, y, z, resultsBuffer, &numElemInUserFuncReturn);

        if(err.status!=0){ /* All we can do is quit. There's no 'nice' way, e.g. passing on the status value, since we're constrained by the function interface. */
          unsetMacros();
          decrefAllUserFuncs();
          if(!silent){
            bail_out(err.message);
          }
exit(1);
        }

        for(i=0;i<numElemInUserFuncReturn;i++)
          B[i] = resultsBuffer[i];
      }else if(funcI<0){ /* There's no library function specified. */
        default_magfield(x, y, z, B);
      }else{
        vectorFunctionSwitch(funcI, funcPar, rVec, B);
      }
  }
}

/*....................................................................*/
void
gasIIdust(double x, double y, double z, double *gas2dust){
  char *ufuncName="gasIIdust";
  double resultsBuffer[UFUNC_BUFFER_SIZE];
  errType err=init_local_err();
  int numElemInUserFuncReturn,i;

  if(pGasIIdust!=NULL){ /* User supplied this function. */
    err = userFuncWrapper(pGasIIdust, ufuncName, 1, x, y, z, resultsBuffer, &numElemInUserFuncReturn);

    if(err.status!=0){ /* All we can do is quit. There's no 'nice' way, e.g. passing on the status value, since we're constrained by the function interface. */
      unsetMacros();
      decrefAllUserFuncs();
      if(!silent){
        bail_out(err.message);
      }
exit(1);
    }

    i=0;
    *gas2dust = resultsBuffer[i];
  }else /* There's no library function specified. */
    default_gasIIdust(x, y, z, gas2dust);
}

/*....................................................................*/
double
gridDensity(configInfo *par, double *rVec){
  char *ufuncName="gridDensity";
  double value;
  double resultsBuffer[UFUNC_BUFFER_SIZE];
  errType err=init_local_err();
  int numElemInUserFuncReturn,i;

  if(pGridDensity!=NULL){ /* User supplied this function. */
    /*  ***** NOTE ***** that we are throwing away the config info! */
    err = userFuncWrapper(pGridDensity, ufuncName, 1, rVec[0], rVec[1], rVec[2], resultsBuffer, &numElemInUserFuncReturn);

    if(err.status!=0){ /* All we can do is quit. There's no 'nice' way, e.g. passing on the status value, since we're constrained by the function interface. */
      unsetMacros();
      decrefAllUserFuncs();
      if(!silent){
        bail_out(err.message);
      }
exit(1);
    }

    i=0;
    value = resultsBuffer[i];
  }else /* There's no library function specified. */
    value = default_gridDensity(par, rVec, density);

  return value;
}



