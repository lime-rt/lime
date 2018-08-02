/*
 *  ml_funcs.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#include "ml_funcs.h"

int funcIs[NUM_MODELS][NUM_RESULTS],funcTypeIs[NUM_MODELS][NUM_RESULTS];
double *funcPars[NUM_MODELS][NUM_RESULTS];

int numScalarFuncPars[] = FUNC_SCALAR_NPARS_LIST;
int numVectorFuncPars[] = FUNC_VECTOR_NPARS_LIST;

/*....................................................................*/
void _freeFuncPars(const int modelI){
  int ri;

  for(ri=0;ri<NUM_RESULTS;ri++){
    free(funcPars[modelI][ri]);
  }
}

/*....................................................................*/
void freeFuncsPars(void){
  int modelI;

  for(modelI=0;modelI<NUM_MODELS;modelI++)
    _freeFuncPars(modelI);
}

/*....................................................................*/
void _setDefaultFuncStuff(const int modelI){
  int ri;

  for(ri=0;ri<NUM_RESULTS;ri++){
    funcIs[    modelI][ri] = -1;
    funcTypeIs[modelI][ri] = -1;
    funcPars[  modelI][ri] = NULL;
  }
}

/*....................................................................*/
void setDefaultFuncStuffs(void){
  int mi;

  for(mi=0;mi<NUM_MODELS;mi++)
    _setDefaultFuncStuff(mi);

}

/*....................................................................*/
int getFuncTypeIFromName(char *funcTypeName){
  int funcTypeI=-1;

  if(     strcmp(funcTypeName,"scalar")==0)
    funcTypeI = FUNC_scalar;
  else if(strcmp(funcTypeName,"vector")==0)
    funcTypeI = FUNC_vector;
  else
    funcTypeI = -1; /* func type name not recognized. */

  return funcTypeI;
}

/*....................................................................*/
int getFuncTypeNameFromI(const int funcTypeI, char *funcTypeName){
  int status=0;

  switch(funcTypeI){
    case FUNC_scalar:
      funcTypeName = "scalar";
      break;
    case FUNC_vector:
      funcTypeName = "vector";
      break;
    default:
      funcTypeName = '\0'; /* funcTypeI not recognized. */
  }

  return status;
}

/*....................................................................*/
int getScalarFuncIFromName(char *funcName){
  int funcI=-1;

  if(strcmp(funcName,"scalarConst")==0)
    funcI = FUNC_scalarConst;
  else if(strcmp(funcName,"scalarPowerR")==0)
    funcI = FUNC_scalarPowerR;
  else if(strcmp(funcName,"scalarPowerRExpZ")==0)
    funcI = FUNC_scalarPowerRExpZ;
  else if(strcmp(funcName,"scalarPowerRTheta")==0)
    funcI = FUNC_scalarPowerRTheta;
  else if(strcmp(funcName,"scalarPowerRZ")==0)
    funcI = FUNC_scalarPowerRZ;
  else
    funcI = -1; /* func name not recognized. */

  return funcI;
}

/*....................................................................*/
int getScalarFuncNameFromI(const int funcI, char *funcName){
  int status=0;

  switch(funcI){
    case FUNC_scalarConst:
      funcName = "scalarConst";
      break;
    case FUNC_scalarPowerR:
      funcName = "scalarPowerR";
      break;
    case FUNC_scalarPowerRExpZ:
      funcName = "scalarPowerRExpZ";
      break;
    case FUNC_scalarPowerRTheta:
      funcName = "scalarPowerRTheta";
      break;
    case FUNC_scalarPowerRZ:
      funcName = "scalarPowerRZ";
      break;
    default:
      funcName = '\0'; /* funcI not recognized. */
  }

  return status;
}

/*....................................................................*/
int getVectorFuncIFromName(char *funcName){
  int funcI=-1;

  if(strcmp(funcName,"vectorConstR")==0)
    funcI = FUNC_vectorConstR;
  else if(strcmp(funcName,"vectorConstXYZ")==0)
    funcI = FUNC_vectorConstXYZ;
  else if(strcmp(funcName,"vectorDipole")==0)
    funcI = FUNC_vectorDipole;
  else if(strcmp(funcName,"vectorRadialPowerR")==0)
    funcI = FUNC_vectorRadialPowerR;
  else if(strcmp(funcName,"vectorRadialPowerRTheta")==0)
    funcI = FUNC_vectorRadialPowerRTheta;
  else if(strcmp(funcName,"vectorToroidalPowerR")==0)
    funcI = FUNC_vectorToroidalPowerR;
  else
    funcI = -1; /* func name not recognized. */

  return funcI;
}

/*....................................................................*/
int getVectorFuncNameFromI(const int funcI, char *funcName){
  int status=0;

  switch(funcI){
    case FUNC_vectorConstR:
      funcName = "vectorConstR";
      break;
    case FUNC_vectorConstXYZ:
      funcName = "vectorConstXYZ";
      break;
    case FUNC_vectorDipole:
      funcName = "vectorDipole";
      break;
    case FUNC_vectorRadialPowerR:
      funcName = "vectorRadialPowerR";
      break;
    case FUNC_vectorRadialPowerRTheta:
      funcName = "vectorRadialPowerRTheta";
      break;
    case FUNC_vectorToroidalPowerR:
      funcName = "vectorToroidalPowerR";
      break;
    default:
      funcName = '\0'; /* funcI not recognized. */
  }

  return status;
}

/*....................................................................*/
int getResultIFromName(char *resultName){
  int resultI=-1;

  if(     strcmp(resultName,"density")==0)
    resultI = RESULT_density;
  else if(strcmp(resultName,"temperature")==0)
    resultI = RESULT_temperature;
  else if(strcmp(resultName,"abundance")==0)
    resultI = RESULT_abundance;
  else if(strcmp(resultName,"doppler")==0)
    resultI = RESULT_doppler;
  else if(strcmp(resultName,"velocity")==0)
    resultI = RESULT_velocity;
  else if(strcmp(resultName,"bmag")==0)
    resultI = RESULT_bmag;
  else if(strcmp(resultName,"tdust")==0)
    resultI = RESULT_tdust;
  else
    resultI = -1; /* result name not recognized. */

  return resultI;
}

/*....................................................................*/
int getResultNameFromI(const int resultI, char *resultName){
  int status=0;

  switch(resultI){
    case RESULT_density:
      resultName = "density";
      break;
    case RESULT_temperature:
      resultName = "temperature";
      break;
    case RESULT_tdust:
      resultName = "tdust";
      break;
    case RESULT_abundance:
      resultName = "abundance";
      break;
    case RESULT_doppler:
      resultName = "doppler";
      break;
    case RESULT_velocity:
      resultName = "velocity";
      break;
    case RESULT_bmag:
      resultName = "bmag";
      break;
    default:
      resultName = '\0'; /* resultI not recognized. */
  }

  return status;
}

/*....................................................................*/
void
_scalarConst_getIndex(char *funcParName, int *parIndex){
  /* **NOTE** you should not change this ordering without doing the same for the argument list of the respective function. */

  if(strcmp(funcParName, "val")==0)
    *parIndex = 0;
  else
    *parIndex = -1; /* parameter name not recognized. */
}

/*....................................................................*/
double
_scalarConst(const double val){
  return val;
}

/*....................................................................*/
void
_scalarPowerR_getIndex(char *funcParName, int *parIndex){
  /* **NOTE** you should not change this ordering without doing the same for the argument list of the respective function. */

  if(     strcmp(funcParName, "lowerR")==0)
    *parIndex = 0;
  else if(strcmp(funcParName, "factor")==0)
    *parIndex = 1;
  else if(strcmp(funcParName, "exponent")==0)
    *parIndex = 2;
  else if(strcmp(funcParName, "offset")==0)
    *parIndex = 3;
  else
    *parIndex = -1; /* parameter name not recognized. */
}

/*....................................................................*/
double
_scalarPowerR(const double lowerR, const double factor, const double exponent\
  , const double offset, double r, const double rVec[ML_NUM_DIMS]){

  double value;

  if(r<0.0)
    r = sqrt(rVec[0]*rVec[0] + rVec[1]*rVec[1] + rVec[2]*rVec[2]);

  if(r >= lowerR){
    value = factor*pow(r, exponent) + offset;
  }else{
    value = factor*pow(lowerR, exponent) + offset;
  }

  return value; 
}

/*....................................................................*/
void
_scalarPowerRExpZ_getIndex(char *funcParName, int *parIndex){
  /* **NOTE** you should not change this ordering without doing the same for the argument list of the respective function. */

  if(     strcmp(funcParName, "lowerR")==0)
    *parIndex = 0;
  else if(strcmp(funcParName, "factR")==0)
    *parIndex = 1;
  else if(strcmp(funcParName, "expR")==0)
    *parIndex = 2;
  else if(strcmp(funcParName, "offsetR")==0)
    *parIndex = 3;
  else if(strcmp(funcParName, "factZ")==0)
    *parIndex = 4;
  else if(strcmp(funcParName, "scaleZ")==0)
    *parIndex = 5;
  else if(strcmp(funcParName, "offsetZ")==0)
    *parIndex = 6;
  else
    *parIndex = -1; /* parameter name not recognized. */
}

/*....................................................................*/
double
_scalarPowerRExpZ(const double lowerR, const double factR, const double expR\
  , const double offsetR, const double factZ, const double scaleZ\
  , const double offsetZ, const double rVec[ML_NUM_DIMS]){

  double r,valR,valZ;

  r = sqrt(rVec[0]*rVec[0] + rVec[1]*rVec[1] + rVec[2]*rVec[2]);

  if(r >= lowerR){
    valR = factR*pow(r, expR) + offsetR;
  }else{
    valR = factR*pow(lowerR, expR) + offsetR;
  }

  valZ = factZ * exp(-abs(rVec[2])/scaleZ) + offsetZ;

  return valR * valZ; 
}

/*....................................................................*/
void
_scalarPowerRTheta_getIndex(char *funcParName, int *parIndex){
  /* **NOTE** you should not change this ordering without doing the same for the argument list of the respective function. */

  if(     strcmp(funcParName, "lowerR")==0)
    *parIndex = 0;
  else if(strcmp(funcParName, "factR")==0)
    *parIndex = 1;
  else if(strcmp(funcParName, "expR")==0)
    *parIndex = 2;
  else if(strcmp(funcParName, "offsetR")==0)
    *parIndex = 3;
  else if(strcmp(funcParName, "lowerTheta")==0)
    *parIndex = 4;
  else if(strcmp(funcParName, "factTheta")==0)
    *parIndex = 5;
  else if(strcmp(funcParName, "expTheta")==0)
    *parIndex = 6;
  else if(strcmp(funcParName, "offsetTheta")==0)
    *parIndex = 7;
  else
    *parIndex = -1; /* parameter name not recognized. */
}

/*....................................................................*/
double
_scalarPowerRTheta(const double lowerR, const double factR, const double expR\
  , const double offsetR, const double lowerTheta, const double factTheta\
  , const double expTheta, const double offsetTheta, double r\
  , const double rVec[ML_NUM_DIMS]){

  double valR,theta,valTheta;

  if(r<0.0)
    r = sqrt(rVec[0]*rVec[0] + rVec[1]*rVec[1] + rVec[2]*rVec[2]);

  theta = acos(abs(rVec[2])/r);

  if(r >= lowerR){
    valR = factR*pow(r, expR) + offsetR;
  }else{
    valR = factR*pow(lowerR, expR) + offsetR;
  }

  if(theta >= lowerTheta){
    valTheta = factTheta*pow(theta,expTheta) + offsetTheta;
  }else{
    valTheta = factTheta*pow(lowerTheta,expTheta) + offsetTheta;
  }

  return valR*valTheta; 
}

/*....................................................................*/
void
_scalarPowerRZ_getIndex(char *funcParName, int *parIndex){
  /* **NOTE** you should not change this ordering without doing the same for the argument list of the respective function. */

  if(     strcmp(funcParName, "lowerR")==0)
    *parIndex = 0;
  else if(strcmp(funcParName, "factR")==0)
    *parIndex = 1;
  else if(strcmp(funcParName, "expR")==0)
    *parIndex = 2;
  else if(strcmp(funcParName, "offsetR")==0)
    *parIndex = 3;
  else if(strcmp(funcParName, "lowerZ")==0)
    *parIndex = 4;
  else if(strcmp(funcParName, "factZ")==0)
    *parIndex = 5;
  else if(strcmp(funcParName, "expZ")==0)
    *parIndex = 6;
  else if(strcmp(funcParName, "offsetZ")==0)
    *parIndex = 7;
  else
    *parIndex = -1; /* parameter name not recognized. */
}

/*....................................................................*/
double
_scalarPowerRZ(const double lowerR, const double factR, const double expR\
  , const double offsetR, const double lowerZ, const double factZ\
  , const double expZ, const double offsetZ, const double rVec[ML_NUM_DIMS]){

  double r,valR,valZ;

  r = sqrt(rVec[0]*rVec[0] + rVec[1]*rVec[1] + rVec[2]*rVec[2]);

  if(r >= lowerR){
    valR = factR*pow(r, expR) + offsetR;
  }else{
    valR = factR*pow(lowerR, expR) + offsetR;
  }

  if(abs(rVec[2]) >= abs(lowerZ)){
    valZ = factZ*pow(abs(rVec[2]),expZ) + offsetZ;
  }else{
    valZ = factZ*pow(abs(lowerZ),expZ) + offsetZ;
  }

  return valR*valZ;
}

/*....................................................................*/
void
_vectorUnscaledConst(double vecResult[ML_NUM_DIMS]){
  int i;

  for(i=0;i<ML_NUM_DIMS;i++)
    vecResult[i] = 1.0;
}

/*....................................................................*/
void
_vectorUnscaledRadial(const double rVec[ML_NUM_DIMS], const double r, double vecResult[ML_NUM_DIMS]){
  int i;

  if(r>0.0){
    for(i=0;i<ML_NUM_DIMS;i++)
      vecResult[i] = rVec[i]*(1.0/r);
  }else{
    for(i=0;i<ML_NUM_DIMS;i++)
      vecResult[i] = 0.0;
  }
}

///*....................................................................*/
//void
//_vectorUnscaledDipole(const double rVec[ML_NUM_DIMS], double vecResult[ML_NUM_DIMS]){
//}

/*....................................................................*/
void
_vectorUnscaledToroidal(const double rVec[ML_NUM_DIMS], double rho\
  , double vecResult[ML_NUM_DIMS]){

  if(rho<0.0)
    rho = sqrt(rVec[0]*rVec[0] + rVec[1]*rVec[1]);

  vecResult[0] = -rVec[1]/rho;
  vecResult[1] =  rVec[0]/rho;
  vecResult[2] = 0.0;
}

/*....................................................................*/
void
_vectorConstR_getIndex(char *funcParName, int *parIndex){
  /* **NOTE** you should not change this ordering without doing the same for the argument list of the respective function. */

  if(     strcmp(funcParName, "val")==0)
    *parIndex = 0;
  else
    *parIndex = -1; /* parameter name not recognized. */
}

/*....................................................................*/
void
_vectorConstR(const double val, const double rVec[ML_NUM_DIMS], double vecResult[ML_NUM_DIMS]){
  double r;
  int i;

  r = sqrt(rVec[0]*rVec[0] + rVec[1]*rVec[1] + rVec[2]*rVec[2]);

  _vectorUnscaledRadial(rVec, r, vecResult);

  for(i=0;i<ML_NUM_DIMS;i++)
    vecResult[i] *= val;
}

/*....................................................................*/
void
_vectorConstXYZ_getIndex(char *funcParName, int *parIndex){
  /* **NOTE** you should not change this ordering without doing the same for the argument list of the respective function. */

  if(     strcmp(funcParName, "valX")==0)
    *parIndex = 0;
  else if(strcmp(funcParName, "valY")==0)
    *parIndex = 1;
  else if(strcmp(funcParName, "valZ")==0)
    *parIndex = 2;
  else
    *parIndex = -1; /* parameter name not recognized. */
}

/*....................................................................*/
void
_vectorConstXYZ(const double valX, const double valY, const double valZ\
  , double vecResult[ML_NUM_DIMS]){

  vecResult[0] = valX;
  vecResult[1] = valY;
  vecResult[2] = valZ;
}

/*....................................................................*/
void
_vectorDipole_getIndex(char *funcParName, int *parIndex){
  /* **NOTE** you should not change this ordering without doing the same for the argument list of the respective function. */

  if(     strcmp(funcParName, "lowerR")==0)
    *parIndex = 0;
  else if(strcmp(funcParName, "mDipole")==0)
    *parIndex = 1;
  else
    *parIndex = -1; /* parameter name not recognized. */
}

/*....................................................................*/
void
_vectorDipole(const double lowerR, const double mDipole\
  , const double rVec[ML_NUM_DIMS], double vecResult[ML_NUM_DIMS]){
  /* Magnitude of mDipole is actually mu_0 * dipole_moment / 4 / pi. The dipole moment is assumed to be pointing in the +ve Z direction. */

  double rSquared,oneOnRSquared,oneOnRCubed;

  rSquared = rVec[0]*rVec[0] + rVec[1]*rVec[1] + rVec[2]*rVec[2];

  if(rSquared >= lowerR*lowerR){
    oneOnRSquared = 1.0/rSquared;
    oneOnRCubed = oneOnRSquared/sqrt(rSquared);
  }else{
    oneOnRSquared = 1.0/(lowerR*lowerR);
    oneOnRCubed = oneOnRSquared/lowerR;
  }

  vecResult[0] = 3.0*mDipole*oneOnRCubed* rVec[2]*rVec[0]*oneOnRSquared;
  vecResult[1] = 3.0*mDipole*oneOnRCubed* rVec[2]*rVec[1]*oneOnRSquared;
  vecResult[2] = 3.0*mDipole*oneOnRCubed*(rVec[2]*rVec[2]*oneOnRSquared - (1.0/3.0));
}

/*....................................................................*/
void
_vectorRadialPowerR_getIndex(char *funcParName, int *parIndex){
  /* **NOTE** you should not change this ordering without doing the same for the argument list of the respective function. */

  if(     strcmp(funcParName, "lowerR")==0)
    *parIndex = 0;
  else if(strcmp(funcParName, "factor")==0)
    *parIndex = 1;
  else if(strcmp(funcParName, "exponent")==0)
    *parIndex = 2;
  else if(strcmp(funcParName, "offset")==0)
    *parIndex = 3;
  else
    *parIndex = -1; /* parameter name not recognized. */
}

/*....................................................................*/
void
_vectorRadialPowerR(const double lowerR, const double factor\
  , const double exponent, const double offset, const double rVec[ML_NUM_DIMS]\
  , double vecResult[ML_NUM_DIMS]){

  double r,scale;
  int i;

  r = sqrt(rVec[0]*rVec[0] + rVec[1]*rVec[1] + rVec[2]*rVec[2]);
  scale = _scalarPowerR(lowerR, factor, exponent, offset, r, rVec);

  _vectorUnscaledRadial(rVec, r, vecResult);

  for(i=0;i<ML_NUM_DIMS;i++)
    vecResult[i] *= scale;
}

/*....................................................................*/
void
_vectorRadialPowerRTheta_getIndex(char *funcParName, int *parIndex){
  /* **NOTE** you should not change this ordering without doing the same for the argument list of the respective function. */

  if(     strcmp(funcParName, "lowerR")==0)
    *parIndex = 0;
  else if(strcmp(funcParName, "factR")==0)
    *parIndex = 1;
  else if(strcmp(funcParName, "expR")==0)
    *parIndex = 2;
  else if(strcmp(funcParName, "offsetR")==0)
    *parIndex = 3;
  else if(strcmp(funcParName, "lowerTheta")==0)
    *parIndex = 4;
  else if(strcmp(funcParName, "factTheta")==0)
    *parIndex = 5;
  else if(strcmp(funcParName, "expTheta")==0)
    *parIndex = 6;
  else if(strcmp(funcParName, "offsetTheta")==0)
    *parIndex = 7;
  else
    *parIndex = -1; /* parameter name not recognized. */
}

/*....................................................................*/
void
_vectorRadialPowerRTheta(const double lowerR, const double factR, const double expR\
  , const double offsetR, const double lowerTheta, const double factTheta\
  , const double expTheta, const double offsetTheta, const double rVec[ML_NUM_DIMS]\
  , double vecResult[ML_NUM_DIMS]){

  double r,scale;
  int i;

  r = sqrt(rVec[0]*rVec[0] + rVec[1]*rVec[1] + rVec[2]*rVec[2]);
  scale = _scalarPowerRTheta(lowerR, factR, expR, offsetR, lowerTheta, factTheta, expTheta, offsetTheta, r, rVec);

  _vectorUnscaledRadial(rVec, r, vecResult);

  for(i=0;i<ML_NUM_DIMS;i++)
    vecResult[i] *= scale;
}

/*....................................................................*/
void
_vectorToroidalPowerR_getIndex(char *funcParName, int *parIndex){
  /* **NOTE** you should not change this ordering without doing the same for the argument list of the respective function. */

  if(     strcmp(funcParName, "lowerRho")==0)
    *parIndex = 0;
  else if(strcmp(funcParName, "factor")==0)
    *parIndex = 1;
  else if(strcmp(funcParName, "exponent")==0)
    *parIndex = 2;
  else if(strcmp(funcParName, "offset")==0)
    *parIndex = 3;
  else
    *parIndex = -1; /* parameter name not recognized. */
}

/*....................................................................*/
void
_vectorToroidalPowerR(const double lowerRho, const double factor\
  , const double exponent, const double offset, const double rVec[ML_NUM_DIMS]\
  , double vecResult[ML_NUM_DIMS]){

  double rho,scale;
  int i;

  rho = sqrt(rVec[0]*rVec[0] + rVec[1]*rVec[1]);
  scale = _scalarPowerR(lowerRho, factor, exponent, offset, rho, rVec);

  _vectorUnscaledToroidal(rVec, rho, vecResult);

  for(i=0;i<ML_NUM_DIMS;i++)
    vecResult[i] *= scale;
}

/*....................................................................*/
int
getFuncParIndexFromName(const int funcI, const int funcTypeI, char *funcParName, int *parIndex){
  int status=0;

  switch(funcTypeI){
  case FUNC_scalar:
    switch(funcI){
    case FUNC_scalarConst:
      _scalarConst_getIndex(funcParName, parIndex);
      break;
    case FUNC_scalarPowerR:
      _scalarPowerR_getIndex(funcParName, parIndex);
      break;
    case FUNC_scalarPowerRExpZ:
      _scalarPowerRExpZ_getIndex(funcParName, parIndex);
      break;
    case FUNC_scalarPowerRTheta:
      _scalarPowerRTheta_getIndex(funcParName, parIndex);
      break;
    case FUNC_scalarPowerRZ:
      _scalarPowerRZ_getIndex(funcParName, parIndex);
      break;
    default:
      status = 1;
      *parIndex = -1; /* nonsensical value. */
    }

    break;

  case FUNC_vector:
    switch(funcI){
    case FUNC_vectorConstR:
      _vectorConstR_getIndex(funcParName, parIndex);
      break;

    case FUNC_vectorConstXYZ:
      _vectorConstXYZ_getIndex(funcParName, parIndex);
      break;

    case FUNC_vectorDipole:
      _vectorDipole_getIndex(funcParName, parIndex);
      break;

    case FUNC_vectorRadialPowerR:
      _vectorRadialPowerR_getIndex(funcParName, parIndex);
      break;

    case FUNC_vectorRadialPowerRTheta:
      _vectorRadialPowerRTheta_getIndex(funcParName, parIndex);
      break;

    case FUNC_vectorToroidalPowerR:
      _vectorToroidalPowerR_getIndex(funcParName, parIndex);
      break;

    default:
      status = 2;
      *parIndex = -1; /* nonsensical value. */
    } 

    break;
  default:
    status = 3;
    *parIndex = -1; /* nonsensical value. */
  }

  if(*parIndex<0 && status==0) /* signals that the *_getIndex() function which was called could not identify the parameter name. */
    status = 4;

return status;
}


/*....................................................................*/
double
scalarFunctionSwitch(const int funcI, double *fpars, const double rVec[ML_NUM_DIMS]){
  double value;

  switch(funcI){
  case FUNC_scalarConst:
    value = _scalarConst(fpars[0]);
    break;

  case FUNC_scalarPowerR:
    value = _scalarPowerR(fpars[0], fpars[1], fpars[2], fpars[3], -1.0, rVec);
    break;

  case FUNC_scalarPowerRExpZ:
    value = _scalarPowerRExpZ(fpars[0], fpars[1], fpars[2], fpars[3], fpars[4], fpars[5], fpars[6], rVec);
    break;

  case FUNC_scalarPowerRTheta:
    value = _scalarPowerRTheta(fpars[0], fpars[1], fpars[2], fpars[3], fpars[4], fpars[5], fpars[6], fpars[7], -1.0, rVec);
    break;

  case FUNC_scalarPowerRZ:
    value = _scalarPowerRZ(fpars[0], fpars[1], fpars[2], fpars[3], fpars[4], fpars[5], fpars[6], fpars[7], rVec);
    break;

  default:
    value = 0.0; //*** error?
  } 

  return value;
}

/*....................................................................*/
void
vectorFunctionSwitch(const int funcI, double *fpars\
  , const double rVec[ML_NUM_DIMS], double vecResult[ML_NUM_DIMS]){

  int i;

  switch(funcI){
    case FUNC_vectorConstR:
      _vectorConstR(fpars[0], rVec, vecResult);
      break;

    case FUNC_vectorConstXYZ:
      _vectorConstXYZ(fpars[0], fpars[1], fpars[2], vecResult);
      break;

    case FUNC_vectorDipole:
      _vectorDipole(fpars[0], fpars[1], rVec, vecResult);
      break;

    case FUNC_vectorRadialPowerR:
      _vectorRadialPowerR(fpars[0], fpars[1], fpars[2], fpars[3], rVec\
        , vecResult);
      break;

    case FUNC_vectorRadialPowerRTheta:
      _vectorRadialPowerRTheta(fpars[0], fpars[1], fpars[2], fpars[3]\
        , fpars[4], fpars[5], fpars[6], fpars[7], rVec, vecResult);
      break;

    case FUNC_vectorToroidalPowerR:
      _vectorToroidalPowerR(fpars[0], fpars[1], fpars[2], fpars[3]\
        , rVec, vecResult);
      break;

    default: //********** should be an error?
      for(i=0;i<ML_NUM_DIMS;i++)
        vecResult[i] = 0.0;
  } 
}

