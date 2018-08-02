/*
 *  aux.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#ifdef IS_PYTHON
#include "lime.h" // actually just for Python.h
#endif

#include "constants.h" // for STR_LEN_0 etc
#include "messages.h" // for warning(), bail_out() etc
#include "aux.h"

/*....................................................................*/
void
reportInfAtOrigin(const double value, const char *funcName){
  char message[STR_LEN_0];

  if((isinf(value) || isnan(value)) && !silent){
    snprintf(message, STR_LEN_0, "You have a singularity at the origin of your %s() function.", funcName);
    warning(message);
  }
}

/*....................................................................*/
void
reportInfsAtOrigin(const int numElements, const double *values, const char *funcName){
  int i;
  char message[STR_LEN_0];

  if(numElements<=0){
    if((isinf(*values) || isnan(*values)) && !silent){
      snprintf(message, STR_LEN_0, "You have a singularity at the origin of your %s() function.", funcName);
      warning(message);
    }
  }else{
    for(i=0;i<numElements;i++){
      if((isinf(values[i]) || isnan(values[i])) && !silent){
        snprintf(message, STR_LEN_0, "You have a singularity at the origin in return %d of your %s() function.", i, funcName);
        warning(message);
      }
    }
  }
}

/*....................................................................*/
void sigintHandler(int sigI){
#ifdef IS_PYTHON
  Py_Finalize();
#endif

//*** write output file?

exit(1);
}

/*....................................................................*/
void checkFgets(char *fgetsResult, char *message){
  char string[STR_LEN_0];

  if(fgetsResult==NULL){
    if(!silent){
      snprintf(string, STR_LEN_0, "fgets() failed to read %s", message);
      bail_out(string);
    }
exit(1);
  }
}

/*....................................................................*/
void checkFscanf(const int fscanfResult, const int expectedNum, char *message){
  char string[STR_LEN_0];

  if(fscanfResult!=expectedNum){
    if(!silent){
      snprintf(string, STR_LEN_0, "fscanf() failed to read %s - read %d bytes when %d expected.", message, fscanfResult, expectedNum);
      bail_out(string);
    }
exit(1);
  }
}

/*....................................................................*/
void checkFread(const size_t freadResult, const size_t expectedNum, char *message){
  char string[STR_LEN_0];

  if(freadResult!=expectedNum){
    if(!silent){
      snprintf(string, STR_LEN_0, "fread() failed to read %s. Expected %d got %d", message, (int)expectedNum, (int)freadResult);
      bail_out(string);
    }
exit(1);
  }
}

/*....................................................................*/
void checkFwrite(const size_t fwriteResult, const size_t expectedNum, char *message){
  char string[STR_LEN_0];

  if(fwriteResult!=expectedNum){
    if(!silent){
      snprintf(string, STR_LEN_0, "fwrite() failed to write %s. Expected %d got %d", message, (int)expectedNum, (int)fwriteResult);
      bail_out(string);
    }
exit(1);
  }
}

/*....................................................................*/
double
gaussline(const double v, const double oneOnSigma){
  double val;
  val = v*v*oneOnSigma*oneOnSigma;
#ifdef FASTEXP
  return FastExp(val);
#else
  return exp(-val);
#endif
}

/*....................................................................*/
double
dotProduct3D(const double *vA, const double *vB){
  return vA[0]*vB[0] + vA[1]*vB[1] + vA[2]*vB[2];
}

/*....................................................................*/
void
copyInparStr(const char *inStr, char **outStr){
  if(inStr==NULL || strlen(inStr)<=0 || strlen(inStr)>STR_LEN_0){
    *outStr = NULL;
  }else{
    *outStr = malloc(sizeof(**outStr)*(STR_LEN_0+1));
    strcpy(*outStr, inStr);
  }
}

/*....................................................................*/
_Bool
charPtrIsNullOrEmpty(const char *inStr){
  if(inStr==NULL || strlen(inStr)<=0)
    return 1;
  else
    return 0;
}

/*....................................................................*/
_Bool allBitsSet(const int flags, const int mask){
  /* Returns true only if all the masked bits of flags are set. */

  if(~flags & mask)
    return 0;
  else
    return 1;
}

/*....................................................................*/
_Bool anyBitSet(const int flags, const int mask){
  /* Returns true if any of the masked bits of flags are set. */

  if(flags & mask)
    return 1;
  else
    return 0;
}

/*....................................................................*/
_Bool bitIsSet(const int flags, const int bitI){
  /* Returns true if the designated bit of flags is set. */

  if(flags & (1 << bitI))
    return 1;
  else
    return 0;
}

/*....................................................................*/
_Bool onlyBitsSet(const int flags, const int mask){
  /* Returns true if flags has no bits set apart from those which are true in mask. */

  if(flags & ~mask)
    return 0;
  else
    return 1;
}

