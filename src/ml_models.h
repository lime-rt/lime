/*
 *  ml_models.h
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#ifndef ML_MODELS_H
#define ML_MODELS_H

#include <Python.h>
#include "ml_types.h"

#define MAX_LEN_PAR_NAME      30
#define MAX_LEN_PAR_TYPE      20

#define MODEL_None        0
#define MODEL_BoEb56      1	/* BonnorEbert56 */
#define MODEL_CG97        2	/* CG97 */
#define MODEL_DDN01       3	/* DDN01 */
#define MODEL_LiSh96      4	/* LiShu96 */
#define MODEL_Ma88        5	/* Mamon88 */
#define MODEL_Me09        6	/* Mendoza09 */
#define MODEL_Shu77       7	/* Shu77 */
#define MODEL_Ul76        8	/* Ulrich76 */
#define MODEL_Al03        9	/* allen03a */


typedef struct{
  char name[MAX_LEN_PAR_NAME+1];
  char dtype[MAX_LEN_PAR_TYPE+1];
  int index;
} modelParamType;

extern int currentModelI,*modelIntPars,numModelParams;
extern double *modelDblPars;
extern char *modelStrPar;
extern modelParamType *modelParams;
extern int copyTemp; /* <0 means tgas should be set to equal tdust, >0 means the other way, ==0 means no action. */


int	getParamI(char *idStr, int *index);
int	getModelIFromName(char *modelName);
int	getModelNameFromI(const int modelI, char *modelName);
int	finalizeModelConfig(const int modelI);

#endif /* ML_MODELS_H */


