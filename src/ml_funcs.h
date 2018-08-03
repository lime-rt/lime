/*
 *  ml_funcs.h
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
************ change sprintfs to snprintfs
 */

#ifndef ML_FUNCS_H
#define ML_FUNCS_H

#include <Python.h>
#include "ml_types.h"

#define FUNC_scalarConst                0
#define FUNC_scalarPowerR               1
#define FUNC_scalarPowerRExpZ           2
#define FUNC_scalarPowerRTheta          3
#define FUNC_scalarPowerRZ              4
#define FUNC_NUM_SCALAR     5

#define FUNC_vectorConstR               0
#define FUNC_vectorConstXYZ             1
#define FUNC_vectorDipole               2
#define FUNC_vectorRadialPowerR         3
#define FUNC_vectorRadialPowerRTheta    4
#define FUNC_vectorToroidalPowerR       5
#define FUNC_NUM_VECTOR     6

#define FUNC_SCALAR_NPARS_LIST   {1,4,7,8,8}
#define FUNC_VECTOR_NPARS_LIST   {1,3,2,4,8,4}

#define FUNC_scalar  0
#define FUNC_vector  1

extern int funcTypeIs[NUM_MODELS][NUM_RESULTS];
extern int numScalarFuncPars[FUNC_NUM_SCALAR];
extern int numVectorFuncPars[FUNC_NUM_VECTOR];

void	freeFuncsPars(void);
void	setDefaultFuncStuffs(void);
int	getFuncTypeIFromName(char *funcTypeName);
int	getFuncTypeNameFromI(const int funcTypeI, char *funcTypeName);
int	getScalarFuncIFromName(char *funcName);
int	getScalarFuncNameFromI(const int funcI, char *funcName);
int	getVectorFuncIFromName(char *funcName);
int	getVectorFuncNameFromI(const int funcI, char *funcName);
int	getResultIFromName(char *resultName);
int	getResultNameFromI(const int resultI, char *resultName);
int	getFuncParIndexFromName(const int funcI, const int funcTypeI, char *funcParName, int *parIndex);
double	scalarFunctionSwitch(const int funcI, double *fpars, const double*);
void	vectorFunctionSwitch(const int funcI, double *fpars, const double*, double*);

#endif /* ML_FUNCS_H */


