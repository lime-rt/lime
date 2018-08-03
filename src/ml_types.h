/*
 *  ml_types.h
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#ifndef ML_TYPES_H
#define ML_TYPES_H

#include <Python.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_erf.h>

#include "local_err.h"

#define MAX_LEN_MODEL_NAME    30
#define MAX_LEN_RESULT_NAME   20

#define ML_NUM_DIMS    3

#define NUM_MODELS         10
#define NUM_RESULTS         7

#define RESULT_density       0
#define RESULT_temperature   1
#define RESULT_abundance     2
#define RESULT_doppler       3
#define RESULT_velocity      4
#define RESULT_bmag          5
#define RESULT_tdust         6

extern int funcIs[NUM_MODELS][NUM_RESULTS];
extern double *funcPars[NUM_MODELS][NUM_RESULTS];

errType	extractParams(PyObject *pModelObj);
errType	extractFuncs(PyObject *pModelObj);
errType	getModelI(PyObject *pModelObj, int *modelI);

#endif /* ML_TYPES_H */


