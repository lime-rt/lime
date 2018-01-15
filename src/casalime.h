/*
 *  casalime.h
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2017 The LIME development team
 *
 */

#ifndef CASALIME_H
#define CASALIME_H

#include "lime_config.h" /* for STR_LEN_0 */

#define STATUS_STR_LEN		STR_LEN_0+11

typedef struct {
  float progressGridBuilding;
//  int statusGridBuilding;
//Progress of LIME's grid building (0 to 1)
  float progressGridSmoothing;
//  int statusGridSmoothing;
//Progress of LIME's grid smoothing (0 to 1)
  int statusGrid;
//1 if LIME's grid is complete, 0 else
  float progressConvergence;
//Progress of LIME's convergence (0 to 1)
  int numberIterations;
  double minsnr;
  double median;
  float progressPhotonPropagation;
//  int statusPhotonPropagation;
//Progress of LIME's photon propagation (0 to 1)
  float progressRayTracing;
//Progress of LIME's ray tracing (0 to 1)
  int statusRayTracing;
//1 if LIME's ray tracing is complete, 0 else
  int statusGlobal;
//1 if LIME's run is complete, 0 else
  int error;
//1 if an internal LIME error occured, 0 else. See messaeg for an error message
  char message[STATUS_STR_LEN];
} statusType;

extern statusType statusObj;
extern _Bool statusObjInitialized;

#endif /* CASALIME_H */

