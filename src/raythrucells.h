/*
 *  raythrucells.h
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2016 The LIME development team
 *
 */

#ifndef RAYTHRUCELLS_H
#define RAYTHRUCELLS_H

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#include "dims.h" /* Should define the macro N_DIMS which specifies the number of spatial dimensions. Note that this is a _maximum_ value, the actual value 'numDims' which is passed via the function interfaces may be <= N_DIMS. */

/* Error codes:
*/
#define RTC_ERR_SVD_FAIL	0
#define RTC_ERR_NON_SPAN	1
#define RTC_ERR_LU_DECOMP_FAIL	2
#define RTC_ERR_LU_SOLVE_FAIL	3
#define RTC_ERR_TOO_MANY_ENTRY	4
#define RTC_ERR_BUG		5
#define RTC_ERR_OLD_NOT_FOUND	6


/* In the following comments, N is short for N_DIMS, the number of spatial dimensions. */

/* NOTE that it is assumed that vertx[i] is opposite the face that abuts with neigh[i] for all i.
*/ 
struct simplex {
  unsigned long id;
  unsigned long vertx[N_DIMS+1];
  double centre[N_DIMS];
  struct simplex *neigh[N_DIMS+1]; /* An entry ==NULL flags an external face. */
};

/* This struct is meant to record all relevant information about the intersection between a ray (defined by a direction unit vector 'dir' and a starting position 'r') and a face of a simplex.
*/
typedef struct {
  int fi;
  /* The index (in the range {0...N}) of the face (and thus of the opposite vertex, i.e. the one 'missing' from the bary[] list of this face).
  */
  int orientation;
  /* >0 means the ray exits, <0 means it enters, ==0 means the face is parallel to ray.
  */
  double bary[N_DIMS], dist, collPar;
  /* 'dist' is defined via r_int = r + dist*dir. 'collPar' is a measure of how close to any edge of the face r_int lies.
  */
} intersectType;

typedef struct {
  double r[N_DIMS][N_DIMS], simplexCentre[N_DIMS];
  /* 'r' is a list of the the N vertices of the face, each of which has N cartesian components. 'simplexCentre' is a convenience pointer which gives the location of the geometric centre of the simplex. */
} faceType;

typedef struct {
  double axes[N_DIMS-1][N_DIMS], r[N_DIMS][N_DIMS-1], origin[N_DIMS];
  /* 'r' expresses the location of the N vertices of a simplicial polytope face in N-space, in terms of components along the N-1 orthogonal axes in the sub-plane of the face. Thus you should malloc r as r[N][N-1]. */
} facePlusBasisType;

typedef struct{
  faceType *faces,*(*facePtrs[N_DIMS+1]);
} faceListType;

faceType *extractFace(const int numDims, double *vertexCoords, struct simplex *dc\
  , const unsigned long dci, const int fi);
int	followRayThroughCells(const int numDims, double *x, double *dir\
  , double *vertexCoords, struct simplex *dc, const unsigned long numCells\
  , const double epsilon, faceType **facePtrs[N_DIMS+1], intersectType *entryIntcpt, unsigned long **chainOfCellIds\
  , intersectType **cellExitIntcpts, int *lenChainPtrs);

#endif /* RAYTHRUCELLS_H */

