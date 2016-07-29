/*
 *  velospline.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2016 The LIME development team
 *
 */

#include "lime.h"

void calcInterpCoeffs(configInfo *par, struct grid *gp){
  /*
The velocity v(d) (or rather the scalar component of vector velocity in the direction of the edge) at a distance d along the line (or 'edge' in triangulation jargon) between a given grid point and its neighbour is approximated in LIME by the polynomial expression

	       __N-1
	       \
	v(d) ~  >     a_j*s^j
	       /_j=0

where
	     d
	s = --- - 0.5.
	     ds

The coefficients a_j are calculated in the present function for each grid point (which uses twice the memory needed). This is done by sampling velocity at N points along the edge with at distances evenly spread between 0 and ds, then solving the system of equations to obtain the interpolation coefficients.

A Chebyshev interpolation would probably be preferred, but we wil leave that for future generations.

Note that, given coefficients calculated for edge AB, the present definition allows easy conversion from the values calculated for grid point A to those for grid point B: coefficient j for B is just (-1)^j times the respective coefficient for A.

The number of coefficients N is currently hard-wired at 5. 
  */

  int i,k,di,ri,ci,dummy;
  const int nCoeffs = NUM_VEL_COEFFS;
  double vel[DIM],x[DIM],dFrac[nCoeffs],velComp,sToPower;
  gsl_matrix *matrix = gsl_matrix_alloc(nCoeffs,nCoeffs);
  gsl_vector *aCoeffs = gsl_vector_alloc(nCoeffs);
  gsl_vector *velVector = gsl_vector_alloc(nCoeffs);
  gsl_permutation *p = gsl_permutation_alloc(nCoeffs);

  for(ri=0;ri<nCoeffs;ri++)
    dFrac[ri] = (double)ri/(double)(nCoeffs-1); /* Should range between 0 and 1 inclusive. */

  for(i=0;i<par->pIntensity;i++){
    gp[i].a0=malloc(gp[i].numNeigh*sizeof(double));
    gp[i].a1=malloc(gp[i].numNeigh*sizeof(double));
    gp[i].a2=malloc(gp[i].numNeigh*sizeof(double));
    gp[i].a3=malloc(gp[i].numNeigh*sizeof(double));
    gp[i].a4=malloc(gp[i].numNeigh*sizeof(double));

    for(k=0;k<gp[i].numNeigh;k++){
      for(ri=0;ri<nCoeffs;ri++){
        for(di=0;di<DIM;di++)
          x[di] = gp[i].x[di] + dFrac[ri]*gp[i].dir[k].xn[di]*gp[i].ds[k];
        velocity(x[0],x[1],x[2],vel);
        velComp = veloproject(gp[i].dir[k].xn,vel); /* Component of velocity in the direction of the neighbour point. */

        sToPower = 1.0;
        for(ci=0;ci<nCoeffs;ci++){
          gsl_matrix_set(matrix,ri,ci,sToPower);
          sToPower *= (dFrac[ri] - 0.5);
        }
        gsl_vector_set(velVector,ri,velComp);
      }

      /*
We obtain the N coefficients by solving the matrix equation

	Ma = v
where
	M_{i,j} = s_i^j
for
	        i
	s_i = ----- - 0.5,
	       N-1

a_i is the ith coefficient, and v_i is the ith sample of the velocity component along the line between point A and point B.
      */
      gsl_linalg_LU_decomp (matrix, p, &dummy);
      gsl_linalg_LU_solve (matrix, p, velVector, aCoeffs);
      gp[i].a0[k]=gsl_vector_get(aCoeffs,0);
      gp[i].a1[k]=gsl_vector_get(aCoeffs,1);
      gp[i].a2[k]=gsl_vector_get(aCoeffs,2);
      gp[i].a3[k]=gsl_vector_get(aCoeffs,3);
      gp[i].a4[k]=gsl_vector_get(aCoeffs,4);
    }
  }

  for(i=par->pIntensity;i<par->ncell;i++){
    gp[i].a0=malloc(gp[i].numNeigh*sizeof(double));
    gp[i].a1=malloc(gp[i].numNeigh*sizeof(double));
    gp[i].a2=malloc(gp[i].numNeigh*sizeof(double));
    gp[i].a3=malloc(gp[i].numNeigh*sizeof(double));
    gp[i].a4=malloc(gp[i].numNeigh*sizeof(double));
    for(k=0;k<gp[i].numNeigh;k++){
      gp[i].a0[k]=0.;
      gp[i].a1[k]=0.;
      gp[i].a2[k]=0.;
      gp[i].a3[k]=0.;
      gp[i].a4[k]=0.;
    }
  }

  gsl_permutation_free (p);
  gsl_vector_free (aCoeffs);
  gsl_vector_free (velVector);
  gsl_matrix_free (matrix);
}

void calcInterpCoeffs_lin(configInfo *par, struct grid *g){
  /*
This is the same as calcInterpCoeffs() except only 2 coefficients are calculated, allowing a linear interpolation of the velocity component along the line between 2 grid points. This can be used if the velocities have been calculated elsewhere at the grid points themselves but the velocity function is not available for further sampling.
  */
  int i,k;
  double v[2];

  for(i=0;i<par->pIntensity;i++){
    g[i].a0=malloc(g[i].numNeigh*sizeof(double));
    g[i].a1=malloc(g[i].numNeigh*sizeof(double));
    for(k=0;k<g[i].numNeigh;k++){
      v[0]=veloproject(g[i].dir[k].xn,g[i].vel);
      v[1]=veloproject(g[i].dir[k].xn,g[i].neigh[k]->vel);
      g[i].a1[k] = v[1] - v[0];
      g[i].a0[k] = (v[0] + v[1])/2.0;
    }
  }

  for(i=par->pIntensity;i<par->ncell;i++){
    g[i].a0=malloc(g[i].numNeigh*sizeof(double));
    g[i].a1=malloc(g[i].numNeigh*sizeof(double));
    for(k=0;k<g[i].numNeigh;k++){
      g[i].a0[k]=0.;
      g[i].a1[k]=0.;
    }
  }
}

