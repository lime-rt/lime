/*
 *  stokesangles.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2016 The LIME development team
 *
 */

#include "lime.h"

void stokesangles(double B[3], double (*rotMat)[3], double *trigFuncs){
  /*
This function rotates the B-field vector from the model frame to the observer frame, then calculates and returns some useful values which will in function sourceFunc_pol() make it easy to obtain the Stokes parameters of polarized submillimetre dust emission. (For an explanation of the reasons for choosing the particular quantities we do, see the comment in that function.)

Whenever one deals with polarized light, it is important to specify the coordinate systems carefully. In LIME the observer frame is defined such that, when the observer looks at the sky, the frame axes appear as follows:

	       ^ Y
	       |
	       |
	       |
	<------+
	X

The Z axis points into the page, away from the observer. Comparing this to normal astronomical coordinates one can see that X is in the direction of +ve right ascension and Y in the direction of +ve declination.

The IAU-recommended coordinate frame for expressing polarized light however is

	       ^ X
	       |
	       |
	       |
	<------O
	Y

with Z now emerging from the page (i.e pointing in the direction of propagation, towards the observer).

A vector defined in the LIME model basis can be converted to the observer basis by post-multiplying it with the image rotation matrix rotMat. (Another way of putting this is that the rows of rotMat are the unit vectors of the model coordinate frame expressed in the observer basis.) For the B field, this is expressed symbolically as

	Bp^T = B^T * rotMat

where ^T denotes transpose.
  */
  const int nDim=3;
  double Bp[nDim];
  int i, j;
  double BxySquared, BSquared;
	
  /* Rotate B to the observer frame.
  */
  for(i=0;i<nDim;i++){
    Bp[i] = 0.0;
    for(j=0;j<nDim;j++)
      Bp[i] += B[j]*rotMat[j][i];
  }

  /* Square of length of B projected into the observer XY plane:
  */
  BxySquared = Bp[0]*Bp[0] + Bp[1]*Bp[1];
  if(BxySquared==0.0){
    trigFuncs[0] = 0.0;
    trigFuncs[1] = 0.0;
    trigFuncs[2] = 0.0;
    return;
  }

  BSquared = BxySquared + Bp[2]*Bp[2];
  trigFuncs[0] = BxySquared / BSquared; /* cos^2 of the angle which Bp makes to the XY plane */

  /* cos(2*phi) = cos^2(phi) - sin^2(phi) */
  trigFuncs[1] = (Bp[0]*Bp[0] - Bp[1]*Bp[1])/BxySquared;

  /* sin(2*phi) = 2*cos(phi)*sin(phi) */
  trigFuncs[2] = 2.0*Bp[0]*Bp[1]/BxySquared;
}

