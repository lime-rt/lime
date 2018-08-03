/*
 *  sourcefunc.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#include "lime.h"

/*....................................................................*/
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

Note that this is called from within a multi-threaded block.
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

/*....................................................................*/
void sourceFunc_line(const molData *md, const double vfac, const struct populations *mol\
  , const int lineI, double *jnu, double *alpha){
  /*
Note that this is called from within a multi-threaded block.
  */

  /* Line part:		j_nu = v*consts*1/b*rho*n_i*A_ij */
  *jnu   += vfac*HPIP*mol->specNumDens[md->lau[lineI]]*md->aeinst[lineI];

  /* Line part: alpha_nu = v*const*1/b*rho*(n_j*B_ij-n_i*B_ji) */
  *alpha += vfac*HPIP*(mol->specNumDens[md->lal[lineI]]*md->beinstl[lineI]
                      -mol->specNumDens[md->lau[lineI]]*md->beinstu[lineI]);

  return;
}

/*....................................................................*/
void sourceFunc_cont(const struct continuumLine cont, double *jnu\
  , double *alpha){
  /*
Note that this is called from within a multi-threaded block.
  */

  /* Emission */
  /* Continuum part:	j_nu = T_dust * kappa_nu */
  *jnu   += cont.dust*cont.knu;

  /* Absorption */
  /* Continuum part: Dust opacity */
  *alpha += cont.knu;

  return;
}

/*....................................................................*/
void sourceFunc_pol(double B[3], const struct continuumLine cont\
  , double (*rotMat)[3], double *snu, double *alpha){
  /*
The theory behind this was originally drawn from

  Padovani, M. et al, A&A 543, A16 (2012)

and references therein. However, as pointed out in Ade, P. A. R. et al, A&A 576, A105 (2015), Padovani's expression for sigma2 is too small by a factor of 2. This correction has been propagated here.

Note that this is called from within a multi-threaded block.
  */

  double jnu, trigFuncs[3];
  const double max_polarization = 0.15;

  stokesangles(B, rotMat, trigFuncs);

  /* Emission */
  /* Continuum part:	j_nu = rho_dust * kappa_nu */
  jnu = cont.dust*cont.knu;
  snu[0] = jnu*(1.0 - max_polarization*(trigFuncs[0] - (2.0/3.0)));
  snu[1] = jnu*max_polarization*trigFuncs[1]*trigFuncs[0];
  snu[2] = jnu*max_polarization*trigFuncs[2]*trigFuncs[0];
  
  /* Absorption */
  /* Continuum part: Dust opacity */
  *alpha = cont.knu;
}

/*....................................................................*/
void calcSourceFn(double dTau, const configInfo *par, double *remnantSnu, double *expDTau){
  /*
  The source function S is defined as j_nu/alpha, which is clearly not
  defined for alpha==0. However S is used in the algorithm only in the
  term (1-exp[-alpha*ds])*S, which is defined for all values of alpha.
  The present function calculates this term and returns it in the
  argument remnantSnu. For values of abs(alpha*ds) less than a pre-
  calculated cutoff supplied in configInfo, a Taylor approximation is
  used.

  Note that the same cutoff condition holds for replacement of
  exp(-dTau) by its Taylor expansion to 3rd order.

  Note that this is called from within the multi-threaded block.
  */

#ifdef FASTEXP
  *expDTau = FastExp(dTau);
  if (fabs(dTau)<par->taylorCutoff){
    *remnantSnu = 1. - dTau*(1. - dTau*(1./3.))*(1./2.);
  } else {
    *remnantSnu = (1.-(*expDTau))/dTau;
  }
#else
  if (fabs(dTau)<par->taylorCutoff){
    *remnantSnu = 1. - dTau*(1. - dTau*(1./3.))*(1./2.);
    *expDTau = 1. - dTau*(*remnantSnu);
  } else {
    *expDTau = exp(-dTau);
    *remnantSnu = (1.-(*expDTau))/dTau;
  }
#endif
}

