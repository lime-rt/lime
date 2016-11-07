/*
 *  sourcefunc.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2016 The LIME development team
 *
 */

#include "lime.h"


/*....................................................................*/
void sourceFunc_line(const molData *md, const double vfac, const struct populations *mol\
  , const int lineI, double *jnu, double *alpha){

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
  */

  double jnu, trigFuncs[3];

  stokesangles(B, rotMat, trigFuncs);

  /* Emission */
  /* Continuum part:	j_nu = rho_dust * kappa_nu */
  jnu = cont.dust*cont.knu;
  snu[0] = jnu*(1.0 - maxp*(trigFuncs[0] - (2.0/3.0)));
  snu[1] = jnu*maxp*trigFuncs[1]*trigFuncs[0];
  snu[2] = jnu*maxp*trigFuncs[2]*trigFuncs[0];
  
  /* Absorption */
  /* Continuum part: Dust opacity */
  *alpha = cont.knu;
}


