/*
 *  sourcefunc.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015 The LIME development team
 *
TODO:
  - Merge sourceFunc_*_raytrace() and sourceFunc_*() after changes to the way grid data is stored makes this possible.
 */

#include "lime.h"


/*....................................................................*/
void
sourceFunc_cont(const struct populations *gm, const int lineI, double *jnu\
  , double *alpha){

  /* Emission */
  /* Continuum part:	j_nu = T_dust * kappa_nu */
  *jnu   += gm->dust[lineI]*gm->knu[lineI];

  /* Absorption */
  /* Continuum part: Dust opacity */
  *alpha += gm->knu[lineI];

  return;
}

/*....................................................................*/
void
sourceFunc_line(const molData *md, const double vfac, const struct populations *gm\
  , const int lineI, double *jnu, double *alpha){

  double factor = vfac*HPIP*gm->binv*gm->nmol;
  /* Line part:		j_nu = v*consts*1/b*rho*n_i*A_ij */
  *jnu   += factor*gm->pops[md->lau[lineI]]*md->aeinst[lineI];

  /* Line part: alpha_nu = v*const*1/b*rho*(n_j*B_ij-n_i*B_ji) */
  *alpha += factor*(gm->pops[md->lal[lineI]]*md->beinstl[lineI]
                                      -gm->pops[md->lau[lineI]]*md->beinstu[lineI]);

  return;
}



/*....................................................................*/
void
sourceFunc_line_raytrace(const molData md, const double vfac\
  , const struct pop2 gm, const int lineI, double *jnu, double *alpha){

  /* Line part:		j_nu = v*consts*1/b*rho*n_i*A_ij */
  *jnu   += vfac*HPIP*gm.specNumDens[md.lau[lineI]]*md.aeinst[lineI];

  /* Line part: alpha_nu = v*const*1/b*rho*(n_j*B_ij-n_i*B_ji) */
  *alpha += vfac*HPIP*(gm.specNumDens[md.lal[lineI]]*md.beinstl[lineI]
                      -gm.specNumDens[md.lau[lineI]]*md.beinstu[lineI]);

  return;
}


/*....................................................................*/
void
sourceFunc_cont_raytrace(const struct pop2 gm, const int lineI, double *jnu\
  , double *alpha){

  /* Emission */
  /* Continuum part:	j_nu = T_dust * kappa_nu */
  *jnu   += gm.dust[lineI]*gm.knu[lineI];

  /* Absorption */
  /* Continuum part: Dust opacity */
  *alpha += gm.knu[lineI];

  return;
}

/*....................................................................*/
void
sourceFunc_pol(double B[3], const struct pop2 gm, int iline\
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
  jnu = gm.dust[iline]*gm.knu[iline];
  snu[0] = jnu*(1.0 - maxp*(trigFuncs[0] - 2.0/3.0));
  snu[1] = jnu*maxp*trigFuncs[1]*trigFuncs[0];
  snu[2] = jnu*maxp*trigFuncs[2]*trigFuncs[0];
  
  /* Absorption */
  /* Continuum part: Dust opacity */
  *alpha = gm.knu[iline];
}

