/*
 *  sourcefunc.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015 The LIME development team
 *
 */

#include "lime.h"

void
sourceFunc(double *snu, double *dtau, double ds, molData *m,double vfac,struct grid *g,int pos,int ispec, int iline, int doline){
  double jnu, alpha;
  
  /* Emission */
  /* Continuum part:	j_nu = T_dust * kappa_nu */
  jnu	 = g[pos].mol[ispec].dust[iline]*g[pos].mol[ispec].knu[iline];
  
  /* Line part:		j_nu = v*consts*1/b*rho*n_i*A_ij */
  if(doline) jnu	+= vfac*HPIP*g[pos].mol[ispec].binv*g[pos].mol[ispec].nmol*g[pos].mol[ispec].pops[m[ispec].lau[iline]]*m[ispec].aeinst[iline];
  
  
  
  /* Absorption */
  /* Continuum part: Dust opacity */
  alpha  = g[pos].mol[ispec].knu[iline];
  
  
  /* Line part: alpha_nu = v*const*1/b*rho*(n_j*B_ij-n_i*B_ji) */
  if(doline) alpha += vfac*HPIP*g[pos].mol[ispec].binv*g[pos].mol[ispec].nmol*(g[pos].mol[ispec].pops[m[ispec].lal[iline]]*m[ispec].beinstl[iline]
                                                                          -g[pos].mol[ispec].pops[m[ispec].lau[iline]]*m[ispec].beinstu[iline]);
  
  
  /* Calculate source function and tau */
  *snu=0.;
  *dtau=0.;
  if(fabs(alpha)>0.){
    *snu=(jnu/alpha)*m[ispec].norminv;
    *dtau= alpha*ds;
  }
  return;
}


void
sourceFunc_line(double *jnu, double *alpha, molData *m,double vfac,struct grid *g,int pos,int ispec, int iline){
  
  /* Line part:		j_nu = v*consts*1/b*rho*n_i*A_ij */
  *jnu   += vfac*HPIP*g[pos].mol[ispec].binv*g[pos].mol[ispec].nmol*g[pos].mol[ispec].pops[m[ispec].lau[iline]]*m[ispec].aeinst[iline];
  
  /* Line part: alpha_nu = v*const*1/b*rho*(n_j*B_ij-n_i*B_ji) */
  *alpha += vfac*HPIP*g[pos].mol[ispec].binv*g[pos].mol[ispec].nmol*(g[pos].mol[ispec].pops[m[ispec].lal[iline]]*m[ispec].beinstl[iline]
                                                                -g[pos].mol[ispec].pops[m[ispec].lau[iline]]*m[ispec].beinstu[iline]);
  
  return;
}

void
sourceFunc_cont(double *jnu, double *alpha,struct grid *g,int pos,int ispec, int iline){
  
  /* Emission */
  /* Continuum part:	j_nu = T_dust * kappa_nu */
  *jnu   += g[pos].mol[ispec].dust[iline]*g[pos].mol[ispec].knu[iline];
  
  /* Absorption */
  /* Continuum part: Dust opacity */
  *alpha += g[pos].mol[ispec].knu[iline];
  
  return;
}

/*....................................................................*/
void
sourceFunc_pol(const double ds, const double B[3], const molData md, const struct pop2 gm\
  , const int lineI, const double incl, double *snu, double *dtau){

  double dSigma, dSigma2, dI, dQ, dU, alpha;
  double trigFuncs[3];

  stokesangles(B, incl, trigFuncs);

  /* Emission */
  /* Continuum part:	j_nu = rho_dust * kappa_nu */
  dSigma  = gm.dust[lineI]*gm.knu[lineI];
  dSigma2 = maxp*gm.dust[lineI]*gm.knu[lineI]*(0.5*trigFuncs[2]*trigFuncs[2]-1./3.);
  dI      = dSigma - dSigma2;
  dQ      = maxp*gm.dust[lineI]*gm.knu[lineI]*(2.*trigFuncs[0]*trigFuncs[0]-1.)*trigFuncs[2]*trigFuncs[2];
  dU      = maxp*gm.dust[lineI]*gm.knu[lineI]*(2.*trigFuncs[0]*trigFuncs[1]*trigFuncs[2]*trigFuncs[2]);

  /* Absorption */
  /* Continuum part: Dust opacity */
  alpha  = gm.knu[lineI];

  /* Calculate source function and tau */
  snu[0] = 0.0;
  snu[1] = 0.0;
  snu[2] = 0.0;
  *dtau=0.;
  if(fabs(alpha)>0.){
    snu[0]= (dI/alpha)*md.norminv;
    snu[1]=-(dQ/alpha)*md.norminv;
    snu[2]=-(dU/alpha)*md.norminv;

    *dtau = alpha*ds;
  }
  return;
}


