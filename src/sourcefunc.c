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
  if(doline) jnu	+= vfac*HPIP*g[pos].mol[ispec].binv*g[pos].nmol[ispec]*g[pos].mol[ispec].pops[m[ispec].lau[iline]]*m[ispec].aeinst[iline];
  
  
  
  /* Absorption */
  /* Continuum part: Dust opacity */
  alpha  = g[pos].mol[ispec].knu[iline];
  
  
  /* Line part: alpha_nu = v*const*1/b*rho*(n_j*B_ij-n_i*B_ji) */
  if(doline) alpha += vfac*HPIP*g[pos].mol[ispec].binv*g[pos].nmol[ispec]*(g[pos].mol[ispec].pops[m[ispec].lal[iline]]*m[ispec].beinstl[iline]
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
sourceFunc_line(double *jnu, double *alpha, molData *md, double vfac, struct grid *gp, int pos, int ispec, int iline){
  
  /* Line part:		j_nu = v*consts*1/b*rho*n_i*A_ij */
  *jnu   += vfac*HPIP*gp[pos].mol[ispec].binv*gp[pos].nmol[ispec]*gp[pos].mol[ispec].pops[md[ispec].lau[iline]]*md[ispec].aeinst[iline];
  
  /* Line part: alpha_nu = v*const*1/b*rho*(n_j*B_ij-n_i*B_ji) */
  *alpha += vfac*HPIP*gp[pos].mol[ispec].binv*gp[pos].nmol[ispec]*(gp[pos].mol[ispec].pops[md[ispec].lal[iline]]*md[ispec].beinstl[iline]
                                                                  -gp[pos].mol[ispec].pops[md[ispec].lau[iline]]*md[ispec].beinstu[iline]);
  
  return;
}

void
sourceFunc_cont(double *jnu, double *alpha, struct grid *gp, int pos, int ispec, int iline){
  
  /* Emission */
  /* Continuum part:	j_nu = T_dust * kappa_nu */
  *jnu   += gp[pos].mol[ispec].dust[iline]*gp[pos].mol[ispec].knu[iline];
  
  /* Absorption */
  /* Continuum part: Dust opacity */
  *alpha += gp[pos].mol[ispec].knu[iline];
  
  return;
}

void
sourceFunc_pol(double *snu, double *dtau, double ds, molData *m, double vfac, struct grid *gp, int pos, int ispec, int iline, double incl){
  double dSigma, dSigma2, dI, dQ, dU, alpha;
  double angle[3];
  
  stokesangles(gp[pos].x[0],gp[pos].x[1],gp[pos].x[2],incl,angle);
  
  /* Emission */
  /* Continuum part:	j_nu = rho_dust * kappa_nu */
  dSigma  = gp[pos].mol[ispec].dust[iline]*gp[pos].mol[ispec].knu[iline];
  dSigma2 = maxp*gp[pos].mol[ispec].dust[iline]*gp[pos].mol[ispec].knu[iline]*(0.5*angle[2]*angle[2]-1./3.);
  dI      = dSigma - dSigma2;
  dQ      = maxp*gp[pos].mol[ispec].dust[iline]*gp[pos].mol[ispec].knu[iline]*(2.*angle[0]*angle[0]-1.)*angle[2]*angle[2];
  dU      = maxp*gp[pos].mol[ispec].dust[iline]*gp[pos].mol[ispec].knu[iline]*(2.*angle[0]*angle[1]*angle[2]*angle[2]);
  
  /* Absorption */
  /* Continuum part: Dust opacity */
  alpha  = gp[pos].mol[ispec].knu[iline];
  
  /* Calculate source function and tau */
  *snu=0.;
  *dtau=0.;
  if(fabs(alpha)>0.){
    snu[0]=(dI/alpha)*m[ispec].norminv;
    snu[1]=-(dQ/alpha)*m[ispec].norminv;
    snu[2]=-(dU/alpha)*m[ispec].norminv;
    
    *dtau= alpha*ds;
  }
  return;
}

