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
sourceFunc_line(double *jnu, double *alpha, molData *m,double vfac,struct grid *g,int pos,int ispec, int iline){
  
  /* Line part:		j_nu = v*consts*1/b*rho*n_i*A_ij */
  *jnu   += vfac*HPIP*g[pos].mol[ispec].binv*g[pos].nmol[ispec]*g[pos].mol[ispec].pops[m[ispec].lau[iline]]*m[ispec].aeinst[iline];
  
  /* Line part: alpha_nu = v*const*1/b*rho*(n_j*B_ij-n_i*B_ji) */
  *alpha += vfac*HPIP*g[pos].mol[ispec].binv*g[pos].nmol[ispec]*(g[pos].mol[ispec].pops[m[ispec].lal[iline]]*m[ispec].beinstl[iline]
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

void sourceFunc_pol(double *snu, double *alpha, struct grid *g, int pos, int ispec, int iline, double (*rotMat)[3]){
  /*
The theory behind this is drawn from

  Padovani, M. et al, A&A 543, A16 (2012)

and references therein.
  */
  double jnu, trigFuncs[3];
  
  stokesangles(g[pos].x[0], g[pos].x[1], g[pos].x[2], rotMat, trigFuncs);
  
  /* Emission */
  /* Continuum part:	j_nu = rho_dust * kappa_nu */
  jnu = g[pos].mol[ispec].dust[iline]*g[pos].mol[ispec].knu[iline];
  snu[0] = jnu*(1.0 - maxp*(0.5*trigFuncs[0] - 1.0/3.0));
  snu[1] = jnu*maxp*trigFuncs[1]*trigFuncs[0];
  snu[2] = jnu*maxp*trigFuncs[2]*trigFuncs[0];
  
  /* Absorption */
  /* Continuum part: Dust opacity */
  *alpha = g[pos].mol[ispec].knu[iline];
}


