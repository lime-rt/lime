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

void sourceFunc_pol(double *snu, double *alpha, struct grid *gp, int pos, int ispec, int iline, double (*rotMat)[3]){
  /*
The theory behind this was originally drawn from

  Padovani, M. et al, A&A 543, A16 (2012)

and references therein. However, as pointed out in Ade, P. A. R. et al, A&A 576, A105 (2015), Padovani's expression for sigma2 is too small by a factor of 2. This correction has been propagated here.
  */
  double jnu, trigFuncs[3];
  
  stokesangles(gp[pos].x[0], gp[pos].x[1], gp[pos].x[2], rotMat, trigFuncs);
  
  /* Emission */
  /* Continuum part:	j_nu = rho_dust * kappa_nu */
  jnu = gp[pos].mol[ispec].dust[iline]*gp[pos].mol[ispec].knu[iline];
  snu[0] = jnu*(1.0 - maxp*(trigFuncs[0] - 2.0/3.0));
  snu[1] = jnu*maxp*trigFuncs[1]*trigFuncs[0];
  snu[2] = jnu*maxp*trigFuncs[2]*trigFuncs[0];
  
  /* Absorption */
  /* Continuum part: Dust opacity */
  *alpha = gp[pos].mol[ispec].knu[iline];
}

