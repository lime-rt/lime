/*
 *  model_BoEb56.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#include "ml_recipes.h"


double m_rhoc,m_T,logQ,or1,cs2;

/*....................................................................*/
int
BoEb56_onFinalizeConfiguration(void){
  int i;

  if(getParamI("rhoc", &i)) return ML_UNRECOG_PARAM;
  m_rhoc   = modelDblPars[i]*AMU*ML_MEAN_MOL_WT*1.0e6; /* Using AMU rather than the proton mass as AJ had it. */
  if(getParamI("T", &i)) return ML_UNRECOG_PARAM;
  m_T      = modelDblPars[i];

  logQ = log(BoEb56_q);
  or1 = 10.0*AU;
  cs2   = KBOLTZ*m_T*ML_oneOnMuMp;

  return 0;
}

/*....................................................................*/
double
BoEb56_density(double x, double y, double z){
  double r0,r1,rho,r,dum,mr05,dr=0.0,ri;
  int nstep,i;

//  if(isFirstCall){
//    isFirstCall = 0;
//    logQ = log(BoEb56_q);
//    oneOnMuMp = 1.0/BoEb56_mu/MPROTON;
//    or1 = 10.0*AU;
//    m_rhoc = modelDblPars[1]*MPROTON*BoEb56_mu*1.0e6;//********************* need to make sure rhoc IS the 2nd param of this model. Pass in a list not a dict!!
//    cs2   = KBOLTZ*modelDblPars[0]*ML_oneOnMuMp;
//  }

//#define pi 3.14159265358979323846264338328
//#define q 1.03 // Relative stepsize in the integration (r_i+1 / r_i)
//#define mp 1.6726000e-27 // Proton mass (SI)
//#define gg 6.674e-11 // gravitational constant (SI)
//#define kk 1.3807e-23 // Boltzman constant (SI)
//#define au 1.496e11 // Astronomical unit (SI)
//#define mu 2.3 // mean molecular weight 



// Set the boundary condition and pre-calculate some variables 

  r0 = 0e0;
//  or1 = 10.*AU;
  r1 = or1;
  rho = m_rhoc;
  r = sqrt(x*x + y*y + z*z);
  
  if (r<=or1){
//return rho/(BoEb56_mu*MPROTON);
return rho*ML_oneOnMuMp;
  }


//  cs2   = KBOLTZ*m_T / (BoEb56_mu*MPROTON);
//  cs2   = KBOLTZ*m_T*ML_oneOnMuMp;
  dum   = log(r/r1) / logQ;
  nstep = (int)ceil(dum);
  mr05  = 0.0;
  for(i=0;i<nstep-1;i++){
    dr    = r1 - r0;
    ri    = r0 + 0.5*dr;
    mr05  = mr05 + 4.0*M_PI*ri*ri*rho*dr;
    rho   = exp(log(rho) - GRAV*mr05/cs2/ri/ri*dr);
    r0    = r1;
    r1    = r0*BoEb56_q;
  }

  r1   = or1 * pow(BoEb56_q, dum);
  ri   = r0 + 0.5*dr;
  dr   = r1 - r0;
  mr05 = mr05 + 4.0*M_PI*ri*ri*rho*dr;
  rho  = exp(log(rho) - GRAV*mr05/cs2/ri/ri*dr) ;

//return rho/(BoEb56_mu*MPROTON);
return rho*ML_oneOnMuMp;
}

/*....................................................................*/
double
BoEb56_temperature(double x, double y, double z){
  return modelDblPars[0];
}

/*....................................................................*/
void
BoEb56_velocity(double x, double y, double z, double *vel){
  vel[0] = 0.0;
  vel[1] = 0.0;
  vel[2] = 0.0;	
}

