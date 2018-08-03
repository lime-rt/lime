/*
 *  model_Al03.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#include "ml_recipes.h"

double m_v0,m_t0,m_rn,m_T,m_r0,m_n0,m_b0;

/*....................................................................*/
int
Al03_onFinalizeConfiguration(void){
  int i;

  if(getParamI("cs", &i)) return ML_UNRECOG_PARAM;
  m_v0 = modelDblPars[i];
  if(getParamI("age", &i)) return ML_UNRECOG_PARAM;
  m_t0 = modelDblPars[i]*YJULIAN;
  if(getParamI("Rn", &i)) return ML_UNRECOG_PARAM;
  m_rn = modelDblPars[i]*AU;
  if(getParamI("T", &i)) return ML_UNRECOG_PARAM;
  m_T  = modelDblPars[i];

  m_r0 = m_v0*m_t0;
  m_n0 = 1.0/(4.0*M_PI*GRAV*m_t0*m_t0*2.0e-20); // mH2
  m_b0 = m_v0/(sqrt(GRAV)*m_t0);

  return 0;
}

/*....................................................................*/
double
Al03_density(const double x, const double y, const double z){
  double r,sinTheta;

  r = sqrt(x*x + y*y + z*z);
  sinTheta = sqrt(x*x + y*y)/r;

  if (r < m_rn)
    r = m_rn;

  return m_n0*pow(1.8*sinTheta/(r/m_r0), 1.5);
}

/*....................................................................*/
double
Al03_temperature(const double x, const double y, const double z){
  return m_T;
}

/*....................................................................*/
void
Al03_velocity(const double x, const double y, const double z, double *vel){
  double r,cosTheta,sinTheta,cosPhi,sinPhi;
  double vr;

  r = sqrt(x*x + y*y + z*z);  

  cosTheta = z / r;
  sinTheta = sqrt(x*x + y*y) / r;
  cosPhi = x / (r * sinTheta);
  sinPhi = y / (r * sinTheta);

  // Spherical components:
  vr = - m_v0 * exp(-(r / m_r0) * (r / m_r0)) / sqrt(r / m_r0);

  // Cartesian components:
  vel[0] = vr * sinTheta * cosPhi;
  vel[1] = vr * sinTheta * sinPhi;
  vel[2] = vr * cosTheta;
}

/*....................................................................*/
void
Al03_bmag(const double x, const double y, const double z, double *b){
  double r, cosTheta, sinTheta, cosPhi, sinPhi;
  double r0975;
  double dphir, dphit;
  double br, btheta; 

  r = sqrt(x*x + y*y + z*z);    
  cosTheta = z / r;
  sinTheta = sqrt(x*x + y*y) / r;
  cosPhi = x / (r * sinTheta);
  sinPhi = y / (r * sinTheta);

  // Spherical components:
  r0975 = (r / m_r0) + 0.975;

  dphir = r0975 * M_PI * M_PI * sinTheta * sinTheta / 2; //  / ( (r / r0 ) * sinTheta ) ? 
  dphit = r0975 * r0975 * M_PI * M_PI * sinTheta * cosTheta / 2; //  / ( (r / r0 ) * sinTheta ) ? 

  br = m_b0 * dphit / r;
  btheta = - m_b0 * dphir;

  // Cartesian components:
  b[0] = br * sinTheta * cosPhi + btheta * cosTheta * cosPhi;
  b[1] = br * sinTheta * sinPhi + btheta * cosTheta * sinPhi;
  b[2] = br * cosTheta - btheta * sinPhi;
}

