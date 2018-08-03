/*
 *  model_Me09.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#include "ml_recipes.h"

#ifndef min
  #define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif
#ifndef max
  #define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

double m_mdot;
double m_mstar;
double m_rc;
double m_mu;
double m_nu;
double m_x;
double m_doppler_b;

double vk,norm;

/*....................................................................*/
int 
Me09_onFinalizeConfiguration(void){
  int i;

  if(getParamI("mdot", &i)) return ML_UNRECOG_PARAM;
  m_mdot      = modelDblPars[i]*MSUN/YJULIAN; // Msun/yr to kg/s
  if(getParamI("mstar", &i)) return ML_UNRECOG_PARAM;
  m_mstar     = modelDblPars[i]*MSUN;
  if(getParamI("rc", &i)) return ML_UNRECOG_PARAM;
  m_rc        = modelDblPars[i]*AU;
  if(getParamI("mu", &i)) return ML_UNRECOG_PARAM;
  m_mu        = modelDblPars[i];
  if(getParamI("nu", &i)) return ML_UNRECOG_PARAM;
  m_nu        = modelDblPars[i];

  vk = sqrt(GRAV*m_mstar/(m_rc));
  norm = 1.0 / (Me09_mmw*AMU) / (4.0 * M_PI * vk * m_rc*m_rc );

  return 0;
}


/*....................................................................*/
double
get_theta0(double r, double theta_inp, double mu, double nu){
  int i;
  double theta0_final;
  double theta0[Me09_ngrid], sint02[Me09_ngrid], diff, rr_old, rr_new, th0_old, th0_new;
  double epsilon, e, phi0, zeta, der, theta, dum, dum2;

  if (theta_inp>M_PI/2.){
    theta = M_PI - theta_inp;
  } else{
    theta = theta_inp;
  }


// Before we would go on with the full grid interpolation
// let's check if theta is zero. If it is the grid interpolation
//  will fail anyway, and we know that if theta=0 then theta0=0 */

  if (theta==0.){
    theta0_final = 0.0;
  }else{

// OK, so theta is not zero and we need to calculate theta0 via interpolation 
// Set up the theta0 grid
//
    for(i=0;i<Me09_ngrid;i++){
      theta0[i] = theta * (i*1.0)/(Me09_ngrid*1.0-1.);
      sint02[i] = sin(theta0[i])*sin(theta0[i]);
    }

    diff     = -1.0;
    rr_old   = 0.;
    rr_new   = 0.;
    th0_new  = 0.;
    th0_old  = 0.;
    i        = 10;

    while(diff<0.0 && i<Me09_ngrid){
      rr_old  = rr_new;
      th0_old = th0_new;
      th0_new = theta0[i];

      epsilon = nu*nu + mu*mu *sint02[i] - 2.0*mu;
      e       = sqrt(1.0 + epsilon*sint02[i]); 
      dum2    = (1.0 - mu*sint02[i])/e; 
      phi0    = acos(min(dum2,1.0)); 
      dum     = min(cos(theta) / cos(th0_new),1.); 
      zeta    = acos(dum) + phi0;

      rr_new  = sint02[i] / (1.0 - e*cos(zeta)) / r;
      diff    = rr_new - 1.0;

      i++;
    }

// OK, now we either bracketed r or reached pi/2 now calculate theta0_final

    if (diff==0.){
      theta0_final = th0_new;
    }else{
      der   = (th0_new - th0_old) / (rr_new - rr_old);
      theta0_final = th0_old + der * (1.0 - rr_old); 
    }
  }

  if (theta_inp>M_PI/2.){
    theta0_final = M_PI - theta0_final;
  }

  return theta0_final;
}

/*....................................................................*/
double 
Me09_density(const double x, const double y, const double z){
//#define pi 3.14159265358979323846264338328
//#define gg 6.674e-11    // gravitational constant [cm^3 g^-1 s^-2]
//#define mmw 2.        // Mean molecular weight
//#define mp 1.6726e-27  // Proton mass	
//#define AU	1.49598e11

  double r, theta, rnorm, theta0, cost0, cost, dum_rho, rho, rho_norm;
  
  r      = sqrt(x*x + y*y + z*z);
  theta  = acos(z / (r+1e-90));

// Rough workaround - should be treated properly later on
//   Problem; if theta=0 -> theta0=0, then in dum3 we divide by zero...
//   so if theta=0 I take the density at theta=0.1)*/

  if (theta==0.){
    theta = theta + 0.1;
  }
  
  rnorm      = r/m_rc;
//  vk         = sqrt(GRAV*m_mstar/(m_rc)); 
  theta0     = get_theta0(rnorm, theta, m_mu, m_nu); 

  cost0  = cos(theta0);
  cost   = cos(theta);
  
  dum_rho = sin(theta0) / (rnorm*rnorm) / ((1.0 + (3. * cost0*cost0- 1.) / rnorm - 2.*m_mu*cost0*cost) *\
            sqrt(1.0 - cost*cost / (cost0*cost0)) + m_nu/sin(theta0)*(1.0 + cost*cost - 2.*cost0*cost0)) ; 

// This should in principle never happen, but to be sure..
// Well, it may happen if we want to calculate the density outside of m_rc / mu
//  TODO: Some rules should be set to give the density if r > m_rc/mu, or the r<m_rc/mu condition should be forced at the input				  

  //dum_rho = fabs(dum_rho);
  if (dum_rho!=dum_rho){
    dum_rho  = 0.0;
  }
//  rho_norm = m_mdot / (Me09_mmw*AMU) / (4.0 * M_PI * vk * m_rc*m_rc ); 
  rho_norm = m_mdot*norm; 

  rho = max(dum_rho, 0.)*rho_norm;

  return rho;
}

/*....................................................................*/
void
Me09_velocity(const double x, const double y, const double z, double* v){
//#define pi 3.14159265358979323846264338328
//#define gg 6.674e-11    // gravitational constant [cm^3 g^-1 s^-2]
//#define AU		1.49598e11


  double r, rnorm, theta, theta0, sint02, epsilon, e, phi0, zeta;//, vk;	
  double rcyl, sinp, cosp, sint, cost;
  double vr, vt, vp;

  r      = sqrt(x*x + y*y + z*z);
  rnorm  = r/m_rc;
  theta  = acos((z) / (r+1e-90));
  
  theta0 = get_theta0(rnorm, theta, m_mu, m_nu);   
  sint02 = sin(theta0)*sin(theta0);

  epsilon = m_nu*m_nu + m_mu*m_mu *sint02 - 2.0*m_mu;
  e       = sqrt(1.0 + epsilon*sint02);
  phi0    = acos(min((1.0 - m_mu*sint02)/e, 1.0));
  zeta    = acos(min(cos(theta) / cos(theta0), 1.0)) + phi0;
//  vk      = sqrt(GRAV*m_mstar/(m_rc));

  if (theta==0.){
    vr = -vk / rnorm ;
    vt = 0.;
    vp = 0.;
  } else{
    vr      = -e*sin(zeta)*sin(theta0) / (rnorm * (1.0 - e*cos(zeta))) * vk ;
    vt      = sin(theta0) / rnorm / sqrt(sin(theta)) * sqrt(cos(theta0)*cos(theta0) - cos(theta)*cos(theta))* vk ;
    vp      = sint02 / (rnorm*sin(theta))* vk ;
  }

  if (z<0.){
    vt = -vt;
  }

  rcyl = sqrt(x*x + y*y);

  sinp = y / (rcyl+1e-90);
  cosp = x / (rcyl+1e-90);
  sint = rcyl / (r+1e-90);
  cost = z    / (r+1e-90);

  v[0] = (sint*cosp*vr + cost*cosp*vt - sinp*vp) ;
  v[1] = (sint*sinp*vr + cost*sinp*vt + cosp*vp) ;
  v[2] = (cost*vr - sint*vt) ;
}


