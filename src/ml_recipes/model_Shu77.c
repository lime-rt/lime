/*
 *  model_Shu77.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#include "ml_recipes.h"

double m_T;
double m_time;

double stab_x[20] = {0.0500000 ,0.100000 ,0.150000 ,0.200000 ,0.250000 ,0.300000 ,0.350000 ,0.400000 ,0.450000 ,0.500000 ,0.550000 ,0.600000 ,0.650000 ,0.700000 ,0.750000 ,0.800000 ,0.850000 ,0.900000 ,0.950000 ,1.00000}; 
double stab_v[20] = {5.44000 ,3.47000 ,2.58000 ,2.05000 ,1.68000 ,1.40000 ,1.18000 ,1.01000 ,0.861000 ,0.735000 ,0.625000 ,0.528000 ,0.442000 ,0.363000 ,0.291000 ,0.225000 ,0.163000 ,0.106000 ,0.0510000 ,0.00000};
double stab_alpha[20] = {71.5000 ,27.8000 ,16.4000 ,11.5000 ,8.76000 ,7.09000 ,5.95000 ,5.14000 ,4.52000 ,4.04000 ,3.66000 ,3.35000 ,3.08000 ,2.86000 ,2.67000 ,2.50000 ,2.35000 ,2.22000 ,2.10000 ,2.00000};

/*....................................................................*/
int 
Shu77_onFinalizeConfiguration(void){
  int i;

  if(getParamI("T", &i)) return ML_UNRECOG_PARAM;
  m_T = modelDblPars[i];
  if(getParamI("time", &i)) return ML_UNRECOG_PARAM;
  m_time = modelDblPars[i]*YJULIAN;

  return 0;
}

/*....................................................................*/
double 
Shu77_density(const double x, const double y, const double z){
  double r, a, rc, mod_x, dlx, dly, alpha;
  unsigned long id ;
  double rho;

  r = sqrt(x*x + y*y + z*z);
  a  = sqrt(KBOLTZ*m_T/(ML_MEAN_MOL_WT*AMU));
  rc = a*m_time;
 
  // Hydrostatic outer envelope
  if (r>=rc){
    mod_x = r / (a*m_time);
    alpha = 2e0 / (mod_x * mod_x);
    rho = alpha * (a*a) / (2e0 * M_PI* GRAV) / (r*r)/(ML_MEAN_MOL_WT*AMU);
  }
 // Collapse solution  
  else if(r>=rc*0.05){

    mod_x = r / (a*m_time);
// 	hunt(stab_x, 20, mod_x, &id);
    id = gsl_interp_bsearch(stab_x, mod_x, 0, 19);
    dlx = log10(stab_x[id+1]/stab_x[id]);
    dly =  log10(stab_alpha[id+1]/stab_alpha[id]);
    alpha = pow(10., dly/dlx * log10(mod_x/stab_x[id])) * stab_alpha[id]; 

    rho=  alpha/ ((4e0 * M_PI * GRAV * m_time*m_time)*(ML_MEAN_MOL_WT*AMU));
  }else{
  // Polynomial extrapolation inwards of 0.05 * Rcrit, which is the innermost
  // point in Shu's table that has a finite alpha and can therefore be used for
  // interpolation/extrapolation
    mod_x = r / (a*m_time);
    dlx = log10(stab_x[1]/stab_x[0]);
    dly =  log10(stab_alpha[1]/stab_alpha[0]);
    alpha = pow(10., dly/dlx * log10(mod_x/stab_x[0])) * stab_alpha[0]; 
    rho   = alpha/ ((4e0 * M_PI * GRAV * m_time*m_time)*(ML_MEAN_MOL_WT*AMU));
  }
  return rho;
}

/*....................................................................*/
double
Shu77_temperature(const double x, const double y, const double z){
  return m_T; 
}

/*....................................................................*/
double
Shu77_t_dust(const double x, const double y, const double z){
  return Shu77_temperature(x, y, z);
}

/*....................................................................*/
void
Shu77_velocity(const double x, const double y, const double z, double* v){
//#define pi 3.14159265358979323846264338328e0
//  double mp = 1.6726000e-27; // Proton mass
//  double kk = 1.3807e-23; // Boltzman constant
  
  double r, a, rc, mod_x, dlx, dly, sinp, sint, cosp, cost,rcyl, vminus , vr;
  unsigned long id; 
  r = sqrt(x*x + y*y + z*z);
// sound speed

  a  = sqrt(KBOLTZ*m_T/(ML_MEAN_MOL_WT*AMU));
  rc = a*m_time;

  // Hydrostatic outer envelope
  if (r>=rc){
    vr = 0.0;
  } 
 // Collapse solution  
  else if(r>=rc*0.05){
    mod_x = r / (a*m_time);
 	
//	hunt(stab_x, 21, mod_x, &id);
    id = gsl_interp_bsearch(stab_x, mod_x, 0, 20);

    dlx = log10(stab_x[id+1]/stab_x[id]);
    dly =  log10(stab_v[id+1]/stab_v[id]);
    vminus = pow(10e0, dly/dlx * log10(mod_x/stab_x[id])) * stab_v[id]; 
    vr =  -a * vminus;
  }
  // Polynomial extrapolation inwards of 0.05 * Rcrit, which is the innermost
  // point in Shu's table that has a finite alpha and can therefore be used for
  // interpolation/extrapolation
  else{
    mod_x = r / (a*m_time);
    dlx = log10(stab_x[1]/stab_x[0]);
    dly =  log10(stab_v[1]/stab_v[0]);
    vminus = pow(10., dly/dlx * log10(mod_x/stab_x[0])) * stab_v[0]; 
    vr  = -a*vminus;
  }
  
  rcyl = sqrt(x*x + y*y);

  sinp = y / (rcyl+1e-90);
  cosp = x / (rcyl+1e-90);
  sint = rcyl / (r+1e-90);
  cost = z    / (r+1e-90);

  v[0] = sint*cosp*vr;
  v[1] = sint*sinp*vr;
  v[2] = cost*vr;
}


