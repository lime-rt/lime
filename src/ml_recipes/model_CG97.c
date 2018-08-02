/*
 *  model_CG97.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#include "ml_recipes.h"

double m_rin,m_rout,m_plsig1,m_hph,m_sig0,m_mstar,m_rstar,m_tstar,m_bgdens,rootTwoPi;

double my_r[CG97_npgrid],my_ri[CG97_npgrid+2],my_hp[CG97_npgrid],my_alpha[CG97_npgrid],my_ti[CG97_npgrid],my_ts[CG97_npgrid];

    
/*....................................................................*/
double
get_kappav(void){
    return 400.0;
}

/*....................................................................*/
double
get_eps(const double temp, const double beta){
    return (8.0*M_PI*ML_gsize*KBOLTZ_cgs*temp/(HPLANCK_cgs*CLIGHT_cgs))*beta;
}

/*....................................................................*/
double
get_ts(const double r, const double rstar, const double tstar){

    double tds_old, tds, diff, cconv, dummy_const, epsilon;

    tds_old = sqrt(rstar/r) * tstar;
    tds = 2000.0;
    diff = (tds_old - tds)/tds;
    cconv = 1e-2;

    dummy_const = sqrt(rstar/r/2.) * tstar;


    while(fabs(diff)>cconv){
        tds_old = tds;
        epsilon = get_eps(tds, 1.0);
        tds = dummy_const / sqrt(sqrt(epsilon));
        diff = (tds_old - tds) / tds;
    }

    return tds;
}

/*....................................................................*/
double
get_ti_thick(const double alpha, const double r, const double rstar, const double tstar){
    
    double ti;

    ti = sqrt(sqrt(alpha/4.0)) * sqrt(rstar/r) * tstar;

    return ti;
    
}

/*....................................................................*/
double
get_ti_thincool(const double alpha, const double r, const double rstar, const double tstar, const double sigma){
//*** doesn't seem to be used?

    double ti_old, ti, diff, cconv, kappa_v, dummy_const, epsilon;

    ti_old = sqrt(rstar/r) * tstar;
    ti = 2000.0;
    diff = (ti_old-ti)/ ti;
    cconv = 1e-2;

    kappa_v = get_kappav();

    dummy_const = sqrt(sqrt(alpha/4./sigma/kappa_v)) *sqrt(rstar/r) * tstar;

    while(fabs(diff)>cconv){
        ti_old = ti;
        epsilon = get_eps(ti,1.0);
        ti = dummy_const / sqrt(sqrt(epsilon));
        diff = (ti_old - ti)/ti;
    }
    return ti;
}

/*....................................................................*/
double
get_ti_thin(double alpha, double r, double rstar, double tstar, double tds){
//*** doesn't seem to be used?

    double ti_old, ti, diff, cconv, eps_s, dummy_const, epsilon;

    ti_old = sqrt(rstar/r) * tstar;
    ti = 2000.0;
    diff = (ti_old-ti)/ ti;
    cconv = 1e-2;

    eps_s = get_eps(tds,1.0);

    dummy_const = sqrt(sqrt(alpha*eps_s*eps_s)) *sqrt(rstar/r) * tstar;

    while(fabs(diff)>cconv){
        ti_old = ti;
        epsilon = get_eps(ti,1.0);
        ti = dummy_const / sqrt(sqrt(epsilon));
        diff = (ti_old - ti)/ti;
    }

    return ti;
}

/*....................................................................*/
void
get_hpcg01(double r[], const double mstar, const double rstar\
  , const double tstar, const double mu, const double hph, double alpha[]\
  , double hpr[]){

    double tc, gamma[CG97_npgrid], ti[CG97_npgrid];
    int niter, ir;
    double ti_old, ti_new, diff, diff_max;

    tc = GRAV_cgs * mstar * mu * AMU_cgs / KBOLTZ_cgs / rstar;
    for(niter=0;niter<2;niter++){
        for(ir=0;ir<CG97_npgrid/2;ir++){
            if(ir==0){
                if (niter==0){
                    gamma[0] = 1.2;
                    gamma[1] = 1.2;
                } else{
                    gamma[2*ir] = 1.5 + 0.5 * log(ti[1]/ti[0]) / log(r[1]/r[0]);
                    gamma[2*ir+1] = gamma[2*ir];
                } 
            } else{
                gamma[2*ir] = 1.5 + 0.5 * log(ti[2*ir-1] / ti[2*ir-2]) / log(r[2*ir-1] / r[2*ir-2]);
                gamma[2*ir+1] = gamma[2*ir];

            }
            
        
            ti_old = 0.;
            ti_new = 10.;
            diff = fabs(ti_new - ti_old);
            diff_max = 0.05;

            while(diff>diff_max){
                ti_old = ti_new;
                ti_new = sqrt(sqrt( ((rstar / (3.*M_PI*r[2*ir])) + \
                        0.25 * (gamma[2*ir]-1.0) * hph * sqrt(ti_old/tc) * sqrt(r[2*ir]/rstar)) ))*\
                         sqrt(rstar/r[2*ir]) * tstar;
                diff = fabs(ti_new-ti_old);
            }
            ti[2*ir] = ti_new;

            ti_old = 0.; 
            ti_new = 10.;
            diff = fabs(ti_new - ti_old);
            diff_max = 0.05;
            
            while(diff>diff_max){
                ti_old = ti_new;
                ti_new = sqrt(sqrt( ((rstar / (3.*M_PI*r[2*ir+1])) + \
                        0.25 * (gamma[2*ir+1]-1.0) * hph * sqrt(ti_old/tc) * sqrt(r[2*ir+1]/rstar)) ))*\
                         sqrt(rstar/r[2*ir+1]) * tstar;
                diff = fabs(ti_new-ti_old);
            }
           ti[2*ir+1] = ti_new ;

           alpha[2*ir] = (4./3./M_PI) * rstar / r[2*ir] + (gamma[2*ir]-1.) * hph * sqrt(ti[2*ir]/tc) * \
                         sqrt(r[2*ir]/rstar);

           alpha[2*ir+1] = (4./3./M_PI) * rstar / r[2*ir+1] + (gamma[2*ir+1]-1.) * hph * sqrt(ti[2*ir+1]/tc) * \
                         sqrt(r[2*ir+1]/rstar);


           hpr[2*ir] = sqrt(ti[2*ir]/tc) * sqrt(r[2*ir] / rstar);
           hpr[2*ir+1] = sqrt(ti[2*ir+1]/tc) * sqrt(r[2*ir+1] / rstar);
        }
    }

}

/*....................................................................*/
int 
CG97_onFinalizeConfiguration(void){
  int ir;

// I expect the input to be cgs and here I convert everything to SI for Lime
//  m_rin    = m_paramDouble["rin"]*AU;
//  m_rout   = m_paramDouble["rout"]*AU;
//  m_plsig1 = m_paramDouble["plsig1"];
//  m_hph    = m_paramDouble["hph"];
//  m_sig0   = m_paramDouble["sig0"]*10.0;
//  m_mstar  = m_paramDouble["Mstar"]*MSUN;
//  m_rstar  = m_paramDouble["Rstar"]*RSUN;
//  m_tstar  = m_paramDouble["Tstar"];
//  m_bgdens = m_paramDouble["bgdens"]*1.0e6;

  int i;

  if(getParamI("rin", &i)) return ML_UNRECOG_PARAM;
  m_rin    = modelDblPars[i]*AU;
  if(getParamI("rout", &i)) return ML_UNRECOG_PARAM;
  m_rout   = modelDblPars[i]*AU;
  if(getParamI("plsig1", &i)) return ML_UNRECOG_PARAM;
  m_plsig1 = modelDblPars[i];
  if(getParamI("hph", &i)) return ML_UNRECOG_PARAM;
  m_hph    = modelDblPars[i];
  if(getParamI("sig0", &i)) return ML_UNRECOG_PARAM;
  m_sig0   = modelDblPars[i]*10.0;
  if(getParamI("Mstar", &i)) return ML_UNRECOG_PARAM;
  m_mstar  = modelDblPars[i]*MSUN;
  if(getParamI("Rstar", &i)) return ML_UNRECOG_PARAM;
  m_rstar  = modelDblPars[i]*RSUN;
  if(getParamI("Tstar", &i)) return ML_UNRECOG_PARAM;
  m_tstar  = modelDblPars[i];
  if(getParamI("bgdens", &i)) return ML_UNRECOG_PARAM;
  m_bgdens = modelDblPars[i]*1.0e6;

  my_ri[0] = 0.;
  for(ir=0; ir<CG97_npgrid+1; ir++){
    my_ri[ir+1] = m_rin*0.9*pow((m_rout*1.1/(m_rin*0.9)), (float)ir / (float)(CG97_npgrid-1))*100.0;
  }
  for (ir=0; ir<CG97_npgrid; ir++){
    my_r[ir] = sqrt(my_ri[ir+1] * my_ri[ir+2]);
  }

// The CG97 functions are written in cgs so use the cgs units for there
//    get_hpcg01(my_r, m_mstar, m_rstar, m_tstar, 2.3, m_hph, my_alpha, my_hp);
  get_hpcg01(my_r, m_mstar*1e3, m_rstar*1e2, m_tstar, ML_MEAN_MOL_WT, m_hph, my_alpha, my_hp);

  for (ir=0; ir<CG97_npgrid+1; ir++){
    my_ri[ir] = my_ri[ir]*1e-2;
  }

  for (ir=0; ir<CG97_npgrid; ir++){
    my_r[ir]  = my_r[ir]*1e-2;
    my_hp[ir] = my_hp[ir] * my_r[ir];
    my_ts[ir] = get_ts(my_r[ir]*1e2, m_rstar*1e2, m_tstar);
    my_ti[ir] = get_ti_thick(my_alpha[ir], my_r[ir]*1e2, m_rstar*1e2, m_tstar);
  }

  rootTwoPi = sqrt(2.0*M_PI);

  return 0;
}

/*....................................................................*/
double 
CG97_density(const double x, const double y, const double z){
//#define pi 3.14159265358979323846264338328
//#define mmw 2.3       // Mean molecular weight
//#define mpsi 1.6726e-27  // Proton mass	
  int id;
  double rc, rho;


  rc    = sqrt(x*x + y*y);

  if (rc > m_rin){
    if (rc < m_rout){
      id = gsl_interp_bsearch(my_ri, rc, 0, CG97_npgrid-1);
//      rho = m_sig0 * pow((rc / m_rout), m_plsig1)/ (my_hp[id] * sqrt(2.*M_PI)) * exp(-0.5* (z / my_hp[id]) * (z/my_hp[id]));
      rho = m_sig0 * pow((rc / m_rout), m_plsig1)/ (my_hp[id] * rootTwoPi) * exp(-0.5* (z / my_hp[id]) * (z/my_hp[id]));
    }
    else {
      rho   = 0.0;
    }    
  }
  else { 
    rho = 0.0;
  }
  
//return rho/(ML_MEAN_MOL_WT*AMU) + m_bgdens;
return rho*ML_oneOnMuMp + m_bgdens;
}

/*....................................................................*/
double
CG97_t_dust(const double x, const double y, const double z){
  int id;
  double rc, temp;

  rc    = sqrt(x*x + y*y);

  if (rc > m_rin){
    if (rc < m_rout){
      id = gsl_interp_bsearch(my_ri, rc, 0, CG97_npgrid-1);

      if (fabs(z)>m_hph*my_hp[id]){
        temp = my_ts[id];
      } else {
        temp = my_ti[id];
      } 
    }
    else {
      temp   = 2.73;
    }    
  }
  else { 
    temp = 2.73;
  }

  return temp;
}

/*....................................................................*/
void
CG97_velocity(const double x, const double y, const double z, double* v){
  /* Re-written to get rid of the singularity at rc==0. */

  double rc,vkep,sinp,cosp,rEff;

  rc   = sqrt(x*x + y*y);
//  vkep = sqrt(GRAV*m_mstar/rc);
//  sinp = y / (rc + 1e-90);
//  cosp = x / (rc + 1e-90);

  if (rc<m_rin)
    rEff = m_rin;
  else
    rEff = rc;

  vkep = sqrt(GRAV*m_mstar/rEff);
  sinp = y / rEff;
  cosp = x / rEff;

  v[0] =  vkep * sinp;
  v[1] = -vkep * cosp;
  v[2] = 0.0e+0;
}

