/*
 *  model_DDN01.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#include "ml_recipes.h"
#include "../aux.h"

  double my_r[DDN01_nrmod];
  double my_ri[DDN01_nrmod+1];
  double my_hp[DDN01_nrmod];
  double my_hs[DDN01_nrmod];
  double my_ti[DDN01_nrmod];
  double my_ts[DDN01_nrmod];
  double my_sig0;

  double m_rin;
  double m_rout;
  double m_mdisk;
  double m_plsig1;
  double m_mstar;
  double m_rstar;
  double m_tstar;
  double m_bgdens;
  const char *m_dustopac_fname;

  double rootTwoPi;

/*....................................................................*/
void
read_kappa(const char *fname, double wav[], double kappa[]){
   
    double wmin=0.01, wmax=1e5; 
    FILE *rfile;
    char str[80],message[80];
    int nlines, i,  j;
    double *wav_tmp, *kappa_tmp, dum;
	gsl_spline *spline;
	gsl_interp_accel *acc=gsl_interp_accel_alloc();

    // Create the wavelength array */
    for (i=0; i<DDN01_nwav; i++){
        wav[i] = wmin * pow(wmax/wmin, (i*1.0)/(1.0*DDN01_nwav-1.));
    }
   
    // Read the opacity file */
    rfile = fopen(fname, "r");
    // Check the line numbers */
    nlines = 0;
    while(fgetc(rfile)!=EOF){
        sprintf(message, "line %d", nlines);
        checkFgets(fgets(str, 80, rfile), message);
        nlines++; 
    }

    rewind(rfile);
    wav_tmp = (double *)malloc(sizeof(double)*nlines);
    kappa_tmp = (double *)malloc(sizeof(double)*nlines);
    for (i=0;i<nlines;i++){
        checkFscanf(fscanf(rfile, "%lf %lf\n", &wav_tmp[i], &kappa_tmp[i]), 2, "wav & kappa");

        // In the Lime dust opacity files the wavelength should be in cm, so stick to that and convert that to micron */
        wav_tmp[i] = log10(wav_tmp[i]);
        kappa_tmp[i] = log10(kappa_tmp[i]);
    }
    fclose(rfile);

// Order the arrays in wavelength if necessary*/
    if(wav_tmp[1]<wav_tmp[0]){
        i=0;
        j=nlines-1;
        while(i<j){

            dum = wav_tmp[i];
            wav_tmp[i] = wav_tmp[j];
            wav_tmp[j] = dum;

            dum = kappa_tmp[i];
            kappa_tmp[i] = kappa_tmp[j];
            kappa_tmp[j] = dum;

            i++;
            j--;
        }
    }


//  OK, the file is read now check if there is wide enough wavelenght coverage for both the 
//   Heating and the cooling radiation fields */

    spline = gsl_spline_alloc(gsl_interp_cspline, nlines);
    gsl_spline_init(spline,wav_tmp, kappa_tmp, nlines);
   
    dum = (kappa_tmp[nlines-1] - kappa_tmp[nlines-2]) / (wav_tmp[nlines-1] - wav_tmp[nlines-2]);

    for(i=0;i<DDN01_nwav;i++){
        if (log10(wav[i])<=wav_tmp[0]){
            kappa[i] = pow(10.,kappa_tmp[0]);
        } else if (log10(wav[i])>=wav_tmp[nlines-1]){
            kappa[i] = pow(10,kappa_tmp[nlines-1]) * pow(10., (log10(wav[i])-wav_tmp[nlines-1]) * dum); 
        } else {
            kappa[i] = pow(10., gsl_spline_eval(spline, log10(wav[i]), acc));
        }
                
    }
    free(wav_tmp);
    free(kappa_tmp);
}

// ---------------------------------------------------------------------------------------------------------------------
//  Simple integrator *
// --------------------------------------------------------------------------------------------------------------------- */
double
integr(int n, double x[], double y[]){
    
    double dum;
    int i;

    // NOTE:
    // I use this function only to integrate in frequency and since all frequency-dependent quantities are ordered in 
    // wavelength instead of frequency the last term in the right-hand side is multiplied by -1
    dum = 0.0;
    for (i=0;i<n-1;i++){
        dum = dum + 0.5 * (y[i+1] + y[i]) * (x[i] - x[i+1]);
    } 

    return dum;
}


// Planck function */
void
bnu_t(double nu[], double T, double bnut[]){
    int i;
    
    for(i=0;i<DDN01_nwav;i++){
        bnut[i] = 2.0*HPLANCK_cgs*pow(nu[i],3)/CLIGHT_cgs/CLIGHT_cgs / (exp(HPLANCK_cgs*nu[i]/KBOLTZ_cgs/T)-1.0);
    }


    return;
}

// ---------------------------------------------------------------------------------------------------------------------
 // Calculates the ratio of the surface height over pressure scale height
 // --------------------------------------------------------------------------------------------------------------------- ///
double
get_chicg(double alpha, double sigma, double kappa_p){

    double chi_old, dum, threshold, chi, diff;

    chi_old = 0.;
    dum = 2.*alpha /sigma/kappa_p;
    threshold = 1e-4;
    chi = 4.0;
    diff = fabs(chi_old - chi)/chi;

    while(diff>threshold){  
        chi_old = chi;
        chi = sqrt(fabs(-2.0 * log(dum * exp(-0.5*chi_old*chi_old) / (1.0 - gsl_sf_erf(chi_old*sqrt(0.5))))));
        diff = fabs(chi_old-chi)/chi;
    }

    return chi;
}

// ---------------------------------------------------------------------------------------------------------------------
 // Calculates the planck mean opacity
 // --------------------------------------------------------------------------------------------------------------------- ///
double
get_kappa_planck(double nu[], double kappa[], double T){

    double bb[DDN01_nwav], dum[DDN01_nwav];
    int i;

    bnu_t(nu, T, bb);

    for (i=0;i<DDN01_nwav;i++){
        dum[i] = bb[i]*kappa[i];
    }
    //return integr(DDN01_nwav, nu, dum) / integr(DDN01_nwav, nu, bb);*/
    return integr(DDN01_nwav, nu, dum) / (STEFANB_cgs*pow(T,4)/M_PI);
}   

// ---------------------------------------------------------------------------------------------------------------------
 // Calculates the planck mean opacity (from an imput frequency depencent flux array instead of the temperature only)
 // --------------------------------------------------------------------------------------------------------------------- ///
double
get_kappa_planck_Firr(double nu[], double kappa[], double Firr[]){

    double dum[DDN01_nwav];
    int i;

    for (i=0;i<DDN01_nwav;i++){
        dum[i] = Firr[i]*kappa[i];
    }

    return integr(DDN01_nwav, nu, dum) / integr(DDN01_nwav, nu, Firr); 
}   

// ---------------------------------------------------------------------------------------------------------------------
 // Calculates the psi parameter for the disk interior
 // --------------------------------------------------------------------------------------------------------------------- ///
double
get_psii(double kappa[], double nu[], double Ti, double sigma){
    
    double bnut[DDN01_nwav], fnu[DDN01_nwav];
    int i;

    bnu_t(nu, Ti, bnut);
    for (i=0;i<DDN01_nwav;i++){
        fnu[i] = bnut[i] * (1.0 - exp(-sigma * kappa[i]));
    }

    return integr(DDN01_nwav, nu, fnu) / (STEFANB_cgs*pow(Ti,4) / M_PI);
    //return integr(DDN01_nwav, nu, fnu) / integr(DDN01_nwav, nu, bnut);*/
}

// ---------------------------------------------------------------------------------------------------------------------
 // Calculates the psi parameter for the disk surface layer
 // --------------------------------------------------------------------------------------------------------------------- ///
double
get_psis(double kappa[], double nu[], double Tsurf, double sigma){
    
    double bnut[DDN01_nwav], fnu[DDN01_nwav];
    int i;

    bnu_t(nu, Tsurf, bnut);
    for (i=0;i<DDN01_nwav;i++){
        bnut[i] = bnut[i] * kappa[i];
        fnu[i] = bnut[i] * (1.0 - exp(-sigma * kappa[i]));
    }

    return integr(DDN01_nwav, nu, fnu) / integr(DDN01_nwav, nu, bnut);
}

// ---------------------------------------------------------------------------------------------------------------------
 // Calculates the surface density at the outer radius from the disk mass
 // --------------------------------------------------------------------------------------------------------------------- ///
double
get_sig0(double rin, double rout, double mdisk, double plsig){

    double dum;

    if (plsig==-2.){
        dum = 2.0 * M_PI * pow(rout, -plsig) * (log(rout) - log(rin));
    } else {
        dum = 2.0 * M_PI * pow(rout, -plsig) / (2.0 + plsig) * (pow(rout, (2.+plsig)) - pow(rin, (2.+plsig)));
    }

    return mdisk/dum;
}

// ---------------------------------------------------------------------------------------------------------------------
 // Calculates the disk interior temperature 
 // --------------------------------------------------------------------------------------------------------------------- ///
double
get_ti(double alpha, double r, double rstar, double tstar, double kappa[], double nu[], double sigma, double tsurf, double ti_init, double Firr){

    double psi_s, psi_i, ti, ti_old, diff, threshold;

    ti = ti_init; 
    psi_s = get_psis(kappa, nu, tsurf, sigma); 
    psi_i = get_psis(kappa, nu, ti, sigma); 

    ti_old = ti;
    ti = sqrt(sqrt(Firr * 0.5 * psi_s / (psi_i * STEFANB_cgs)));
    diff = fabs(ti_old - ti) / ti;
    threshold = 0.01;

    while(diff>threshold){
        ti_old = ti;
        psi_i = get_psis(kappa, nu, ti, sigma); 
        ti = sqrt(sqrt(Firr * 0.5 * psi_s / (psi_i * STEFANB_cgs)));
        diff = fabs(ti_old - ti) / ti;
    }

    return ti;
}

// ---------------------------------------------------------------------------------------------------------------------
 // Calculates the rim temperature including self-irradiation
 // --------------------------------------------------------------------------------------------------------------------- ///
double
get_trim(double rstar, double rin, double mstar, double sigma, double kappa_p, double tstar){

    double chi, trim, trim_old, h, threshold, diff;

    chi = get_chicg(1.0, sigma, kappa_p*8.0);
    trim_old = 2000.;
    h = sqrt(KBOLTZ_cgs*trim_old*pow(rin,3) / (ML_MEAN_MOL_WT*AMU_cgs*GRAV_cgs*mstar));
    trim = sqrt(rstar/rin) * tstar * sqrt(sqrt((1. + h*chi/rin)));

    diff = fabs(trim_old -trim)/trim;
    threshold = 0.01;

    while(diff>threshold){
        trim_old = trim;
        h = sqrt(KBOLTZ_cgs*trim_old*pow(rin,3) / (ML_MEAN_MOL_WT*AMU_cgs*GRAV_cgs*mstar));
        trim = sqrt(rstar/rin) * tstar * sqrt(sqrt((1. + h*chi/rin)));
        diff = fabs(trim_old -trim)/trim;
    }


    return trim;
}

// ---------------------------------------------------------------------------------------------------------------------
 // Calculates the temperature in the surface layer
 // --------------------------------------------------------------------------------------------------------------------- ///
double
get_tsurf(double r, double rstar, double tstar, double kappa[], double nu[], double Firr[]){

    double ts, ts_old, diff, threshold, c1, eps;

    ts_old = sqrt(rstar/r) * tstar;
    ts = 2000.;
    diff = fabs(ts_old - ts)/ts;
    threshold = 0.01;

    c1 = sqrt(rstar/r/2.) * tstar;


    while(diff>threshold){
        ts_old = ts;
        eps = get_kappa_planck(nu, kappa, ts) / get_kappa_planck_Firr(nu, kappa, Firr);
        ts = c1 /sqrt(sqrt(eps));
        diff = fabs(ts_old - ts) / ts;
        
    }
   return ts; 
}

// ---------------------------------------------------------------------------------------------------------------------
 // Function to calculate the disk structure  (r, Hs, Hp, ti, ts)
 // --------------------------------------------------------------------------------------------------------------------- ///
void
get_disk_struct(double rin, double rdisk, double mdisk, double plsig1, const char *opac_fname, double mstar, double rstar, double tstar, 
    double r[], double ri[], double Hp[], double Hs[], double ti[], double ts[]){

    double wav[DDN01_nwav], nu[DDN01_nwav], kappa[DDN01_nwav], Firr[DDN01_nwav], Firr_int, Firr_star[DDN01_nwav][DDN01_nrmod], Firr_star_int[DDN01_nrmod], bnut[DDN01_nwav], dum[DDN01_nwav];
    double alpha[DDN01_nrmod], gamma[DDN01_nrmod], chi[DDN01_nrmod], sigma[DDN01_nrmod];
    double chirim[DDN01_nrmod], Hsrim[DDN01_nrmod], Hprim[DDN01_nrmod], alpha_iter[DDN01_nrmod], dum_arwav[DDN01_nwav], bnut_rim[DDN01_nwav];
//, trim_diff[DDN01_nrmod]
    double sig0, kappa_planck_star, trim, tc, alpha_old, delta, theta, diff1, diff2;
    int ir,inu,iiter, iish=0, iicg=0;
    int cg97_converged, converged;

    // Read the opacities */ 
    read_kappa(opac_fname, wav, kappa);
    for (inu=0; inu<DDN01_nwav; inu++){
        // The wavelength should be in micron */
        nu[inu] = CLIGHT_cgs/wav[inu]*1e4;
    }
    
    // Set up radial grid */
    sig0 = get_sig0(rin, rdisk, mdisk, plsig1)*0.01;

    for (ir=0; ir<DDN01_nrmod+1; ir++){
        ri[ir] = rin*0.9 * pow((rdisk*1.1/(rin*0.9)), (float)ir / (float)(DDN01_nrmod-1));
    }
    for (ir=0; ir<DDN01_nrmod; ir++){
        r[ir] = sqrt(ri[ir] * ri[ir+1]);
    }
    for (ir=0; ir<DDN01_nrmod+1; ir++){
        ri[ir] = ri[ir] * 1e-2; 
    }

    for (ir=0; ir<DDN01_nrmod; ir++){
        r[ir] = rin * pow(rdisk/rin, (1.0*ir) / (1.0*DDN01_nrmod -1.));
        sigma[ir] = sig0 * pow(r[ir]/rdisk, plsig1);
        alpha[ir] = 0.05;
    }
    

    kappa_planck_star = get_kappa_planck(nu, kappa, tstar);

    trim = get_trim(rstar, rin, mstar, sigma[0], kappa_planck_star, tstar);
    Hprim[0] = sqrt(KBOLTZ_cgs*trim*pow(rin,3) / (GRAV_cgs*mstar*ML_MEAN_MOL_WT*AMU_cgs));

    for (ir=0;ir<DDN01_nrmod;ir++){
        chirim[ir]    = get_chicg(1.0, sigma[ir], kappa_planck_star*8.0);
        Hsrim[ir]     = chirim[0] * Hprim[0] - (r[ir] - r[0])/8.0;
        Hsrim[ir]     = (Hsrim[ir] + fabs(Hsrim[ir]))*0.5;
        Hprim[ir]     = Hsrim[ir]/chirim[ir];
//        trim_diff[ir] = trim * (Hsrim[ir]/Hsrim[0]) * (Hsrim[ir]/Hsrim[0]);
        if (Hsrim[ir]>0)
            iish = ir+1;
        ti[ir] = 100.;
    }

    
// First go for a CG97 disk */
    tc = GRAV_cgs*mstar*ML_MEAN_MOL_WT*AMU_cgs/KBOLTZ_cgs/rstar;

    gamma[0] = 1.28;
    gamma[1] = 1.28;
    alpha[0] = 0.05;


    for (ir=0; ir<DDN01_nrmod; ir++){
        bnu_t(nu, tstar, bnut);
        for (inu=0; inu<DDN01_nwav; inu++){
            Firr_star[inu][ir] = bnut[inu] * M_PI * pow((rstar/r[ir]), 2);
            dum[inu] = Firr_star[inu][ir] * alpha[ir];
        }
        Firr_star_int[ir] = integr(DDN01_nwav, nu, dum);
    }


    iiter=0;
    cg97_converged = 0;

    while(cg97_converged!=1){
        iiter++;
        diff2 = 0.;
        for (ir=0; ir<DDN01_nrmod; ir++){

            alpha_iter[ir] = alpha[ir];
            converged = 0;

            if ((ir%2)==0){
                if (ir>1){
                    gamma[ir] = (r[ir-1] + r[ir-2]) / (Hs[ir-1] + Hs[ir-2]) * (Hs[ir-2] - Hs[ir-1]) / (r[ir-2] - r[ir-1]);
                    gamma[ir+1] = gamma[ir];
                }
            }
            if (ir==0){
                if (iiter>1){
                    gamma[ir]   = (r[ir+2] + r[ir+3]) / (Hs[ir+2] + Hs[ir+3]) * (Hs[ir+3] - Hs[ir+2]) / (r[ir+3] - r[ir+2]);
                    gamma[ir+1] = gamma[ir];
                }
            }

            for (inu=0;inu<DDN01_nwav;inu++){
                dum_arwav[inu] = Firr_star[inu][ir];
            }

            while (converged!=1){
                alpha_old = alpha[ir];

                ts[ir] = get_tsurf(r[ir], rstar, tstar, kappa, nu, dum_arwav);
                ti[ir] = sqrt(sqrt(0.5*Firr_star_int[ir])) * trim;
                ti[ir] = get_ti(alpha[ir], r[ir], rstar, tstar, kappa, nu, sigma[ir], ts[ir], ti[ir], Firr_star_int[ir]);
                chi[ir] = get_chicg(alpha[ir], sigma[ir], kappa_planck_star);
                Hp[ir] = sqrt(ti[ir]/tc) * sqrt(r[ir]/rstar) * r[ir];
                Hs[ir] = Hp[ir] * chi[ir];
                

                alpha[ir] = (4./3./M_PI) * rstar / r[ir] + (gamma[ir]-1.) * chi[ir] * sqrt(ti[ir]/tc) * sqrt(r[ir]/rstar);
                //alpha[ir] = (gamma[ir]-1.) * Hs[ir]/r[ir];

                // Check the convergence criterion for alpha at r[ir] */
                diff1 = fabs((alpha[ir] - alpha_old) / alpha[ir]);
                if (diff1<0.001){
                    converged = 1;
                }
                // Check the deviation of alpha[ir] in this iteration compared to the previous big iteration */
                diff1 = fabs((alpha[ir] - alpha_iter[ir])/alpha[ir]);
                if (diff1>diff2){
                    diff2 = diff1;
                }

            }
        if (diff2<0.001)
            cg97_converged=1;
        }
    }
    
    // Now we have a CG97 disk, let's include the rim and the self-irradiation */
    
    bnu_t(nu, trim, bnut_rim);

    for (ir=iish;ir<DDN01_nrmod;ir++){
        converged = 0;
        delta = -1.;
        theta = -1.;

        while (converged!=1){
            delta = (r[ir] / rin - 1.0) / (Hs[ir] / Hsrim[0] - 1.);
            theta = atan( (r[ir]-rin) / (Hs[ir] - Hsrim[0]) );
                
            for (inu=0;inu<DDN01_nwav;inu++){
                 Firr[inu] = 0.;
            }

            if ((Hs[ir]/r[ir])<(Hsrim[0]/r[0])){
                iicg = ir+1;
            } else {
                for (inu=0;inu<DDN01_nwav;inu++){
                    Firr[inu] = Firr_star[inu][ir];
                }
            }

            alpha_old = alpha[ir];
            if(delta>=0){
                if (delta<1.){
                    for (inu=0;inu<DDN01_nwav;inu++){
                        Firr[inu] = Firr[inu] + 2. /M_PI * bnut_rim[inu] * pow(rin/r[ir], 2) * \
                                    cos(theta) * (delta * sqrt(1.0 - delta*delta) + asin(delta));
                        dum[inu] = Firr[inu] * alpha[ir];
                    }
                } else {
                    for (inu=0;inu<DDN01_nwav;inu++){
                        Firr[inu] = Firr[inu] + bnut_rim[inu] * pow(rin/r[ir], 2) * cos(theta);
                        dum[inu] = Firr[inu] * alpha[ir];
                    }
                }
                    
                Firr_int = integr(DDN01_nwav, nu, dum);
                
                ts[ir] = get_tsurf(r[ir], rstar, tstar, kappa, nu, Firr);
                ti[ir] = sqrt(sqrt(0.5*Firr_int)) * trim;
                ti[ir] = get_ti(alpha[ir], r[ir], rstar, tstar, kappa, nu, sigma[ir], ts[ir], ti[ir], Firr_int);
                chi[ir] = get_chicg(alpha[ir], sigma[ir], kappa_planck_star);
                Hp[ir] = sqrt(ti[ir]/tc) * sqrt(r[ir]/rstar) * r[ir];
                Hs[ir] = Hp[ir] * chi[ir];


                alpha[ir] = (gamma[ir]-1.) * Hs[ir]/r[ir];
                diff1 = fabs( (alpha[ir]-alpha_old) / alpha[ir]);
                if (diff1<1e-4)
                    converged = 1;
            } else {
                Hs[ir] = 0.;
                Hp[ir] = 0.;
                ti[ir] = 0.;
                ts[ir] = 0.;
                chi[ir] = 0.;
                iish = ir-1;
                converged = 1;
            }
        }

    }

    // Now we have a CG97 disk + puffed-up rim + self-irradiation, let's glue them together */

    for (ir=0;ir<iish+1;ir++){
        Hs[ir] = Hsrim[ir];
        Hp[ir] = Hprim[ir];
        chi[ir] = chirim[ir];
        ti[ir] = Hprim[ir]*Hprim[ir] * tc * rstar/pow(r[ir], 3);
        ts[ir] = ti[ir];
    }

    if (iish<iicg){
        if (iicg>0.){
            for (ir=iish;ir<iicg;ir++){
                diff1 = ( (Hs[ir+1]/r[ir+1]) - ( log(r[ir+1]) - log(r[ir]) ) / 8.) * r[ir];
                Hs[ir] = pow((pow(diff1, 8) + pow(Hs[ir],8)), 1./8.);
                chi[ir] = get_chicg(alpha[ir], sigma[ir], kappa_planck_star);
                Hp[ir] = Hs[ir] / chi[ir];
                ti[ir] = Hp[ir]*Hp[ir] * (GRAV_cgs*mstar*ML_MEAN_MOL_WT*AMU_cgs) / KBOLTZ_cgs / pow(r[ir],3);
                if (ti[ir]>ts[ir])
                    ts[ir] = ti[ir];
            }
         }
    }


    for (ir=0; ir<DDN01_nrmod; ir++){
        if (Hs[ir]/r[ir]<=0.04){
            Hs[ir] = 0.04*r[ir];
            Hp[ir] = 0.04*r[ir]/4.;
            ti[ir] = Hp[ir]*Hp[ir] * (GRAV_cgs*mstar*ML_MEAN_MOL_WT*AMU_cgs) / KBOLTZ_cgs / pow(r[ir],3);
            ts[ir] = ti[ir];
        }
        Hs[ir] = Hs[ir] * 1e-2;
        Hp[ir] = Hp[ir] * 1e-2;

    }

} 

/*....................................................................*/
int 
DDN01_onFinalizeConfiguration(void){

// I expect the input to be cgs and here I convert everything to SI for Lime
//  m_rin    = m_paramDouble["rin"]*1.496e+11;
//  m_rout   = m_paramDouble["rout"]*1.496e+11;
//  m_mdisk  = m_paramDouble["mdisk"]*1.99e+30;
//  m_plsig1 = m_paramDouble["plsig1"];
//  m_mstar  = m_paramDouble["Mstar"]*1.99e+30;
//  m_rstar  = m_paramDouble["Rstar"]*6.96e+8;
//  m_tstar  = m_paramDouble["Tstar"];
//  m_bgdens = m_paramDouble["bgdens"]*1e6;
//  m_dustopac_fname = (m_paramString["dustopac_fname"]).c_str();
  
  int i;

  if(getParamI("rin", &i)) return ML_UNRECOG_PARAM;
  m_rin    = modelDblPars[i]*1.496e+11;
  if(getParamI("rout", &i)) return ML_UNRECOG_PARAM;
  m_rout   = modelDblPars[i]*1.496e+11;
  if(getParamI("mdisk", &i)) return ML_UNRECOG_PARAM;
  m_mdisk  = modelDblPars[i]*1.99e+30;
  if(getParamI("plsig1", &i)) return ML_UNRECOG_PARAM;
  m_plsig1 = modelDblPars[i];
  if(getParamI("Mstar", &i)) return ML_UNRECOG_PARAM;
  m_mstar  = modelDblPars[i]*1.99e+30;
  if(getParamI("Rstar", &i)) return ML_UNRECOG_PARAM;
  m_rstar  = modelDblPars[i]*6.96e+8;
  if(getParamI("Tstar", &i)) return ML_UNRECOG_PARAM;
  m_tstar  = modelDblPars[i];
  if(getParamI("bgdens", &i)) return ML_UNRECOG_PARAM;
  m_bgdens = modelDblPars[i]*1e6;
  m_dustopac_fname = modelStrPar;

//  kappa_fname = "testopac_micron.dat";
//  read_kappa(m_dustopac_fname, wav, kappa);

  my_sig0 = get_sig0(m_rin, m_rout, m_mdisk, m_plsig1);
  get_disk_struct(m_rin*1e2, m_rout*1e2, m_mdisk*1e3, m_plsig1, m_dustopac_fname, m_mstar*1e3, m_rstar*1e2, m_tstar, my_r, my_ri, my_hp, my_hs, my_ti, my_ts);

  rootTwoPi = sqrt(2.*M_PI);

  return 0;
}

/*....................................................................*/
double
DDN01_t_dust(const double x, const double y, const double z){
  int id;
  double rc, temp;

  rc    = sqrt(x*x + y*y);

  if (rc > m_rin){
    if (rc < m_rout){
      id = gsl_interp_bsearch(my_ri, rc, 0, DDN01_nrmod-1);

      if (fabs(z)>my_hs[id]){
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
double
DDN01_density(const double x, const double y, const double z){
//#define pi 3.14159265358979323846264338328
//#define mmw 2.3       // Mean molecular weight
//#define mpsi 1.6726e-27  // Proton mass	
  int id;
  double rc, rho;


  rc    = sqrt(x*x + y*y);

  if (rc > m_rin){
    if (rc < m_rout){
      id = gsl_interp_bsearch(my_ri, rc, 0, DDN01_nrmod-1);
      rho = my_sig0 * pow((rc/m_rout), m_plsig1)/ (my_hp[id] * rootTwoPi) * exp(-0.5* (z/my_hp[id]) * (z/my_hp[id])) ;
    }
    else {
      rho   = 0.0;
    }    
  }
  else { 
    rho = 0.0;
  }

//  return (rho)/(ML_MEAN_MOL_WT*AMU) + m_bgdens;
  return rho*ML_oneOnMuMp + m_bgdens;
}

/*....................................................................*/
void
DDN01_velocity(const double x, const double y, const double z, double* v){
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

