/*
 *  model_Ul76.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#include "ml_recipes.h"

double m_mdot;
double m_mstar;
double m_rc;
double m_x;
double m_doppler_b;

#ifndef WITH_C99
typedef struct{
  double r,i;
} complexType; 

/*....................................................................*/
complexType
dcomp(const double real, const double imag){
  complexType result;

  result.r = real;
  result.i = imag;

  return result;
}

/*....................................................................*/
complexType
cx_conjugate(complexType a){
  return dcomp(a.r, -a.i);
}

/*....................................................................*/
complexType
cx_multiply(complexType a, complexType b){
  /* (ar+i*ai)*(br+i*bi) = (ar*br-ai*bi) + i*(ai*br+ar*bi) */
  complexType result;

  result.r = a.r*b.r - a.i*b.i;
  result.i = a.i*b.r + a.r*b.i;

  return result;
}

/*....................................................................*/
complexType
cx_divide(complexType a, complexType b){
  complexType bConj,bMagSquared,result;

  bConj = cx_conjugate(b);
  bMagSquared = cx_multiply(b, bConj);
  result = cx_multiply(a, bConj);
  result.r /= bMagSquared.r;
  result.i /= bMagSquared.r;

  return result;
}

/*....................................................................*/
complexType
cx_add(complexType a, complexType b){
  complexType result;

  result.r = a.r + b.r;
  result.i = a.i + b.i;

  return result;
}

/*....................................................................*/
complexType
cx_subtract(complexType a, complexType b){
  complexType result;

  result.r = a.r - b.r;
  result.i = a.i - b.i;

  return result;
}

/*....................................................................*/
complexType
cx_power(complexType a, complexType b){
  /*
	(ar+i*ai)^(br+i*bi) = exp([ln{r}+i*theta]*[br+i*bi])
	                    = exp([ln{r}*br-theta*bi]+i*[ln{r}*bi+theta*br])
where
	r^2 = ar*ar+ai*ai
and
	theta = atan2(ai,ar).
  */
  complexType result;
  double lnR,theta;

  lnR = log(a.r*a.r + a.i*a.i)*0.5;
  theta = atan2(a.i,a.r);
  result.r = cos(lnR*b.r - theta*b.i);
  result.i = sin(lnR*b.i + theta*b.r);

  return result;
}

/*....................................................................*/
// A function to calculate the mu0 = cos(theta0) parameter in the model
// of Ulrich 1976
double
getmu0(double r, double rc, double mu){
  complexType a, b, c, d, cc, qq;
  complexType dum1, dum2, dum3, isq3,dum1a,dum1b,dum1c,t1,t2,t3,ipos,ineg;
  complexType x1, x2, x3;
  double mu0;
//  double theta, stheta;

  a = dcomp(rc,0.0);//rc + 0.*I;
  b = dcomp(0.0, 0.0);//0.0 + 0.*I;
  c = dcomp((r-rc),0.0); 
  d = dcomp((r*mu), 0.0);
//  theta = acos(mu);
//  stheta = sin(theta);

//**** since b==0 the first 2 terms are pointless:
////  dum1 = 2.0 * b*b*b - 9.0 * a*b*c + 27.0 * a*a*d;
//  dum1 = cx_multiply(cx_multiply(dcomp(2.0,0.0),b),cx_multiply(b,b))
//       - cx_multiply(cx_multiply(dcomp(9.0,0.0),a),cx_multiply(b,c))
//       + cx_multiply(cx_multiply(dcomp(27.0,0.0),a),cx_multiply(a,d));
  dum1a = cx_multiply(cx_multiply(dcomp(2.0,0.0),b),cx_multiply(b,b));
  dum1b = cx_multiply(cx_multiply(dcomp(9.0,0.0),a),cx_multiply(b,c));
  dum1c = cx_multiply(cx_multiply(dcomp(27.0,0.0),a),cx_multiply(a,d));
  dum1 = cx_subtract(dum1a,dum1b);
  dum1 = cx_add(dum1,dum1c);

//**** since b==0 the first term is pointless:
//  dum2 = b*b - 3.0*a*c;
  dum2 = cx_subtract(cx_multiply(b,b),cx_multiply(dcomp(3.0,0.0),cx_multiply(a,c)));

//  dum3 = 3.0*a;
  dum3 = cx_multiply(dcomp(3.0,0.0),a);

  isq3 = dcomp(0.0, sqrt(3.0));

//  qq = sqrt(dum1*dum1 - 4.0 * dum2*dum2*dum2);
//  qq = cx_power(cx_subtract(cx_multiply(dum1,dum1),cx_multiply(cx_multiply(dcomp(4.0,0.0),dum2),cx_multiply(dum2,dum2)))), dcomp(0.5,0.0));
  complexType qqarg1,qqarg2;
  qqarg1 = cx_subtract(cx_multiply(dum1,dum1),cx_multiply(cx_multiply(dcomp(4.0,0.0),dum2),cx_multiply(dum2,dum2)));
  qqarg2 = dcomp(0.5,0.0);
  qq = cx_power(qqarg1, qqarg2);

//  cc = pow((0.5 * (qq + dum1)), 1.0/3.0);
  cc = cx_power((cx_multiply(dcomp(0.5,0.0),cx_add(qq,dum1))),dcomp(1.0/3.0,0.0));

////  x1 = -b / dum3 - cc / dum3 - (b*b - dum3*c) / (dum3*cc);
//  x1 = cx_divide(cx_multiply(dcomp(-1.0,0.0),b),dum3)
//     - cx_divide(cc,dum3)
//     - cx_divide(cx_add(cx_multiply(b,b),cx_multiply(dum3,c)),cx_multiply(dum3,cc));
  x1 = cx_subtract(cx_divide(cx_multiply(dcomp(-1.0,0.0),b),dum3),cx_divide(cc,dum3));
  x1 = cx_subtract(x1,cx_divide(cx_subtract(cx_multiply(b,b),cx_multiply(dum3,c)),cx_multiply(dum3,cc)));

  t1 = cx_divide(cx_multiply(dcomp(-1.0,0.0),b),dum3);
  t2 = cx_divide(cc,cx_multiply(dum3,dcomp(2.0,0.0)));
  t3 = cx_divide(cx_subtract(cx_multiply(b,b),cx_multiply(dum3,c)),cx_multiply(dcomp(2.0,0.0),cx_multiply(dum3,cc)));
  ipos = cx_add(dcomp(1.0,0.0),isq3);
  ineg = cx_subtract(dcomp(1.0,0.0),isq3);

////  x2 = -b / dum3 + cc * (1.+isq3) / (dum3*2.0) + (1.0 - isq3) * (b*b - dum3*c) / (2.*dum3*cc);
//  x2 = cx_divide(cx_multiply(dcomp(-1.0,0.0),b),dum3)
//     + cx_divide(cx_multiply(cc,cx_add(dcomp(1.0,0.0),isq3)),cx_multiply(dum3,dcomp(2.0,0.0)))
//     + cx_divide(cx_multiply(cx_subtract(dcomp(1.0,0.0),isq3),cx_subtract(cx_multiply(b,b),cx_multiply(dum3,c))),cx_multiply(dcomp(2.0,0.0),cx_multiply(dum3,cc)));
  x2 = cx_add(cx_add(t1,cx_multiply(t2,ipos)),cx_multiply(t3,ineg));

////  x3 = -b / dum3 + cc * (1.-isq3) / (dum3*2.0) + (1.0 + isq3) * (b*b - dum3*c) / (2.*dum3*cc);
//  x3 = cx_divide(cx_multiply(dcomp(-1.0,0.0),b),dum3)
//     + cx_divide(cx_multiply(cc,cx_subtract(dcomp(1.0,0.0),isq3)),cx_multiply(dum3,dcomp(2.0,0.0)))
//     + cx_divide(cx_multiply(cx_add(dcomp(1.0,0.0),isq3),cx_subtract(cx_multiply(b,b),cx_multiply(dum3,c))),cx_multiply(dcomp(2.0,0.0),cx_multiply(dum3,cc)));
  x3 = cx_add(cx_add(t1,cx_multiply(t2,ineg)),cx_multiply(t3,ipos));


/*	if (abs(imag(x1))==0.){
		if (sin(acos(x1))*stheta>0.){
			mu0 = creal(x1); 
		}
	}

	if (abs(imag(x2))==0.){
		if (sin(acos(x2))*stheta>0.){
			mu0 =creal(x2); 
		}
	}	
	if (abs(imag(x3))==0.){
		if (sin(acos(x3))*stheta>0.){
			mu0= creal(x3); 
		}
	}

	return mu0;
*/
	
  mu0 = 0;
//  if (imag(x1)<1e-5){
  if (x1.i<1e-5){
//    if (abs(real(x1))>mu0){
    if (abs(x1.r)>mu0){
//      mu0 = fabs(real(x1));
      mu0 = fabs(x1.r);
    }
  }

  if (x2.i<1e-5){
    if (abs(x2.r)>mu0){
      mu0 = fabs(x2.r);
    }
  }	
  if (x3.i<1e-5){
    if (abs(x3.r)>mu0){
      mu0 = fabs(x3.r);
    }
  }	
	
  return mu0;
}

#else
const double getmu0(double r, double rc, double mu){

  dcomp a, b, c, d, cc, qq;
  dcomp dum1, dum2, dum3, isq3;
  dcomp x1, x2, x3;
  double mu0;
  double theta, stheta;

  a = dcomp(rc,0.0);//rc + 0.*I;
  b = dcomp(0.0, 0.0);//0.0 + 0.*I;
  c = dcomp((r-rc),0.0); 
  d = dcomp((r*mu), 0.0);
  theta = acos(mu);
  stheta = sin(theta);	

//  dum1 = 2.0 * pow(b,3.0) - 9.0 * a * b * c + 27.0 * pow(a,2.0) * d;
  dum1 = 2.0 * b*b*b - 9.0 * a * b * c + 27.0 * a*a * d;
//  dum2 = pow(b,2.0) - 3.0*a*c;
  dum2 = b*b - 3.0*a*c;
  dum3 = 3.*a;
  isq3 = dcomp(0.0, sqrt(3.0));

//  qq = sqrt(pow(dum1,2.0) - 4.0 * pow(dum2,3.0));
  qq = sqrt(dum1*dum1 - 4.0 * dum2*dum2*dum2);
  cc = pow((0.5 * (qq + dum1)), (1.0/3.0));

//  x1 = -b / dum3 - cc / dum3 - (pow(b, 2.0) - dum3*c) / (dum3*cc);
  x1 = -b / dum3 - cc / dum3 - (b*b - dum3*c) / (dum3*cc);
//  x2 = -b / dum3 + cc * (1.+isq3) / (dum3*2.0) + (1.0 - isq3) * (pow(b,2.0) - dum3*c) / (2.*dum3*cc);
  x2 = -b / dum3 + cc * (1.+isq3) / (dum3*2.0) + (1.0 - isq3) * (b*b - dum3*c) / (2.*dum3*cc);
//  x3 = -b / dum3 + cc * (1.-isq3) / (dum3*2.0) + (1.0 + isq3) * (pow(b,2.0) - dum3*c) / (2.*dum3*cc);
  x3 = -b / dum3 + cc * (1.-isq3) / (dum3*2.0) + (1.0 + isq3) * (b*b - dum3*c) / (2.*dum3*cc);


/*	if (abs(imag(x1))==0.){
		if (sin(acos(x1))*stheta>0.){
			mu0 = creal(x1); 
		}
	}

	if (abs(imag(x2))==0.){
		if (sin(acos(x2))*stheta>0.){
			mu0 =creal(x2); 
		}
	}	
	if (abs(imag(x3))==0.){
		if (sin(acos(x3))*stheta>0.){
			mu0= creal(x3); 
		}
	}

	return mu0;
*/

  mu0 = 0;
  if (imag(x1)<1e-5){
    if (abs(real(x1))>mu0){
      mu0 = fabs(real(x1));
    }
  }

  if (imag(x2)<1e-5){
    if (abs(real(x2))>mu0){
      mu0 = fabs(real(x2));
    }
  }

  if (imag(x3)<1e-5){
    if (abs(real(x3))>mu0){
      mu0 = fabs(real(x3));
    }
  }
	
  return mu0;
}
#endif

/*....................................................................*/
int 
Ul76_onFinalizeConfiguration(void){
////  m_mdot      = m_paramDouble["mdot"]*6.3102486e+22;
//  m_mdot      = m_paramDouble["mdot"]*MSUN/YJULIAN; // Msun/yr to kg/s
//  m_mstar     = m_paramDouble["mstar"]*MSUN;
//  m_rc        = m_paramDouble["rc"]*AU;

  int i;

  if(getParamI("mdot", &i)) return ML_UNRECOG_PARAM;
  m_mdot      = modelDblPars[i]*MSUN/YJULIAN; // Msun/yr to kg/s
  if(getParamI("mstar", &i)) return ML_UNRECOG_PARAM;
  m_mstar     = modelDblPars[i]*MSUN;
  if(getParamI("rc", &i)) return ML_UNRECOG_PARAM;
  m_rc        = modelDblPars[i]*AU;

  return 0;
}


/*....................................................................*/
double 
Ul76_density(const double x, const double y, const double z){
//#define pi 3.14159265358979323846264338328
//#define gg 6.674e-11    // gravitational constant [cm^3 g^-1 s^-2]
//#define mmw 2.0        // Mean molecular weight
//#define mp 1.6726e-27  // Proton mass	
//#define AU 1.49598e11

  double r, mu, mu0, rho,oneOnRoot;
  double dum2, dum3;//dum1, 

  r   = sqrt(x*x + y*y + z*z) ;
  // For z=0 we have a divergence in the equations so we'll cheat here a bit
  // and add 1.0 to the z value. Since the usual scale is several AU  this doesn't do any harm
  // to the solution
  mu  = fabs(z+1.) / (r+1e-90);
  mu0 = getmu0(r/m_rc, 1.0, mu);

  oneOnRoot = 1.0/(8.0 * M_PI * sqrt(GRAV*m_mstar*r*r*r));
//  dum1 =  m_mdot / (8.0 * M_PI * sqrt(GRAV*m_mstar*r*r*r)) ;
//  dum1 =  m_mdot*oneOnRoot;
  dum2 = sqrt(1. + mu/mu0) ;
  dum3 = (0.5*mu/mu0 + mu0*mu0 * m_rc/r) ;

//  rho = m_mdot / (8.0 * M_PI * sqrt(GRAV*m_mstar*pow(r,3.))) / sqrt(1. + mu/mu0) / (0.5*mu/mu0 + mu0*mu0 * m_rc/r) / (Ul76_mmw*AMU);
  rho = m_mdot*oneOnRoot / dum2 / dum3 / (Ul76_mmw*AMU);

  return rho;
}

/*....................................................................*/
void
Ul76_velocity(const double x, const double y, const double z, double* v){
//#define pi 3.14159265358979323846264338328
//#define gg 6.674e-11    // gravitational constant [cm^3 g^-1 s^-2]
//#define AU		1.49598e11

  double r, rcyl, sinp, cosp, sint, cost, mu, mu0,root1,root2;
  double vr, vt, vp, dum_cost;

  r   = sqrt(x*x + y*y + z*z);
  rcyl = sqrt(x*x + y*y);

  sinp = y / (rcyl+1e-90);
  cosp = x / (rcyl+1e-90);
  sint = rcyl / (r+1e-90);
  cost = z    / (r+1e-90);
  dum_cost = (fabs(z)+1.0) / (r+1e-90); 

  mu  = cost;
  mu0 = getmu0(r/m_rc, 1.0, dum_cost);
  if (cost<0.){
    mu0 = -mu0;
  }

  root1 = sqrt(GRAV*m_mstar/r);
  root2 = sqrt(1.0 - mu*mu);
//  vr  = -sqrt(GRAV*m_mstar/r) * sqrt(1.0 + mu/mu0);
//  vt  =  sqrt(GRAV*m_mstar/r) * (mu0 - mu) * sqrt( (mu0+mu) / (mu0 * sqrt(1.0 - mu*mu)));
//  vp  =  sqrt(GRAV*m_mstar/r) * sqrt(1.0 - mu0*mu0) / sqrt(1.0 - mu*mu) * sqrt(1.0 - mu/mu0);
  vr  = -root1 * sqrt(1.0 + mu/mu0);
  vt  =  root1 * (mu0 - mu) * sqrt( (mu0+mu) / (mu0 * root2));
  vp  =  root1 * sqrt(1.0 - mu0*mu0) / root2 * sqrt(1.0 - mu/mu0);

  v[0] = (sint*cosp*vr + cost*cosp*vt - sinp*vp) ;
  v[1] = (sint*sinp*vr + cost*sinp*vt + cosp*vp) ;
  v[2] = (cost*vr - sint*vt) ;
}

