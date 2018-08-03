/*
 *  model_LiSh96.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

/* The c++ code was so monstrously wasteful it just has to be CB's work. Have speeded up some of the grosser time-wasters. */

#include "ml_recipes.h"
#include "../messages.h"

double m_cs,m_cs_cgs,rootGravCgs,degToRn;
unsigned m_h0Index;

/*....................................................................*/
int
LiSh96_onFinalizeConfiguration(void){
  int i;

  if(getParamI("cs", &i)) return ML_UNRECOG_PARAM;
  m_cs = modelDblPars[i];
  if(getParamI("h0", &i)) return ML_UNRECOG_PARAM;
  m_h0Index = (unsigned)modelIntPars[i];

  m_cs_cgs = m_cs*100.0;
  rootGravCgs = sqrt(GRAV_cgs);
  degToRn = M_PI/180.0;

  return 0;
}


// Helper functions:

/*....................................................................*/
double
calcProduct(double a[11], double par){
  double product;

  product = a[0]+par*(a[1]+par*(a[2]+par*(a[3]+par*(a[4]+par*(a[5]+par*(a[6]+par*(a[7]+par*(a[8]+par*(a[9]+par*a[10])))))))));

  return product;
}

/*....................................................................*/
void
FIT_FI(double theta, unsigned h0Index, double *val_FI){
  double a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, par,power;

  if (h0Index==0){
    // h0 = 0.125
    if (theta<=(5.*degToRn) || (theta>=175.*degToRn)) {
                                a0 = -5.7890867;
                                a1 = 205.39736;
                                a2 = -369.88535;
                                a3 = -1179228.9;
                                a4 = 1.0046149e+08;
                                a5 = -4.2354367e+09;
                                a6 = 1.0577018e+11;
                                a7 = -1.6309162e+12;
                                a8 = 1.5258949e+13;
                                a9 = -7.9475118e+13;
                                a10 = 1.7687771e+14;
    }else {
                                a0 = -3.5670285;
                                a1 = 17.991643;
                                a2 = -76.407211;
                                a3 = 223.00690;
                                a4 = -427.81213;
                                a5 = 537.19752;
                                a6 = -436.24204;
                                a7 = 220.25914;
                                a8 = -62.770856;
                                a9 = 7.7068707;
                                a10 = 0.;
    }
  }else if (h0Index==1) {
    // h0 = 0.5
    if (theta<=(5.*degToRn) || (theta>=175.*degToRn)) {
                                a0 = -5.4055998;
                                a1 = 116.98580;
                                a2 = 12582.815;
                                a3 = -2220227.1;
                                a4 = 1.5158678e+08;
                                a5 = -5.8474702e+09;
                                a6 = 1.3904732e+11;
                                a7 = -2.0779845e+12;
                                a8 = 1.9022528e+13;
                                a9 = -9.7496539e+13;
                                a10 = 2.1431890e+14;
    }else {
                                a0 = -3.4408879;
                                a1 = 18.526773;
                                a2 = -77.668760;
                                a3 = 226.27281;
                                a4 = -433.73861;
                                a5 = 544.38973;
                                a6 = -441.95413;
                                a7 = 223.09314;
                                a8 = -63.567150;
                                a9 = 7.8036012;
                                a10 = 0.;
    }
  }else if (h0Index==2) {
    // h0 = 0.75
    if (theta<=(5.*degToRn) || (theta>=175.*degToRn)) {
                                a0 = -4.9212098;
                                a1 = 31.161348;
                                a2 = 19056.104;
                                a3 = -2350688.8;
                                a4 = 1.4209935e+08;
                                a5 = -5.1251353e+09;
                                a6 = 1.1666610e+11;
                                a7 = -1.6905691e+12;
                                a8 = 1.5124313e+13;
                                a9 = -7.6149159e+13;
                                a10 = 1.6503404e+14;
    }else{
                                a0 = -3.3707918;
                                a1 = 18.500380;
                                a2 = -77.447786;
                                a3 = 225.67114;
                                a4 = -432.48958;
                                a5 = 542.71313;
                                a6 = -440.52205;
                                a7 = 222.34270;
                                a8 = -63.349783;
                                a9 = 7.7771164;
                                a10 = 0.;
    }
  }else if (h0Index==3) {
    // h0 = 1.0
    if (theta<=(5.*degToRn) || (theta>=175.*degToRn)) {
                                a0 = -4.5116379;
                                a1 = 1.2452873;
                                a2 = 15309.445;
                                a3 = -1576926.2;
                                a4 = 86017008.;
                                a5 = -2.8944098e+09;
                                a6 = 6.2668488e+10;
                                a7 = -8.7450044e+11;
                                a8 = 7.5975161e+12;
                                a9 = -3.7368812e+13;
                                a10 = 7.9460580e+13;
    }else{
                                a0 = -3.3057712;
                                a1 = 18.449954;
                                a2 = -77.137220;
                                a3 = 224.62869;
                                a4 = -430.27816;
                                a5 = 539.80764;
                                a6 = -438.12405;
                                a7 = 221.14260;
                                a8 = -63.020393;
                                a9 = 7.7392590;
                                a10 = 0.;
    }
  }else if (h0Index==4) {
    // h0 = 1.25
    if (theta<=(5.*degToRn) || (theta>=175.*degToRn)) {
                                a0 = -5.0316989;
                                a1 = 70.411691;
                                a2 = 17370.700;
                                a3 = -2472974.7;
                                a4 = 1.5853121e+08;
                                a5 = -5.9193516e+09;
                                a6 = 1.3791697e+11;
                                a7 = -2.0323001e+12;
                                a8 = 1.8413217e+13;
                                a9 = -9.3630244e+13;
                                a10 = 2.0453706e+14;
    }else {
                                a0 = -3.2558545;
                                a1 = 18.520204;
                                a2 = -77.581352;
                                a3 = 226.20610;
                                a4 = -433.70647;
                                a5 = 544.57657;
                                a6 = -442.35889;
                                a7 = 223.47259;
                                a8 = -63.743803;
                                a9 = 7.8357126;
                                a10 = 0.;
    }
  }else{
    bail_out("Unrecognized value of h0Index");
exit(1);
  }

  if (theta>(M_PI/2.)) {
                par=M_PI-theta;
  }else{
                par=theta;
  }

//        *val_FI = pow(10.,a0+a1*par+a2*pow(par,2.)+a3*pow(par,3.)+a4*pow(par,4.)+a5*pow(par,5.)
//                          +a6*pow(par,6.)+a7*pow(par,7.)+a8*pow(par,8.)+a9*pow(par,9.)+a10*pow(par,10.));
  power = a0+par*(a1+par*(a2+par*(a3+par*(a4+par*(a5+par*(a6+par*(a7+par*(a8+par*(a9+par*a10)))))))));
  *val_FI = pow(10.,power);

  if (par<1.57080e-6) {
    *val_FI = 0;
  }
}

/*....................................................................*/
void
FIT_D1FI(double theta, unsigned h0Index, double *val_D1FI){

  double b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, par,power;

  if (h0Index==0) {
    // h0 = 0.125
    if (theta<=(2.*degToRn) || (theta>=178.*degToRn)) {
                                b0 = -5.7978819;
                                b1 = 3420.8656;
                                b2 = -1447319.6;
                                b3 = 3.3474335e+08;
                                b4 = -4.6303963e+10;
                                b5 = 4.0455798e+12;
                                b6 = -2.2870203e+14;
                                b7 = 8.3419323e+15;
                                b8 = -1.8942095e+17;
                                b9 = 2.4341291e+18;
                                b10 = -1.3516589e+19;
    }else{
                                b0 = 0.0010432216;
                                b1 = 0.94998080;
                                b2 = -1.1636316;
                                b3 = 1.6589067;
                                b4 = -3.1972786;
                                b5 = 4.0764528;
                                b6 = -3.3552368;
                                b7 = 1.7329784;
                                b8 = -0.50606902;
                                b9 = 0.063502876;
                                b10 = 0.;
    }
  }else if (h0Index==1) {
    // h0 = 0.5
    if (theta<=(2.*degToRn) || (theta>=178.*degToRn)) {
                                b0 = -5.6841594;
                                b1 = 3421.0258;
                                b2 = -1446535.1;
                                b3 = 3.3455240e+08;
                                b4 = -4.6278010e+10;
                                b5 = 4.0433633e+12;
                                b6 = -2.2857917e+14;
                                b7 = 8.3375219e+15;
                                b8 = -1.8932207e+17;
                                b9 = 2.4328716e+18;
                                b10 = -1.3509665e+19;
    }else{
                                b0 = 1.7283135e-05;
                                b1 = 1.3289417;
                                b2 = 0.010731257;
                                b3 = -0.63967496;
                                b4 = 0.30219781;
                                b5 = -0.73875582;
                                b6 = 1.0052509;
                                b7 = -0.84865141;
                                b8 = 0.39641631;
                                b9 = -0.071769539;
                                b10 = 0.;
    }
  }else if (h0Index==2) {
    // h0 = 0.75
    if (theta<=(2.*degToRn) || (theta>=178.*degToRn)) {
                                b0 = -5.6159326;
                                b1 = 3421.0256;
                                b2 = -1446534.8;
                                b3 = 3.3455235e+08;
                                b4 = -4.6278001e+10;
                                b5 = 4.0433623e+12;
                                b6 = -2.2857911e+14;
                                b7 = 8.3375191e+15;
                                b8 = -1.8932200e+17;
                                b9 = 2.4328706e+18;
                                b10 = -1.3509659e+19;
    }else {
                                b0 = 0.00064503765;
                                b1 = 1.5327040;
                                b2 = 0.28541206;
                                b3 = -1.9447386;
                                b4 = 5.2895923;
                                b5 = -10.590988;
                                b6 = 12.311262;
                                b7 = -8.5059899;
                                b8 = 3.1409106;
                                b9 = -0.46616486;
                                b10 = 0.;
    }
  }else if (h0Index==3) {
    // h0 = 1.0
    if (theta<=(2.*degToRn) || (theta>=178.*degToRn)) {
                                b0 = -5.5545308;
                                b1 = 3421.0264;
                                b2 = -1446535.5;
                                b3 = 3.3455259e+08;
                                b4 = -4.6278041e+10;
                                b5 = 4.0433664e+12;
                                b6 = -2.2857936e+14;
                                b7 = 8.3375290e+15;
                                b8 = -1.8932224e+17;
                                b9 = 2.4328738e+18;
                                b10 = -1.3509677e+19;
    }else{
                                b0 = 0.0020473146;
                                b1 = 1.7204718;
                                b2 = 0.85401373;
                                b3 = -5.1359032;
                                b4 = 15.012772;
                                b5 = -27.400935;
                                b6 = 29.635387;
                                b7 = -18.744299;
                                b8 = 6.2062771;
                                b9 = -0.81039269;
                                b10 = 0.;
    }
  }else if (h0Index==4) {
                    // h0 = 1.25
    if (theta<=(2.*degToRn) || (theta>=178.*degToRn)) {
                                b0 = -5.4996702;
                                b1 = 3421.0261;
                                b2 = -1446535.3;
                                b3 = 3.3455251e+08;
                                b4 = -4.6278030e+10;
                                b5 = 4.0433654e+12;
                                b6 = -2.2857932e+14;
                                b7 = 8.3375278e+15;
                                b8 = -1.8932222e+17;
                                b9 = 2.4328737e+18;
                                b10 = -1.3509677e+19;
    }else {
                                b0 = 0.0027084211;
                                b1 = 1.9422695;
                                b2 = 1.0291320;
                                b3 = -5.7758578;
                                b4 = 15.435726;
                                b5 = -24.933720;
                                b6 = 22.899781;
                                b7 = -11.189803;
                                b8 = 2.1936396;
                            b9 = b10 = 0.;
    }
  }else{
    bail_out("Unrecognized value of h0Index");
exit(1);
  }

  if(theta > (M_PI/2.)) {
    par=M_PI-theta;
  }else{
    par=theta;
  }

  power = b0+par*(b1+par*(b2+par*(b3+par*(b4+par*(b5+par*(b6+par*(b7+par*(b8+par*(b9+par*b10)))))))));

  if (par<=(2.*degToRn)) {
//    *val_D1FI = pow(10.,b0+b1*par+b2*pow(par,2.)+b3*pow(par,3.)+b4*pow(par,4.)+b5*pow(par,5.)
//                                                +b6*pow(par,6.)+b7*pow(par,7.)+b8*pow(par,8.)+b9*pow(par,9.)+b10*pow(par,10.));
    *val_D1FI = pow(10.,power);
  }else{
//    *val_D1FI = b0+b1*par+b2*pow(par,2.)+b3*pow(par,3.)+b4*pow(par,4.)+b5*pow(par,5.)
//                +b6*pow(par,6.)+b7*pow(par,7.)+b8*pow(par,8.)+b9*pow(par,9.)+b10*pow(par,10.);
    *val_D1FI = power;
  }

  if (theta > (M_PI/2.)) {
    *val_D1FI = -(*val_D1FI);
  }

  if (par==M_PI/2. || par<1.57080e-6) {
    *val_D1FI = 0.;
  }
}

/*....................................................................*/
void
FIT_RR(double theta, unsigned h0Index, double *val_RR){

  double c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, par,power;

  if (h0Index==0) {
    // h0 = 0.125
    if (theta<=(5.*degToRn) || (theta>=175.*degToRn)) {
                                c0 = -1.3062647;
                                c1 = 51.351209;
                                c2 = -92.607231;
                                c3 = -294770.22;
                                c4 = 25113156.;
                                c5 = -1.0587786e+09;
                                c6 = 2.6440686e+10;
                                c7 = -4.0770175e+11;
                                c8 = 3.8144900e+12;
                                c9 = -1.9867518e+13;
                                c10 = 4.4216664e+13;
    }else{
                                c0 = -0.77212340;
                                c1 = 5.0361044;
                                c2 = -24.210442;
                                c3 = 83.164993;
                                c4 = -191.74154;
                                c5 = 298.26990;
                                c6 = -313.11445;
                                c7 = 218.19150;
                                c8 = -96.580345;
                                c9 = 24.555825;
                                c10 = -2.7286711;
    }
  }else if (h0Index==1) {
    // h0 = 0.5
    if (theta<=(5.*degToRn) || (theta>=175.*degToRn)) {
                                c0 = -5.2701139;
                                c1 = 116.98355;
                                c2 = 12583.584;
                                c3 = -2220272.5;
                                c4 = 1.5158944e+08;
                                c5 = -5.8475666e+09;
                                c6 = 1.3904956e+11;
                                c7 = -2.0780177e+12;
                                c8 = 1.9022830e+13;
                                c9 = -9.7498099e+13;
                                c10 = 2.1432234e+14     ;
    }else {
                                c0 = -3.3053204;
                                c1 = 18.524631;
                                c2 = -77.322550;
                                c3 = 226.17076;
                                c4 = -433.47829;
                                c5 = 543.86354;
                                c6 = -441.38772;
                                c7 = 222.71068;
                                c8 = -63.431596;
                                c9 = 7.7849548;
                                c10 = 0.;
    }
  }else if (h0Index==2) {
    // h0 = 0.75
    if (theta<=(5.*degToRn) || (theta>=175.*degToRn)) {
                                c0 = -7.4070974;
                                c1 = 46.768169;
                                c2 = 28579.319;
                                c3 = -3525500.9;
                                c4 = 2.1311782e+08;
                                c5 = -7.6865772e+09;
                                c6 = 1.7497333e+11;
                                c7 = -2.5354761e+12;
                                c8 = 2.2683060e+13;
                                c9 = -1.1420644e+14;
                                c10 = 2.4751332e+14;
    }else {
                                c0 = -5.0810940;
                                c1 = 27.743121;
                                c2 = -115.72260;
                                c3 = 338.18057;
                                c4 = -647.83734;
                                c5 = 812.61889;
                                c6 = -659.38952;
                                c7 = 332.72674;
                                c8 = -94.809204;
                                c9 = 11.645955;
                                c10 = 0.;
    }
  }else if (h0Index==3) {
    // h0 = 1.0
    if (theta<=(5.*degToRn) || (theta>=175.*degToRn)) {
                                c0 = -9.2899541;
                                c1 = 2.5580199;
                                c2 = 30605.924;
                                c3 = -3152489.6;
                                c4 = 1.7195433e+08;
                                c5 = -5.7859487e+09;
                                c6 = 1.2527122e+11;
                                c7 = -1.7480406e+12;
                                c8 = 1.5186368e+13;
                                c9 = -7.4693720e+13;
                                c10 = 1.5882537e+14;
    }else {
                                c0 = -7.0495662;
                                c1 = 41.213965;
                                c2 = -196.52842;
                                c3 = 671.34639;
                                c4 = -1544.3087;
                                c5 = 2400.9562;
                                c6 = -2520.1512;
                                c7 = 1756.2736;
                                c8 = -777.48434;
                                c9 = 197.65936;
                                c10 = -21.953171;
    }
  }else if (h0Index==4) {
    // h0 = 1.25
    if (theta<=(5.*degToRn) || (theta>=175.*degToRn)) {
                                c0 = -13.149963;
                                c1 = 176.19380;
                                c2 = 43394.895;
                                c3 = -6179185.7;
                                c4 = 3.9613946e+08;
                                c5 = -1.4791622e+10;
                                c6 = 3.4463814e+11;
                                c7 = -5.0785019e+12;
                                c8 = 4.6012797e+13;
                                c9 = -2.3397308e+14;
                                c10 = 5.1111916e+14;
    }else {
                                c0 = -8.7113248;
                                c1 = 46.333987;
                                c2 = -193.78955;
                                c3 = 567.16615;
                                c4 = -1089.1455;
                                c5 = 1370.3674;
                                c6 = -1115.9978;
                                c7 = 565.59125;
                                c8 = -161.97286;
                                c9 = 20.001390;
                                c10 = 0.;
    }
  }else{
    bail_out("Unrecognized value of h0Index");
exit(1);
  }

  if (theta>(M_PI/2.)) {
    par=M_PI-theta;
  } else {
    par=theta;
  }

  power = c0+par*(c1+par*(c2+par*(c3+par*(c4+par*(c5+par*(c6+par*(c7+par*(c8+par*(c9+par*c10)))))))));

//  *val_RR = pow(10.,c0+c1*par+c2*pow(par,2.)+c3*pow(par,3.)+c4*pow(par,4.)+c5*pow(par,5.)
//                          +c6*pow(par,6.)+c7*pow(par,7.)+c8*pow(par,8.)+c9*pow(par,9.)+c10*pow(par,10.));
  *val_RR = pow(10.,power);

  if (par<1.57080e-6) {
    *val_RR = 0;
  }
}


/*....................................................................*/
double
LiSh96_density(const double x, const double y, const double z){
  double r, theta;
  double val_RR;

  r=sqrt(x*x+y*y+z*z);
  theta = acos(z/r);

  FIT_RR(theta,m_h0Index,&val_RR);

//  return pow(m_cs,2.)*val_RR/(2.*M_PI*GRAV*pow(r,2.))/(2.3*AMU);
  return m_cs*m_cs*val_RR/(2.0*M_PI*GRAV*r*r)/(ML_MEAN_MOL_WT*AMU);
}

/*....................................................................*/
void
LiSh96_bmag(const double x, const double y, const double z, double* B){
  double r, theta, phi,tempVal;
  double val_FI,val_D1FI;
  double Br, Bt, Bf,sinTheta,cosTheta,sinPhi,cosPhi;

  r=sqrt(x*x+y*y+z*z);
  theta = acos(z/r);
  phi   = atan2(y,x);

  FIT_FI(theta,m_h0Index, &val_FI);
  FIT_D1FI(theta,m_h0Index, &val_D1FI);

//  Br = 2.*pow(m_cs*1e2,2.)*val_D1FI/(sqrt(GRAV*1e3)*r*1e2*sin(theta))/1e4;
//  Bt = -2.*pow(m_cs*1e2,2.)*val_FI/(sqrt(GRAV*1e3)*r*1e2*sin(theta))/1e4;
//  Br =  2.0*m_cs_cgs*m_cs_cgs*val_D1FI/(rootGravCgs*r*1e2*sin(theta))/1e4;
//  Bt = -2.0*m_cs_cgs*m_cs_cgs*val_FI  /(rootGravCgs*r*1e2*sin(theta))/1e4;
  sinTheta = sin(theta);
  cosTheta = z/r;//cos(theta);
  sinPhi = sin(phi);
  cosPhi = cos(phi);
  tempVal = 2.0*m_cs_cgs*m_cs_cgs/(rootGravCgs*r*1e2*sinTheta)/1e4;
  Br =  tempVal*val_D1FI;
  Bt = -tempVal*val_FI;
  Bf = 0.;

//  B[0] = Br*sin(theta)*cos(phi) + Bt*cos(theta)*cos(phi) - Bf*sin(phi);
//  B[1] = Br*sin(theta)*sin(phi) + Bt*cos(theta)*sin(phi) + Bf*cos(phi);
//  B[2] = Br*cos(theta)          - Bt*sin(theta);
  B[0] = Br*sinTheta*cosPhi + Bt*cosTheta*cosPhi - Bf*sinPhi;
  B[1] = Br*sinTheta*sinPhi + Bt*cosTheta*sinPhi + Bf*cosPhi;
  B[2] = Br*cosTheta        - Bt*sinTheta;
}

/*....................................................................*/
void
LiSh96_velocity(const double x, const double y, const double z, double* v){
  v[0] = 0.0;
  v[1] = 0.0;
  v[2] = 0.0;	
}

/*....................................................................*/
double
LiSh96_temperature(const double x, const double y, const double z){
return m_cs*m_cs*ML_MEAN_MOL_WT*AMU/KBOLTZ;
}

