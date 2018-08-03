/*
 *  model_Ma88.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#include "ml_recipes.h"

double m_mdot;
double m_rin;
double m_ve;
double m_tin;
double m_ab0;
double m_alpha;
double m_rhalf;



double
bilinearInterpol(const double x, const double y, const int xiSize, const int yiSize, double xgrid[], double ygrid[], double zgrid[][yiSize]){
	// Do a bi-linear interpolation to determine the value z at (x, y)
	// with zgrid given at (xgrid, ygrid)
	//
	// For values of (x, y) outside the grid, an extrapolation is
	// performed with the interpolation parameters of the related
	// grid cell.

  // Find i with xgrid[i] <= x < xgrid[i+1]:
  int i = -1;
  int imax = xiSize - 1;
  while (i < imax && x >= xgrid[i+1])
    i++;

  // Same with j and y:
  int j = -1;
//  int jmax = ygrid.size() - 1;
  int jmax = yiSize - 1;
  while (j < jmax && y >= ygrid[j+1])
    j++;

  // Define grid rectangle:
//  double il, iu;
  int il, iu;
  if (i < 0) {
		// Extrapolation:
		il = 0;
		iu = 1;
  } else if (i == imax) {
		// Extrapolation:
		il = imax - 1;
		iu = imax;
  } else {
		// Interpolation:
		il = i;
		iu = i + 1;
  }

//  double jl, ju;
  int jl, ju;
  if (j < 0) {
		// Extrapolation:
		jl = 0;
		ju = 1;
  } else if (j == jmax) {
		// Extrapolation:
		jl = jmax - 1;
		ju = jmax;
  } else {
		// Interpolation:
		jl = j;
		ju = j + 1;
  }

  // Interpolate or extrapolate:
  double t = (x - xgrid[il]) / (xgrid[iu] - xgrid[il]);
  double u = (y - ygrid[jl]) / (ygrid[ju] - ygrid[jl]);

  double z =   (1. - t) * (1. - u) * zgrid[il][jl]\
	           + t * (1. - u) * zgrid[iu][jl]\
	           + t * u * zgrid[iu][ju]\
	           + (1. - t) * u * zgrid[il][ju];
  return z;
}

//vector<double>
//logMdotTable(){
void
logMdotTable(double logMdots[13]){
  int i;

	// Values of Mdot as in table 3 of Mamon et al. 1988
	// log(Mdot / (Msun/yr)) is returned

	double mdots[13] =\
	{1.e-8, 2.e-8, 5.e-8, 1.e-7, 2.e-7, 5.e-7, 1.e-6, 2.e-6, 5.e-6, 1.e-5, 2.e-5, 5.e-5, 1.e-4};
	// Unit is Msun/year

//	vector<double> logMdots(13, 0.);
	for (i=0; i<13; i++)
		logMdots[i] = log10(mdots[i]);

//	return logMdots;
}

//vector<double>
//veTable(){
void
veTable(double ves[3]){
	// Values of Mve as in table 3 of Mamon et al. 1988
	// ve / (km/s) is returned

//	vector<double> ves(3, 0.);
	ves[0] = 7.5;  // km/s
	ves[1] = 15.;
	ves[2] = 30.;

//	return ves;
}

double
getAlpha(const double mdot, const double ve){
  const int xiSize=13,yiSize=3;
  double logMdots[xiSize],ves[yiSize],result;

	// mdot in Msun/yr
	// ve in km/s !!!

	// Interpolation is done for alpha
	// on logMdot-ve grid

//  double alphaTable[xiSize][yiSize] = {
  double alphaTable[13][3] = {\
    {1.71, 1.39, 1.20},\
    {1.81, 1.46, 1.23},\
    {1.96, 1.60, 1.31},\
    {2.09, 1.74, 1.39},\
    {2.22, 1.89, 1.51},\
    {2.38, 2.09, 1.71},\
    {2.51, 2.24, 1.88},\
    {2.66, 2.39, 2.05},\
    {2.90, 2.61, 2.29},\
    {3.07, 2.79, 2.47},\
    {3.26, 2.96, 2.66},\
    {3.51, 3.20, 2.89},\
    {3.71, 3.39, 3.07}\
  };

//  vector< vector<double> > vAlphaTable;
//  vAlphaTable.resize(13, vector<double>(3, 0.));

//  for (int i=0; i<13; i++)
//		for (int j=0; j<3; j++)
//			vAlphaTable[i][j] = alphaTable[i][j];

//	return bilinearInterpol(log10(mdot), ve , logMdotTable(), veTable(), vAlphaTable);

  logMdotTable(logMdots);
  veTable(ves);
//  result = bilinearInterpol(log10(mdot), ve, logMdots, ves, alphaTable, 13, 3);
  result = bilinearInterpol(log10(mdot), ve, xiSize, yiSize, logMdots, ves, alphaTable);

  return result;
}

double
getRhalf(const double mdot, const double ve){
  const int xiSize=13,yiSize=3;
  int i,j;
  double logMdots[xiSize],ves[yiSize],logRhalf;

	// mdot in Msun/yr
	// ve in km/s !!!

	// Unit of output is m

	// Interpolation is done for log(r1/2)
	// on logMdot-ve grid

//  double rhalfTable[xiSize][yiSize] = {
  double rhalfTable[13][3] = {\
    {7.50e15, 9.01e15, 1.39e16},\
    {9.79e15, 1.05e16, 1.48e16},\
    {1.49e16, 1.40e16, 1.71e16},\
    {2.12e16, 1.85e16, 2.01e16},\
    {3.10e16, 2.54e16, 2.49e16},\
    {5.23e16, 4.05e16, 3.55e16},\
    {7.91e16, 5.95e16, 4.88e16},\
    {1.21e17, 8.88e16, 6.94e16},\
    {2.14e17, 1.54e17, 1.15e17},\
    {3.35e17, 2.35e17, 1.72e17},\
    {5.31e17, 3.65e17, 2.61e17},\
    {9.99e17, 6.67e17, 4.63e17},\
    {1.64e18, 1.07e18, 7.26e17}\
  };
  // Unit is cm

  double logRhalfTable[xiSize][yiSize];

//	vector< vector<double> > logRhalfTable;
//	logRhalfTable.resize(13, vector<double>(3, 0.));

  for (i=0; i<xiSize; i++){
    for (j=0; j<yiSize; j++) {
      logRhalfTable[i][j] = log10(rhalfTable[i][j] / 100.); // Convert to m
    }
  }

//	double logRhalf = bilinearInterpol(log10(mdot), ve , logMdotTable(), veTable(), logRhalfTable);

  logMdotTable(logMdots);
  veTable(ves);
//  logRhalf = bilinearInterpol(log10(mdot), ve, logMdots, ves, logRhalfTable, 13, 3);
  logRhalf = bilinearInterpol(log10(mdot), ve, xiSize, yiSize, logMdots, ves, logRhalfTable);

  return pow(10., logRhalf);
}



/*....................................................................*/
int
Ma88_onFinalizeConfiguration(void){
////  m_mdot = m_paramDouble["mdot"] * 6.3102486e+22; // Msun/yr to kg/s
//  m_mdot = m_paramDouble["mdot"] * MSUN/YJULIAN; // Msun/yr to kg/s
//  m_rin  = m_paramDouble["rin"]  * AU;
//  m_ve   = m_paramDouble["ve"]  * 1.e3;
//  m_tin  = m_paramDouble["tin"];
//  m_ab0  = m_paramDouble["ab0"];

  int i;
  double mDotRaw,veRaw;

  if(getParamI("mdot", &i)) return ML_UNRECOG_PARAM;
  mDotRaw = modelDblPars[i];
  m_mdot = mDotRaw * MSUN/YJULIAN; // Msun/yr to kg/s
  if(getParamI("rin", &i)) return ML_UNRECOG_PARAM;
  m_rin  = modelDblPars[i] * AU;
  if(getParamI("ve", &i)) return ML_UNRECOG_PARAM;
  veRaw = modelDblPars[i];
  m_ve   = veRaw * 1.e3;
  if(getParamI("tin", &i)) return ML_UNRECOG_PARAM;
  m_tin  = modelDblPars[i];
  if(getParamI("ab0", &i)) return ML_UNRECOG_PARAM;
  m_ab0  = modelDblPars[i];

  // Alpha and r1/2 are interpolated from table 3 of Mamon et al. 1988
  // The interpolation routines take mdot in Msun/yr and ve in km/s
//  m_alpha = getAlpha(m_paramDouble["mdot"], m_paramDouble["ve"]);
//  m_rhalf= getRhalf(m_paramDouble["mdot"], m_paramDouble["ve"]);

  m_alpha = getAlpha(mDotRaw, veRaw);
  m_rhalf = getRhalf(mDotRaw, veRaw);

  // cout << "RS_DEBUG: alpha = " << m_alpha << " r1/2 = " << m_rhalf << endl;

  return 0;
}

/*....................................................................*/
double
Ma88_density(const double x, const double y, const double z){
  double r = sqrt(x*x + y*y + z*z);

  double density = 0.;
  if(r > m_rin){
    density = m_mdot / (4.0 * M_PI * r * r * m_ve);
  }

  // convert to number density:
  // TODO: set correct molar mass for each collision partner
  density /= 2.0 * AMU;

  return density;
}

/*....................................................................*/
double
Ma88_temperature(const double x, const double y, const double z){
  double r =sqrt(x*x+y*y+z*z);

  double temperature = 0.;
  if(r > m_rin && r<6000 * AU){
    temperature = m_tin * pow(r / m_rin, -0.54);
  }

  if(r > 6000 * AU){
    temperature = m_tin * pow((6000. * AU) / m_rin, -0.54) * pow(r / (6000 * AU), -0.72);
  }

  return temperature;
}

/*....................................................................*/
double
Ma88_abundance(const double x, const double y, const double z){
  double r = sqrt(x*x+y*y+z*z);

  double t1 = pow((r / m_rhalf), m_alpha);
  return m_ab0 * exp(-1 * log(2) * t1);
}

/*....................................................................*/
void
Ma88_velocity(const double x, const double y, const double z, double* v){
  double theta=atan2(sqrt(x*x+y*y),z);
  double phi=atan2(y,x);

  //Vector transformation back into Cartesian basis
  v[0] = m_ve*sin(theta)*cos(phi);
  v[1] = m_ve*sin(theta)*sin(phi);
  v[2] = m_ve*cos(theta);
}

