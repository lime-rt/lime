/*
 *  lime.h
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2016 The LIME development team
 *
 */

#ifndef LIME_H
#define LIME_H

#include "inpars.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_linalg.h>

#ifdef OLD_QHULL
#include <qhull/qhull_a.h>
#else
#include <libqhull/qhull_a.h>
#endif

#ifdef OLD_FITSIO
#include <cfitsio/fitsio.h>
#else
#include <fitsio.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_num_threads() 0
#define omp_get_thread_num() 0
#define omp_set_dynamic(int) 0
#endif

#define DIM 3
#define VERSION	"1.5"
#define DEFAULT_NTHREADS 1
#ifndef NTHREADS /* Value passed from the LIME script */
#define NTHREADS DEFAULT_NTHREADS
#endif

/* Physical constants */
// - NIST values as of 23 Sept 2015:
#define AMU             1.66053904e-27		// atomic mass unit             [kg]
#define CLIGHT          2.99792458e8		// speed of light in vacuum     [m / s]
#define HPLANCK         6.626070040e-34		// Planck constant              [J * s]
#define KBOLTZ          1.38064852e-23		// Boltzmann constant           [J / K]

// From IAU 2009:
#define GRAV            6.67428e-11		// gravitational constant       [m^3 / kg / s^2]
#define AU              1.495978707e11		// astronomical unit            [m]

// Derived:
#define PC              3.08567758e16		// parsec (~3600*180*AU/PI)     [m]
#define HPIP            8.918502221e-27		// HPLANCK*CLIGHT/4.0/PI/SPI
#define HCKB            1.43877735		// 100.*HPLANCK*CLIGHT/KBOLTZ

/* Other constants */
#define PI                      3.14159265358979323846	// pi
#define SPI                     1.77245385091		// sqrt(pi)
#define maxp                    0.15
#define OtoP                    3.
#define NITERATIONS             16
#define max_phot                10000		/* don't set this value higher unless you have enough memory. */
#define ininphot                9
#define minpop                  1.e-6
#define eps                     1.0e-30
#define TOL                     1e-6
#define MAXITER                 50
#define goal                    50
#define fixset                  1e-6
#define maxBlendDeltaV		1.e4		/* m/s */
#define MAX_NSPECIES            100
#define MAX_NIMAGES             100
#define N_RAN_PER_SEGMENT       3
#define FAST_EXP_MAX_TAYLOR     3
#define FAST_EXP_NUM_BITS       8
#define N_SMOOTH_ITERS          20
#define TYPICAL_ISM_DENS        1000.0

/* Collision partner ID numbers from LAMDA */
#define CP_H2			1
#define CP_p_H2			2
#define CP_o_H2			3
#define CP_e			4
#define CP_H			5
#define CP_He			6
#define CP_Hplus		7


typedef struct {
  double radius,minScale,tcmb;
  double radiusSqu,minScaleSqu,taylorCutoff;
  int sinkPoints,pIntensity,blend;
  int ncell,nImages,nSpecies,collPart,doPregrid;
  char *outputfile, *binoutputfile;
  char *gridfile;
  char *pregrid;
  char *restart;
  char *dust;
  int sampling,lte_only,init_lte,antialias,polarization,nThreads;
  char **moldatfile;
} configInfo;

/* Molecular data: shared attributes */
typedef struct {
  int nlev,nline,*ntrans,*ntemp,npart;
  int *lal,*lau,*lcl,*lcu;
  double *aeinst,*freq,*beinstu,*beinstl,*eterm,*gstat;
  double **down;
  double norm,norminv,*cmb,*local_cmb;
} molData;

/* Data concerning a single grid vertex which is passed from photon() to stateq(). This data needs to be thread-safe. */
typedef struct {
  double *jbar,*phot,*vfac;
} gridPointData;

typedef struct {
  double *intensity;
} surfRad;

/* Point coordinate */
typedef struct {
  double x[3];
  double xn[3];
} point;

struct rates {
  int t_binlow;
  double interp_coeff;
};


struct populations {
  double * pops, *knu, *dust;
  double dopb, binv;
  struct rates *partner;
};

/* Grid properties */
struct grid {
  int id;
  double x[3];
  double vel[3];
  double *a0,*a1,*a2,*a3,*a4;
  int numNeigh;
  point *dir;
  struct grid **neigh;
  double *w;
  int sink;
  int nphot;
  int conv;
  double *dens,t[2],*nmol,*abun, dopb;
  double *ds;
  struct populations* mol;
};

typedef struct {
  double *intense;
  double *tau;
  double stokes[3];
} spec;

/* Image information */
typedef struct {
  int doline;
  int nchan,trans;
  spec *pixel;
  double velres;
  double imgres;
  int pxls;
  int unit;
  double freq,bandwidth;
  char *filename;
  double source_vel;
  double theta,phi;
  double distance;
  double rotMat[3][3];
} image;

typedef struct {double x,y, *intensity, *tau;} rayData;

struct blend{
  int molJ, lineJ;
  double deltaV;
};

struct lineWithBlends{
  int lineI, numBlends;
  struct blend *blends;
};

struct molWithBlends{
  int molI, numLinesWithBlends;
  struct lineWithBlends *lines;
};

struct blendInfo{
  int numMolsWithBlends;
  struct molWithBlends *mols;
};

/* Some global variables */
int silent;

/* User-specifiable functions */
void density(double,double,double,double *);
void temperature(double,double,double,double *);
void abundance(double,double,double,double *);
void doppler(double,double,double, double *);
void velocity(double,double,double,double *);
void magfield(double,double,double,double *);
void gasIIdust(double,double,double,double *);
void gridDensity(configInfo,double,double,double,double*);

/* More functions */
void	run(inputPars, image *);

void	binpopsout(configInfo *, struct grid *, molData *);
void	buildGrid(configInfo *, struct grid *);
void	calcFastExpRange(const int, const int, int*, int*, int*);
void	calcSourceFn(double, const configInfo*, double*, double*);
void	calcTableEntries(const int, const int);
void	checkGridDensities(configInfo*, struct grid*);
void	continuumSetup(int, image*, molData*, configInfo*, struct grid*);
void	distCalc(configInfo*, struct grid*);
int	factorial(const int);
double	FastExp(const float);
void	fit_d1fi(double, double, double*);
void	fit_fi(double, double, double*);
void	fit_rr(double, double, double*);
void	freeGrid(const configInfo*, const molData*, struct grid*);
void	freeMoldata(const int, molData*);
void	freeMolsWithBlends(struct molWithBlends*, const int);
void	freePopulation(const configInfo*, const molData*, struct populations*);
double	gaussline(double, double);
void	getArea(configInfo *, struct grid *, const gsl_rng *);
void	getclosest(double, double, double, long *, long *, double *, double *, double *);
void	getjbar(int, molData*, struct grid*, const int, configInfo*, struct blendInfo, int, gridPointData*, double*);
void	getMass(configInfo *, struct grid *, const gsl_rng *);
void	getmatrix(int, gsl_matrix *, molData *, struct grid *, int, gridPointData *);
int	getNextEdge(double*, int, struct grid*, const gsl_rng*);
void	getVelosplines(configInfo *, struct grid *);
void	getVelosplines_lin(configInfo *, struct grid *);
void	gridAlloc(configInfo *, struct grid **);
void	input(inputPars *, image *);
float	invSqrt(float);
void	kappa(molData *, struct grid *, configInfo *,int);
void	levelPops(molData *, configInfo *, struct grid *, int *);
void	line_plane_intersect(struct grid *, double *, int , int *, double *, double *, double);
void	lineBlend(molData*, configInfo*, struct blendInfo*);
void	LTE(configInfo *, struct grid *, molData *);
void	lteOnePoint(configInfo*, molData*, const int, const double, double*);
void	molinit(molData *, configInfo *, struct grid *,int);
void	openSocket(char*);
void	parseInput(inputPars, configInfo*, image**, molData**);
void	photon(int, struct grid*, molData*, int, const gsl_rng*, configInfo*, const int, struct blendInfo, gridPointData*, double*);
double	planckfunc(int, double, molData *, int);
int	pointEvaluation(configInfo*, double, double, double, double);
void	popsin(configInfo *, struct grid **, molData **, int *);
void	popsout(configInfo *, struct grid *, molData *);
void	predefinedGrid(configInfo *, struct grid *);
void	qhull(configInfo *, struct grid *);
double	ratranInput(char *, char *, double, double, double);
void	raytrace(int, configInfo *, struct grid *, molData *, image *);
void	readUserInput(inputPars *, image **, int *, int *);
void	report(int, configInfo *, struct grid *);
void	setUpConfig(configInfo *, image **, molData **);
void	smooth(configInfo *, struct grid *);
void	sourceFunc(double*, double*, double, molData*, double, struct grid*, int, int, int, int);
void	sourceFunc_cont(double*, double*, struct grid*, int, int, int);
void	sourceFunc_line(double*, double*, molData*, double, struct grid*, int, int, int);
void	sourceFunc_pol(double*, double*, struct grid*, int, int, int, double (*rotMat)[3]);
void	stateq(int, struct grid*, molData*, const int, configInfo*, struct blendInfo, int, gridPointData*, double*, _Bool*);
void	statistics(int, molData *, struct grid *, int *, double *, double *, int *);
void	stokesangles(double, double, double, double (*rotMat)[3], double*);
double	taylor(const int, const float);
void	traceray(rayData, int, int, configInfo*, struct grid*, molData*, image*, double);
void	velocityspline(struct grid *, int, int, double, double, double*);
void	velocityspline2(double *, double *, double, double, double, double*);
double	veloproject(double *, double *);
void	writefits(int, configInfo *, molData *, image *);
void	write_VTK_unstructured_Points(configInfo *, struct grid *);


/* Curses functions */

void 	greetings();
void 	greetings_parallel(int);
void	screenInfo();
void 	done(int);
void 	progressbar(double,int);
void 	progressbar2(int,int,double,double,double);
void	casaStyleProgressBar(const int,int);
void 	goodnight(int, char *);
void	quotemass(double);
void 	warning(char *);
void	bail_out(char *);
void    collpartmesg(char *, int);
void    collpartmesg2(char *, int);
void    collpartmesg3(int, int);

#endif /* LIME_H */

