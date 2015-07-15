/*
 *  lime.h
 *  LIME, The versatile 3D line modeling environment
 *
 *  Created by Christian Brinch on 13/11/06.
 *  Copyright 2006-2014, Christian Brinch,
 *  <brinch@nbi.dk>
 *  Niels Bohr institutet
 *  University of Copenhagen
 *	All rights reserved.
 *
 */

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

#define silent 0
#define DIM 3
#define VERSION	"1.5"
#define DEFAULT_NTHREADS 1
#ifndef NTHREADS /* Value passed from the LIME script */
#define NTHREADS DEFAULT_NTHREADS
#endif

/* Physical constants */
#define PI			3.14159265358979323846
#define SPI			1.77245385
#define CLIGHT	    2.997924562e8
#define HPLANCK	    6.626196e-34
#define KBOLTZ	    1.380622e-23
#define AMU			1.6605402e-27
#define HPIP		HPLANCK*CLIGHT/4.0/PI/SPI
#define HCKB		100.*HPLANCK*CLIGHT/KBOLTZ
#define PC			3.08568025e16
#define AU			1.49598e11
#define maxp		0.15
#define OtoP		3.
#define GRAV		6.67428e-11

/* Other constants */
#define NITERATIONS 	16
#define max_phot		10000		/* don't set this value higher unless you have enough memory. */
#define ininphot		9
#define minpop			1.e-6
#define eps				1.0e-30
#define TOL				1e-6
#define MAXITER			50
#define goal			50
#define fixset			1e-6
#define blendmask		1.e4
#define MAX_NSPECIES            100
#define N_RAN_PER_SEGMENT       3
#define FAST_EXP_MAX_TAYLOR	3
#define FAST_EXP_NUM_BITS	8


/* input parameters */
typedef struct {
  double radius,radiusSqu,minScale,minScaleSqu,tcmb,taylorCutoff;
  int ncell,sinkPoints,pIntensity,nImages,nSpecies,blend;
  char *outputfile, *binoutputfile, *inputfile;
  char *gridfile;
  char *pregrid;
  char *restart;
  char *dust;
  int sampling,collPart,lte_only,antialias,polarization,doPregrid,nThreads;
  char **moldatfile;
} inputPars;

/* Molecular data: shared attributes */
typedef struct {
  int nlev,nline,*ntrans,npart;
  int *lal,*lau,*lcl,*lcu;
  double *aeinst,*freq,*beinstu,*beinstl,*up,*down,*eterm,*gstat;
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
  double *up, *down;
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

typedef struct {
  int line1, line2;
  double deltav;
} blend;

typedef struct {double x,y, *intensity, *tau;} rayData;



/* Some functions */
void density(double,double,double,double *);
void temperature(double,double,double,double *);
void abundance(double,double,double,double *);
void doppler(double,double,double, double *);
void velocity(double,double,double,double *);
void magfield(double,double,double,double *);
void gasIIdust(double,double,double,double *);

/* More functions */

void   	binpopsout(inputPars *, struct grid *, molData *);
void   	buildGrid(inputPars *, struct grid *);
void    calcSourceFn(double dTau, const inputPars *par, double *remnantSnu, double *expDTau);
void	continuumSetup(int, image *, molData *, inputPars *, struct grid *);
void	distCalc(inputPars *, struct grid *);
void	fit_d1fi(double, double, double*);
void    fit_fi(double, double, double*);
void    fit_rr(double, double, double*);
void   	input(inputPars *, image *);
float  	invSqrt(float);
void    freeInput(inputPars *, image*, molData* m );
void   	freeGrid(const inputPars * par, const molData* m, struct grid * g);
void   	freePopulation(const inputPars * par, const molData* m, struct populations * pop);
double 	gaussline(double, double);
void    getArea(inputPars *, struct grid *, const gsl_rng *);
void    getjbar(int, molData *, struct grid *, inputPars *,gridPointData *,double *);
void    getMass(inputPars *, struct grid *, const gsl_rng *);
void   	getmatrix(int, gsl_matrix *, molData *, struct grid *, int, gridPointData *);
void	getclosest(double, double, double, long *, long *, double *, double *, double *);
void	getVelosplines(inputPars *, struct grid *);
void	getVelosplines_lin(inputPars *, struct grid *);
void	gridAlloc(inputPars *, struct grid **);
void   	kappa(molData *, struct grid *, inputPars *,int);
void	levelPops(molData *, inputPars *, struct grid *, int *);
void	line_plane_intersect(struct grid *, double *, int , int *, double *, double *, double);
void	lineBlend(molData *, inputPars *, blend **);
void    lineCount(int,molData *,int **, int **, int *);
void	LTE(inputPars *, struct grid *, molData *);
void   	molinit(molData *, inputPars *, struct grid *,int);
void    openSocket(inputPars *par, int);
void	qhull(inputPars *, struct grid *);
void  	photon(int, struct grid *, molData *, int, const gsl_rng *,inputPars *,blend *,gridPointData *,double *);
void	parseInput(inputPars *, image **, molData **);
double 	planckfunc(int, double, molData *, int);
int     pointEvaluation(inputPars *,double, double, double, double);
void   	popsin(inputPars *, struct grid **, molData **, int *);
void   	popsout(inputPars *, struct grid *, molData *);
void	predefinedGrid(inputPars *, struct grid *);
double 	ratranInput(char *, char *, double, double, double);
void   	raytrace(int, inputPars *, struct grid *, molData *, image *);
void	report(int, inputPars *, struct grid *);
void	smooth(inputPars *, struct grid *);
int     sortangles(double *, int, struct grid *, const gsl_rng *);
void	sourceFunc(double *, double *, double, molData *,double,struct grid *,int,int, int,int);
void    sourceFunc_line(double *,double *,molData *, double, struct grid *, int, int,int);
void    sourceFunc_cont(double *,double *, struct grid *, int, int,int);
void    sourceFunc_pol(double *, double *, double, molData *, double, struct grid *, int, int, int, double);
void   	stateq(int, struct grid *, molData *, int, inputPars *,gridPointData *,double *);
void	statistics(int, molData *, struct grid *, int *, double *, double *, int *);
void    stokesangles(double, double, double, double, double *);
void    traceray(rayData, int, int, inputPars *, struct grid *, molData *, image *, int, int *, int *, double);
void   	velocityspline(struct grid *, int, int, double, double, double*);
void   	velocityspline2(double *, double *, double, double, double, double*);
double 	veloproject(double *, double *);
void	writefits(int, inputPars *, molData *, image *);
void    write_VTK_unstructured_Points(inputPars *, struct grid *);
int	factorial(const int n);
double	taylor(const int maxOrder, const float x);
void	calcFastExpRange(const int maxTaylorOrder, const int maxNumBitsPerMantField, int *numMantissaFields, int *lowestExponent, int *numExponentsUsed);
void	calcTableEntries(const int maxTaylorOrder, const int maxNumBitsPerMantField);
inline double	FastExp(const float negarg);


/* Curses functions */

void 	greetings();
void 	greetings_parallel(int);
void	screenInfo();
void 	done(int);
void 	progressbar(double,int);
void 	progressbar2(int,double,double,double);
void	casaStyleProgressBar(const int,int);
void 	goodnight(int, char *);
void	quotemass(double);
void 	warning(char *);
void	bail_out(char *);
void    collpartmesg(char *, int);
void    collpartmesg2(char *, int);
void    collpartmesg3(int, int);



