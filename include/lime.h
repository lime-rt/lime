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

#define DIM 3
#define VERSION "0.1"

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
#ifndef DEBUG
#define NITERATIONS 	16
#else
#define NITERATIONS 	2
#endif

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

#define MAX_LINE 512      /*!< Maximum number of characters in each input file line */

int silent;
typedef char filename_t [MAX_LINE];
typedef char pythonname_t [MAX_LINE];

/* input parameters */
typedef struct {
  double radius,minScale,tcmb;
  int ncell,sinkPoints,pIntensity,nImages,nSpecies,blend;
  filename_t outputfile;
  filename_t binoutputfile;
  filename_t gridfile;
  filename_t pregrid;
  filename_t restart;
  filename_t dust;
  int sampling,collPart,lte_only,antialias,polarization;
  filename_t* moldatfile;
  pythonname_t python_module_name;
  filename_t   python_module_path;
  pythonname_t density_func_name;
  pythonname_t velocity_func_name;
  pythonname_t temperature_func_name;
  pythonname_t doppler_func_name;
  pythonname_t abundance_func_name;
} inputPars;

/* Molecular data and radiation field */
typedef struct {
  int nlev,nline,*ntrans,npart;
  int *lal,*lau,*lcl,*lcu;
  double *aeinst,*freq,*beinstu,*beinstl,*up,*down,*eterm,*gstat;
  double *jbar,norm,norminv,*cmb,*local_cmb;
  double *phot, *ds, *vfac;
} molData;

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

/* Grid proporties */
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
  filename_t filename;
  double source_vel;
  double theta,phi;
  double distance;
} image;

typedef struct {
  int line1, line2;
  double deltav;
} blend;


extern void* density_eval;

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
void	continuumSetup(int, image *, molData *, inputPars *, struct grid *);
void	distCalc(inputPars *, struct grid *);
void	fit_d1fi(double, double, double*);
void    fit_fi(double, double, double*);
void    fit_rr(double, double, double*);
int   	input(char* input_file, inputPars *, image *);
float  	invSqrt(float);
void    freeInput(inputPars *, image*, molData* m );
void   	freeGrid(const inputPars * par, const molData* m, struct grid * g);
void   	freePopulation(const inputPars * par, const molData* m, struct populations * pop);
double 	gaussline(double, double);
void    getArea(inputPars *, struct grid *, const gsl_rng *);
void    getjbar(int, molData *, struct grid *, inputPars *);
void    getMass(inputPars *, struct grid *, const gsl_rng *);
void   	getmatrix(int, gsl_matrix *, molData *, struct grid *, int);
void	getclosest(double, double, double, long *, long *, double *, double *, double *);
void	getVelosplines(inputPars *, struct grid *);
void	getVelosplines_lin(inputPars *, struct grid *);
void	gridAlloc(inputPars *, struct grid **);
void   	kappa(molData *, struct grid *, inputPars *,int);
void	levelPops(molData *, inputPars *, struct grid *, int *);
void	line_plane_intersect(struct grid *, double *, int , int *, double *, double *);
void	lineBlend(molData *, inputPars *, blend **);
void    lineCount(int,molData *,int **, int **, int *);
void	LTE(inputPars *, struct grid *, molData *);
void   	molinit(molData *, inputPars *, struct grid *,int);
void    openSocket(inputPars *par, int);
void	qhull(inputPars *, struct grid *);
void  	photon(int, struct grid *, molData *, int, const gsl_rng *,inputPars *,blend *);
void	parseInput(char* input_file, inputPars *, image **, molData **);
double 	planckfunc(int, double, molData *, int);
int		pointEvaluation(inputPars *,double, double, double, double);
void   	popsin(inputPars *, struct grid **, molData **, int *);
void   	popsout(inputPars *, struct grid *, molData *);
void	predefinedGrid(inputPars *, struct grid *);
int     python_call_initialize( const inputPars *par );
void    python_call_finalize();
double 	ratranInput(char *, char *, double, double, double);
void   	raytrace(int, inputPars *, struct grid *, molData *, image *);
void	report(int, inputPars *, struct grid *);
void	smooth(inputPars *, struct grid *);
int		sortangles(double *, int, struct grid *, const gsl_rng *);
void	sourceFunc(double *, double *, double, molData *,double,struct grid *,int,int, int,int);
void    sourceFunc_line(double *,double *,molData *, double, struct grid *, int, int,int);
void    sourceFunc_cont(double *,double *, struct grid *, int, int,int);
void    sourceFunc_pol(double *, double *, double, molData *, double, struct grid *, int, int, int, double);
void   	stateq(int, struct grid *, molData *, double *, int, inputPars *);
void	statistics(int, molData *, struct grid *, int *, double *, double *, int *);
void    stokesangles(double, double, double, double, double *);
void   	velocityspline(struct grid *, int, int, double, double, double*);
void   	velocityspline2(double *, double *, double, double, double, double*);
double 	veloproject(double *, double *);
void	writefits(int, inputPars *, molData *, image *);
void    write_VTK_unstructured_Points(inputPars *, struct grid *);


/* Curses functions */

void 	greetings();
void	screenInfo();
void 	done(int);
void 	progressbar(double,int);
void 	progressbar2(int,double,double,double);
void 	goodnight(int, char *);
void	quotemass(double);
void 	warning(char *);
void	bail_out(char *);
void    collpartmesg(char *, int);
void    collpartmesg2(char *, int);
void    collpartmesg3(int, int);



