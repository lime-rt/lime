/*
 *  lime.h
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2017 The LIME development team
 *
 */

#ifndef LIME_H
#define LIME_H

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

#include "dims.h"

#define VERSION	"1.7.3"
#define DEFAULT_NTHREADS 1
#ifndef NTHREADS /* Value passed from the LIME script */
#define NTHREADS DEFAULT_NTHREADS
#endif

/* Physical constants */
/* - NIST values as of 23 Sept 2015: */
#define AMU             1.66053904e-27		/* atomic mass unit             [kg]	*/
#define CLIGHT          2.99792458e8		/* speed of light in vacuum     [m / s]	*/
#define HPLANCK         6.626070040e-34		/* Planck constant              [J * s]	*/
#define KBOLTZ          1.38064852e-23		/* Boltzmann constant           [J / K]	*/

/* From IAU 2009: */
#define GRAV            6.67428e-11		/* gravitational constant       [m^3 / kg / s^2]	*/
#define AU              1.495978707e11		/* astronomical unit            [m]	*/

/* Derived: */
#define PC              3.08567758e16		/* parsec (~3600*180*AU/PI)     [m]	*/
#define HPIP            8.918502221e-27		/* HPLANCK*CLIGHT/4.0/PI/SPI	*/
#define HCKB            1.43877735		/* 100.*HPLANCK*CLIGHT/KBOLTZ	*/

/* Other constants */
#define PI                      3.14159265358979323846	/* pi	*/
#define SPI                     1.77245385091		/* sqrt(pi)	*/
#define maxp                    0.15
#define NITERATIONS             16
#define MAX_RAYS_PER_POINT      10000
#define RAYS_PER_POINT          200
#define minpop                  1.e-6
#define eps                     1.0e-30
#define IMG_MIN_ALLOWED         1.0e-30
#define TOL                     1e-6
#define MAXITER                 50
#define goal                    50
#define fixset                  1e-6
#define maxBlendDeltaV          1.e4		/* m/s */
#define MAX_NSPECIES            100
#define MAX_NIMAGES             100
#define N_RAN_PER_SEGMENT       3
#define FAST_EXP_MAX_TAYLOR     3
#define FAST_EXP_NUM_BITS       8
#define NUM_GRID_STAGES         4
#define MAX_N_COLL_PART         20
#define N_SMOOTH_ITERS          20
#define TYPICAL_ISM_DENS        1000.0
#define STR_LEN_0               80
#define DENSITY_POWER           0.2
#define MAX_N_HIGH              10
#define TREE_POWER              2.0
#define ERF_TABLE_LIMIT         6.0             /* For x>6 erf(x)-1<double precision machine epsilon, so no need to store the values for larger x. */
#define ERF_TABLE_SIZE          6145
#define BIN_WIDTH               (ERF_TABLE_LIMIT/(ERF_TABLE_SIZE-1.))
#define IBIN_WIDTH              (1./BIN_WIDTH)
#define N_VEL_SEG_PER_HALF      1
#define NUM_VEL_COEFFS          (1+2*N_VEL_SEG_PER_HALF) /* This is the number of velocity samples per edge (not including the grid vertices at each end of the edge). Currently this is elsewhere hard-wired at 3, the macro just being used in the file I/O modules. Note that we want an odd number of velocity samples per edge if we want to have the ability to do 2nd-order interpolation of velocity within Delaunay tetrahedra. */

/* Bit locations for the grid data-stage mask, that records the information which is present in the grid struct: */
#define DS_bit_x             0	/* id, x, sink */
#define DS_bit_neighbours    1	/* neigh, dir, ds, numNeigh */
#define DS_bit_velocity      2	/* vel */
#define DS_bit_density       3	/* dens */
#define DS_bit_abundance     4	/* abun, nmol */
#define DS_bit_turb_doppler  5	/* dopb */
#define DS_bit_temperatures  6	/* t */
#define DS_bit_magfield      7	/* B */
#define DS_bit_ACOEFF        8	/* a0, a1, a2, a3, a4 */
#define DS_bit_populations   9	/* mol */

#define DS_mask_x             1<<DS_bit_x
#define DS_mask_neighbours   (1<<DS_bit_neighbours   | DS_mask_x)
#define DS_mask_velocity     (1<<DS_bit_velocity     | DS_mask_x)
#define DS_mask_density      (1<<DS_bit_density      | DS_mask_x)
#define DS_mask_abundance    (1<<DS_bit_abundance    | DS_mask_x)
#define DS_mask_turb_doppler (1<<DS_bit_turb_doppler | DS_mask_x)
#define DS_mask_temperatures (1<<DS_bit_temperatures | DS_mask_x)
#define DS_mask_magfield     (1<<DS_bit_magfield     | DS_mask_x)
#define DS_mask_ACOEFF       (1<<DS_bit_ACOEFF       | DS_mask_neighbours | DS_mask_velocity)

#define DS_mask_1            DS_mask_x
#define DS_mask_2            DS_mask_neighbours
#define DS_mask_3            (DS_mask_2|DS_mask_density|DS_mask_abundance|DS_mask_turb_doppler|DS_mask_temperatures|DS_mask_ACOEFF)
#define DS_mask_populations  (1<<DS_bit_populations | DS_mask_3)
#define DS_mask_4            DS_mask_populations
#define DS_mask_all          (DS_mask_populations | DS_mask_magfield)
#define DS_mask_all_but_mag  DS_mask_all & ~(1<<DS_bit_magfield)

#include "collparts.h"
#include "inpars.h"

typedef struct {
  double radius,minScale,tcmb,*nMolWeights,*dustWeights;
  double radiusSqu,minScaleSqu,taylorCutoff,gridDensGlobalMax;
  double (*gridDensMaxLoc)[DIM],*gridDensMaxValues,*collPartMolWeights;
  int sinkPoints,pIntensity,blend,*collPartIds,traceRayAlgorithm,samplingAlgorithm;
  int ncell,nImages,nSpecies,numDensities,doPregrid,numGridDensMaxima;
  int sampling,lte_only,init_lte,antialias,polarization,nThreads,numDims;
  int nLineImages,nContImages;
  int dataFlags,nSolveIters;
  char *outputfile,*binoutputfile,*gridfile,*pregrid,*restart,*dust;
  char *gridInFile,**gridOutFiles;
  char **girdatfile,**moldatfile,**collPartNames;
  _Bool writeGridAtStage[NUM_GRID_STAGES],resetRNG,doInterpolateVels,useAbun;
} configInfo;

struct cpData {
  double *down,*temp;
  int collPartId,ntemp,ntrans,*lcl,*lcu,densityIndex;
  char *name;
};

/* Molecular data: shared attributes */
typedef struct {
  int nlev,nline,npart;
  int *lal,*lau;
  double *aeinst,*freq,*beinstu,*beinstl,*eterm,*gstat,*gir;
  double *cmb,amass;
  struct cpData *part;
  char molName[80];
} molData;

/* Data concerning a single grid vertex which is passed from photon() to stateq(). This data needs to be thread-safe. */
typedef struct {
  double *jbar,*phot,*vfac,*vfac_loc;
} gridPointData;

/* Point coordinate */
typedef struct {
  double x[DIM];
  double xn[DIM];
} point;

struct rates {
  int t_binlow;
  double interp_coeff;
};

struct continuumLine{
  double dust, knu;
};

struct populations {
  double *pops,*specNumDens;
  double dopb,binv,nmol,abun;
  struct rates *partner;
  struct continuumLine *cont;
};

/* Grid properties */
struct grid {
  int id;
  double x[DIM], vel[DIM], B[3]; /* B field only makes physical sense in 3 dimensions. */
  double *v1,*v2,*v3;
  int numNeigh;
  point *dir;
  struct grid **neigh;
  double *w;
  int sink;
  int nphot;
  int conv;
  double *dens,t[2],dopb_turb;
  double *ds;
  struct populations *mol;
  struct continuumLine cont;
};

typedef struct{
  double x[DIM], xCmpntRay, B[3];
  struct populations *mol;
  struct continuumLine cont;
} gridInterp;

typedef struct {
  double *intense;
  double *tau;
  double stokes[3];
  int numRays;
} spec;

/* Image information */
typedef struct {
  int doline;
  int nchan,trans,molI;
  spec *pixel;
  double velres;
  double imgres;
  int pxls;
  char *units;
  int *imgunits;
  int numunits;
  double freq,bandwidth;
  char *filename;
  double source_vel;
  double theta,phi,incl,posang,azimuth;
  double distance;
  double rotMat[3][3];
  _Bool doInterpolateVels;
} imageInfo;

typedef struct {
  double x,y, *intensity, *tau;
  unsigned int ppi;
} rayData;

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

/* NOTE that it is assumed that vertx[i] is opposite the face that abuts with neigh[i] for all i.
*/ 
struct cell {
  struct grid *vertx[DIM+1];
  struct cell *neigh[DIM+1]; /* ==NULL flags an external face. */
  unsigned long id;
  double centre[DIM];
};

/* Some global variables */
extern int silent;

/* User-specifiable functions */
void density(double,double,double,double *);
void temperature(double,double,double,double *);
void abundance(double,double,double,double *);
void molNumDensity(double,double,double,double *);
void doppler(double,double,double, double *);
void velocity(double,double,double,double *);
void magfield(double,double,double,double *);
void gasIIdust(double,double,double,double *);
double gridDensity(configInfo*, double*);

/* More functions */
void	run(inputPars, image*, const int);

_Bool	allBitsSet(const int flags, const int mask);
_Bool	anyBitSet(const int flags, const int mask);
_Bool	bitIsSet(const int flags, const int bitI);
_Bool	onlyBitsSet(const int flags, const int mask);

void	binpopsout(configInfo*, struct grid*, molData*);
void	calcFastExpRange(const int, const int, int*, int*, int*);
void	calcGridCollRates(configInfo*, molData*, struct grid*);
void	calcGridContDustOpacity(configInfo*, const double, double*, double*, const int, struct grid*);
void	calcGridLinesDustOpacity(configInfo*, molData*, double*, double*, const int, struct grid*);
void	calcGridMolDensities(configInfo*, struct grid**);
void	calcGridMolDoppler(configInfo*, molData*, struct grid*);
void	calcGridMolSpecNumDens(configInfo*, molData*, struct grid*);
void	calcInterpCoeffs(configInfo*, struct grid*);
void	calcInterpCoeffs_lin(configInfo*, struct grid*);
void	calcSourceFn(double, const configInfo*, double*, double*);
void	calcTableEntries(const int, const int);
void	checkFirstLineMolDat(FILE *fp, char *moldatfile);
void	checkGridDensities(configInfo*, struct grid*);
void	checkUserDensWeights(configInfo*);
void    copyInparStr(const char*, char**);
void	delaunay(const int, struct grid*, const unsigned long, const _Bool, const _Bool, struct cell**, unsigned long*);
void	distCalc(configInfo*, struct grid*);
double	dotProduct3D(const double*, const double*);
int	factorial(const int);
double	FastExp(const float);
void    fillErfTable();
void	fit_d1fi(double, double, double*);
void	fit_fi(double, double, double*);
void	fit_rr(double, double, double*);
void	freeConfigInfo(configInfo par);
void	freeGrid(const unsigned int, const unsigned short, struct grid*);
void	freeGridPointData(const int, gridPointData*);
void	freeImgInfo(const int, imageInfo*);
void	freeInputPars(inputPars par);
void	freeMolData(const int, molData*);
void	freeMolsWithBlends(struct molWithBlends*, const int);
void	freePopulation(const unsigned short, struct populations*);
void	freeSomeGridFields(const unsigned int, const unsigned short, struct grid*);
double  gaussline(const double, const double);
void	getArea(configInfo*, struct grid*, const gsl_rng*);
void	getclosest(double, double, double, long*, long*, double*, double*, double*);
double	geterf(const double, const double);
void	getjbar(int, molData*, struct grid*, const int, configInfo*, struct blendInfo, int, gridPointData*, double*);
void	getMass(configInfo*, struct grid*, const gsl_rng*);
void	getVelocities(configInfo *, struct grid *);
void	getVelocities_pregrid(configInfo *, struct grid *);
void	gridPopsInit(configInfo*, molData*, struct grid*);
void	input(inputPars*, image*);
float	invSqrt(float);
void	levelPops(molData*, configInfo*, struct grid*, int*, double*, double*, const int);
void	lineBlend(molData*, configInfo*, struct blendInfo*);
void	LTE(configInfo*, struct grid*, molData*);
void	lteOnePoint(molData*, const int, const double, double*);
void	mallocAndSetDefaultGrid(struct grid**, const size_t, const size_t);
void	molInit(configInfo*, molData*);
void	openSocket(char*);
void	parseInput(inputPars, image*, const int, configInfo*, imageInfo**, molData**);
void	photon(int, struct grid*, molData*, const gsl_rng*, configInfo*, const int, struct blendInfo, gridPointData*, double*);
double	planckfunc(const double, const double);
int	pointEvaluation(configInfo*, const double, double*);
void	popsin(configInfo*, struct grid**, molData**, int*);
void	popsout(configInfo*, struct grid*, molData*);
void	predefinedGrid(configInfo*, struct grid*);
void	processFitsError(int);
double	ratranInput(char*, char*, double, double, double);
void	raytrace(int, configInfo*, struct grid*, molData*, imageInfo*, double*, double*, const int);
void	readDustFile(char*, double**, double**, int*);
void	readOrBuildGrid(configInfo*, struct grid**);
void	readUserInput(inputPars*, imageInfo**, int*, int*);
unsigned long reorderGrid(const unsigned long, struct grid*);
void	report(int, configInfo*, struct grid*);
void	reportInfsAtOrigin(const int, const double*, const char*);
void	setUpConfig(configInfo*, imageInfo**, molData**);
void	setUpDensityAux(configInfo*, int*, const int);
void	smooth(configInfo*, struct grid*);
void    sourceFunc_line(const molData*, const double, const struct populations*, const int, double*, double*);
void    sourceFunc_cont(const struct continuumLine, double*, double*);
void	sourceFunc_pol(double*, const struct continuumLine, double (*rotMat)[3], double*, double*);
void	stateq(int, struct grid*, molData*, const int, configInfo*, struct blendInfo, int, gridPointData*, double*, _Bool*);
void	statistics(int, molData*, struct grid*, int*, double*, double*, int*);
void	stokesangles(double*, double (*rotMat)[3], double*);
double	taylor(const int, const float);
void	writeFitsAllUnits(const int, configInfo*, imageInfo*);
void	writeGridIfRequired(configInfo*, struct grid*, molData*);
void	write_VTK_unstructured_Points(configInfo*, struct grid*);


/* Curses functions */

void	bail_out(char*);
void	casaStyleProgressBar(const int, int);
void	collpartmesg(char*, int);
void	collpartmesg2(char*);
void	collpartmesg3(int, int);
void	goodnight(int, char*);
void	greetings(void);
void	greetings_parallel(int);
void	printDone(int);
void	printMessage(char *);
void	progressbar(double, int);
void	progressbar2(configInfo*, int, int, double, double, double);
void	quotemass(double);
void	screenInfo(void);
void	warning(char*);

#ifdef FASTEXP
extern double EXP_TABLE_2D[128][10];
extern double EXP_TABLE_3D[256][2][10];
#else
extern double EXP_TABLE_2D[1][1]; /* nominal definitions so the fastexp.c module will compile. */
extern double EXP_TABLE_3D[1][1][1];
#endif

extern double ERF_TABLE[ERF_TABLE_SIZE];
extern double oneOver_i[FAST_EXP_MAX_TAYLOR+1];

#endif /* LIME_H */

