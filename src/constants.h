/*
 *  constants.h
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#ifndef CONSTANTS_H
#define CONSTANTS_H

/* Physical constants */
/* - NIST values as of 23 Sept 2015: */
#define AMU             1.66053904e-27       /* atomic mass unit             [kg]		*/
#define CLIGHT          2.99792458e8         /* speed of light in vacuum     [m/s]		*/
#define HPLANCK         6.626070040e-34      /* Planck constant              [J s]		*/
#define KBOLTZ          1.38064852e-23       /* Boltzmann constant           [J/K]		*/
#define YJULIAN         (365.25*24.0*3600.0) /* Length of the Julian year    [s]		*/
#define STEFANB         5.670367e-8          /* Stefan-Boltzmann constant    [W/m^2/K^4]	*/

/* From IAU 2009: */
#define GRAV            6.67428e-11          /* gravitational constant       [m^3/kg/s^2]	*/
#define AU              1.495978707e11       /* astronomical unit            [m]		*/

#define LOCAL_CMB_TEMP  2.72548              /* local mean CMB temp. from Fixsen (2009) [K]	*/

/* Derived: */
#define PC              3.08567758e16        /* parsec (~3600*180*AU/PI)     [m]		*/
#define MSUN            1.9891e30            /* Solar mass                   [kg]		*/
#define RSUN            6.957e8              /* Solar radius                 [m]		*/

/* CGS: */
#define AMU_cgs         (AMU*1000.0)         /* atomic mass unit             [g]		*/
#define CLIGHT_cgs      (CLIGHT*100.0)       /* speed of light in vacuum     [cm/s]		*/
#define HPLANCK_cgs     (HPLANCK*1.0e7)      /* Planck constant              [erg s]		*/
#define KBOLTZ_cgs      (KBOLTZ*1.0e7)       /* Boltzmann constant           [erg/K]		*/
#define STEFANB_cgs     (STEFANB*1000.0)     /* Stefan-Boltzmann constant    [erg/cm^2/K^4/s]	*/
#define GRAV_cgs        (GRAV*1000.0)        /* gravitational constant       [cm^3/g/s^2]	*/
#define AU_cgs          (AU*100.0)           /* astronomical unit            [cm]		*/
#define MSUN_cgs        (MSUN*1000.0)        /* Solar mass                   [g]		*/
#define RSUN_cgs        (RSUN*100.0)         /* Solar radius                 [cm]		*/

#define SQRT_PI                 (sqrt(M_PI))           /* sqrt(pi)	*/
#define ARCSEC_TO_RAD           (M_PI/180.0/3600.0)
#define EPS                     1.0e-30                /* general use small number */

#define TYPICAL_ISM_DENS        1000.0
#define STR_LEN_0               80
#define STR_LEN_1               127
#define STR_LEN_2               255

#ifndef TRUE
#define TRUE                   (_Bool)1
#endif

#ifndef FALSE
#define FALSE                  (_Bool)0
#endif

extern int silent;

#endif /* CONSTANTS_H */

