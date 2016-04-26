/*
 *  model.c
 *  LIME, The versatile 3D line modeling tool 
 *
 *  Created by Christian Brinch on 11/05/07.
 *  Copyright 2006-2008, Christian Brinch, 
 *  <cbrinch@astro.uni-bonn.de>
 *      Argelander-Institut fur Astronomie, 
 *      Universitat Bonn. 
 *      All rights reserved.
 *
 */

#include "lime.h"
#include <math.h>

/*****************************************************************/
void
input(inputPars *par, image *img){
  par->radius         = 16100*AU;
  par->minScale       = 34.6*AU;
  par->tcmb           = 2.725;
  par->pIntensity     = 30000;
  par->sinkPoints     = 6000;
  par->dust           = "jena_thin_e6.tab";
  par->moldatfile[0]  = "co.dat";
  par->outputfile     = "hh46_CO_cavitywall.pop";
  par->sampling       = 23;
  par->species        = "CO";
  par->lte_only       = 0;

// Image #0    
  img[0].nchan      = 61;               // Number of channels
  img[0].velres     = 400.;             // Channel resolution in m/s
  img[0].trans      = 0;                // zero-indexed J quantum number
  img[0].pxls       = 211;              // Pixels per dimension
  img[0].imgres     = 0.1;              // Resolution in arc seconds
  img[0].theta      = 0.5*(53./90.)*PI; // 0: face-on, pi/2: edge-on
  img[0].phi        = 0.0;              // azimuthal direction
  img[0].distance   = 450*PC;           // source distance in pc
  img[0].source_vel = 0;                // source velocity in m/s
  img[0].unit       = 0;                // 0:Kelvin 1:Jy/px 2:SI 3:Lsun/pixel
  img[0].filename   = "hh46_CO_cavitywall_001.fits"; // Output filename

// Image #1    
  img[1].nchan      = 61;               // Number of channels
  img[1].velres     = 400.;             // Channel resolution in m/s
  img[1].trans      = 1;                // zero-indexed J quantum number
  img[1].pxls       = 573;              // Pixels per dimension
  img[1].imgres     = 0.1;              // Resolution in arc seconds
  img[1].theta      = 0.5*(53./90.)*PI; // 0: face-on, pi/2: edge-on
  img[1].phi        = 0.0;              // azimuthal direction
  img[1].distance   = 450*PC;           // source distance in pc
  img[1].source_vel = 0;                // source velocity in m/s
  img[1].unit       = 0;                // 0:Kelvin 1:Jy/px 2:SI 3:Lsun/pixel
  img[1].filename   = "hh46_CO_cavitywall_002.fits"; // Output filename

// Image #2    
  img[2].nchan      = 61;               // Number of channels
  img[2].velres     = 400.;             // Channel resolution in m/s
  img[2].trans      = 2;                // zero-indexed J quantum number
  img[2].pxls       = 381;              // Pixels per dimension
  img[2].imgres     = 0.1;              // Resolution in arc seconds
  img[2].theta      = 0.5*(53./90.)*PI; // 0: face-on, pi/2: edge-on
  img[2].phi        = 0.0;              // azimuthal direction
  img[2].distance   = 450*PC;           // source distance in pc
  img[2].source_vel = 0;                // source velocity in m/s
  img[2].unit       = 0;                // 0:Kelvin 1:Jy/px 2:SI 3:Lsun/pixel
  img[2].filename   = "hh46_CO_cavitywall_003.fits"; // Output filename

// Image #3    
  img[3].nchan      = 61;               // Number of channels
  img[3].velres     = 400.;             // Channel resolution in m/s
  img[3].trans      = 3;                // zero-indexed J quantum number
  img[3].pxls       = 287;              // Pixels per dimension
  img[3].imgres     = 0.1;              // Resolution in arc seconds
  img[3].theta      = 0.5*(53./90.)*PI; // 0: face-on, pi/2: edge-on
  img[3].phi        = 0.0;              // azimuthal direction
  img[3].distance   = 450*PC;           // source distance in pc
  img[3].source_vel = 0;                // source velocity in m/s
  img[3].unit       = 0;                // 0:Kelvin 1:Jy/px 2:SI 3:Lsun/pixel
  img[3].filename   = "hh46_CO_cavitywall_004.fits"; // Output filename

// Image #4    
  img[4].nchan      = 61;               // Number of channels
  img[4].velres     = 400.;             // Channel resolution in m/s
  img[4].trans      = 4;                // zero-indexed J quantum number
  img[4].pxls       = 229;               // Pixels per dimension
  img[4].imgres     = 0.1;              // Resolution in arc seconds
  img[4].theta      = 0.5*(53./90.)*PI; // 0: face-on, pi/2: edge-on
  img[4].phi        = 0.0;              // azimuthal direction
  img[4].distance   = 450*PC;           // source distance in pc
  img[4].source_vel = 0;                // source velocity in m/s
  img[4].unit       = 0;                // 0:Kelvin 1:Jy/px 2:SI 3:Lsun/pixel
  img[4].filename   = "hh46_CO_cavitywall_005.fits"; // Output filename

// Image #5    
  img[5].nchan      = 61;               // Number of channels
  img[5].velres     = 400.;             // Channel resolution in m/s
  img[5].trans      = 5;                // zero-indexed J quantum number
  img[5].pxls       = 191;              // Pixels per dimension
  img[5].imgres     = 0.1;              // Resolution in arc seconds
  img[5].theta      = 0.5*(53./90.)*PI; // 0: face-on, pi/2: edge-on
  img[5].phi        = 0.0;              // azimuthal direction
  img[5].distance   = 450*PC;           // source distance in pc
  img[5].source_vel = 0;                // source velocity in m/s
  img[5].unit       = 0;                // 0:Kelvin 1:Jy/px 2:SI 3:Lsun/pixel
  img[5].filename   = "hh46_CO_cavitywall_006.fits"; // Output filename

// Image #6    
  img[6].nchan      = 61;               // Number of channels
  img[6].velres     = 400.;             // Channel resolution in m/s
  img[6].trans      = 6;                // zero-indexed J quantum number
  img[6].pxls       = 163;              // Pixels per dimension
  img[6].imgres     = 0.1;              // Resolution in arc seconds
  img[6].theta      = 0.5*(53./90.)*PI; // 0: face-on, pi/2: edge-on
  img[6].phi        = 0.0;              // azimuthal direction
  img[6].distance   = 450*PC;           // source distance in pc
  img[6].source_vel = 0;                // source velocity in m/s
  img[6].unit       = 0;                // 0:Kelvin 1:Jy/px 2:SI 3:Lsun/pixel
  img[6].filename   = "hh46_CO_cavitywall_007.fits"; // Output filename

// Image #7    
  img[7].nchan      = 61;               // Number of channels
  img[7].velres     = 400.;             // Channel resolution in m/s
  img[7].trans      = 7;                // zero-indexed J quantum number
  img[7].pxls       = 231;              // Pixels per dimension
  img[7].imgres     = 0.1;              // Resolution in arc seconds
  img[7].theta      = 0.5*(53./90.)*PI; // 0: face-on, pi/2: edge-on
  img[7].phi        = 0.0;              // azimuthal direction
  img[7].distance   = 450*PC;           // source distance in pc
  img[7].source_vel = 0;                // source velocity in m/s
  img[7].unit       = 0;                // 0:Kelvin 1:Jy/px 2:SI 3:Lsun/pixel
  img[7].filename   = "hh46_CO_cavitywall_008.fits"; // Output filename

// Image #8    
  img[8].nchan      = 61;               // Number of channels
  img[8].velres     = 400.;             // Channel resolution in m/s
  img[8].trans      = 8;                // zero-indexed J quantum number
  img[8].pxls       = 281;              // Pixels per dimension
  img[8].imgres     = 0.1;              // Resolution in arc seconds
  img[8].theta      = 0.5*(53./90.)*PI; // 0: face-on, pi/2: edge-on
  img[8].phi        = 0.0;              // azimuthal direction
  img[8].distance   = 450*PC;           // source distance in pc
  img[8].source_vel = 0;                // source velocity in m/s
  img[8].unit       = 0;                // 0:Kelvin 1:Jy/px 2:SI 3:Lsun/pixel
  img[8].filename   = "hh46_CO_cavitywall_009.fits"; // Output filename

// Image #9    
  img[9].nchan      = 61;               // Number of channels
  img[9].velres     = 400.;             // Channel resolution in m/s
  img[9].trans      = 9;                // zero-indexed J quantum number
  img[9].pxls       = 323;              // Pixels per dimension
  img[9].imgres     = 0.1;              // Resolution in arc seconds
  img[9].theta      = 0.5*(53./90.)*PI; // 0: face-on, pi/2: edge-on
  img[9].phi        = 0.0;              // azimuthal direction
  img[9].distance   = 450*PC;           // source distance in pc
  img[9].source_vel = 0;                // source velocity in m/s
  img[9].unit       = 0;                // 0:Kelvin 1:Jy/px 2:SI 3:Lsun/pixel
  img[9].filename   = "hh46_CO_cavitywall_010.fits"; // Output filename

// Image #10    
  img[10].nchan      = 61;               // Number of channels
  img[10].velres     = 400.;             // Channel resolution in m/s
  img[10].trans      = 10;               // zero-indexed J quantum number
  img[10].pxls       = 357;              // Pixels per dimension
  img[10].imgres     = 0.1;              // Resolution in arc seconds
  img[10].theta      = 0.5*(53./90.)*PI; // 0: face-on, pi/2: edge-on
  img[10].phi        = 0.0;              // azimuthal direction
  img[10].distance   = 450*PC;           // source distance in pc
  img[10].source_vel = 0;                // source velocity in m/s
  img[10].unit       = 0;                // 0:Kelvin 1:Jy/px 2:SI 3:Lsun/pixel
  img[10].filename   = "hh46_CO_cavitywall_011.fits"; // Output filename

// Image #11    
  img[11].nchan      = 61;               // Number of channels
  img[11].velres     = 400.;             // Channel resolution in m/s
  img[11].trans      = 11;               // zero-indexed J quantum number
  img[11].pxls       = 327;              // Pixels per dimension
  img[11].imgres     = 0.1;              // Resolution in arc seconds
  img[11].theta      = 0.5*(53./90.)*PI; // 0: face-on, pi/2: edge-on
  img[11].phi        = 0.0;              // azimuthal direction
  img[11].distance   = 450*PC;           // source distance in pc
  img[11].source_vel = 0;                // source velocity in m/s
  img[11].unit       = 0;                // 0:Kelvin 1:Jy/px 2:SI 3:Lsun/pixel
  img[11].filename   = "hh46_CO_cavitywall_012.fits"; // Output filename

// Image #12    
  img[12].nchan      = 61;               // Number of channels
  img[12].velres     = 400.;             // Channel resolution in m/s
  img[12].trans      = 12;               // zero-indexed J quantum number
  img[12].pxls       = 303;              // Pixels per dimension
  img[12].imgres     = 0.1;              // Resolution in arc seconds
  img[12].theta      = 0.5*(53./90.)*PI; // 0: face-on, pi/2: edge-on
  img[12].phi        = 0.0;              // azimuthal direction
  img[12].distance   = 450*PC;           // source distance in pc
  img[12].source_vel = 0;                // source velocity in m/s
  img[12].unit       = 0;                // 0:Kelvin 1:Jy/px 2:SI 3:Lsun/pixel
  img[12].filename   = "hh46_CO_cavitywall_013.fits"; // Output filename

// Image #13    
  img[13].nchan      = 61;               // Number of channels
  img[13].velres     = 400.;             // Channel resolution in m/s
  img[13].trans      = 13;               // zero-indexed J quantum number
  img[13].pxls       = 301;              // Pixels per dimension
  img[13].imgres     = 0.1;              // Resolution in arc seconds
  img[13].theta      = 0.5*(53./90.)*PI; // 0: face-on, pi/2: edge-on
  img[13].phi        = 0.0;              // azimuthal direction
  img[13].distance   = 450*PC;           // source distance in pc
  img[13].source_vel = 0;                // source velocity in m/s
  img[13].unit       = 0;                // 0:Kelvin 1:Jy/px 2:SI 3:Lsun/pixel
  img[13].filename   = "hh46_CO_cavitywall_014.fits"; // Output filename

// Image #14    
  img[14].nchan      = 61;               // Number of channels
  img[14].velres     = 400.;             // Channel resolution in m/s
  img[14].trans      = 14;               // zero-indexed J quantum number
  img[14].pxls       = 301;              // Pixels per dimension
  img[14].imgres     = 0.1;              // Resolution in arc seconds
  img[14].theta      = 0.5*(53./90.)*PI; // 0: face-on, pi/2: edge-on
  img[14].phi        = 0.0;              // azimuthal direction
  img[14].distance   = 450*PC;           // source distance in pc
  img[14].source_vel = 0;                // source velocity in m/s
  img[14].unit       = 0;                // 0:Kelvin 1:Jy/px 2:SI 3:Lsun/pixel
  img[14].filename   = "hh46_CO_cavitywall_015.fits"; // Output filename

// Image #15    
  img[15].nchan      = 61;               // Number of channels
  img[15].velres     = 400.;             // Channel resolution in m/s
  img[15].trans      = 15;               // zero-indexed J quantum number
  img[15].pxls       = 301;              // Pixels per dimension
  img[15].imgres     = 0.1;              // Resolution in arc seconds
  img[15].theta      = 0.5*(53./90.)*PI; // 0: face-on, pi/2: edge-on
  img[15].phi        = 0.0;              // azimuthal direction
  img[15].distance   = 450*PC;           // source distance in pc
  img[15].source_vel = 0;                // source velocity in m/s
  img[15].unit       = 0;                // 0:Kelvin 1:Jy/px 2:SI 3:Lsun/pixel
  img[15].filename   = "hh46_CO_cavitywall_016.fits"; // Output filename

// Image #16    
  img[16].nchan      = 61;               // Number of channels
  img[16].velres     = 400.;             // Channel resolution in m/s
  img[16].trans      = 16;               // zero-indexed J quantum number
  img[16].pxls       = 301;              // Pixels per dimension
  img[16].imgres     = 0.1;              // Resolution in arc seconds
  img[16].theta      = 0.5*(53./90.)*PI; // 0: face-on, pi/2: edge-on
  img[16].phi        = 0.0;              // azimuthal direction
  img[16].distance   = 450*PC;           // source distance in pc
  img[16].source_vel = 0;                // source velocity in m/s
  img[16].unit       = 0;                // 0:Kelvin 1:Jy/px 2:SI 3:Lsun/pixel
  img[16].filename   = "hh46_CO_cavitywall_017.fits"; // Output filename

// Image #17    
  img[17].nchan      = 61;               // Number of channels
  img[17].velres     = 400.;             // Channel resolution in m/s
  img[17].trans      = 17;               // zero-indexed J quantum number
  img[17].pxls       = 301;              // Pixels per dimension
  img[17].imgres     = 0.1;              // Resolution in arc seconds
  img[17].theta      = 0.5*(53./90.)*PI; // 0: face-on, pi/2: edge-on
  img[17].phi        = 0.0;              // azimuthal direction
  img[17].distance   = 450*PC;           // source distance in pc
  img[17].source_vel = 0;                // source velocity in m/s
  img[17].unit       = 0;                // 0:Kelvin 1:Jy/px 2:SI 3:Lsun/pixel
  img[17].filename   = "hh46_CO_cavitywall_018.fits"; // Output filename

// Image #18    
  img[18].nchan      = 61;               // Number of channels
  img[18].velres     = 400.;             // Channel resolution in m/s
  img[18].trans      = 18;               // zero-indexed J quantum number
  img[18].pxls       = 301;              // Pixels per dimension
  img[18].imgres     = 0.1;              // Resolution in arc seconds
  img[18].theta      = 0.5*(53./90.)*PI; // 0: face-on, pi/2: edge-on
  img[18].phi        = 0.0;              // azimuthal direction
  img[18].distance   = 450*PC;           // source distance in pc
  img[18].source_vel = 0;                // source velocity in m/s
  img[18].unit       = 0;                // 0:Kelvin 1:Jy/px 2:SI 3:Lsun/pixel
  img[18].filename   = "hh46_CO_cavitywall_019.fits"; // Output filename

// Image #19    
  img[19].nchan      = 61;               // Number of channels
  img[19].velres     = 400.;             // Channel resolution in m/s
  img[19].trans      = 19;               // zero-indexed J quantum number
  img[19].pxls       = 301;              // Pixels per dimension
  img[19].imgres     = 0.1;              // Resolution in arc seconds
  img[19].theta      = 0.5*(53./90.)*PI; // 0: face-on, pi/2: edge-on
  img[19].phi        = 0.0;              // azimuthal direction
  img[19].distance   = 450*PC;           // source distance in pc
  img[19].source_vel = 0;                // source velocity in m/s
  img[19].unit       = 0;                // 0:Kelvin 1:Jy/px 2:SI 3:Lsun/pixel
  img[19].filename   = "hh46_CO_cavitywall_020.fits"; // Output filename

// Image #20    
  img[20].nchan      = 61;               // Number of channels
  img[20].velres     = 400.;             // Channel resolution in m/s
  img[20].trans      = 20;               // zero-indexed J quantum number
  img[20].pxls       = 301;              // Pixels per dimension
  img[20].imgres     = 0.1;              // Resolution in arc seconds
  img[20].theta      = 0.5*(53./90.)*PI; // 0: face-on, pi/2: edge-on
  img[20].phi        = 0.0;              // azimuthal direction
  img[20].distance   = 450*PC;           // source distance in pc
  img[20].source_vel = 0;                // source velocity in m/s
  img[20].unit       = 0;                // 0:Kelvin 1:Jy/px 2:SI 3:Lsun/pixel
  img[20].filename   = "hh46_CO_cavitywall_021.fits"; // Output filename

// Image #21    
  img[21].nchan      = 61;               // Number of channels
  img[21].velres     = 400.;             // Channel resolution in m/s
  img[21].trans      = 21;               // zero-indexed J quantum number
  img[21].pxls       = 301;              // Pixels per dimension
  img[21].imgres     = 0.1;              // Resolution in arc seconds
  img[21].theta      = 0.5*(53./90.)*PI; // 0: face-on, pi/2: edge-on
  img[21].phi        = 0.0;              // azimuthal direction
  img[21].distance   = 450*PC;           // source distance in pc
  img[21].source_vel = 0;                // source velocity in m/s
  img[21].unit       = 0;                // 0:Kelvin 1:Jy/px 2:SI 3:Lsun/pixel
  img[21].filename   = "hh46_CO_cavitywall_022.fits"; // Output filename

// Image #22    
  img[22].nchan      = 61;               // Number of channels
  img[22].velres     = 400.;             // Channel resolution in m/s
  img[22].trans      = 22;               // zero-indexed J quantum number
  img[22].pxls       = 301;              // Pixels per dimension
  img[22].imgres     = 0.1;              // Resolution in arc seconds
  img[22].theta      = 0.5*(53./90.)*PI; // 0: face-on, pi/2: edge-on
  img[22].phi        = 0.0;              // azimuthal direction
  img[22].distance   = 450*PC;           // source distance in pc
  img[22].source_vel = 0;                // source velocity in m/s
  img[22].unit       = 0;                // 0:Kelvin 1:Jy/px 2:SI 3:Lsun/pixel
  img[22].filename   = "hh46_CO_cavitywall_023.fits"; // Output filename

// Image #23    
  img[23].nchan      = 61;               // Number of channels
  img[23].velres     = 400.;             // Channel resolution in m/s
  img[23].trans      = 23;               // zero-indexed J quantum number
  img[23].pxls       = 301;              // Pixels per dimension
  img[23].imgres     = 0.1;              // Resolution in arc seconds
  img[23].theta      = 0.5*(53./90.)*PI; // 0: face-on, pi/2: edge-on
  img[23].phi        = 0.0;              // azimuthal direction
  img[23].distance   = 450*PC;           // source distance in pc
  img[23].source_vel = 0;                // source velocity in m/s
  img[23].unit       = 0;                // 0:Kelvin 1:Jy/px 2:SI 3:Lsun/pixel
  img[23].filename   = "hh46_CO_cavitywall_024.fits"; // Output filename

// Image #24    
  img[24].nchan      = 61;               // Number of channels
  img[24].velres     = 400.;             // Channel resolution in m/s
  img[24].trans      = 24;               // zero-indexed J quantum number
  img[24].pxls       = 301;              // Pixels per dimension
  img[24].imgres     = 0.1;              // Resolution in arc seconds
  img[24].theta      = 0.5*(53./90.)*PI; // 0: face-on, pi/2: edge-on
  img[24].phi        = 0.0;              // azimuthal direction
  img[24].distance   = 450*PC;           // source distance in pc
  img[24].source_vel = 0;                // source velocity in m/s
  img[24].unit       = 0;                // 0:Kelvin 1:Jy/px 2:SI 3:Lsun/pixel
  img[24].filename   = "hh46_CO_cavitywall_025.fits"; // Output filename

// Image #25    
  img[25].nchan      = 61;               // Number of channels
  img[25].velres     = 400.;             // Channel resolution in m/s
  img[25].trans      = 25;               // zero-indexed J quantum number
  img[25].pxls       = 301;              // Pixels per dimension
  img[25].imgres     = 0.1;              // Resolution in arc seconds
  img[25].theta      = 0.5*(53./90.)*PI; // 0: face-on, pi/2: edge-on
  img[25].phi        = 0.0;              // azimuthal direction
  img[25].distance   = 450*PC;           // source distance in pc
  img[25].source_vel = 0;                // source velocity in m/s
  img[25].unit       = 0;                // 0:Kelvin 1:Jy/px 2:SI 3:Lsun/pixel
  img[25].filename   = "hh46_CO_cavitywall_026.fits"; // Output filename

// Image #26    
  img[26].nchan      = 61;               // Number of channels
  img[26].velres     = 400.;             // Channel resolution in m/s
  img[26].trans      = 26;               // zero-indexed J quantum number
  img[26].pxls       = 301;              // Pixels per dimension
  img[26].imgres     = 0.1;              // Resolution in arc seconds
  img[26].theta      = 0.5*(53./90.)*PI; // 0: face-on, pi/2: edge-on
  img[26].phi        = 0.0;              // azimuthal direction
  img[26].distance   = 450*PC;           // source distance in pc
  img[26].source_vel = 0;                // source velocity in m/s
  img[26].unit       = 0;                // 0:Kelvin 1:Jy/px 2:SI 3:Lsun/pixel
  img[26].filename   = "hh46_CO_cavitywall_027.fits"; // Output filename

// Image #27    
  img[27].nchan      = 61;               // Number of channels
  img[27].velres     = 400.;             // Channel resolution in m/s
  img[27].trans      = 27;               // zero-indexed J quantum number
  img[27].pxls       = 301;              // Pixels per dimension
  img[27].imgres     = 0.1;              // Resolution in arc seconds
  img[27].theta      = 0.5*(53./90.)*PI; // 0: face-on, pi/2: edge-on
  img[27].phi        = 0.0;              // azimuthal direction
  img[27].distance   = 450*PC;           // source distance in pc
  img[27].source_vel = 0;                // source velocity in m/s
  img[27].unit       = 0;                // 0:Kelvin 1:Jy/px 2:SI 3:Lsun/pixel
  img[27].filename   = "hh46_CO_cavitywall_028.fits"; // Output filename

// Image #28    
  img[28].nchan      = 61;               // Number of channels
  img[28].velres     = 400.;             // Channel resolution in m/s
  img[28].trans      = 28;               // zero-indexed J quantum number
  img[28].pxls       = 301;              // Pixels per dimension
  img[28].imgres     = 0.1;              // Resolution in arc seconds
  img[28].theta      = 0.5*(53./90.)*PI; // 0: face-on, pi/2: edge-on
  img[28].phi        = 0.0;              // azimuthal direction
  img[28].distance   = 450*PC;           // source distance in pc
  img[28].source_vel = 0;                // source velocity in m/s
  img[28].unit       = 0;                // 0:Kelvin 1:Jy/px 2:SI 3:Lsun/pixel
  img[28].filename   = "hh46_CO_cavitywall_029.fits"; // Output filename

// Image #29    
  img[29].nchan      = 61;               // Number of channels
  img[29].velres     = 400.;             // Channel resolution in m/s
  img[29].trans      = 29;               // zero-indexed J quantum number
  img[29].pxls       = 301;              // Pixels per dimension
  img[29].imgres     = 0.1;              // Resolution in arc seconds
  img[29].theta      = 0.5*(53./90.)*PI; // 0: face-on, pi/2: edge-on
  img[29].phi        = 0.0;              // azimuthal direction
  img[29].distance   = 450*PC;           // source distance in pc
  img[29].source_vel = 0;                // source velocity in m/s
  img[29].unit       = 0;                // 0:Kelvin 1:Jy/px 2:SI 3:Lsun/pixel
  img[29].filename   = "hh46_CO_cavitywall_030.fits"; // Output filename

// Image #30    
  img[30].nchan      = 61;               // Number of channels
  img[30].velres     = 400.;             // Channel resolution in m/s
  img[30].trans      = 30;               // zero-indexed J quantum number
  img[30].pxls       = 301;              // Pixels per dimension
  img[30].imgres     = 0.1;              // Resolution in arc seconds
  img[30].theta      = 0.5*(53./90.)*PI; // 0: face-on, pi/2: edge-on
  img[30].phi        = 0.0;              // azimuthal direction
  img[30].distance   = 450*PC;           // source distance in pc
  img[30].source_vel = 0;                // source velocity in m/s
  img[30].unit       = 0;                // 0:Kelvin 1:Jy/px 2:SI 3:Lsun/pixel
  img[30].filename   = "hh46_CO_cavitywall_031.fits"; // Output filename

// Image #31    
  img[31].nchan      = 61;               // Number of channels
  img[31].velres     = 400.;             // Channel resolution in m/s
  img[31].trans      = 31;               // zero-indexed J quantum number
  img[31].pxls       = 301;              // Pixels per dimension
  img[31].imgres     = 0.1;              // Resolution in arc seconds
  img[31].theta      = 0.5*(53./90.)*PI; // 0: face-on, pi/2: edge-on
  img[31].phi        = 0.0;              // azimuthal direction
  img[31].distance   = 450*PC;           // source distance in pc
  img[31].source_vel = 0;                // source velocity in m/s
  img[31].unit       = 0;                // 0:Kelvin 1:Jy/px 2:SI 3:Lsun/pixel
  img[31].filename   = "hh46_CO_cavitywall_032.fits"; // Output filename

// Image #32    
  img[32].nchan      = 61;               // Number of channels
  img[32].velres     = 400.;             // Channel resolution in m/s
  img[32].trans      = 32;               // zero-indexed J quantum number
  img[32].pxls       = 301;              // Pixels per dimension
  img[32].imgres     = 0.1;              // Resolution in arc seconds
  img[32].theta      = 0.5*(53./90.)*PI; // 0: face-on, pi/2: edge-on
  img[32].phi        = 0.0;              // azimuthal direction
  img[32].distance   = 450*PC;           // source distance in pc
  img[32].source_vel = 0;                // source velocity in m/s
  img[32].unit       = 0;                // 0:Kelvin 1:Jy/px 2:SI 3:Lsun/pixel
  img[32].filename   = "hh46_CO_cavitywall_033.fits"; // Output filename

// Image #33    
  img[33].nchan      = 61;               // Number of channels
  img[33].velres     = 400.;             // Channel resolution in m/s
  img[33].trans      = 33;               // zero-indexed J quantum number
  img[33].pxls       = 301;              // Pixels per dimension
  img[33].imgres     = 0.1;              // Resolution in arc seconds
  img[33].theta      = 0.5*(53./90.)*PI; // 0: face-on, pi/2: edge-on
  img[33].phi        = 0.0;              // azimuthal direction
  img[33].distance   = 450*PC;           // source distance in pc
  img[33].source_vel = 0;                // source velocity in m/s
  img[33].unit       = 0;                // 0:Kelvin 1:Jy/px 2:SI 3:Lsun/pixel
  img[33].filename   = "hh46_CO_cavitywall_034.fits"; // Output filename

// Image #34    
  img[34].nchan      = 61;               // Number of channels
  img[34].velres     = 400.;             // Channel resolution in m/s
  img[34].trans      = 34;               // zero-indexed J quantum number
  img[34].pxls       = 301;              // Pixels per dimension
  img[34].imgres     = 0.1;              // Resolution in arc seconds
  img[34].theta      = 0.5*(53./90.)*PI; // 0: face-on, pi/2: edge-on
  img[34].phi        = 0.0;              // azimuthal direction
  img[34].distance   = 450*PC;           // source distance in pc
  img[34].source_vel = 0;                // source velocity in m/s
  img[34].unit       = 0;                // 0:Kelvin 1:Jy/px 2:SI 3:Lsun/pixel
  img[34].filename   = "hh46_CO_cavitywall_035.fits"; // Output filename

// Image #35    
  img[35].nchan      = 61;               // Number of channels
  img[35].velres     = 400.;             // Channel resolution in m/s
  img[35].trans      = 35;               // zero-indexed J quantum number
  img[35].pxls       = 301;              // Pixels per dimension
  img[35].imgres     = 0.1;              // Resolution in arc seconds
  img[35].theta      = 0.5*(53./90.)*PI; // 0: face-on, pi/2: edge-on
  img[35].phi        = 0.0;              // azimuthal direction
  img[35].distance   = 450*PC;           // source distance in pc
  img[35].source_vel = 0;                // source velocity in m/s
  img[35].unit       = 0;                // 0:Kelvin 1:Jy/px 2:SI 3:Lsun/pixel
  img[35].filename   = "hh46_CO_cavitywall_036.fits"; // Output filename

}

/*****************************************************************/

void
density (double x, double y, double z, double gastemp, double *density){ 
  double dist,dens,cav_a,cav_b,xycav,xy,opr,dum[2];

  dist=sqrt(x*x+y*y+z*z);
  if(dist<34.6*AU) {
    //inside inner edge: set density to zero
    dens = 1e-5;
  } else {
    dens = 1e6 * 1.1e9 * pow(dist/(34.6*AU),-2.0);
  }

  /*
    Define the cavity wall.
  */
  cav_a = 34e3 * AU;
  cav_b = 8e3 * AU;
  if(fabs(z)>=2.*cav_a) {
    xycav = 0.0;
  }
  else {
    xycav = cav_b * sqrt(1.0-(fabs(z)/cav_a-1.0)*(fabs(z)/cav_a-1.0));
  }
  xy = sqrt(x*x+y*y);
  if(xy<xycav) {
    // inside the cavity
    dens = 1e-5;
  }

  /*
    Calculate the H2 ortho/para ratio
  */
  if(gastemp<0.) {
    opr = 3.0;
  } else {
    opr = 9.0*exp(-170.6/gastemp);
    if(opr>3.0) opr = 3.0;
  }

  density[0] = (1.0/(1.0+opr))*dens; //para-H2
  density[1] = (opr/(1.0+opr))*dens; //ortho-H2
  //density[0] = 0.25*dens; //para-H2
  //density[1] = 0.75*dens; //ortho-H2
}


void
temperature(double x, double y, double z, double *temperature){
#include "tgasgrid.k99.h"
#include "dusttemp.h"
#include "uvflux.h"
  int ir,ith,jr,jth,idens,i,ig0;
  double r,theta,g0,g0att,g0geom,av,dens,dens0[2];
  double epsdens,epsg0,temp1,temp2,gastemp,dusttemp,cav_a,cav_b,xycav,xy;
  
  r=sqrt(x*x+y*y+z*z);
  theta=atan2(sqrt(x*x+y*y),fabs(z));

  /*
    Find the dust temperature from the RADMC output.
  */
  if(r<=gr[0]) {
    jr=0;
  }
  else if(r>=gr[nr-1]) {
    jr=nr-1;
  }
  else {
    for(ir=0;ir<=nr-2;ir++) {
      if(r>=gr[ir] && r<gr[ir+1]) {
        jr=ir;
      }
    }
  }

  if(theta<=gtheta[0]) {
    jth=0;
  }
  else if(theta>=gtheta[nth-1]) {
    jth=nth-1;
  }
  else {
    for(ith=0;ith<=nth-2;ith++) {
      if(theta>=gtheta[ith] && theta<gtheta[ith+1]) {
        jth=ith;
      }
    }
  }

  dusttemp = gdusttemp[jr][jth];

  /*
    Find the UV flux from the RADMC output and calculate the attenuation.
  */
  g0 = gg0[jr][jth];
  g0att = g0;
  g0geom = gg0ref * (AU/r)*(AU/r);
  if(g0>g0geom) {
    av = 0.;
  } else {
    if(g0<1e-20*g0geom) {
      g0 = 1e-20*g0geom;
    }
    av = -log(g0/g0geom) / 3.02;
    g0 = g0geom;
  }
  
  /*
    Calculate the total gas density.
  */
  if(r<34.6*AU) {
    //inside inner edge: set density to zero
    dens = 1e-5;
  } else {
    dens = 1e6 * 1.1e9 * pow(r/(34.6*AU),-2.0);
  }
  cav_a = 34e3 * AU;
  cav_b = 8e3 * AU;
  if(fabs(z)>=2.*cav_a) {
    xycav = 0.0;
  }
  else {
    xycav = cav_b * sqrt(1.0-(fabs(z)/cav_a-1.0)*(fabs(z)/cav_a-1.0));
  }
  xy = sqrt(x*x+y*y);
  if(xy<xycav) {
    // inside the cavity
    dens = 1e-5;
  }
  
  /*
    Find the index and the weight for the density in the PDR temperature grid.
  */
  if(dens<=gtgdens[0]) {
    idens=0;
    epsdens=0.0;
  }
  else if(dens>=gtgdens[ntgn-1]) {
    idens=ntgn-2;
    epsdens=1.0;
  }
  else {
    for(i=0;i<=ntgn-2;i++) {
      if(dens>=gtgdens[i] && dens<gtgdens[i+1]) {
        idens=i;
        break;
      }
    }
    epsdens = (log(dens)-log(gtgdens[idens])) / (log(gtgdens[idens+1])-log(gtgdens[idens]));
    if(epsdens<0.0 || epsdens>1.0) {
      printf("Error with epsdens in function 'temperature'\n%i %e %e %e %e\n",idens,epsdens,dens,gtgdens[idens],gtgdens[idens+1]);
      exit(0);
    }
  }
  
  /*
    Find the index and the weight for G0 in the PDR temperature grid.
  */
  if(g0<=gtgg0[0]) {
    ig0=0;
    epsg0=0.0;
  }
  else if(g0>=gtgg0[ntgg-1]) {
    ig0=ntgg-2;
    epsg0=1.0;
  }
  else {
    for(i=0;i<=ntgg-2;i++) {
  if(g0>=gtgg0[i] && g0<gtgg0[i+1]) {
    ig0=i;
    break;
  }
    }
    epsg0 = (log(g0)-log(gtgg0[ig0])) / (log(gtgg0[ig0+1])-log(gtgg0[ig0]));
    if(epsg0<0.0 || epsg0>1.0) {
      printf("Error with epsg0 in function 'temperature'\n%i %e %e %e %e\n",ig0,epsg0,g0,gtgg0[ig0],gtgg0[ig0+1]);
      exit(0);
    }
  }
  
  /*
    Interpolate for G0 and density.
  */
  temp1 = (1.0-epsg0) * log(ggastemp[idens][ig0]) + epsg0 * log(ggastemp[idens][ig0+1]);
  temp2 = (1.0-epsg0) * log(ggastemp[idens+1][ig0]) + epsg0 * log(ggastemp[idens+1][ig0+1]);
  gastemp = (1.0-epsdens) * temp1 + epsdens * temp2;
  gastemp = exp(gastemp);

  /*
    Correct the temperature for extinction.
  */
  gastemp = gastemp * exp(-0.6*av);

  /*
    Make sure the gas temperature is not below the dust temperature.
  */
  if(gastemp<dusttemp) {
    gastemp = dusttemp;
  }

  temperature[1] = dusttemp;
  temperature[0] = gastemp;
}


void
abundance(double x, double y, double z, double *abundance){
#include "dusttemp.h"
#include "uvflux.h"
  FILE *in,*out;
  int ir,ith,jr,jth;
  double r,theta,dens,temp,temp0[2],g0,av,abun,g0geom,dens0[2];
  double cav_a,cav_b,xycav,xy;

  r=sqrt(x*x+y*y+z*z);
  theta=atan2(sqrt(x*x+y*y),fabs(z));

  /*
    Find the dust temperature from the RADMC grid and limit it according to
    values allowed in the chemical grid.
  */
  if(r<=gr[0]) {
    jr=0;
  }
  else if(r>=gr[nr-1]) {
    jr=nr-1;
  }
  else {
    for(ir=0;ir<=nr-2;ir++) {
      if(r>=gr[ir] && r<gr[ir+1]) {
        jr=ir;
      }
    }
  }

  if(theta<=gtheta[0]) {
    jth=0;
  }
  else if(theta>=gtheta[nth-1]) {
    jth=nth-1;
  }
  else {
    for(ith=0;ith<=nth-2;ith++) {
      if(theta>=gtheta[ith] && theta<gtheta[ith+1]) {
        jth=ith;
      }
    }
  }

  temp = gdusttemp[jr][jth];

  if(temp<10.) temp=10.;
  if(temp>1000.) temp=1000.;

  /*
    Calculate the visual extinction.
  */
  g0 = gg0[jr][jth];
  g0geom = gg0ref * (AU/r)*(AU/r);
  if(g0>g0geom) {
    av = 1e-3;
  } else {
    av = -log(g0/g0geom) / 3.02;
    g0 = g0geom;
  }
  if(g0<1e-1) g0=1e-1;
  if(g0>1e6) g0=1e6;
  if(av<1e-3) av=1e-3;
  if(av>1e2) av=1e2;

  /*
    Calculate the density in cm^-3 and limit it according to the chemical grid.
  */
  if(r<34.6*AU) {
    //inside inner edge: set density to zero
    dens = 1e-5;
  } else {
    dens = 1e6 * 1.1e9 * pow(r/(34.6*AU),-2.0);
  }
  cav_a = 34e3 * AU;
  cav_b = 8e3 * AU;
  if(fabs(z)>=2.*cav_a) {
    xycav = 0.0;
  }
  else {
    xycav = cav_b * sqrt(1.0-(fabs(z)/cav_a-1.0)*(fabs(z)/cav_a-1.0));
  }
  xy = sqrt(x*x+y*y);
  if(xy<xycav) {
    // inside the cavity
    dens = 1e-5;
  }
  dens = 2.0 * 1e-6 * dens; //from m^-3 to cm^-3 and from n(H2) to nH
  if(dens<1e2) dens=1e2;
  if(dens>1e10) dens=1e10;

  /*
    Write the input file for the interpolation.
  */
  out = fopen("grid_in.dat","a");
  fprintf(out,"1  %9.3E    %9.3E    %9.3E    %9.3E    %9.3E\n",dens,temp,g0,av,5e-17);
  fclose(out);

  /*
    Set the abundance to -1 as a flag to grid.c to get the actual abundances
    by an external program call.
  */
  abun = -1.;
  abundance[0] = abun;
}


void
doppler(double x, double y, double z, double *doppler){
  *doppler = 800.;
}


void
velocity(double x, double y, double z, double *vel){
  double r,theta,phi,vr;
  
  r = sqrt(x*x+y*y+z*z);
  theta = acos(z/r);
  phi = atan2(y,x);
  
  vr = -2000.0 * pow(r/(34.6*AU),-0.5);
  if(vr<-2000.0) vr=-2000.0;

  vel[0] = vr*sin(theta)*cos(phi);
  vel[1] = vr*sin(theta)*sin(phi);
  vel[2] = vr*cos(theta);
}


int
exclude (double x, double y, double z){ 
/*
  Function returns 1 if a point is inside the outflow cavity or 0 otherwise.
*/
  int exclude;
  double dist,cav_a,cav_b,xycav,xy;

  exclude = 0;
  dist = sqrt(x*x+y*y+z*z);
  if(dist<34.6*AU) {
    //inside inner edge: set density to zero
    exclude = 1;
  }

  /*
    Define the cavity wall.
  */
  cav_a = 34e3 * AU;
  cav_b = 8e3 * AU;
  if(fabs(z)>=2.*cav_a) {
    xycav = 0.0;
  }
  else {
    xycav = cav_b * sqrt(1.0-(fabs(z)/cav_a-1.0)*(fabs(z)/cav_a-1.0));
  }
  xy = sqrt(x*x+y*y);
  if(xy<xycav) {
    // inside the cavity
    exclude = 1;
  }

  return exclude;
}


double
angletocavity (double x, double y, double z){ 
  double cav_a,cav_b,angle,rsq,qa,qb,qc,qd,z0,xy0,angle0;

  cav_a = 34e3 * AU;
  cav_b = 8e3 * AU;

  rsq = x*x + y*y + z*z;
  angle = asin(fabs(z)/sqrt(rsq));

  qa = 1. - (cav_b*cav_b)/(cav_a*cav_a);
  qb = 2.*cav_b*cav_b/cav_a;
  qc = -rsq;
  qd = qb*qb - 4.*qa*qc;
  z0 = (-qb+sqrt(qd)) / (2.*qa);
  angle0 = asin(z0/sqrt(rsq));

  return angle/angle0;
}


double
cavity_phi (double r){ 
  double cav_a,cav_b,angle,rsq,qa,qb,qc,qd,z;

  cav_a = 34e3 * AU;
  cav_b = 8e3 * AU;

  rsq = r*r;

  qa = 1. - (cav_b*cav_b)/(cav_a*cav_a);
  qb = 2.*cav_b*cav_b/cav_a;
  qc = -rsq;
  qd = qb*qb - 4.*qa*qc;
  z = (-qb+sqrt(qd)) / (2.*qa);
  angle = acos(z/sqrt(rsq));

  return angle;
}


/******************************************************************/
