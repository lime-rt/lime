/*
 *  casa_messages.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2017 The LIME development team
 *
TODOs:
 */

#include "lime.h"
#include "casalime.h" /* For the status object and STATUS_STR_LEN */
#include <time.h>

/*....................................................................*/
void
greetings(void){
  snprintf(statusObj.message, STATUS_STR_LEN, "*** LIME, The versatile line modeling engine, version %s\n", VERSION);
}

/*....................................................................*/
void
greetings_parallel(int numThreads){
  if (numThreads>1){
    snprintf(statusObj.message, STATUS_STR_LEN, "*** LIME, The versatile line modeling engine, Ver. %s (parallel running, %d threads)\n", VERSION, numThreads);
  } else {
    snprintf(statusObj.message, STATUS_STR_LEN, "*** LIME, The versatile line modeling engine, Ver. %s\n", VERSION);
  }
}

/*....................................................................*/
void screenInfo(void){
  // do nothing
}

/*....................................................................*/
void
printDone(int line){
  if (line == 4){
    snprintf(statusObj.message, STATUS_STR_LEN, "   Building grid: DONE                               \n\n");
  }else if (line == 5){
    snprintf(statusObj.message, STATUS_STR_LEN, "   Smoothing grid: DONE                              \n\n");
//    statusObj.statusGrid = 1;
  }else if (line == 10){
    snprintf(statusObj.message, STATUS_STR_LEN, "   Propagating photons: DONE                         \n\n");
  }else if (line == 13){
    snprintf(statusObj.message, STATUS_STR_LEN, "   Raytracing model: DONE                            \n\n");
//    statusObj.statusRayTracing = 1;
  }else if (line == 15){
    snprintf(statusObj.message, STATUS_STR_LEN, "\n   Writing fits file: DONE                           \n\n");
  }
}

/*....................................................................*/
void
progressbar(double fraction, int line){
  if (line == 4){
    snprintf(statusObj.message, STATUS_STR_LEN, "   Building grid: %.2f percent done\r", fraction * 100.);
//    statusObj.progressGridBuilding = fraction;
  }else if (line == 5){
    snprintf(statusObj.message, STATUS_STR_LEN, "   Smoothing grid: %.2f percent done\r", fraction * 100.);
//    statusObj.progressGridSmoothing = fraction;
  }else if (line == 10){
    snprintf(statusObj.message, STATUS_STR_LEN, "   Propagating photons: %.2f percent done\r", fraction * 100.);
//    statusObj.progressPhotonPropagation = fraction;
  }else if (line == 13){
    snprintf(statusObj.message, STATUS_STR_LEN, "   Raytracing model: %.2f percent done\r", fraction * 100.);
//    statusObj.progressRayTracing = fraction;
  }else if (line == 15){
    snprintf(statusObj.message, STATUS_STR_LEN, "   Writing fits file\n");
  }
}

/*....................................................................*/
void
progressbar2(configInfo *par, int flag, int prog, double percent, double minsnr, double median){
  if (flag == 0) {
    snprintf(statusObj.message, STATUS_STR_LEN, "  Iteration %i / max %i: Starting\n", prog + 1, par->nSolveIters);
//    statusObj.minsnr = 0.0;
//    statusObj.median = 0.0;

  } else if (flag == 1){
    if (minsnr < 1000)
      snprintf(statusObj.message, STATUS_STR_LEN, "      Statistics: Min(SNR)    %3.3f                     \n", minsnr); 
    else 
      snprintf(statusObj.message, STATUS_STR_LEN, "      Statistics: Min(SNR)    %.3e                      \n", minsnr);

    if (median < 1000)
      snprintf(statusObj.message, STATUS_STR_LEN, "      Statistics: Median(SNR) %3.3f                     \n", median);
    else 
      snprintf(statusObj.message, STATUS_STR_LEN, "      Statistics: Median(SNR) %.3e                      \n", median);

//    statusObj.minsnr = minsnr;
//    statusObj.median = median;

    snprintf(statusObj.message, STATUS_STR_LEN, "  Iteration %i / max %i: DONE\n\n", prog+1, par->nSolveIters);
  }

//  statusObj.numberIterations = prog;
}

/*....................................................................*/
void casaStyleProgressBar(const int maxI, int i){
  // do nothing
}

/*....................................................................*/
void
reportOutput(char *filename){
  snprintf(statusObj.message, STATUS_STR_LEN, "Output written to %s\n", filename);
}

/*....................................................................*/
void
goodnight(int initime){
  int runtime=time(0)-initime;
  snprintf(statusObj.message, STATUS_STR_LEN, "*** Runtime: %3dh %2dm %2ds\n\n", runtime / 3600, runtime / 60 % 60, runtime % 60);
//  statusObj.statusGlobal = 1;
}

/*....................................................................*/
void
printMessage(char *message){
  if(strlen(message)>0)
    snprintf(statusObj.message, STATUS_STR_LEN, "%s\n", message );
}

/*....................................................................*/
void
warning(char *message){
  if(strlen(message)>0)
    snprintf(statusObj.message, STATUS_STR_LEN, "Warning : %s\n", message );
}

/*....................................................................*/
void error(char *message){
  if(!silent) bail_out(message);
  exit(1);
}

/*....................................................................*/
void
bail_out(char *message){
  snprintf(statusObj.message, STATUS_STR_LEN, "Error: %s\n", message );
}

/*....................................................................*/
void
collpartmesg(char *molecule, int partners){//, int specnumber){
  if (partners==1)
    snprintf(statusObj.message, STATUS_STR_LEN, "   Molecule: %.25s\n   %d collision partner:\n", molecule, partners);
  else
    snprintf(statusObj.message, STATUS_STR_LEN, "   Molecule: %.25s\n   %d collision partners:\n", molecule, partners);
}

/*....................................................................*/
void
collpartmesg2(char name[10]){
  snprintf(statusObj.message, STATUS_STR_LEN, "      %s\n ", name);
}

/*....................................................................*/
void
collpartmesg3(int number, int flag){
  if (number==1){
    if(flag==1)
      snprintf(statusObj.message, STATUS_STR_LEN, "   Model provides: %d density profile\n\n*** Warning! ***: Too few density profiles", number);
    else
      snprintf(statusObj.message, STATUS_STR_LEN, "   Model provides: %d density profile\n\n", number);
  }else{
    if(flag==1)
      snprintf(statusObj.message, STATUS_STR_LEN, "   Model provides: %d density profiles\n\n*** Warning! ***: Too few density profiles", number);
    else
      snprintf(statusObj.message, STATUS_STR_LEN, "   Model provides: %d density profiles\n\n", number);
  }
}

/*....................................................................*/
void
processFitsError(int status){
  /*****************************************************/
  /* Print out cfitsio error messages and exit program */
  /*****************************************************/

  if(status){
    if(!silent){
      fits_report_error(stderr, status); /* print error report */
    }
    exit( status );    /* terminate the program, returning error status */
  }
  return;
}

