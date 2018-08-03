/*
 *  messages.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
TODOs:
	- Are all these char arguments better as pointers, or arrays of unspecified length? Do we need to specify the length?
	- Print a blank line of this len before each curses-style warning or message?
	- Wouldn't it be better if 'silent' was tested inside these functions rather than at every single point in the rest of the code where they are called?
 */

#include "constants.h" // for STR_LEN_0 etc
#include "messages.h"
#include "local_err.h"

/*....................................................................*/
errType
init_local_err(void){
  errType err;

  err.status = 0;
  err.message[0] = '\0';

return err;
}

/*....................................................................*/
errType
write_local_err(int status, char *message){
  errType err;

  err.status = status;
  strncpy(err.message, message, ERR_STR_LEN);

return err;
}

/*....................................................................*/
void
error(int exitStatus, char *message){
  printf("ERROR: %s\n", message);
exit(exitStatus);
}

/*....................................................................*/
void
greetings(char *version){
  printf("*** LIME, The versatile line modeling engine, version %s\n", version);
#ifdef TEST
  printf(">>> NOTE! Test flag is set in the Makefile. <<<\n");
#endif
#ifdef FASTEXP
  printf(">>> NOTE! Fast-exponential routine is enabled. <<<\n");
#endif
}

/*....................................................................*/
void
greetings_parallel(int numThreads, char *version){
  if (numThreads>1){
    printf("*** LIME, The versatile line modeling engine, Ver. %s (parallel running, %d threads)\n", version, numThreads);
    fflush(stdout);
  } else {
    printf("*** LIME, The versatile line modeling engine, Ver. %s\n", version);
    fflush(stdout);
  }
#ifdef TEST
  printf(">>> NOTE! Test flag is set in the Makefile. <<<\n");
#endif
#ifdef FASTEXP
  printf(">>> NOTE! Fast-exponential routine is enabled. <<<\n");
#endif
}

/*....................................................................*/
void
screenInfo(){
  /* NOP */
}

/*....................................................................*/
void
printDone(int line){
  if (line == 4){
    printf(  "   Building grid: DONE                               \n\n"); 
    fflush(stdout);
  }else if (line == 5){
    printf(  "   Smoothing grid: DONE                              \n\n");
    fflush(stdout);
  }else if (line == 10){
    printf(  "   Propagating photons: DONE                         \n\n");
    fflush(stdout);
  }else if (line == 13){
    printf(  "   Raytracing model: DONE                            \n\n");
    fflush(stdout);
  }else if (line == 15){
    printf("\n   Writing fits file: DONE                           \n\n");
    fflush(stdout);
  }
}

/*....................................................................*/
void
progressbar(double percent, int line){
  if (line == 4)
    printf("   Building grid: %.2f percent done\r", percent * 100.);
  else if (line == 5){
    printf("   Smoothing grid: %.2f percent done\r", percent * 100.);
    fflush(stdout);
  }else if (line == 10){
    printf("   Propagating photons: %.2f percent done\r", percent * 100.);
    fflush(stdout);
  }else if (line == 13)
    printf("   Raytracing model: %.2f percent done\r", percent * 100.);
  else if (line == 15)
    printf("   Writing fits file\n");
}

/*....................................................................*/
void
progressbar2(int nSolveIters, int flag, int prog, double percent, double minsnr, double median){
  if (flag == 0) {
    printf("  Iteration %i / max %i: Starting\n", prog + 1, nSolveIters);
    fflush(stdout);
  } else if (flag == 1){
    if (minsnr < 1000)
      printf("      Statistics: Min(SNR)    %3.3f                     \n", minsnr); 
    else 
      printf("      Statistics: Min(SNR)    %.3e                      \n", minsnr);
    fflush(stdout);

    if (median < 1000)
      printf("      Statistics: Median(SNR) %3.3f                     \n", median);
    else 
      printf("      Statistics: Median(SNR) %.3e                      \n", median);
    fflush(stdout);

    printf("  Iteration %i / max %i: DONE\n\n", prog+1, nSolveIters);
    fflush(stdout);
  }
}

/*....................................................................*/
void casaStyleProgressBar(const int maxI, int i){
  static int minorCounter=0,majorCounter=0;
  static float counter=0.0;
  const int minorsPerMajor=5, maxMajor=10, minorInterval=2; /* product must be 100 */
  int percentI;

  while (counter<=i){
    /* decide whether to print minor or major symbol. */
    if (minorCounter==majorCounter){
      percentI = minorCounter*minorInterval;
      if (percentI==100){
        printf("%d%%\n", percentI);
      }else{
        printf("%d%%", percentI);
        fflush(stdout);
      }

      majorCounter += minorsPerMajor;
    } else { /* assume minorCounter<majorCounter, because I can't see how it could be >! */
      printf(".");
      fflush(stdout);
    }

    minorCounter++;
    counter = (maxI-1)*(minorCounter/(float)(minorsPerMajor*maxMajor));
  }
}

/*....................................................................*/
void
reportOutput(char filename[STR_LEN_0]){
  printf("Output written to %s\n", filename);
}

/*....................................................................*/
void
goodnight(int initime){
  int runtime=time(0)-initime;
  printf("*** Program ended successfully               \n");
  printf("    Runtime: %3dh %2dm %2ds\n\n", runtime / 3600, runtime / 60 % 60, runtime % 60);
}

/*....................................................................*/
void
printMessage(char *message){
  if(strlen(message)>0)
    {
      printf("%s\n", message );
    }
}

/*....................................................................*/
void
warning(char *message){
  if(strlen(message)>0)
    {
      printf("Warning : %s\n", message );
    }
}

/*....................................................................*/
//void error(char *message){
//  if(!silent) bail_out(message);
//  exit(1);
//}

/*....................................................................*/
void
bail_out(char *message){
  printf("Error: %s\n", message );
}

/*....................................................................*/
void
collpartmesg(char molecule[STR_LEN_0], int partners){//, int specnumber){
  printf("   Molecule: %.25s\n", molecule);
  if (partners==1)
    printf("   %d collision partner:\n", partners);
  else
    printf("   %d collision partners:\n", partners);
}

/*....................................................................*/
void
collpartmesg2(char name[10]){
  printf("      %s\n ", name);
}

/*....................................................................*/
void
collpartmesg3(int number, int flag){
  if (number==1)
    printf("   Model provides: %d density profile\n\n", number);
  else
    printf("   Model provides: %d density profiles\n\n", number);

  if(flag==1) printf("*** Warning! ***: Too few density profiles");  
}

