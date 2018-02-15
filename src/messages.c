/*
 *  messages.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2017 The LIME development team
 *
TODOs:
	- Are all these char arguments better as pointers, or arrays of unspecified length? Do we need to specify the length?
	- Print a blank line of this len before each curses-style warning or message?
	- Wouldn't it be better if 'silent' was tested inside these functions rather than at every single point in the rest of the code where they are called?
 */

#include "constants.h" // for STR_LEN_0 etc
#include "messages.h"

#ifdef NO_STDOUT
#include "pyshared_io.h"
#endif

/*....................................................................*/
void
greetings(char *version){
#ifdef NO_NCURSES

  printf("*** LIME, The versatile line modeling engine, version %s\n", version);
#ifdef TEST
  printf(">>> NOTE! Test flag is set in the Makefile. <<<\n");
#endif
#ifdef FASTEXP
  printf(">>> NOTE! Fast-exponential routine is enabled. <<<\n");
#endif

#elif defined (NO_STDOUT)

  snprintf(statusObj.message, STATUS_STR_LEN, "*** LIME, The versatile line modeling engine, version %s\n", version);

#else

  initscr();
  printw("*** LIME, The versatile line modeling engine, version %s\n", version);
#ifdef TEST
  printw(">>> NOTE! Test flag is set in the Makefile. <<<\n");
#endif
#ifdef FASTEXP
  printw(">>> NOTE! Fast-exponential routine is enabled. <<<\n");
#endif
  refresh();

#endif
}

/*....................................................................*/
void
greetings_parallel(int numThreads, char *version){
#ifdef NO_NCURSES

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

#elif defined (NO_STDOUT)

  if (numThreads>1){
    snprintf(statusObj.message, STATUS_STR_LEN, "*** LIME, The versatile line modeling engine, Ver. %s (parallel running, %d threads)\n", version, numThreads);
  } else {
    snprintf(statusObj.message, STATUS_STR_LEN, "*** LIME, The versatile line modeling engine, Ver. %s\n", version);
  }

#else

  initscr();
  if (numThreads>1){
    printw("*** LIME, The versatile line modeling engine, Ver. %s (parallel running, %d threads)\n", version, numThreads);
  } else {
    printw("*** LIME, The versatile line modeling engine, Ver. %s\n", version);
  }
#ifdef TEST
  printw(">>> NOTE! Test flag is set in the Makefile. <<<\n");
#endif
#ifdef FASTEXP
  printw(">>> NOTE! Fast-exponential routine is enabled. <<<\n");
#endif
  refresh();

#endif
}

/*....................................................................*/
void
screenInfo(){
#ifndef NO_NCURSES
  move(4,4);  printw("Building grid      :");
  move(4,51); printw("|");
  move(5,4);  printw("Smoothing grid     :");
  move(5,51); printw("|");
  move(7,4);  printw("Statistics         :");
  move(9,4);  printw("Iterations         :");
  move(10,4); printw("Propagating photons:");
  move(10,51);printw("|");
  move(13,4); printw("Ray-tracing model  :");
  move(13,51);printw("|");
  move(4,60); printw("|      Molecular data");
  move(5,60); printw("|");
  move(6,60); printw("|");
  move(7,60); printw("|");
  move(8,60); printw("|");
  move(9,60); printw("|");
  move(10,60); printw("|");
  move(11,60); printw("|");
  move(12,60); printw("|");
  move(13,60); printw("|");
  move(14,60); printw("|");
  refresh();	
#endif
}

/*....................................................................*/
void
printDone(int line){
#ifdef NO_NCURSES
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

#elif defined (NO_STDOUT)

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

#else
  move(line,52); printw(" [ok]");
  refresh();
#endif
}

/*....................................................................*/
void
progressbar(double percent, int line){
#ifdef NO_NCURSES
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

#elif defined (NO_STDOUT)

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

#else
  int i;
  for(i=0;i<(int)(percent*25.);i++){
    move(line,25+i);
    printw("#");
  }
  refresh();
#endif
}

/*....................................................................*/
void
progressbar2(int nSolveIters, int flag, int prog, double percent, double minsnr, double median){
#ifdef NO_NCURSES
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

#elif defined (NO_STDOUT)

  if (flag == 0) {
    snprintf(statusObj.message, STATUS_STR_LEN, "  Iteration %i / max %i: Starting\n", prog + 1, nSolveIters);
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

    snprintf(statusObj.message, STATUS_STR_LEN, "  Iteration %i / max %i: DONE\n\n", prog+1, nSolveIters);
  }

//  statusObj.numberIterations = prog;

#else
  if (flag == 0) {
    move(9,25+prog); printw("#");
    if(percent<100) {
      move(10,25); printw("                         ");
    }
    refresh();
  } else if (flag == 1){
    move(7,38); printw("                    ");            
    move(8,38); printw("                    ");
    if(minsnr<1000){
      move(7,25); printw("Min(SNR)    %3.3f", minsnr);
    } else {
      move(7,25); printw("Min(SNR)    %.3e", minsnr);
    }
    if(median<1000){
      move(8,25); printw("Median(SNR) %3.3f", median);
    } else {
      move(8,25); printw("Median(SNR) %.3e", median);
    }
    refresh();
  }
#endif
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
#ifdef NO_NCURSES
  printf("Output written to %s\n", filename);
#elif defined (NO_STDOUT)
  snprintf(statusObj.message, STATUS_STR_LEN, "Output written to %s\n", filename);
#else
  move(14,4); printw("Output written to %s", filename);
  refresh();
#endif
}

/*....................................................................*/
void
goodnight(int initime){
  int runtime=time(0)-initime;
#ifdef NO_NCURSES
  printf("*** Program ended successfully               \n");
  printf("    Runtime: %3dh %2dm %2ds\n\n", runtime / 3600, runtime / 60 % 60, runtime % 60);

#elif defined (NO_STDOUT)
  snprintf(statusObj.message, STATUS_STR_LEN, "*** Runtime: %3dh %2dm %2ds\n\n", runtime / 3600, runtime / 60 % 60, runtime % 60);
//  statusObj.statusGlobal = 1;

#else
  move(22,0); printw("*** Program ended successfully               ");
  move(22,58); printw("runtime: %3dh %2dm %2ds", runtime/3600, runtime/60%60, runtime%60);
  move(23,0); printw("*** [Press any key to quit]");
  refresh();
  getch();
  endwin();
#endif
}

/*....................................................................*/
void
printMessage(char *message){
#ifdef NO_NCURSES
  if(strlen(message)>0)
    {
      printf("%s\n", message );
    }
#elif defined (NO_STDOUT)
  if(strlen(message)>0)
    snprintf(statusObj.message, STATUS_STR_LEN, "%s\n", message );
#else
  move(22,0); printw("*** %s\n",message);
  refresh();
#endif
}

/*....................................................................*/
void
warning(char *message){
#ifdef NO_NCURSES
  if(strlen(message)>0)
    {
      printf("Warning : %s\n", message );
    }
#elif defined (NO_STDOUT)
  if(strlen(message)>0)
    snprintf(statusObj.message, STATUS_STR_LEN, "Warning : %s\n", message );
#else
  move(22,0); printw("*** %s\n",message);
  refresh();
#endif
}

/*....................................................................*/
void error(char *message){
  if(!silent) bail_out(message);
  exit(1);
}

/*....................................................................*/
void
bail_out(char *message){
#ifdef NO_NCURSES
  printf("Error: %s\n", message );
#elif defined (NO_STDOUT)
  snprintf(statusObj.message, STATUS_STR_LEN, "Error: %s\n", message );
#else
  move(22,0); printw("*** %s",message);
  move(23,0); printw("*** [Press any key to quit]");
  refresh();
  getch();
  endwin();
#endif
}

/*....................................................................*/
void
collpartmesg(char molecule[STR_LEN_0], int partners){//, int specnumber){
#ifdef NO_NCURSES
  printf("   Molecule: %.25s\n", molecule);
  if (partners==1)
    printf("   %d collision partner:\n", partners);
  else
    printf("   %d collision partners:\n", partners);

#elif defined (NO_STDOUT)
  if (partners==1)
    snprintf(statusObj.message, STATUS_STR_LEN, "   Molecule: %.25s\n   %d collision partner:\n", molecule, partners);
  else
    snprintf(statusObj.message, STATUS_STR_LEN, "   Molecule: %.25s\n   %d collision partners:\n", molecule, partners);

#else
  move(6,63); printw("%.25s", molecule);
  move(7,63);
  if (partners==1)
    printw("%d collision partner:", partners);
  else
    printw("%d collision partners:", partners);

  refresh();
#endif
}

/*....................................................................*/
void
collpartmesg2(char name[10]){
#ifdef NO_NCURSES
  printf("      %s\n ", name);
#elif defined (NO_STDOUT)
  snprintf(statusObj.message, STATUS_STR_LEN, "      %s\n ", name);
#else
  move(8,63); printw("%s ",name);
  refresh();
#endif
}

/*....................................................................*/
void
collpartmesg3(int number, int flag){
#ifdef NO_NCURSES
  if (number==1)
    printf("   Model provides: %d density profile\n\n", number);
  else
    printf("   Model provides: %d density profiles\n\n", number);

  if(flag==1) printf("*** Warning! ***: Too few density profiles");  
#elif defined (NO_STDOUT)
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
#else
  move(10,63); printw("Model provides:");
  move(11,63);
  if (number==1)
    printw("%d density profile", number);
  else
    printw("%d density profiles", number);

  if(flag==1) {
    move(13,63); printw("*** Warning! ***");
    move(14,63); printw("Too few density profiles");
  }
  refresh();
#endif
}

