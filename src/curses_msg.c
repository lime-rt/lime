/*
 *  curses_msg.c
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

/*....................................................................*/
void
greetings(char *version){
#ifdef TEST
  printf(">>> NOTE! Test flag is set in the Makefile. <<<\n");
#endif
#ifdef FASTEXP
  printf(">>> NOTE! Fast-exponential routine is enabled. <<<\n");
#endif

  initscr();
  printw("*** LIME, The versatile line modeling engine, version %s\n", version);
#ifdef TEST
  printw(">>> NOTE! Test flag is set in the Makefile. <<<\n");
#endif
#ifdef FASTEXP
  printw(">>> NOTE! Fast-exponential routine is enabled. <<<\n");
#endif
  refresh();
}

/*....................................................................*/
void
greetings_parallel(int numThreads, char *version){
#ifdef TEST
  printf(">>> NOTE! Test flag is set in the Makefile. <<<\n");
#endif
#ifdef FASTEXP
  printf(">>> NOTE! Fast-exponential routine is enabled. <<<\n");
#endif

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
}

/*....................................................................*/
void
screenInfo(){
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
}

/*....................................................................*/
void
printDone(int line){
  move(line,52); printw(" [ok]");
  refresh();
}

/*....................................................................*/
void
progressbar(double percent, int line){
  int i;
  for(i=0;i<(int)(percent*25.);i++){
    move(line,25+i);
    printw("#");
  }
  refresh();
}

/*....................................................................*/
void
progressbar2(int nSolveIters, int flag, int prog, double percent, double minsnr, double median){
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
}

/*....................................................................*/
void casaStyleProgressBar(const int maxI, int i){
  /* NOP */
}

/*....................................................................*/
void
reportOutput(char filename[STR_LEN_0]){
  move(14,4); printw("Output written to %s", filename);
  refresh();
}

/*....................................................................*/
void
goodnight(int initime){
  int runtime=time(0)-initime;
  move(22,0); printw("*** Program ended successfully               ");
  move(22,58); printw("runtime: %3dh %2dm %2ds", runtime/3600, runtime/60%60, runtime%60);
  move(23,0); printw("*** [Press any key to quit]");
  refresh();
  getch();
  endwin();
}

/*....................................................................*/
void
printMessage(char *message){
  move(22,0); printw("*** %s\n",message);
  refresh();
}

/*....................................................................*/
void
warning(char *message){
  move(22,0); printw("*** %s\n",message);
  refresh();
}

/*....................................................................*/
void error(char *message){
  if(!silent) bail_out(message);
  exit(1);
}

/*....................................................................*/
void
bail_out(char *message){
  move(22,0); printw("*** %s",message);
  move(23,0); printw("*** [Press any key to quit]");
  refresh();
  getch();
  endwin();
}

/*....................................................................*/
void
collpartmesg(char molecule[STR_LEN_0], int partners){//, int specnumber){
  move(6,63); printw("%.25s", molecule);
  move(7,63);
  if (partners==1)
    printw("%d collision partner:", partners);
  else
    printw("%d collision partners:", partners);

  refresh();
}

/*....................................................................*/
void
collpartmesg2(char name[10]){
  move(8,63); printw("%s ",name);
  refresh();
}

/*....................................................................*/
void
collpartmesg3(int number, int flag){
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
}

