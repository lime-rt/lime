/*
 *  curses.c
 *  LIME, The versatile 3D line modeling environment 
 *
 *  Created by Christian Brinch on 29/10/08.
 *  Copyright 2006-2012, Christian Brinch, 
 *  <brinch@nbi.dk>
 *  Niels Bohr institutet
 *  University of Copenhagen
 *	All rights reserved.
 *
 */

#include <curses.h>
#include <time.h>

void
greetings(){
	initscr();
	printw("*** LIME, The versatile line modeling engine, Ver.1.21 (Devel)\n*** Copyright 2006--2012, Christian Brinch <brinch@nbi.dk>\n");
	refresh();	
}

void
screenInfo(){
	move(4,6);  printw("Building Delaunay grid:");
	move(4,56); printw("|");
	move(5,6);  printw("Smoothing grid        :");
	move(5,56); printw("|");
	move(7,6);  printw("Statistics            :");
	move(9,6);  printw("Iterations            :"); 
	move(10,6); printw("Propagating photons   :");
	move(10,56);printw("|");
	move(13,6); printw("Ray-tracing model     :");
	move(13,56);printw("|");	
	refresh();	
}

void
done(int line){
	move(line,57); printw(" [done]");
    refresh();
}

void
progressbar(double percent, int line){
  int i;
//if((int)(100*percent)%4==0){
    for(i=0;i<(int)(percent*25.);i++){
      move(line,30+i);
      printw("#");
    }
    refresh();
//  }	
}

void 
progressbar2(int prog, double percent, double minsnr, double median){
	move(7,30);		 printw("Min(SNR)    %3.3f", minsnr);
	move(8,30);		 printw("Median(SNR) %4.3f",median);
	move(9,30+prog); printw("#");
	if(percent<100) {
		move(10,30);	 printw("                         ");
	}
	refresh();	
}

void 
goodnight(int initime, char filename[80]){
	int runtime=time(0)-initime;
	move(14,6); printw("Output written to %s", filename);
	move(22,0); printw("*** Program ended succesfully               ");
	move(22,58); printw("runtime: %3dh %2dm %2ds", runtime/3600, runtime/60%60, runtime%60);
	move(23,0); printw("*** [Press any key to quit]");
    refresh();
	getch();
	endwin();
}

void
quotemass(double mass){
  move(21,6); printw("Total mass contained in model: %3.2e solar masses", mass);
  refresh();
}



void 
warning(char message[80]){
	move(22,0); printw("*** %s\n",message);
	refresh();
}

void 
bail_out(char message[80]){
	move(22,0); printw("*** %s",message);
	move(23,0); printw("*** [Press any key to quit]");
    refresh();
	getch();
	endwin();
}
