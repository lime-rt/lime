/*
 *  messages.h
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#ifndef MESSAGES_H
#define MESSAGES_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifdef NCURSES
#include <curses.h>
#endif
#include <time.h>

void	bail_out(char*);
void	casaStyleProgressBar(const int, int);
void	collpartmesg(char*, int);
void	collpartmesg2(char*);
void	collpartmesg3(int, int);
void	goodnight(int);
void	greetings(char*);
void	greetings_parallel(int, char*);
void	printDone(int);
void	printMessage(char *);
void	progressbar(double, int);
//void	progressbar2(configInfo*, int, int, double, double, double);
void	progressbar2(int, int, int, double, double, double);
void	reportOutput(char*);
void	screenInfo(void);
void	warning(char*);

#endif /* MESSAGES_H */

