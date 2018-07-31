/*
 *  messages.h
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2017 The LIME development team
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
void	error(char*);
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

