/*
 *  stokesangles.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015 The LIME development team
 *
 */

#include "lime.h"

void
stokesangles(double x, double y, double z, double incl, double *trigFuncs){

	double B[3],Bp[3];
	double sPsi, cPsi;
	double cGam,cosIncl,sinIncl;
	
	magfield(x,y,z,B);

        cosIncl=cos(incl);
        sinIncl=sin(incl);
	Bp[0] =  B[0];
	Bp[1] =  B[1]*cosIncl - B[2]*sinIncl;
	Bp[2] =  B[1]*sinIncl + B[2]*cosIncl;
	
	if (!(Bp[0]*Bp[0]+Bp[1]*Bp[1])==0) {
		cPsi = Bp[1]/sqrt(Bp[0]*Bp[0]+Bp[1]*Bp[1]);
	}
	else {
		cPsi = 0.;
	}
	
	sPsi = sqrt(1-cPsi*cPsi);
	if (B[0]>0) {
		sPsi=-sPsi;
	}
	
	if (!(Bp[0]*Bp[0]+Bp[1]*Bp[1]+Bp[2]*Bp[2])==0) {
		cGam = sqrt((Bp[0]*Bp[0]+Bp[1]*Bp[1])/(Bp[0]*Bp[0]+Bp[1]*Bp[1]+Bp[2]*Bp[2]));
	}
	else {
		cGam = 0.;
	}
	trigFuncs[0] = cPsi;	//cosinus of Psi
	trigFuncs[1] = sPsi;	//sinus of Psi
	trigFuncs[2] = cGam;	//cosinus of Gamma
	
}
