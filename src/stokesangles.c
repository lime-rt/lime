/*
 *  stokesangles.c
 *  LIME, The versatile 3D line modeling environment 
 *
 *  Created by Marco Padovani on 06/04/10.
 *  Copyright 2006-2012, Marco Padovani 
 *
 *  Institut de Ciències de l'Espai (IEEC-CSIC) - 
 *  UAB - Catalunya, Spain. 
 *  All rights reserved.
 *
 */

#include "lime.h"

void
stokesangles(double x, double y, double z, double incl, double *angle){

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
	angle[0] = cPsi;	//cosinus of Psi
	angle[1] = sPsi;	//sinus of Psi
	angle[2] = cGam;	//cosinus of Gamma
	
}
