/*
 *  writefits.c
 *  LIME, The versatile 3D line modeling environment 
 *
 *  Created by Christian Brinch on 24/04/08.
 *  Copyright 2006-2011, Christian Brinch, 
 *  <brinch@nbi.dk>
 *  Niels Bohr institutet
 *  University of Copenhagen
 *	All rights reserved.
 *
 */

#include <cfitsio/fitsio.h>
#include "lime.h"

void 
writefits(int im, inputPars *par, molData *m, image *img){
	double keyval;
	float *row;
	int px,py,ichan;
	fitsfile *fptr;
	int status = 0;
	int naxis=3, bitpix=-32;
	long naxes[3];
	long int fpixels[3],lpixels[3];
	char negfile[100]="! ";

	
	row = malloc(sizeof(float)*img[im].pxls);
	
	naxes[0]=img[im].pxls;
	naxes[1]=img[im].pxls;
	if(img[im].doline==1) naxes[2]=img[im].nchan;
	else if(img[im].doline==0 && par->polarization) naxes[2]=3;
	else naxes[2]=1;

	fits_create_file(&fptr, img[im].filename, &status);
	
	if(status!=0){
		if(!silent) warning("Overwriting existing fits file                   ");
		status=0;
		strcat(negfile,img[im].filename);
		fits_create_file(&fptr, negfile, &status);
	}

	/* Write FITS header */ 
	fits_create_img(fptr, bitpix, naxis, naxes, &status);
	
	keyval=1.e0;
	fits_write_key(fptr, TDOUBLE, "BSCALE  ", &keyval, "", &status);
	keyval=0.e0;
	fits_write_key(fptr, TDOUBLE, "BZERO   ", &keyval, "", &status);
	if(img[im].unit==0) fits_write_key(fptr, TSTRING, "BUNIT", &"K       ", "", &status);
	if(img[im].unit==1) fits_write_key(fptr, TSTRING, "BUNIT", &"JY/PIXEL", "", &status);
	if(img[im].unit==2) fits_write_key(fptr, TSTRING, "BUNIT", &"WM2HZSR ", "", &status);
	if(img[im].unit==3) fits_write_key(fptr, TSTRING, "BUNIT", &"Lsun/PX ", "", &status);
	if(img[im].unit==3) fits_write_key(fptr, TSTRING, "BUNIT", &"        ", "", &status);
	fits_write_key(fptr, TSTRING, "OBJECT  ", &"LIMEMDL ", "", &status);
	keyval=2.e3;
	fits_write_key(fptr, TDOUBLE, "EPOCH   ", &keyval, "", &status);
	keyval=1.8e2;
	fits_write_key(fptr, TDOUBLE, "LONPOLE ", &keyval, "", &status);
	keyval=img[im].freq;
	fits_write_key(fptr, TDOUBLE, "RESTFREQ", &keyval, "", &status);
		
	fits_write_key(fptr, TSTRING, "CTYPE1", &"RA---SIN", "", &status);
  	keyval=-1.8e2*img[im].imgres/PI;
	fits_write_key(fptr, TDOUBLE, "CDELT1", &keyval, "", &status);
  	keyval=(double) img[im].pxls/2+0.5;
	fits_write_key(fptr, TDOUBLE, "CRPIX1", &keyval, "", &status);
  	keyval=0.e0;
	fits_write_key(fptr, TDOUBLE, "CRVAL1", &keyval, "", &status);
	fits_write_key(fptr, TSTRING, "CUNIT1", &"DEG     ", "", &status);
		
	fits_write_key(fptr, TSTRING, "CTYPE2", &"DEC--SIN", "", &status);
  	keyval=1.8e2*img[im].imgres/PI;
	fits_write_key(fptr, TDOUBLE, "CDELT2", &keyval, "", &status);
  	keyval=(double) img[im].pxls/2+0.5;
	fits_write_key(fptr, TDOUBLE, "CRPIX2", &keyval, "", &status);  	
  	keyval=0.e0;
	fits_write_key(fptr, TDOUBLE, "CRVAL2", &keyval, "", &status);
	fits_write_key(fptr, TSTRING, "CUNIT2", &"DEG     ", "", &status);

	fits_write_key(fptr, TSTRING, "CTYPE3", &"VELO-LSR", "", &status);  
  	keyval=img[im].velres;
	fits_write_key(fptr, TDOUBLE, "CDELT3", &keyval, "", &status);
 	keyval=(double) img[im].nchan/2.+1;
	fits_write_key(fptr, TDOUBLE, "CRPIX3", &keyval, "", &status);
  	keyval=0.e0;
	fits_write_key(fptr, TDOUBLE, "CRVAL3", &keyval, "", &status);
	fits_write_key(fptr, TSTRING, "CUNIT3", &"M/S     ", "", &status);


	/* Write FITS data */
	for(py=0;py<img[im].pxls;py++){
		for(ichan=0;ichan<img[im].nchan;ichan++){
			for(px=0;px<img[im].pxls;px++){
			  if(img[im].unit==0) row[px]=(float) img[im].pixel[px+py*img[im].pxls].intense[ichan]*pow(CLIGHT/img[im].freq,2)/2./KBOLTZ*m[0].norm; 
			  if(img[im].unit==1) row[px]=(float) img[im].pixel[px+py*img[im].pxls].intense[ichan]*1e26*img[im].imgres*img[im].imgres*m[0].norm; 
			  if(img[im].unit==2) row[px]=(float) img[im].pixel[px+py*img[im].pxls].intense[ichan]*m[0].norm; 
			  if(img[im].unit==3) row[px]=(float) img[im].pixel[px+py*img[im].pxls].intense[ichan]*4.*PI*pow(img[im].distance/1.975e13,2)*img[im].freq*img[im].imgres*img[im].imgres*m[0].norm; 
			  if(img[im].unit==4) row[px]=(float) img[im].pixel[px+py*img[im].pxls].tau[ichan]; 
			  if (fabs(row[px])<(float) eps) row[px]=(float)eps;
			}
			fpixels[0]=1;
			fpixels[1]=py+1;
			fpixels[2]=ichan+1;
			lpixels[0]=img[im].pxls;
			lpixels[1]=py+1;
			lpixels[2]=ichan+1;
			fits_write_subset(fptr, TFLOAT, fpixels, lpixels, row, &status);
		}
	}
	fits_close_file(fptr, &status);
	
	free(row);
	
	if(!silent) done(13);
}
