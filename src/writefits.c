/*
 *  writefits.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015 The LIME development team
 *
 */

#include "lime.h"

void 
writefits(int im, inputPars *par, molData *m, image *img){
  double bscale,bzero,epoch,lonpole,equinox,restfreq;
  double cdelt1,crpix1,crval1,cdelt2,crpix2,crval2;
  double cdelt3,crpix3,crval3,ru3;
  int velref;
  float *row;
  int px,py,ichan;
  fitsfile *fptr;
  int status = 0;
  int naxis=3, bitpix=-32;
  long naxes[3];
  long int fpixels[3],lpixels[3];
  char negfile[100]="! ";

  row = malloc(sizeof(*row)*img[im].pxls);

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
  epoch   =2.0e3;
  lonpole =1.8e2;
  equinox =2.0e3;
  restfreq=img[im].freq;
  velref  =257;
  cdelt1  =-1.8e2*img[im].imgres/PI;
  crpix1  =(double) img[im].pxls/2+0.5;
  crval1  =0.0e0;
  cdelt2  =1.8e2*img[im].imgres/PI;
  crpix2  =(double) img[im].pxls/2+0.5;
  crval2  =0.0e0;
  cdelt3  =img[im].velres;
  crpix3  =(double) (img[im].nchan-1)/2.+1;
  crval3  =0.0e0;
  bscale  =1.0e0;
  bzero   =0.0e0;

  fits_write_key(fptr, TSTRING, "OBJECT  ", &"LIMEMDL ",    "", &status);
  fits_write_key(fptr, TDOUBLE, "EPOCH   ", &epoch,         "", &status);
  fits_write_key(fptr, TDOUBLE, "LONPOLE ", &lonpole,       "", &status);
  fits_write_key(fptr, TDOUBLE, "EQUINOX ", &equinox,       "", &status);
  fits_write_key(fptr, TSTRING, "SPECSYS ", &"LSRK    ",    "", &status);
  fits_write_key(fptr, TDOUBLE, "RESTFREQ", &restfreq,      "", &status);
  fits_write_key(fptr, TINT,    "VELREF  ", &velref,        "", &status);
  fits_write_key(fptr, TSTRING, "CTYPE1  ", &"RA---SIN",    "", &status);
  fits_write_key(fptr, TDOUBLE, "CDELT1  ", &cdelt1,        "", &status);
  fits_write_key(fptr, TDOUBLE, "CRPIX1  ", &crpix1,        "", &status);
  fits_write_key(fptr, TDOUBLE, "CRVAL1  ", &crval1,        "", &status);
  fits_write_key(fptr, TSTRING, "CUNIT1  ", &"DEG     ",    "", &status);	
  fits_write_key(fptr, TSTRING, "CTYPE2  ", &"DEC--SIN",    "", &status);
  fits_write_key(fptr, TDOUBLE, "CDELT2  ", &cdelt2,        "", &status);
  fits_write_key(fptr, TDOUBLE, "CRPIX2  ", &crpix2,        "", &status);  	
  fits_write_key(fptr, TDOUBLE, "CRVAL2  ", &crval2,        "", &status);
  fits_write_key(fptr, TSTRING, "CUNIT2  ", &"DEG     ",    "", &status);
  fits_write_key(fptr, TSTRING, "CTYPE3  ", &"VELO-LSR",    "", &status);  
  fits_write_key(fptr, TDOUBLE, "CDELT3  ", &cdelt3,        "", &status);
  fits_write_key(fptr, TDOUBLE, "CRPIX3  ", &crpix3,        "", &status);
  fits_write_key(fptr, TDOUBLE, "CRVAL3  ", &crval3,        "", &status);
  fits_write_key(fptr, TSTRING, "CUNIT3  ", &"M/S     ",    "", &status);
  fits_write_key(fptr, TDOUBLE, "BSCALE  ", &bscale,        "", &status);
  fits_write_key(fptr, TDOUBLE, "BZERO   ", &bzero,         "", &status);
  if(img[im].unit==0) fits_write_key(fptr, TSTRING, "BUNIT", &"K       ", "", &status);
  if(img[im].unit==1) fits_write_key(fptr, TSTRING, "BUNIT", &"JY/PIXEL", "", &status);
  if(img[im].unit==2) fits_write_key(fptr, TSTRING, "BUNIT", &"WM2HZSR ", "", &status);
  if(img[im].unit==3) fits_write_key(fptr, TSTRING, "BUNIT", &"Lsun/PX ", "", &status);
  if(img[im].unit==3) fits_write_key(fptr, TSTRING, "BUNIT", &"        ", "", &status);

  /* Write FITS data */
  for(py=0;py<img[im].pxls;py++){
    for(ichan=0;ichan<img[im].nchan;ichan++){
      for(px=0;px<img[im].pxls;px++){
        if(img[im].unit==0) row[px]=(float) img[im].pixel[px+py*img[im].pxls].intense[ichan]*(CLIGHT/img[im].freq)*(CLIGHT/img[im].freq)/2./KBOLTZ*m[0].norm; 
        else if(img[im].unit==1) row[px]=(float) img[im].pixel[px+py*img[im].pxls].intense[ichan]*1e26*img[im].imgres*img[im].imgres*m[0].norm;
        else if(img[im].unit==2) row[px]=(float) img[im].pixel[px+py*img[im].pxls].intense[ichan]*m[0].norm;
        else if(img[im].unit==3) {
          ru3 = img[im].distance/1.975e13;
          row[px]=(float) img[im].pixel[px+py*img[im].pxls].intense[ichan]*4.*PI*ru3*ru3*img[im].freq*img[im].imgres*img[im].imgres*m[0].norm;
        }
        else if(img[im].unit==4) row[px]=(float) img[im].pixel[px+py*img[im].pxls].tau[ichan];
        else {
          if(!silent) bail_out("Image unit number invalid");
          exit(0);
        }
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
