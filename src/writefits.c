/*
 *  writefits.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#include "lime.h"

void
writeWCS(fitsfile *fptr, const int i, int axesOrder[4], float cdelt[4], double crpix[4], double crval[4], char ctype[4][9], char cunit[4][9]){
  char myStr[9];
  int status = 0;

  sprintf(myStr, "CTYPE%d  ", i+1);
  fits_write_key(fptr, TSTRING, myStr, &ctype[axesOrder[i]], "", &status);
  sprintf(myStr, "CDELT%d  ", i+1);
  fits_write_key(fptr, TFLOAT, myStr, &cdelt[axesOrder[i]], "", &status);
  sprintf(myStr, "CRPIX%d  ", i+1);
  fits_write_key(fptr, TDOUBLE, myStr, &crpix[axesOrder[i]], "", &status);
  sprintf(myStr, "CRVAL%d  ", i+1);
  fits_write_key(fptr, TDOUBLE, myStr, &crval[axesOrder[i]], "", &status);
  sprintf(myStr, "CUNIT%d  ", i+1);
  fits_write_key(fptr, TSTRING, myStr, &cunit[axesOrder[i]], "", &status);	
}

void 
write4Dfits(int im, int unit_index, configInfo *par, imageInfo *img){
  /*
Users have complained that downstream packages (produced by lazy coders >:8) will not deal with FITS cubes having less that 4 axes. Thus all LIME output images are now sent to the present function.
  */
  const int numAxes=4;
  double bscale,bzero,epoch,lonpole,equinox,restfreq;
  int axesOrder[] = {0,1,2,3};
  char ctype[numAxes][9],cunit[numAxes][9];
  double crpix[numAxes],crval[numAxes];
  float cdelt[numAxes];
  double ru3,scale=1.0;
  int velref,unitI,i;
  float *row;
  int px,py,ichan;
  fitsfile *fptr;
  int status = 0;
  int naxis=numAxes, bitpix=-32;
  long naxes[numAxes];
  long int fpixels[numAxes],lpixels[numAxes];
  char negfile[100]="! ",message[STR_LEN_0];
  unsigned long ppi;

  unitI = img[im].imgunits[unit_index];
  row = malloc(sizeof(*row)*img[im].pxls);

  naxes[axesOrder[0]] = img[im].pxls;
  naxes[axesOrder[1]] = img[im].pxls;

  if(img[im].doline)
    naxes[axesOrder[2]] = img[im].nchan;
  else
    naxes[axesOrder[2]] = 1; /* In this case nchan can =3, the number of active Stokes parameters, if the dust is polarized. */

  if(!img[im].doline && par->polarization)
    naxes[axesOrder[3]]=4;
  else
    naxes[axesOrder[3]]=1;

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

  sprintf(ctype[axesOrder[0]], "RA---SIN");
  cdelt[axesOrder[0]] = -1.8e2*(float)(img[im].imgres/M_PI);
  crpix[axesOrder[0]] = ((double)img[im].pxls)/2.0 + 0.5;
  crval[axesOrder[0]] = 0.0e0;
  sprintf(cunit[axesOrder[0]], "DEG    ");

  sprintf(ctype[axesOrder[1]], "DEC--SIN");
  cdelt[axesOrder[1]] = 1.8e2*(float)(img[im].imgres/M_PI);
  crpix[axesOrder[1]] = ((double)img[im].pxls)/2.0 + 0.5;
  crval[axesOrder[1]] = 0.0e0;
  sprintf(cunit[axesOrder[1]], "DEG    ");

  sprintf(ctype[axesOrder[2]], "VELO-LSR");
  if(img[im].doline)
    cdelt[axesOrder[2]] = (float)img[im].velres;
  else
    cdelt[axesOrder[2]] = 1.0;
  crpix[axesOrder[2]] = (double) (naxes[axesOrder[2]]-1)/2.+1;
  crval[axesOrder[2]] = 0.0e0;
  sprintf(cunit[axesOrder[2]], "M/S    ");

  sprintf(ctype[axesOrder[3]], "STOKES  ");
  cdelt[axesOrder[3]] = 1.0;
  crpix[axesOrder[3]] = (double) (naxes[axesOrder[3]]-1)/2.+1;
  crval[axesOrder[3]] = 1.0e0;
  sprintf(cunit[axesOrder[3]], "       ");

  bscale  =1.0e0;
  bzero   =0.0e0;

  fits_write_key(fptr, TSTRING, "OBJECT  ", &"LIMEMDL ",    "", &status);
  fits_write_key(fptr, TDOUBLE, "EPOCH   ", &epoch,         "", &status);
  fits_write_key(fptr, TDOUBLE, "LONPOLE ", &lonpole,       "", &status);
  fits_write_key(fptr, TDOUBLE, "EQUINOX ", &equinox,       "", &status);
  fits_write_key(fptr, TSTRING, "SPECSYS ", &"LSRK    ",    "", &status);
  fits_write_key(fptr, TDOUBLE, "RESTFREQ", &restfreq,      "", &status);
  fits_write_key(fptr, TINT,    "VELREF  ", &velref,        "", &status);

  for(i=0;i<numAxes;i++)
    writeWCS(fptr, i, axesOrder, cdelt, crpix, crval, ctype, cunit);

  fits_write_key(fptr, TDOUBLE, "BSCALE  ", &bscale,        "", &status);
  fits_write_key(fptr, TDOUBLE, "BZERO   ", &bzero,         "", &status);

  if(unitI==0) fits_write_key(fptr, TSTRING, "BUNIT", &"K       ", "", &status);
  if(unitI==1) fits_write_key(fptr, TSTRING, "BUNIT", &"JY/PIXEL", "", &status);
  if(unitI==2) fits_write_key(fptr, TSTRING, "BUNIT", &"WM2HZSR ", "", &status);
  if(unitI==3) fits_write_key(fptr, TSTRING, "BUNIT", &"Lsun/PX ", "", &status);
  if(unitI==4) fits_write_key(fptr, TSTRING, "BUNIT", &"        ", "", &status);

  if(     unitI==0)
    scale=0.5*(CLIGHT/img[im].freq)*(CLIGHT/img[im].freq)/KBOLTZ;
  else if(unitI==1)
    scale=1e26*img[im].imgres*img[im].imgres;
  else if(unitI==2)
    scale=1.0;
  else if(unitI==3) {
    ru3 = img[im].distance/1.975e13;
    scale=4.*M_PI*ru3*ru3*img[im].freq*img[im].imgres*img[im].imgres;
  }
  else if(unitI!=4) {
    if(!silent) bail_out("Image unit number invalid");
    exit(0);
  }

  /* Write FITS data */
  if(img[im].doline){
    for(ichan=0;ichan<img[im].nchan;ichan++){
      for(py=0;py<img[im].pxls;py++){
        for(px=0;px<img[im].pxls;px++){
          ppi = py*img[im].pxls + px;

          if(unitI>-1 && unitI<4)
            row[px]=(float) img[im].pixel[ppi].intense[ichan]*scale;
          else if(unitI==4)
            row[px]=(float) img[im].pixel[ppi].tau[ichan];
          else {
            if(!silent) bail_out("Image unit number invalid");
            exit(0);
          }
          if (fabs(row[px])<IMG_MIN_ALLOWED) row[px]=IMG_MIN_ALLOWED;
        }
        fpixels[axesOrder[0]] = 1;
        fpixels[axesOrder[1]] = py+1;
        fpixels[axesOrder[2]] = ichan+1;
        fpixels[axesOrder[3]] = 1;
        lpixels[axesOrder[0]] = img[im].pxls;
        lpixels[axesOrder[1]] = py+1;
        lpixels[axesOrder[2]] = ichan+1;
        lpixels[axesOrder[3]] = 1;
        fits_write_subset(fptr, TFLOAT, fpixels, lpixels, row, &status);
      }
    }
  }else{
    for(ichan=0;ichan<img[im].nchan;ichan++){
      for(py=0;py<img[im].pxls;py++){
        for(px=0;px<img[im].pxls;px++){
          ppi = py*img[im].pxls + px;

          if(unitI>-1 && unitI<4)
            row[px]=(float) img[im].pixel[ppi].intense[ichan]*scale;
          else if(unitI==4)
            row[px]=(float) img[im].pixel[ppi].tau[ichan];
          else {
            if(!silent) bail_out("Image unit number invalid");
            exit(0);
          }
          if (fabs(row[px])<IMG_MIN_ALLOWED) row[px]=IMG_MIN_ALLOWED;
        }
        fpixels[axesOrder[0]] = 1;
        fpixels[axesOrder[1]] = py+1;
        fpixels[axesOrder[3]] = ichan+1;
        fpixels[axesOrder[2]] = 1;
        lpixels[axesOrder[0]] = img[im].pxls;
        lpixels[axesOrder[1]] = py+1;
        lpixels[axesOrder[3]] = ichan+1;
        lpixels[axesOrder[2]] = 1;
        fits_write_subset(fptr, TFLOAT, fpixels, lpixels, row, &status);
      }
    }

    if(par->polarization){ /* ichan should have run from 0 to 2 in this case. Stokes I, Q and U but no V. Load zeros into the last pol channel: */
      if(img[im].nchan!=3){
        if(!silent){
          sprintf(message, "%d pol channels found but %d expected.", img[im].nchan, 3);
          bail_out(message);
        }
exit(1);
      }
      ichan = 3;
      for(px=0;px<img[im].pxls;px++)
        row[px] = IMG_MIN_ALLOWED;
      for(py=0;py<img[im].pxls;py++){
        fpixels[axesOrder[0]] = 1;
        fpixels[axesOrder[1]] = py+1;
        fpixels[axesOrder[3]] = ichan+1;
        fpixels[axesOrder[2]] = 1;
        lpixels[axesOrder[0]] = img[im].pxls;
        lpixels[axesOrder[1]] = py+1;
        lpixels[axesOrder[3]] = ichan+1;
        lpixels[axesOrder[2]] = 1;
        fits_write_subset(fptr, TFLOAT, fpixels, lpixels, row, &status);
      }
    }
  }

  fits_close_file(fptr, &status);

  free(row);
}

void writeFits(const int i, const int unit_index, configInfo *par, imageInfo *img){
  int unitI = img[i].imgunits[unit_index];

  if(unitI>5){
    if(!silent) bail_out("Image unit number invalid");
exit(1);
  }
  write4Dfits(i, unit_index, par, img);
}

char *removeFilenameExtension(char* inStr, char extensionChar, char pathSeparator) {
    char *outStr, *lastDotInFilename, *lastPathSeparatorInFilename;

    if (inStr == NULL)
        return NULL;

    outStr = malloc(strlen(inStr) + 1);
    if(!outStr){
        if(!silent) bail_out("Error allocating memory for filename extension removal");
        exit(0);
    }
    strcpy(outStr, inStr);
    /* Find last occurrences of extension character and path separator character */
    lastDotInFilename = strrchr(outStr, extensionChar);
    lastPathSeparatorInFilename = (pathSeparator == 0) ? NULL : strrchr(outStr, pathSeparator);

    /* Truncate filename at occurrence of last extension character assuming it comes after the last path separator character */
    if (lastDotInFilename != NULL) {
        if (lastPathSeparatorInFilename != NULL) {
            if (lastPathSeparatorInFilename < lastDotInFilename) {
                *lastDotInFilename = '\0';
            }
        } else {
            *lastDotInFilename = '\0';
        }
    }
    return outStr;
}

void insertUnitStrInFilename(char *img_filename_root, configInfo *par, imageInfo *img, const int im, const int unit_index){
  char *temp_filename, *temp_extensionless_filename, message[STR_LEN_0];
  static char* unit_names[] = {"Kelvin", "Jansky-per-px", "SI", "LSun-per-px", "Tau", "#Rays"};
  char *ext;

  /* Check if unit index falls outside range of possible unit names */
  if(unit_index < 0 || unit_index > sizeof(unit_names)/sizeof(*unit_names) - 1){
    sprintf(message, "Image unit index '%d' does not have a corresponding unit name", unit_index);
    if(!silent) bail_out(message);
    exit(0);
  }

  copyInparStr(img_filename_root, &(temp_filename));
  /* Extract filename extension */
  ext = strrchr(img_filename_root, '.');
  if (!ext) {
    /* Set to blank string if no filename extension was extracted */
    ext = "";
  } else {
    /* Remove extension from temporary filename */
      temp_extensionless_filename = removeFilenameExtension(temp_filename, '.', '/');
      strcpy(temp_filename, temp_extensionless_filename);
      free(temp_extensionless_filename);
  }
  /* Append unit name to temporary filename */
  strcat(temp_filename, "_");
  strcat(temp_filename, unit_names[img[im].imgunits[unit_index]]);
  strcat(temp_filename, ext);

  /* Update image filename from temporary filename */
  copyInparStr(temp_filename, &(img[im].filename));
  free(temp_filename);
}

void writeFitsAllUnits(const int i, configInfo *par, imageInfo *img){
  int j;
  char *img_filename_root;

  if(img[i].numunits == 1){
    writeFits(i,0,par,img);
  }else{
    copyInparStr(img[i].filename, &(img_filename_root));
    for(j=0;j<img[i].numunits;j++) {
      insertUnitStrInFilename(img_filename_root, par, img, i, j);
      writeFits(i,j,par,img);
    }
  }
}

