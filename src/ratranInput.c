/*
 *  ratranInput.c
 *  LIME, The versatile 3D line modeling environment
 *
 *  Created by Christian Brinch on 03/06/09.
 *  Copyright 2006-2014, Christian Brinch,
 *  <brinch@nbi.dk>
 *  Niels Bohr institutet
 *  University of Copenhagen
 *	All rights reserved.
 *
 */

#include "lime.h"

double
ratranInput(char * modelfile, char * variable, double x, double y, double z){
  FILE *fp;
  char string[80];
  char *p,*q;
  double model[13];
  int i=0,dummy,ra=0,rb=0,za=0,zb=0,twode=0,count=0;

  if((fp=fopen(modelfile, "r"))==NULL){
    printf("Error opening RATRAN model file\n");
    exit(1);
  }

  while(strncmp(fgets(string,9,fp),"@",1) != 0){
    if(strncmp(string,"zmax=",5)==0) twode=1;
    if(strcmp(string,"columns=")==0) {
      fgets(string,80,fp);
      p=string;
      q=variable;
      while(*p!='\0'){
        if(*p==*q  && *(p+1)==*(q+1)) i=count;
        if(*p=='r' && *(p+1)=='a') ra=count;
        if(*p=='r' && *(p+1)=='b') rb=count;
        if(*p=='z' && *(p+1)=='a') za=count;
        if(*p=='z' && *(p+1)=='b') zb=count;
        count++;
        p++;
      }
    }
  }
  i=(i+1)/3-1;
  ra=(ra+1)/3-1;
  rb=(rb+1)/3-1;
  za=(za+1)/3-1;
  zb=(zb+1)/3-1;


  while(!feof(fp)) {
    if(twode) fscanf(fp,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                     &dummy, &model[0], &model[1], &model[2], &model[3], &model[4],
                     &model[5], &model[6], &model[7], &model[8], &model[9],
                     &model[10], &model[11], &model[12]);
    else fscanf(fp,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                &dummy, &model[0], &model[1], &model[2], &model[3], &model[4],
                &model[5], &model[6], &model[7], &model[8]);

    // if(twode && sqrt(x*x+y*y)>model[ra] && sqrt(x*x+y*y)<model[rb] && fabs(z)>model[za] && fabs(z)<model[zb]){
    if(twode && (x*x+y*y)>model[ra]*model[ra] && (x*x+y*y)<model[rb]*model[rb] && fabs(z)>model[za] && fabs(z)<model[zb]){
      fclose(fp);
      return model[i];
    }
    // if(!twode && sqrt(x*x+z*z+y*y)>model[ra] && sqrt(x*x+z*z+y*y)<model[rb]){
    if(!twode && (x*x+z*z+y*y)>model[ra]*model[ra] && (x*x+z*z+y*y)<model[rb]*model[rb]){
      fclose(fp);
      return model[i];
    }
  }

  fclose(fp);
  return 0.;
}

