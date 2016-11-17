/*
 *  report.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2016 The LIME development team
 *
 */

#include "lime.h"
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_fit.h>


void
report(int i, configInfo *par, struct grid *g){
  FILE *fp;
  int j,k,p,q,min_q,max_q;
  const int bins = 100;
  double x,min_l,max_l,min_p,max_p,min_p2,max_p2,r2;
  size_t n=bins;
  double hx[bins],hy[bins],hw[bins];
  double c0, c1, cov00, cov01, cov11, chisq1,b0, b1, chisq2;

  gsl_histogram * h = gsl_histogram_alloc (n);
  gsl_histogram * f = gsl_histogram_alloc (n);
  
  if(i==1){
    fp=fopen("LimeReport","w");
    fprintf(fp,"*** LIME report\n***\n");
    fprintf(fp,"*** Grid statistics\n\n");

    j=0;
    p=g[j].numNeigh;
    q=g[j].nphot;
    r2=g[j].x[0]*g[j].x[0]+g[j].x[1]*g[j].x[1]+g[j].x[2]*g[j].x[2];
    min_p2=r2;
    max_p2=r2;
    k=0;
    min_l=g[j].ds[k];
    max_l=g[j].ds[k];
    min_q=q;
    max_q=q;
    for(k=1;k<g[j].numNeigh;k++){
      if(g[j].ds[k]<min_l) min_l=g[j].ds[k];
      if(g[j].ds[k]>max_l) max_l=g[j].ds[k];
    }
    for(j=1;j<par->ncell;j++) {
      p+=g[j].numNeigh;
      if (!g[j].sink) {
        q+=g[j].nphot;
        if (g[j].nphot>max_q) max_q=g[j].nphot;
        if (g[j].nphot<min_q) min_q=g[j].nphot;
      }
      r2=g[j].x[0]*g[j].x[0]+g[j].x[1]*g[j].x[1]+g[j].x[2]*g[j].x[2];
      if(r2<min_p2) min_p2=r2;
      if(r2>max_p2) max_p2=r2;
      for(k=0;k<g[j].numNeigh;k++){
        if(g[j].ds[k]<min_l) min_l=g[j].ds[k];
        if(g[j].ds[k]>max_l) max_l=g[j].ds[k];
      }
    }
    min_p=sqrt(min_p2);
    max_p=sqrt(max_p2);
    x=(double) p/par->ncell;
    p=p/2;

    gsl_histogram_set_ranges_uniform (h, par->minScale, par->radius);
    gsl_histogram_set_ranges_uniform (f, min_l, max_l);

    for(j=0;j<par->ncell-par->sinkPoints;j++) {
      gsl_histogram_increment (h, sqrt(g[j].x[0]*g[j].x[0]+g[j].x[1]*g[j].x[1]+g[j].x[2]*g[j].x[2])); 
      for(k=0;k<g[j].numNeigh;k++){
        gsl_histogram_increment (f, g[j].ds[k]);
      }
    }

    for(j=0;j<bins;j++){
      hx[j]=log10(j*(par->radius-par->minScale)/100+par->minScale);
      if(gsl_histogram_get(h,j) >0) hy[j]=log10(gsl_histogram_get(h,j));
      else hy[j]=0;
      hw[j]=j;
    }
    gsl_fit_wlinear(hx,1,hw,1,hy,1,bins,&c0, &c1, &cov00, &cov01, &cov11, &chisq1);

    for(j=0;j<bins;j++){
      hx[j]=log10(j*(max_l-min_l)/100+min_l);
      if(gsl_histogram_get(f,j) >0) hy[j]=log10(gsl_histogram_get(f,j));
      else hy[j]=0;
      hw[j]=j;
    }
    gsl_fit_wlinear(hx,1,hw,1,hy,1,bins,&b0, &b1, &cov00, &cov01, &cov11, &chisq2);

    fprintf(fp,"    Grid points    Surface points    Delaunay lines    Average neighbors\n");
    fprintf(fp,"      %6d           %5d           %7d               %2.2f\n\n", par->ncell-par->sinkPoints, par->sinkPoints, p,x);
    fprintf(fp,"    Min. and max. line length         Min. and max. point position\n");
    fprintf(fp,"      %1.2e     %1.2e              %1.2e      %1.2e\n\n",min_l,max_l,min_p,max_p);
    fprintf(fp,"    Line distribution slope (Chi^2)   Point distribution slope (Chi^2)\n");
    fprintf(fp,"	  %2.1f (%2.1f)                            %2.1f (%2.1f)\n", b1,chisq2/bins, c1, chisq1/bins);


    fprintf(fp,"\n***\n*** Photon statistic\n\n");
    fprintf(fp,"    Total photons   Average photons    Minimum photons  Maximum photons\n");
    fprintf(fp,"      %6d           %5d                %7d              %7d\n", q,q/(par->ncell-par->sinkPoints),min_q,max_q);
    fclose(fp);
  }
  /*    gsl_histogram_fprintf (stdout, h, "%g", "%g");	 */

  gsl_histogram_free (h);
}
