/*
 *  raytrace.c
 *  LIME, The versatile 3D line modeling environment
 *
 *  Created by Christian Brinch on 16/12/06.
 *  Copyright 2006-2014, Christian Brinch,
 *  <brinch@nbi.dk>
 *  Niels Bohr institutet
 *  University of Copenhagen
 *	All rights reserved.
 *
 */

#include "lime.h"


void
velocityspline2(double x[3], double dx[3], double ds, double binv, double deltav, double *vfac){
  int i,steps=10;
  double v,d,val,vel[3];

  *vfac=0.;
  for(i=0;i<steps;i++){
    d=i*ds/steps;
    velocity(x[0]+(dx[0]*d),x[1]+(dx[1]*d),x[2]+(dx[2]*d),vel);
    v=deltav+veloproject(dx,vel);
    val=fabs(v)*binv;
    if(val <=  2500.){
      *vfac+= exp(-(val*val));
    }
  }
  *vfac=*vfac/steps;
  return;
}


void
line_plane_intersect(struct grid *g, double *ds, int posn, int *nposn, double *dx, double *x, double cutoff){
  double newdist, numerator, denominator ;
  int i;

  for(i=0;i<g[posn].numNeigh;i++) {
    /* Find the shortest distance between (x,y,z) and any of the posn Voronoi faces */
    /* ds=(p0-l0) dot n / l dot n */

    numerator=((g[posn].x[0]+g[posn].dir[i].x[0]/2. - x[0]) * g[posn].dir[i].x[0]+
               (g[posn].x[1]+g[posn].dir[i].x[1]/2. - x[1]) * g[posn].dir[i].x[1]+
               (g[posn].x[2]+g[posn].dir[i].x[2]/2. - x[2]) * g[posn].dir[i].x[2]);

    denominator=(dx[0]*g[posn].dir[i].x[0]+dx[1]*g[posn].dir[i].x[1]+dx[2]*g[posn].dir[i].x[2]);
    
    if(fabs(denominator) > 0){
      newdist=numerator/denominator;
      if(newdist<*ds && newdist > cutoff){
        *ds=newdist;
        *nposn=g[posn].neigh[i]->id;
      }
    }
  }
  if(*nposn==-1) *nposn=posn;
}


void
traceray(rayData ray, int tmptrans, int im, inputPars *par, struct grid *g, molData *m, image *img, int nlinetot, int *counta, int *countb, double cutoff){
  int ichan,posn,nposn,i,iline;
  double vfac=0.,x[3],dx[3],vThisChan;
  double deltav,ds,dist2,ndist2,xp,yp,zp,col,shift,jnu,alpha,remnantSnu,dtau,expDTau,snu_pol[3];

  for(ichan=0;ichan<img[im].nchan;ichan++){
    ray.tau[ichan]=0.0;
    ray.intensity[ichan]=0.0;
  }

  xp=ray.x;
  yp=ray.y;

  if((xp*xp+yp*yp)/par->radiusSqu <= 1 ) {
    zp=sqrt(par->radiusSqu-(xp*xp+yp*yp));

    for(i=0;i<3;i++){
      x[i]=xp*img[im].rotMat[i][0] + yp*img[im].rotMat[i][1] + zp*img[im].rotMat[i][2];
      dx[i]= -img[im].rotMat[i][2];
    }

    i=0;
    dist2=(x[0]-g[i].x[0])*(x[0]-g[i].x[0]) + (x[1]-g[i].x[1])*(x[1]-g[i].x[1]) + (x[2]-g[i].x[2])*(x[2]-g[i].x[2]);
    posn=i;
    for(i=1;i<par->ncell;i++){
      ndist2=(x[0]-g[i].x[0])*(x[0]-g[i].x[0]) + (x[1]-g[i].x[1])*(x[1]-g[i].x[1]) + (x[2]-g[i].x[2])*(x[2]-g[i].x[2]);
      if(ndist2<dist2){
        posn=i;
        dist2=ndist2;
      }
    }

    col=0;
    do{
      ds=2.*zp-col;
      nposn=-1;
      line_plane_intersect(g,&ds,posn,&nposn,dx,x,cutoff);
      if(par->polarization){
        for(ichan=0;ichan<img[im].nchan;ichan++){
          sourceFunc_pol(snu_pol,&dtau,ds,m,vfac,g,posn,0,0,img[im].theta);
          ray.intensity[ichan]+=exp(-ray.tau[ichan])*(1.-exp(-dtau))*snu_pol[ichan];
          ray.tau[ichan]+=dtau;
        }
      } else {
        for(ichan=0;ichan<img[im].nchan;ichan++){
          jnu=.0;
          alpha=0.;

          for(iline=0;iline<nlinetot;iline++){
            if(img[im].doline && m[counta[iline]].freq[countb[iline]] > img[im].freq-img[im].bandwidth/2.
            && m[counta[iline]].freq[countb[iline]] < img[im].freq+img[im].bandwidth/2.){
              if(img[im].trans > -1){
                shift=(m[counta[iline]].freq[countb[iline]]-m[counta[iline]].freq[img[im].trans])/m[counta[iline]].freq[img[im].trans]*CLIGHT;
              } else {
                shift=(m[counta[iline]].freq[countb[iline]]-img[im].freq)/img[im].freq*CLIGHT;
              }
              vThisChan=(ichan-(img[im].nchan-1)/2.)*img[im].velres; // Consistent with the WCS definition in writefits().
              deltav = vThisChan - img[im].source_vel + shift;

              if(!par->pregrid) velocityspline2(x,dx,ds,g[posn].mol[counta[iline]].binv,deltav,&vfac);
              else vfac=gaussline(deltav+veloproject(dx,g[posn].vel),g[posn].mol[counta[iline]].binv);

              sourceFunc_line(&jnu,&alpha,m,vfac,g,posn,counta[iline],countb[iline]);
            }
          }

          if(img[im].doline && img[im].trans > -1) sourceFunc_cont(&jnu,&alpha,g,posn,0,img[im].trans);
          else if(img[im].doline && img[im].trans == -1) sourceFunc_cont(&jnu,&alpha,g,posn,0,tmptrans);
          else sourceFunc_cont(&jnu,&alpha,g,posn,0,0);
          dtau=alpha*ds;
//???          if(dtau < -30) dtau = -30; // as in photon()?
          calcSourceFn(dtau, par, &remnantSnu, &expDTau);
          remnantSnu *= jnu*m[0].norminv*ds;
          ray.intensity[ichan]+=exp(-ray.tau[ichan])*remnantSnu;
          ray.tau[ichan]+=dtau;
        }
      }

      /* new coordinates */
      for(i=0;i<3;i++) x[i]+=ds*dx[i];
      col+=ds;
      posn=nposn;
    } while(col < 2*zp);

    /* add or subtract cmb */
    for(ichan=0;ichan<img[im].nchan;ichan++){
      ray.intensity[ichan]+=(exp(-ray.tau[ichan])-1.)*m[0].local_cmb[tmptrans];
    }
  }
}


void
raytrace(int im, inputPars *par, struct grid *g, molData *m, image *img){
  int *counta, *countb,nlinetot,aa;
  int ichan,px,iline,tmptrans;
  double size,minfreq,absDeltaFreq;
  double cutoff;
  rayData ray;

  gsl_rng *ran = gsl_rng_alloc(gsl_rng_ranlxs2);	/* Random number generator */
#ifdef TEST
  gsl_rng_set(ran,178490);
#else
  gsl_rng_set(ran,time(0));
#endif

  size=img[im].distance*img[im].imgres;

  /* Determine whether there are blended lines or not */
  lineCount(par->nSpecies, m, &counta, &countb, &nlinetot);
  if(img[im].doline==0) nlinetot=1;

  /* Fix the image parameters */
  if(img[im].freq < 0) img[im].freq=m[0].freq[img[im].trans];
  if(img[im].nchan == 0 && img[im].bandwidth>0){
    img[im].nchan=(int) (img[im].bandwidth/(img[im].velres/CLIGHT*img[im].freq));
  } else if (img[im].velres<0 && img[im].bandwidth>0){
    img[im].velres = img[im].bandwidth*CLIGHT/img[im].freq/img[im].nchan;
  } else img[im].bandwidth = img[im].nchan*img[im].velres/CLIGHT * img[im].freq;

  if(img[im].trans<0){
    iline=0;
    minfreq=fabs(img[im].freq-m[0].freq[iline]);
    tmptrans=iline;
    for(iline=1;iline<m[0].nline;iline++){
      absDeltaFreq=fabs(img[im].freq-m[0].freq[iline]);
      if(absDeltaFreq<minfreq){
        minfreq=absDeltaFreq;
        tmptrans=iline;
      }
    }
  } else tmptrans=img[im].trans;

  ray.intensity=malloc(sizeof(double) * img[im].nchan);
  ray.tau=malloc(sizeof(double) * img[im].nchan);

  cutoff = par->minScale*1.0e-7;

  /* Main loop through pixel grid */
  for(px=0;px<(img[im].pxls*img[im].pxls);px++){
    for(ichan=0;ichan<img[im].nchan;ichan++){
      img[im].pixel[px].intense[ichan]=0.0;
      img[im].pixel[px].tau[ichan]=0.0;
    }
    for(aa=0;aa<par->antialias;aa++){
      ray.x = size*(gsl_rng_uniform(ran)+px%img[im].pxls)-size*img[im].pxls/2.;
      ray.y = size*(gsl_rng_uniform(ran)+px/img[im].pxls)-size*img[im].pxls/2.;

      traceray(ray, tmptrans, im, par, g, m, img, nlinetot, counta, countb, cutoff);

      for(ichan=0;ichan<img[im].nchan;ichan++){
        img[im].pixel[px].intense[ichan] += ray.intensity[ichan]/(double) par->antialias;
        img[im].pixel[px].tau[ichan] += ray.tau[ichan]/(double) par->antialias;
      }
    }

    if(!silent) progressbar((double)(px)/(double)(img[im].pxls*img[im].pxls-1), 13);
  }

  img[im].trans=tmptrans;

  free(ray.tau);
  free(ray.intensity);
  free(counta);
  free(countb);
  gsl_rng_free(ran);
}


void
raytrace_1_4(int im, inputPars *par, struct grid *g, molData *m, image *img){
  /*
This is an alternative raytracing algorithm which was implemented by
C Brinch in version 1.4 (the original parallelized version) of LIME.
I've adapted it slightly so it makes use of the function traceray(),
which was modified from the function tracerays() in v1.4. This algorithm
is not currently used, but may be useful as an option; that's why I have
kept it.
  */

  int *counta, *countb,nlinetot;
  int ichan,i,px,iline,tmptrans,count;
  double size,xp,yp,minfreq,absDeltaFreq;
  double cutoff;

  gsl_rng *ran = gsl_rng_alloc(gsl_rng_ranlxs2);	/* Random number generator */
#ifdef TEST
  gsl_rng_set(ran,178490);
#else
  gsl_rng_set(ran,time(0));
#endif
  rayData *rays;
  
  int sg,n;
  double cx,cy;

  double x1,x2,x3,y1,y2,y3,z1,z2,z3,xt[3],yt[3],di,p,d1,d2,d3,temp1;
  int zt[3];
  int c;
  
  char flags[255];
  boolT ismalloc = False;
  facetT *facet, *neighbor, **neighborp;;
  vertexT *vertex,**vertexp;
  coordT *pt_array;

  int id;
  coordT point[3];
  boolT isoutside;
  realT bestdist;

  size=img[im].distance*img[im].imgres;

  /* Determine whether there are blended lines or not */
  lineCount(par->nSpecies, m, &counta, &countb, &nlinetot);
  if(img[im].doline==0) nlinetot=1;

  /* Fix the image parameters */
  if(img[im].freq < 0) img[im].freq=m[0].freq[img[im].trans];
  if(img[im].nchan == 0 && img[im].bandwidth>0){
    img[im].nchan=(int) (img[im].bandwidth/(img[im].velres/CLIGHT*img[im].freq));
  } else if (img[im].velres<0 && img[im].bandwidth>0){
    img[im].velres = img[im].bandwidth*CLIGHT/img[im].freq/img[im].nchan;
  } else img[im].bandwidth = img[im].nchan*img[im].velres/CLIGHT * img[im].freq;

  if(img[im].trans<0){
    iline=0;
    minfreq=fabs(img[im].freq-m[0].freq[iline]);
    tmptrans=iline;
    for(iline=1;iline<m[0].nline;iline++){
      absDeltaFreq=fabs(img[im].freq-m[0].freq[iline]);
      if(absDeltaFreq<minfreq){
        minfreq=absDeltaFreq;
        tmptrans=iline;
      }
    }
  } else tmptrans=img[im].trans;

  /* Allocate dynamical arrays */
  rays = malloc(sizeof(rayData) * (par->pIntensity));
  
  for(i=0;i<par->pIntensity;i++){
    rays[i].x=g[i].x[0];
    rays[i].y=g[i].x[1];
    rays[i].tau=malloc(sizeof(double) * img[im].nchan);
    rays[i].intensity=malloc(sizeof(double) * img[im].nchan);
    for(ichan=0;ichan<img[im].nchan;ichan++) {
      rays[i].tau[ichan]=0.0;
      rays[i].intensity[ichan]=0.0;
    }
  }
  

  /* Smooth out the distribution of rays */
  for(sg=0;sg<20;sg++){
    pt_array=malloc(2*sizeof(coordT)*par->pIntensity);
          
    for(i=0;i<par->pIntensity;i++) {
      pt_array[i*2+0]=rays[i].x;
      pt_array[i*2+1]=rays[i].y;
    }
    
    sprintf(flags,"qhull v s Qbb T0");
    if (!qh_new_qhull(2, par->pIntensity, pt_array, ismalloc, flags, NULL, NULL)) {

      qh_setvoronoi_all();
      
      FORALLvertices {
        i=qh_pointid(vertex->point);
        
        cx=0.;
        cy=0.;
        n=0;
        FOREACHneighbor_(vertex) {
          if (!neighbor->upperdelaunay) n++;
        }
        if(n>0){
        
          
        } else {
          if(!silent) bail_out("Qhull error");
          exit(0);
        }
        
        FOREACHneighbor_(vertex) {
          if (!neighbor->upperdelaunay) {
            cx+=neighbor->center[0];
            cy+=neighbor->center[1];
          }
        }

        rays[i].x = rays[i].x - (rays[i].x-cx/ (double) n)*0.1;
        rays[i].y = rays[i].y - (rays[i].y-cy/ (double) n)*0.1;
      }
    } else {
      printf("qhull error\n");
    }
    
    qh_freeqhull(!qh_ALL);
    free(pt_array);  
  }

  cutoff = par->minScale*1.0e-7;

  /* Main loop through rays */
  count=0;
  for(px=0;px<par->pIntensity;px++){
    traceray(rays[px], tmptrans, im, par, g, m, img, nlinetot, counta, countb, cutoff);
    ++count;
    if(!silent) progressbar((double)(count)/(double)(par->pIntensity-1), 13);
  }

  /* Remap rays onto pixel grid */
  pt_array=malloc(2*sizeof(coordT)*par->pIntensity);

  for(i=0;i<par->pIntensity;i++) {
    pt_array[i*2+0]=rays[i].x;
    pt_array[i*2+1]=rays[i].y;
  }

/* This allocation belongs to "Shepard's method" below
  d=malloc(sizeof(double)*par->pIntensity);
*/
  size=img[im].distance*img[im].imgres;

  sprintf(flags,"qhull d Qbb");
  if (!qh_new_qhull(2, par->pIntensity, pt_array, ismalloc, flags, NULL, NULL)) {
    for(px=0;px<img[im].pxls*img[im].pxls;px++){
      for(ichan=0;ichan<img[im].nchan;ichan++){
        img[im].pixel[px].intense[ichan]=0.0;
        img[im].pixel[px].tau[ichan]=0.0;
      }
      xp=size*(0.5+px%img[im].pxls)-size*img[im].pxls/2.;
      yp=size*(0.5+px/img[im].pxls)-size*img[im].pxls/2.;

/*
This part works great! This is "Shepard's method" with a weight of 8. Slow, unfortunately. Could it be parallelized?

      for(ichan=0;ichan<img[im].nchan;ichan++){
        img[im].pixel[px].intense[ichan] = 0.;
        di=0;
        for(i=0;i<par->pIntensity;i++){
          // d[i]=1./pow(sqrt(pow(xp-rays[i].x,2)+ pow(yp-rays[i].y,2)),8.);
          temp1 = (xp-rays[i].x)*(xp-rays[i].x)+ (yp-rays[i].y)*(yp-rays[i].y)
          d[i]=1./(temp1*temp1*temp1*temp1);
          img[im].pixel[px].intense[ichan] += rays[i].intensity[ichan]*d[i];
**** how to handle img[im].pixel[px].tau[ichan]?
          di+=d[i];
        }
        img[im].pixel[px].intense[ichan] /= di;
      }
*/
      
      
      point[0]=xp;
      point[1]=yp;
      point[2]=0.;
     
      qh_setdelaunay (3, 1, point);
      facet= qh_findbestfacet (point, qh_ALL, &bestdist, &isoutside);
      
      c=0;
      FOREACHvertex_( facet->vertices ) {
        id=qh_pointid(vertex->point);
        xt[c]=rays[id].x; yt[c]=rays[id].y; zt[c]=id;
        c++;
      }
      
      
      x1=xt[0];x2=xt[1];x3=xt[2];
      y1=yt[0];y2=yt[1];y3=yt[2];

      for(ichan=0;ichan<img[im].nchan;ichan++){
        z1=rays[zt[2]].intensity[ichan];z2=rays[zt[1]].intensity[ichan];z3=rays[zt[2]].intensity[ichan];
      
        
        p=1.;
        // d1=1./pow(sqrt(pow(xp-x1,2)+ pow(yp-y1,2)),p);
        // d2=1./pow(sqrt(pow(xp-x2,2)+ pow(yp-y2,2)),p);
        // d3=1./pow(sqrt(pow(xp-x3,2)+ pow(yp-y3,2)),p);

        d1=1./sqrt((xp-x1)*(xp-x1) + (yp-y1)*(yp-y1));
        d2=1./sqrt((xp-x2)*(xp-x2) + (yp-y2)*(yp-y2));
        d3=1./sqrt((xp-x3)*(xp-x3) + (yp-y3)*(yp-y3));

        di=d1+d2+d3;
        img[im].pixel[px].intense[ichan] = 1./di * (z1*d1 + z2*d2 + z3*d3);
//**** how to handle img[im].pixel[px].tau[ichan]?
      }

      
    }
  } else {
	if(!silent) bail_out("Qhull failed to triangulate");
	exit(1);
  }

  img[im].trans=tmptrans;

  free(pt_array);
  for(i=0;i<par->pIntensity;i++){
    free(rays[i].tau);
    free(rays[i].intensity);
  }
  free(rays);
  free(counta);
  free(countb);
}


