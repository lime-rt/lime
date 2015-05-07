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
    v=deltav-veloproject(dx,vel);
    val=fabs(v)*binv;
    if(val <=  2500.){
      *vfac+= exp(-(val*val));
    }
  }
  *vfac=*vfac/steps;
  return;
}


void
line_plane_intersect(struct grid *g, double *ds, int posn, int *nposn, double *dx, double *x){
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
      if(newdist<*ds && newdist > 1e4){
        *ds=newdist;
        *nposn=g[posn].neigh[i]->id;
      }
    }
  }
}

void
raytrace(int im, inputPars *par, struct grid *g, molData *m, image *img){
  int *counta, *countb,nlinetot,aa;
  int ichan, posn,nposn,i,px,iline,tmptrans;
  double *tau, *subintens;
  double vfac=0.,x[3],dx[3];
  double deltav,ds,dist2,ndist2,size,xp,yp,zp,col,shift,minfreq,absDeltaFreq,jnu,alpha,snu,dtau,snu_pol[3];
  gsl_rng *ran = gsl_rng_alloc(gsl_rng_ranlxs2);	/* Random number generator */
#ifdef TEST
  gsl_rng_set(ran,178490);
#else
  gsl_rng_set(ran,time(0));
#endif
  
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
  tau = malloc(sizeof(double)*img[im].nchan);
  subintens = malloc(sizeof(double)*img[im].nchan);
  
  /* Main loop through pixel grid */
  for(px=0;px<(img[im].pxls*img[im].pxls);px++){
    for(ichan=0;ichan<img[im].nchan;ichan++){
      img[im].pixel[px].intense[ichan]=0.0;
      img[im].pixel[px].tau[ichan]=0.0;
    }
    for(aa=0;aa<par->antialias;aa++){
      for(ichan=0;ichan<img[im].nchan;ichan++){
        tau[ichan]=0.0;
        subintens[ichan]=0.0;
      }
      size=img[im].distance*img[im].imgres;

      xp=size*(gsl_rng_uniform(ran)+px%img[im].pxls)-size*img[im].pxls/2.;
      yp=size*(gsl_rng_uniform(ran)+px/img[im].pxls)-size*img[im].pxls/2.;

      if((xp*xp+yp*yp)/par->radiusSqu <= 1 ) {
        zp=sqrt(par->radiusSqu-(xp*xp+yp*yp));

        x[0]=xp*img[im].rotMat[0][0] + yp*img[im].rotMat[0][1] + zp*img[im].rotMat[0][2];
        x[1]=xp*img[im].rotMat[1][0] + yp*img[im].rotMat[1][1] + zp*img[im].rotMat[1][2];
        x[2]=xp*img[im].rotMat[2][0] + yp*img[im].rotMat[2][1] + zp*img[im].rotMat[2][2];

        dx[0]= -img[im].rotMat[0][2];
        dx[1]= -img[im].rotMat[1][2];
        dx[2]= -img[im].rotMat[2][2];
        
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
          line_plane_intersect(g,&ds,posn,&nposn,dx,x);
          if(par->polarization){
            for(ichan=0;ichan<img[im].nchan;ichan++){
              sourceFunc_pol(snu_pol,&dtau,ds,m,vfac,g,posn,0,0,img[im].theta);
              subintens[ichan]+=exp(-tau[ichan])*(1.-exp(-dtau))*snu_pol[ichan];
              tau[ichan]+=dtau;
            }
          } else {
            for(ichan=0;ichan<img[im].nchan;ichan++){
              jnu=.0;
              alpha=0.;
              snu=0.;
              dtau=0.;
              for(iline=0;iline<nlinetot;iline++){
                if(img[im].doline && m[counta[iline]].freq[countb[iline]] > img[im].freq-img[im].bandwidth/2. && m[counta[iline]].freq[countb[iline]] < img[im].freq+img[im].bandwidth/2.){
                  if(img[im].trans > -1){
                    shift=(m[counta[iline]].freq[countb[iline]]-m[counta[iline]].freq[img[im].trans])/m[counta[iline]].freq[img[im].trans]*CLIGHT;
                  } else {
                    shift=(m[counta[iline]].freq[countb[iline]]-img[im].freq)/img[im].freq*CLIGHT;
                  }
                  deltav=(ichan-(int)(img[im].nchan/2.))*img[im].velres-img[im].source_vel + shift;

                  if(!par->pregrid) velocityspline2(x,dx,ds,g[posn].mol[counta[iline]].binv,deltav,&vfac);
                  else vfac=gaussline(deltav-veloproject(dx,g[posn].vel),g[posn].mol[counta[iline]].binv);

                  sourceFunc_line(&jnu,&alpha,m,vfac,g,posn,counta[iline],countb[iline]);
                }
              }

              if(img[im].doline && img[im].trans > -1) sourceFunc_cont(&jnu,&alpha,g,posn,0,img[im].trans);
              else if(img[im].doline && img[im].trans == -1) sourceFunc_cont(&jnu,&alpha,g,posn,0,tmptrans);
              else sourceFunc_cont(&jnu,&alpha,g,posn,0,0);
              if(fabs(alpha)>0.){
                snu=(jnu/alpha)*m[0].norminv;
                dtau=alpha*ds;
              }
              subintens[ichan]+=exp(-tau[ichan])*(1.-exp(-dtau))*snu;
              tau[ichan]+=dtau;

            }
          }

          /* new coordinates */
          for(i=0;i<3;i++) x[i]+=ds*dx[i];
          col+=ds;
          posn=nposn;
        } while(col < 2*zp);

        /* add or subtract cmb */
        for(ichan=0;ichan<img[im].nchan;ichan++){
          subintens[ichan]+=(exp(-tau[ichan])-1.)*m[0].local_cmb[tmptrans];
        }

        for(ichan=0;ichan<img[im].nchan;ichan++){
          img[im].pixel[px].intense[ichan]+=subintens[ichan]/(double) par->antialias;
          img[im].pixel[px].tau[ichan]+=tau[ichan]/(double) par->antialias;
        }
      }
    }
    if(!silent) progressbar((double)(px)/(double)(img[im].pxls*img[im].pxls-1), 13);
  }

  img[im].trans=tmptrans;
  free(tau);
  free(subintens);
  free(counta);
  free(countb);
  gsl_rng_free(ran);
}



