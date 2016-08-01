/*
 *  raytrace.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2016 The LIME development team
 *
 */

#include "lime.h"


void
velocityspline2(double x[3], double dx[3], double ds, double binv, double deltav, double *vfac){
  /*
The bulk velocity of the model material can vary significantly with position, thus so can the value of the line-shape function at a given frequency and direction. The present function calculates 'vfac', an approximate average of the line-shape function along a path of length ds in the direction of the line of sight.
  */
  int i,steps=10;
  double v,d,val,vel[3];

  *vfac=0.;
  for(i=0;i<steps;i++){
    d=i*ds/steps;
    velocity(x[0]+(dx[0]*d),x[1]+(dx[1]*d),x[2]+(dx[2]*d),vel);
    v=deltav-veloproject(dx,vel); /* veloproject returns the component of the local bulk velocity in the direction of dx, whereas deltav is the recession velocity of the channel we are interested in (corrected for bulk source velocity and line displacement from the nominal frequency). Remember also that, since dx points away from the observer, positive values of the projected velocity also represent recessions. Line centre occurs when v==0, i.e. when deltav==veloproject(dx,vel). That is the reason for the subtraction here. */
    val=fabs(v)*binv;
    if(val <=  2500.){
#ifdef FASTEXP
      *vfac+= FastExp(val*val);
#else
      *vfac+=   exp(-(val*val));
#endif
    }
  }
  *vfac=*vfac/steps;
  return;
}


void
line_plane_intersect(struct grid *g, double *ds, int posn, int *nposn, double *dx, double *x, double cutoff){
  /*
This function returns ds as the (always positive-valued) distance between the present value of x and the next Voronoi face in the direction of vector dx, and nposn as the id of the grid cell that abuts that face. 
  */
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
traceray(rayData ray, int cmbMolI, int cmbLineI, int im, configInfo *par, struct grid *gp, molData *md, image *img, double cutoff){
  /*
For a given image pixel position, this function evaluates the intensity of the total light emitted/absorbed along that line of sight through the (possibly rotated) model. The calculation is performed for several frequencies, one per channel of the output image.

Note that the algorithm employed here is similar to that employed in the function photon() which calculates the average radiant flux impinging on a grid cell: namely the notional photon is started at the side of the model near the observer and 'propagated' in the receding direction until it 'reaches' the far side. This is rather non-physical in conception but it makes the calculation easier.
  */
  const int stokesIi=0;
  int ichan,posn,nposn,i,molI,lineI;
  double vfac=0.,x[3],dx[3],vThisChan;
  double deltav,ds,dist2,ndist2,xp,yp,zp,col,lineRedShift,jnu,alpha,remnantSnu,dtau,expDTau,snu_pol[3];

  for(ichan=0;ichan<img[im].nchan;ichan++){
    ray.tau[ichan]=0.0;
    ray.intensity[ichan]=0.0;
  }

  xp=ray.x;
  yp=ray.y;

  if((xp*xp+yp*yp)/par->radiusSqu <= 1 ) {
    zp=-sqrt(par->radiusSqu-(xp*xp+yp*yp)); /* There are two points of intersection between the line of sight and the spherical model surface; this is the Z coordinate (in the unrotated frame) of the one nearer to the observer. */

    /* Rotate the line of sight as desired. */
    for(i=0;i<3;i++){
      x[i]=xp*img[im].rotMat[i][0] + yp*img[im].rotMat[i][1] + zp*img[im].rotMat[i][2];
      dx[i]= img[im].rotMat[i][2]; /* This points away from the observer. */
    }

    /* Find the grid point nearest to the starting x. */
    i=0;
    dist2=(x[0]-gp[i].x[0])*(x[0]-gp[i].x[0]) + (x[1]-gp[i].x[1])*(x[1]-gp[i].x[1]) + (x[2]-gp[i].x[2])*(x[2]-gp[i].x[2]);
    posn=i;
    for(i=1;i<par->ncell;i++){
      ndist2=(x[0]-gp[i].x[0])*(x[0]-gp[i].x[0]) + (x[1]-gp[i].x[1])*(x[1]-gp[i].x[1]) + (x[2]-gp[i].x[2])*(x[2]-gp[i].x[2]);
      if(ndist2<dist2){
        posn=i;
        dist2=ndist2;
      }
    }

    col=0;
    do{
      ds=-2.*zp-col; /* This default value is chosen to be as large as possible given the spherical model boundary. */
      nposn=-1;
      line_plane_intersect(gp,&ds,posn,&nposn,dx,x,cutoff); /* Returns a new ds equal to the distance to the next Voronoi face, and nposn, the ID of the grid cell that abuts that face. */ 
      if(par->polarization){
        sourceFunc_pol(snu_pol,&alpha,gp,posn,0,0,img[im].rotMat);
        dtau=alpha*ds;
        calcSourceFn(dtau, par, &remnantSnu, &expDTau);
        remnantSnu *= md[0].norminv*ds;

        for(ichan=0;ichan<img[im].nchan;ichan++){ /* Loop over I, Q and U */
#ifdef FASTEXP
          ray.intensity[ichan]+=FastExp(ray.tau[ichan])*remnantSnu*snu_pol[ichan];
#else
          ray.intensity[ichan]+=   exp(-ray.tau[ichan])*remnantSnu*snu_pol[ichan];
#endif
          ray.tau[ichan]+=dtau; //**** But this will be the same for I, Q or U.
        }
      } else {
        for(ichan=0;ichan<img[im].nchan;ichan++){
          jnu=.0;
          alpha=0.;

          if(img[im].doline){
            for(molI=0;molI<par->nSpecies;molI++){
              for(lineI=0;lineI<md[molI].nline;lineI++){
                if(md[molI].freq[lineI] > img[im].freq-img[im].bandwidth*0.5
                && md[molI].freq[lineI] < img[im].freq+img[im].bandwidth*0.5){
                  /* Calculate the red shift of the transition wrt to the frequency specified for the image. */
                  if(img[im].trans > -1){
                    lineRedShift=(md[molI].freq[img[im].trans]-md[molI].freq[lineI])/md[molI].freq[img[im].trans]*CLIGHT;
                  } else {
                    lineRedShift=(img[im].freq-md[molI].freq[lineI])/img[im].freq*CLIGHT;
                  }

                  vThisChan=(ichan-(img[im].nchan-1)/2.)*img[im].velres; /* Consistent with the WCS definition in writefits(). */
                  deltav = vThisChan - img[im].source_vel - lineRedShift;
                  /* Line centre occurs when deltav = the recession velocity of the radiating material. Explanation of the signs of the 2nd and 3rd terms on the RHS: (i) A bulk source velocity (which is defined as >0 for the receding direction) should be added to the material velocity field; this is equivalent to subtracting it from deltav, as here. (ii) A positive value of lineRedShift means the line is red-shifted wrt to the frequency specified for the image. The effect is the same as if the line and image frequencies were the same, but the bulk recession velocity were higher. lineRedShift should thus be added to the recession velocity, which is equivalent to subtracting it from deltav, as here. */

                  /* Calculate an approximate average line-shape function at deltav within the Voronoi cell. */
                  if(!par->pregrid) velocityspline2(x,dx,ds,gp[posn].mol[molI].binv,deltav,&vfac);
                  else vfac=gaussline(deltav-veloproject(dx,gp[posn].vel),gp[posn].mol[molI].binv);

                  /* Increment jnu and alpha for this Voronoi cell by the amounts appropriate to the spectral line. */
                  sourceFunc_line(&jnu,&alpha,md,vfac,gp,posn,molI,lineI);
                }
              }
            }
          }

          sourceFunc_cont(&jnu,&alpha,gp,posn,cmbMolI,cmbLineI);
          dtau=alpha*ds;
          calcSourceFn(dtau, par, &remnantSnu, &expDTau);
          remnantSnu *= jnu*md[0].norminv*ds;
#ifdef FASTEXP
          ray.intensity[ichan]+=FastExp(ray.tau[ichan])*remnantSnu;
#else
          ray.intensity[ichan]+=   exp(-ray.tau[ichan])*remnantSnu;
#endif
          ray.tau[ichan]+=dtau;
        }
      }

      /* Move the working point to the edge of the next Voronoi cell. */
      for(i=0;i<3;i++) x[i]+=ds*dx[i];
      col+=ds;
      posn=nposn;
    } while(col < 2.0*fabs(zp));

    /* Add or subtract cmb. */
    if(par->polarization){ /* just add it to Stokes I */
#ifdef FASTEXP
      ray.intensity[stokesIi]+=FastExp(ray.tau[stokesIi])*md[cmbMolI].local_cmb[cmbLineI];
#else
      ray.intensity[stokesIi]+=exp(   -ray.tau[stokesIi])*md[cmbMolI].local_cmb[cmbLineI];
#endif

    }else{
#ifdef FASTEXP
      for(ichan=0;ichan<img[im].nchan;ichan++){
        ray.intensity[ichan]+=FastExp(ray.tau[ichan])*md[cmbMolI].local_cmb[cmbLineI];
      }
#else
      for(ichan=0;ichan<img[im].nchan;ichan++){
        ray.intensity[ichan]+=exp(-ray.tau[ichan])*md[cmbMolI].local_cmb[cmbLineI];
      }
#endif
    }
  }
}


void
raytrace(int im, configInfo *par, struct grid *gp, molData *md, image *img){
  int aa,ichan,px,cmbMolI,cmbLineI,i,threadI,nRaysDone,molI,lineI;
  double size,minfreq,absDeltaFreq,totalNumPixelsMinus1=(double)(img[im].pxls*img[im].pxls-1);
  double cutoff;
  const gsl_rng_type *ranNumGenType = gsl_rng_ranlxs2;
  gsl_error_handler_t *defaultErrorHandler=NULL;

  gsl_rng *ran = gsl_rng_alloc(ranNumGenType);	/* Random number generator */
#ifdef TEST
  gsl_rng_set(ran,178490);
#else
  gsl_rng_set(ran,time(0));
#endif

  gsl_rng **threadRans;
  threadRans = malloc(sizeof(gsl_rng *)*par->nThreads);

  for (i=0;i<par->nThreads;i++){
    threadRans[i] = gsl_rng_alloc(ranNumGenType);
    gsl_rng_set(threadRans[i],(int)(gsl_rng_uniform(ran)*1e6));
  }

  size=img[im].distance*img[im].imgres;

  if(img[im].doline){
    /* The user may have set img.trans/img.molI but not img.freq. If so, we calculate freq now.
    */
    if(img[im].trans>-1)
      img[im].freq = md[img[im].molI].freq[img[im].trans];

    /* Fill in the missing one of the triplet nchan/velres/bandwidth.
    */
    if(img[im].bandwidth > 0 && img[im].velres > 0){
      img[im].nchan = (int)(img[im].bandwidth/(img[im].velres/CLIGHT*img[im].freq));

    }else if(img[im].bandwidth > 0 && img[im].nchan > 0){
      img[im].velres = img[im].bandwidth*CLIGHT/img[im].freq/img[im].nchan;

    }else{ /*(img[im].velres > 0 && img[im].nchan > 0 */
      img[im].bandwidth = img[im].nchan*img[im].velres/CLIGHT*img[im].freq;
    }
  } /* If not doline, we already have img.freq and nchan by now anyway. */

  /*
For both line and continuum images we have to access (currently in traceray() and sourceFunc_cont()) array elements m[i].local_cmb[j], g[id].mol[i].dust[j] and g[id].mol[i].knu[j]. For a continuum image, the molData object and the grid.mol object are simply convenient (if misleadingly named in this instance) repositories of appropriate cmb, dust and knu information; there is no actual molecule involved, and i and j are both simply 0.

For a line image however, we may have a problem, because of the pair of quantities img.freq and img.trans, the user is allowed to specify freq and not trans. Since the cmb, dust and knu are assumed to be slowly-varying quantities, we can for 'trans' in this case use the line closest to the image frequency. (There must be some reasonably close line, else we would see no line emission in the image.)
  */
  if(img[im].doline){
    if (img[im].trans>=0){
      cmbMolI  = img[im].molI;
      cmbLineI = img[im].trans;

    }else{ /* User didn't set trans. Find the nearest line to the image frequency. */
      for(molI=0;molI<par->nSpecies;molI++){
        for(lineI=0;lineI<md[molI].nline;lineI++){
          absDeltaFreq = fabs(img[im].freq - md[molI].freq[lineI]);
          if((molI==0 && lineI==0) || absDeltaFreq < minfreq){
            minfreq = absDeltaFreq;
            cmbMolI = molI;
            cmbLineI = lineI;
          }
        }
      }
    }
  }else{ /* continuum image */
    cmbMolI = 0;
    cmbLineI = 0;
  }

  cutoff = par->minScale*1.0e-7;

  for(px=0;px<(img[im].pxls*img[im].pxls);px++){
    for(ichan=0;ichan<img[im].nchan;ichan++){
      img[im].pixel[px].intense[ichan]=0.0;
      img[im].pixel[px].tau[ichan]=0.0;
    }
  }

  defaultErrorHandler = gsl_set_error_handler_off();
  /*
The GSL documentation does not recommend leaving the error handler at the default within multi-threaded code.

While this is off however, gsl_* calls will not exit if they encounter a problem. We may need to pay some attention to trapping their errors.
  */

  nRaysDone=0;
  omp_set_dynamic(0);
  #pragma omp parallel private(px,aa,threadI) num_threads(par->nThreads)
  {
    threadI = omp_get_thread_num();

    /* Declaration of thread-private pointers. */
    rayData ray;
    ray.intensity=malloc(sizeof(double) * img[im].nchan);
    ray.tau=malloc(sizeof(double) * img[im].nchan);

    #pragma omp for
    /* Main loop through pixel grid. */
    for(px=0;px<(img[im].pxls*img[im].pxls);px++){
      #pragma omp atomic
      ++nRaysDone;

      for(aa=0;aa<par->antialias;aa++){
        ray.x = -size*(gsl_rng_uniform(threadRans[threadI]) + px%img[im].pxls - 0.5*img[im].pxls);
        ray.y =  size*(gsl_rng_uniform(threadRans[threadI]) + px/img[im].pxls - 0.5*img[im].pxls);

        traceray(ray, cmbMolI, cmbLineI, im, par, gp, md, img, cutoff);

        #pragma omp critical
        {
          for(ichan=0;ichan<img[im].nchan;ichan++){
            img[im].pixel[px].intense[ichan] += ray.intensity[ichan]/(double) par->antialias;
            img[im].pixel[px].tau[ichan] += ray.tau[ichan]/(double) par->antialias;
          }
        }
      }
      if (threadI == 0){ /* i.e., is master thread */
        if(!silent) progressbar((double)(nRaysDone)/totalNumPixelsMinus1, 13);
      }
    }

    free(ray.tau);
    free(ray.intensity);
  } /* End of parallel block. */

  gsl_set_error_handler(defaultErrorHandler);

  for (i=0;i<par->nThreads;i++){
    gsl_rng_free(threadRans[i]);
  }
  free(threadRans);
  gsl_rng_free(ran);
}

