/*
 *  raytrace.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2016 The LIME development team
 *
TODO:
  - In raytrace(), look at rearranging the code to do the qhull step before choosing the rays. This would allow cells with all vertices outside the image boundaries to be excluded. If the image is much smaller than the model, this could lead to significant savings in time. The only downside might be memory useage...
  - We should not be multiplying remnantSnu by md[0].norminv, but to make it correct would be pretty fiddly.
 */

#include "lime.h"


/*....................................................................*/
void calcLineAmpSample(const double x[3], const double dx[3], const double ds\
  , const double binv, double *projVels, const int nSteps\
  , const double oneOnNSteps, const double deltav, double *vfac){
  /*
The bulk velocity of the model material can vary significantly with position, thus so can the value of the line-shape function at a given frequency and direction. The present function calculates 'vfac', an approximate average of the line-shape function along a path of length ds in the direction of the line of sight.
  */
  int i;
  double v,val;

  *vfac=0.;
  for(i=0;i<nSteps;i++){
    v = deltav - projVels[i]; /* projVels contains the component of the local bulk velocity in the direction of dx, whereas deltav is the recession velocity of the channel we are interested in (corrected for bulk source velocity and line displacement from the nominal frequency). Remember also that, since dx points away from the observer, positive values of the projected velocity also represent recessions. Line centre occurs when v==0, i.e. when deltav==veloproject(dx,vel). That is the reason for the subtraction here. */
    val=fabs(v)*binv;
    if(val <=  2500.){
#ifdef FASTEXP
      *vfac+= FastExp(val*val);
#else
      *vfac+=   exp(-(val*val));
#endif
    }
  }
  *vfac *= oneOnNSteps;
  return;
}

/*....................................................................*/
void calcLineAmpInterp(const double projVelRay, const double binv\
  , const double deltav, double *vfac){
  /*
The bulk velocity of the model material can vary significantly with position, thus so can the value of the line-shape function at a given frequency and direction. The present function calculates 'vfac', an approximate average of the line-shape function along a path of length ds in the direction of the line of sight.
  */
  double v,val;

  v = deltav - projVelRay; /* projVelRay is the component of the local bulk velocity in the direction of the ray, whereas deltav is the recession velocity of the channel we are interested in (corrected for bulk source velocity and line displacement from the nominal frequency). Remember also that, since the ray points away from the observer, positive values of the projected velocity also represent recessions. Line centre occurs when v==0, i.e. when deltav==projVelRay. That is the reason for the subtraction here. */
  val = fabs(v)*binv;
  if(val <=  2500.){
#ifdef FASTEXP
    *vfac = FastExp(val*val);
#else
    *vfac =   exp(-(val*val));
#endif
  }else
    *vfac = 0.0;
}

/*....................................................................*/
void
line_plane_intersect(struct grid *gp, double *ds, int posn, int *nposn, double *dx, double *x, double cutoff){
  /*
This function returns ds as the (always positive-valued) distance between the present value of x and the next Voronoi face in the direction of vector dx, and nposn as the id of the grid cell that abuts that face. 
  */
  double newdist, numerator, denominator ;
  int i;

  for(i=0;i<gp[posn].numNeigh;i++) {
    /* Find the shortest distance between (x,y,z) and any of the posn Voronoi faces */
    /* ds=(p0-l0) dot n / l dot n */

    numerator=((gp[posn].x[0]+gp[posn].dir[i].x[0]/2. - x[0]) * gp[posn].dir[i].x[0]+
               (gp[posn].x[1]+gp[posn].dir[i].x[1]/2. - x[1]) * gp[posn].dir[i].x[1]+
               (gp[posn].x[2]+gp[posn].dir[i].x[2]/2. - x[2]) * gp[posn].dir[i].x[2]);

    denominator=(dx[0]*gp[posn].dir[i].x[0]+dx[1]*gp[posn].dir[i].x[1]+dx[2]*gp[posn].dir[i].x[2]);
    
    if(fabs(denominator) > 0){
      newdist=numerator/denominator;
      if(newdist<*ds && newdist > cutoff){
        *ds=newdist;
        *nposn=gp[posn].neigh[i]->id;
      }
    }
  }
  if(*nposn==-1) *nposn=posn;
}

/*....................................................................*/
void
traceray(rayData ray, const int cmbMolI, const int cmbLineI, const int im\
  , configInfo *par, struct grid *gp, molData *md, image *img, struct gAuxType *gAux\
  , const int nlinetot, int *allLineMolIs, int *allLineLineIs, const double cutoff\
  , const int nSteps, const double oneOnNSteps){
  /*
For a given image pixel position, this function evaluates the intensity of the total light emitted/absorbed along that line of sight through the (possibly rotated) model. The calculation is performed for several frequencies, one per channel of the output image.

Note that the algorithm employed here is similar to that employed in the function photon() which calculates the average radiant flux impinging on a grid cell: namely the notional photon is started at the side of the model near the observer and 'propagated' in the receding direction until it 'reaches' the far side. This is rather non-physical in conception but it makes the calculation easier.
  */
  const int stokesIi=0;
  int ichan,stokesId,di,i,posn,nposn,polMolI,polLineI,molI,lineI;
  double xp,yp,zp,x[DIM],dx[DIM],dist2,ndist2,col,ds,snu_pol[3],dtau;
  double contJnu,contAlpha,jnu,alpha,lineRedShift,vThisChan,deltav,vfac=0.;
  double remnantSnu,expDTau,brightnessIncrement;
  double projVels[nSteps],d,vel[DIM];

  for(ichan=0;ichan<img[im].nchan;ichan++){
    ray.tau[ichan]=0.0;
    ray.intensity[ichan]=0.0;
  }

  xp=ray.x;
  yp=ray.y;

  /* The model is circular in projection. We only follow the ray if it will intersect the model.
  */
  if((xp*xp+yp*yp)>par->radiusSqu)
    return;

  zp=-sqrt(par->radiusSqu-(xp*xp+yp*yp)); /* There are two points of intersection between the line of sight and the spherical model surface; this is the Z coordinate (in the unrotated frame) of the one nearer to the observer. */

  /* Rotate the line of sight as desired. */
  for(di=0;di<DIM;di++){
    x[di]=xp*img[im].rotMat[di][0] + yp*img[im].rotMat[di][1] + zp*img[im].rotMat[di][2];
    dx[di]= img[im].rotMat[di][2]; /* This points away from the observer. */
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
      polMolI = 0; /****** Always?? */
      polLineI = 0; /****** Always?? */
      sourceFunc_pol(gp[posn].B, gAux[posn].mol[polMolI], polLineI, img[im].rotMat, snu_pol, &alpha);
      dtau=alpha*ds;
      calcSourceFn(dtau, par, &remnantSnu, &expDTau);
      remnantSnu *= md[polMolI].norminv*ds;

      for(stokesId=0;stokesId<img[im].nchan;stokesId++){ /* Loop over I, Q and U */
#ifdef FASTEXP
        brightnessIncrement = FastExp(ray.tau[stokesId])*remnantSnu*snu_pol[stokesId];
#else
        brightnessIncrement =    exp(-ray.tau[stokesId])*remnantSnu*snu_pol[stokesId];
#endif
        ray.intensity[stokesId] += brightnessIncrement;
        ray.tau[stokesId]+=dtau; //**** But this will be the same for I, Q or U.
      }
    } else {
      if(!par->doPregrid){
        for(i=0;i<nSteps;i++){
          d = i*ds*oneOnNSteps;
          velocity(x[0]+(dx[0]*d),x[1]+(dx[1]*d),x[2]+(dx[2]*d),vel);
          projVels[i] = veloproject(dx,vel);
        }
      }

      /* Calculate first the continuum stuff because it is the same for all channels:
      */
      contJnu = 0.0;
      contAlpha = 0.0;
      sourceFunc_cont_raytrace(gAux[posn].mol[cmbMolI], cmbLineI, &contJnu, &contAlpha);

      for(ichan=0;ichan<img[im].nchan;ichan++){
        jnu = contJnu;
        alpha = contAlpha;
        vThisChan = (ichan-(img[im].nchan-1)*0.5)*img[im].velres; /* Consistent with the WCS definition in writefits(). */

        if(img[im].doline){
          for(molI=0;molI<par->nSpecies;molI++){
            for(lineI=0;lineI<md[molI].nline;lineI++){
              if(md[molI].freq[lineI] > img[im].freq-img[im].bandwidth*0.5\
              && md[molI].freq[lineI] < img[im].freq+img[im].bandwidth*0.5){
                /* Calculate the red shift of the transition wrt to the frequency specified for the image.
                */
                if(img[im].trans > -1){
                  lineRedShift=(md[molI].freq[img[im].trans]-md[molI].freq[lineI])/md[molI].freq[img[im].trans]*CLIGHT;
                } else {
                  lineRedShift=(img[im].freq-md[molI].freq[lineI])/img[im].freq*CLIGHT;
                }

                deltav = vThisChan - img[im].source_vel - lineRedShift;
                /* Line centre occurs when deltav = the recession velocity of the radiating material. Explanation of the signs of the 2nd and 3rd terms on the RHS: (i) A bulk source velocity (which is defined as >0 for the receding direction) should be added to the material velocity field; this is equivalent to subtracting it from deltav, as here. (ii) A positive value of lineRedShift means the line is red-shifted wrt to the frequency specified for the image. The effect is the same as if the line and image frequencies were the same, but the bulk recession velocity were higher. lineRedShift should thus be added to the recession velocity, which is equivalent to subtracting it from deltav, as here. */

                /* Calculate an approximate average line-shape function at deltav within the Voronoi cell. */
                if(!par->doPregrid)
                  calcLineAmpSample(x,dx,ds,gp[posn].mol[molI].binv,projVels,nSteps,oneOnNSteps,deltav,&vfac);
                else
                  vfac = gaussline(deltav-veloproject(dx,gp[posn].vel),gp[posn].mol[molI].binv);

                /* Increment jnu and alpha for this Voronoi cell by the amounts appropriate to the spectral line. */
                sourceFunc_line_raytrace(md[molI],vfac,gAux[posn].mol[molI],lineI,&jnu,&alpha);
              }
            }
          }
        }

        dtau=alpha*ds;
//???          if(dtau < -30) dtau = -30; // as in photon()?
        calcSourceFn(dtau, par, &remnantSnu, &expDTau);
        remnantSnu *= jnu*md[0].norminv*ds;
#ifdef FASTEXP
        brightnessIncrement = FastExp(ray.tau[ichan])*remnantSnu;
#else
        brightnessIncrement =    exp(-ray.tau[ichan])*remnantSnu;
#endif
        ray.intensity[ichan] += brightnessIncrement;
        ray.tau[ichan]+=dtau;
      }
    }

    /* Move the working point to the edge of the next Voronoi cell. */
    for(di=0;di<DIM;di++) x[di]+=ds*dx[di];
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

/*....................................................................*/
void traceray_smooth(rayData ray, const int cmbMolI, const int cmbLineI, const int im\
  , configInfo *par, struct grid *gp, molData *md, image *img, struct gAuxType *gAux\
  , const int nlinetot, int *allLineMolIs, int *allLineLineIs, struct cell *dc\
  , const unsigned long numCells, const double epsilon, gridInterp gips[3]\
  , const int numSegments, const double oneOnNumSegments, const int nSteps, const double oneOnNSteps){
  /*
For a given image pixel position, this function evaluates the intensity of the total light emitted/absorbed along that line of sight through the (possibly rotated) model. The calculation is performed for several frequencies, one per channel of the output image.

Note that the algorithm employed here to solve the RTE is similar to that employed in the function photon() which calculates the average radiant flux impinging on a grid cell: namely the notional photon is started at the side of the model near the observer and 'propagated' in the receding direction until it 'reaches' the far side. This is rather non-physical in conception but it makes the calculation easier.

This version of traceray implements a new algorithm in which the population values are interpolated linearly from those at the vertices of the Delaunay cell which the working point falls within.
  */
  const int stokesIi=0;
  const int numFaces = DIM+1, nVertPerFace=3;
  int ichan,stokesId,di,status,lenChainPtrs,entryI,exitI,vi,vvi,ci;
  int si,polMolI,polLineI,molI,lineI;
  double xp,yp,zp,x[DIM],dir[DIM],projVelRay,vel[DIM];
  double xCmpntsRay[nVertPerFace], ds, snu_pol[3], dtau, contJnu, contAlpha;
  double jnu, alpha, lineRedShift, vThisChan, deltav, vfac, remnantSnu, expDTau;
  double brightnessIncrement;
  intersectType entryIntcptFirstCell, *cellExitIntcpts=NULL;
  unsigned long *chainOfCellIds=NULL, dci;
  unsigned long gis[2][nVertPerFace];

  for(ichan=0;ichan<img[im].nchan;ichan++){
    ray.tau[ichan]=0.0;
    ray.intensity[ichan]=0.0;
  }

  xp=ray.x;
  yp=ray.y;

  /* The model is circular in projection. We only follow the ray if it will intersect the model.
  */
  if((xp*xp+yp*yp)>par->radiusSqu)
    return;

  zp=-sqrt(par->radiusSqu-(xp*xp+yp*yp)); /* There are two points of intersection between the line of sight and the spherical model surface; this is the Z coordinate (in the unrotated frame) of the one nearer to the observer. */

  /* Rotate the line of sight as desired. */
  for(di=0;di<DIM;di++){
    x[di]=xp*img[im].rotMat[di][0] + yp*img[im].rotMat[di][1] + zp*img[im].rotMat[di][2];
    dir[di]= img[im].rotMat[di][2]; /* This points away from the observer. */
  }

  /* Find the chain of cells the ray passes through.
  */
  status = followRayThroughDelCells(x, dir, gp, dc, numCells, epsilon\
    , &entryIntcptFirstCell, &chainOfCellIds, &cellExitIntcpts, &lenChainPtrs);

  if(status!=0){
    free(chainOfCellIds);
    free(cellExitIntcpts);
    return;
  }

  entryI = 0;
  exitI  = 1;
  dci = chainOfCellIds[0];

  /* Obtain the indices of the grid points on the vertices of the entry face.
  */
  vvi = 0;
  for(vi=0;vi<numFaces;vi++){
    if(vi!=entryIntcptFirstCell.fi){
      gis[entryI][vvi++] = dc[dci].vertx[vi]->id;
    }
  }

  /* Calculate, for each of the 3 vertices of the entry face, the displacement components in the direction of 'dir'. *** NOTE *** that if all the rays are parallel, we could precalculate these for all the vertices.
  */
  for(vi=0;vi<nVertPerFace;vi++)
    xCmpntsRay[vi] = veloproject(dir, gp[gis[entryI][vi]].x);

  doBaryInterp(entryIntcptFirstCell, gp, gAux, xCmpntsRay, gis[entryI]\
    , md, par->nSpecies, &gips[entryI]);

  for(ci=0;ci<lenChainPtrs;ci++){
    /* For each cell we have 2 data structures which give information about respectively the entry and exit points of the ray, including the barycentric coordinates of the intersection point between the ray and the appropriate face of the cell. (If we follow rays in 3D space then the cells will be tetrahedra and the faces triangles.) If we know the value of a quantity Q for each of the vertices, then the linear interpolation of the Q values for any face is (for a 3D space) bary[0]*Q[0] + bary[1]*Q[1] + bary[2]*Q[2], where the indices are taken to run over the vertices of that face. Thus we can calculate the interpolated values Q_entry and Q_exit. Further linear interpolation along the path between entry and exit is straightforward.
    */

    dci = chainOfCellIds[ci];

    /* Obtain the indices of the grid points on the vertices of the exit face. */
    vvi = 0;
    for(vi=0;vi<numFaces;vi++){
      if(vi!=cellExitIntcpts[ci].fi){
        gis[exitI][vvi++] = dc[dci].vertx[vi]->id;
      }
    }

    /* Calculate, for each of the 3 vertices of the exit face, the displacement components in the direction of 'dir'. *** NOTE *** that if all the rays are parallel, we could precalculate these for all the vertices.
    */
    for(vi=0;vi<nVertPerFace;vi++)
      xCmpntsRay[vi] = veloproject(dir, gp[gis[exitI][vi]].x);

    doBaryInterp(cellExitIntcpts[ci], gp, gAux, xCmpntsRay, gis[exitI]\
      , md, par->nSpecies, &gips[exitI]);

    /* At this point we have interpolated all the values of interest to both the entry and exit points of the cell. Now we break the path between entry and exit into several segments and calculate all these values at the midpoint of each segment.

At the moment I will fix the number of segments, but it might possibly be faster to rather have a fixed segment length (in barycentric coordinates) and vary the number depending on how many of these lengths fit in the path between entry and exit.
    */
    ds = (gips[exitI].xCmpntRay - gips[entryI].xCmpntRay)*oneOnNumSegments;

    for(si=0;si<numSegments;si++){
      doSegmentInterp(gips, entryI, md, par->nSpecies, oneOnNumSegments, si);

      if(par->polarization){
        polMolI = 0; /****** Always?? */
        polLineI = 0; /****** Always?? */
        sourceFunc_pol(gips[2].B, gips[2].mol[polMolI], polLineI, img[im].rotMat, snu_pol, &alpha);
        dtau=alpha*ds;
        calcSourceFn(dtau, par, &remnantSnu, &expDTau);
        remnantSnu *= md[polMolI].norminv*ds;

        for(stokesId=0;stokesId<img[im].nchan;stokesId++){ /* Loop over I, Q and U */
#ifdef FASTEXP
          brightnessIncrement = FastExp(ray.tau[stokesId])*remnantSnu*snu_pol[stokesId];
#else
          brightnessIncrement =    exp(-ray.tau[stokesId])*remnantSnu*snu_pol[stokesId];
#endif
          ray.intensity[stokesId] += brightnessIncrement;
          ray.tau[stokesId]+=dtau; //**** But this will be the same for I, Q or U.
        }
      } else {
        /* It appears to be necessary to sample the velocity function in the following way rather than interpolating it from the vertices of the Delaunay cell in the same way as with all the other quantities of interest. Velocity varies too much across the cells, and in a nonlinear way, for linear interpolation to yield a totally satisfactory result.
        */
        velocity(gips[2].x[0], gips[2].x[1], gips[2].x[2], vel);
        projVelRay = veloproject(dir, vel);

        /* Calculate first the continuum stuff because it is the same for all channels:
        */
        contJnu = 0.0;
        contAlpha = 0.0;
        sourceFunc_cont_raytrace(gips[2].mol[cmbMolI], cmbLineI, &contJnu, &contAlpha);

        for(ichan=0;ichan<img[im].nchan;ichan++){
          jnu = contJnu;
          alpha = contAlpha;
          vThisChan=(ichan-(img[im].nchan-1)*0.5)*img[im].velres; /* Consistent with the WCS definition in writefits(). */

          if(img[im].doline){
            for(molI=0;molI<par->nSpecies;molI++){
              for(lineI=0;lineI<md[molI].nline;lineI++){
                if(md[molI].freq[lineI] > img[im].freq-img[im].bandwidth*0.5\
                && md[molI].freq[lineI] < img[im].freq+img[im].bandwidth*0.5){
                  /* Calculate the red shift of the transition wrt to the frequency specified for the image.
                  */
                  if(img[im].trans > -1){
                    lineRedShift=(md[molI].freq[img[im].trans]-md[molI].freq[lineI])/md[molI].freq[img[im].trans]*CLIGHT;
                  } else {
                    lineRedShift=(img[im].freq-md[molI].freq[lineI])/img[im].freq*CLIGHT;
                  }

                  deltav = vThisChan - img[im].source_vel - lineRedShift;
                  /* Line centre occurs when deltav = the recession velocity of the radiating material. Explanation of the signs of the 2nd and 3rd terms on the RHS: (i) A bulk source velocity (which is defined as >0 for the receding direction) should be added to the material velocity field; this is equivalent to subtracting it from deltav, as here. (ii) A positive value of lineRedShift means the line is red-shifted wrt to the frequency specified for the image. The effect is the same as if the line and image frequencies were the same, but the bulk recession velocity were higher. lineRedShift should thus be added to the recession velocity, which is equivalent to subtracting it from deltav, as here. */

                  calcLineAmpInterp(projVelRay, gips[2].mol[molI].binv, deltav, &vfac);

                  /* Increment jnu and alpha for this Voronoi cell by the amounts appropriate to the spectral line.
                  */
                  sourceFunc_line_raytrace(md[molI], vfac, gips[2].mol[molI], lineI, &jnu, &alpha);
                } /* end if within freq range. */
              } /* end loop over lines this mol. */
            } /* end loop over all mols. */
          } /* end if doLine. */

          dtau = alpha*ds;
//???          if(dtau < -30) dtau = -30; // as in photon()?
          calcSourceFn(dtau, par, &remnantSnu, &expDTau);
          remnantSnu *= jnu*md[0].norminv*ds;
#ifdef FASTEXP
          brightnessIncrement = FastExp(ray.tau[ichan])*remnantSnu;
#else
          brightnessIncrement =    exp(-ray.tau[ichan])*remnantSnu;
#endif
          ray.intensity[ichan] += brightnessIncrement;
          ray.tau[ichan] += dtau;
        } /* End loop over channels. */
      } /* End if(par->polarization). */
    } /* End loop over segments within cell. */

    entryI = exitI;
    exitI = 1 - exitI;
  } /* End loop over cells in the chain traversed by the ray. */

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

  free(chainOfCellIds);
  free(cellExitIntcpts);
}

/*....................................................................*/
void
raytrace(int im, configInfo *par, struct grid *gp, molData *md, image *img){
  /*
This function constructs an image cube by following sets of rays (at least 1 per image pixel) through the model, solving the radiative transfer equations as appropriate for each ray. The ray locations within each pixel are chosen randomly within the pixel, but the number of rays per pixel is set equal to the number of projected model grid points falling within that pixel, down to a minimum equal to par->alias.
  */
  const int maxNumRaysPerPixel=20; /**** Arbitrary - could make this a global, or an argument. Set it to zero to indicate there is no maximum. */
  const double cutoff = par->minScale*1.0e-7;
  const int numFaces=1+DIM, numInterpPoints=3, numSegments=5, minNumRaysForAverage=2;
  const double oneOnNFaces=1.0/(double)numFaces, oneOnNumSegments = 1.0/(double)numSegments;
  const double epsilon = 1.0e-6; // Needs thinking about. Double precision is much smaller than this.
  const int nStepsThruCell=10;
  const double oneOnNSteps=1.0/(double)nStepsThruCell;

  double size,oneOnNumActiveRaysMinus1,imgCentreXPixels,imgCentreYPixels,minfreq,absDeltaFreq,x[2],sum,oneOnNumRays;
  unsigned int totalNumImagePixels,ppi,numPixelsForInterp;
  int gi,molI,lineI,ei,li,nlinetot,tmptrans,ichan,numActiveRays,i,di,xi,yi,ri,c,id,ids[3],vi;
  int *allLineMolIs,*allLineLineIs,cmbMolI,cmbLineI;
  rayData *rays;
  struct cell *dc=NULL;
  unsigned long numCells,dci;
  struct gAuxType *gAux=NULL; /* This will hold some precalculated values for the grid points. */
  coordT *pt_array, point[3];
  char flags[255];
  boolT ismalloc = False,isoutside;
  realT bestdist;
  facetT *facet;
  vertexT *vertex,**vertexp;
  int curlong, totlong;
  double triangle[3][2],barys[3];
  _Bool isOutsideImage;
  gsl_error_handler_t *defaultErrorHandler=NULL;

  size = img[im].distance*img[im].imgres;
  totalNumImagePixels = img[im].pxls*img[im].pxls;
  imgCentreXPixels = img[im].pxls/2.0;
  imgCentreYPixels = img[im].pxls/2.0;

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

  for(ppi=0;ppi<totalNumImagePixels;ppi++){
    for(ichan=0;ichan<img[im].nchan;ichan++){
      img[im].pixel[ppi].intense[ichan] = 0.0;
      img[im].pixel[ppi].tau[    ichan] = 0.0;
    }
  }

  for(ppi=0;ppi<totalNumImagePixels;ppi++)
    img[im].pixel[ppi].numRays = 0;

  /* The following is the first of the 3 main loops in raytrace. Here we loop over the (internal or non-sink) grid points. We're doing 2 things: loading the rotated, projected coordinates into the rays list, and counting the rays per image pixel.
  */
  rays = malloc(sizeof(rayData)*par->pIntensity); /* We may need to reallocate this later. */
  numActiveRays = 0;
  for(gi=0;gi<par->pIntensity;gi++){
    /* Apply the inverse (i.e. transpose) rotation matrix. (We use the inverse matrix here because here we want to rotate grid coordinates to the observer frame, whereas inside traceray() we rotate observer coordinates to the grid frame.)
    */
    for(i=0;i<2;i++){
      x[i]=0.0;
      for(di=0;di<DIM;di++){
        x[i] += gp[gi].x[di]*img[im].rotMat[di][i];
      }
    }

    /* Calculate which pixel the projected position (x[0],x[1]) falls within.
    */
    xi = floor(x[0]/size + imgCentreXPixels);
    yi = floor(x[1]/size + imgCentreYPixels);
    if(xi<0 || xi>=img[im].pxls || yi<0 || yi>=img[im].pxls){
      isOutsideImage = 1;
      ppi = 0; /* Under these circumstances it ought never to be accessed, but it is not good to leave it without a value at all. */
    }else{
      isOutsideImage = 0;
      ppi = (unsigned int)yi*(unsigned int)img[im].pxls + (unsigned int)xi;
    }

    /* See if we want to keep the ray. For the time being we will include those outside the image bounds, but a cleverer algorithm would exclude some of them. Note that maxNumRaysPerPixel<1 is used to flag that there is no upper limit to the number of rays per pixel.
    */
    if(isOutsideImage || maxNumRaysPerPixel<1 || img[im].pixel[ppi].numRays<maxNumRaysPerPixel){
      if(!isOutsideImage)
        img[im].pixel[ppi].numRays++;

      rays[numActiveRays].ppi = ppi;
      rays[numActiveRays].x = x[0];
      rays[numActiveRays].y = x[1];
      rays[numActiveRays].tau       = malloc(sizeof(double)*img[im].nchan);
      rays[numActiveRays].intensity = malloc(sizeof(double)*img[im].nchan);
      for(ichan=0;ichan<img[im].nchan;ichan++) {
        rays[numActiveRays].tau[ichan] = 0.0;
        rays[numActiveRays].intensity[ichan] = 0.0;
      }

      numActiveRays++;
    }
  } /* End loop 1, over grid points. */

  oneOnNumActiveRaysMinus1 = 1.0/(double)(numActiveRays-1);

  if(numActiveRays<par->pIntensity)
    rays = realloc(rays, sizeof(rayData)*numActiveRays);

  if(par->traceRayAlgorithm==1){
    delaunay(DIM, gp, (unsigned long)par->ncell, 1, &dc, &numCells); /* mallocs dc if getCells==T */

    /* We need to process the list of cells a bit further - calculate their centres, and reset the id values to be the same as the index of the cell in the list. (This last because we are going to construct other lists to indicate which cells have been visited etc.)
    */
    for(dci=0;dci<numCells;dci++){
      for(di=0;di<DIM;di++){
        sum = 0.0;
        for(vi=0;vi<numFaces;vi++){
          sum += dc[dci].vertx[vi]->x[di];
        }
        dc[dci].centre[di] = sum*oneOnNFaces;
      }

      dc[dci].id = dci;
    }

  }else if(par->traceRayAlgorithm!=0){
    if(!silent) bail_out("Unrecognized value of par.traceRayAlgorithm");
    exit(1);
  }

  /* Precalculate binv*nmol*pops for all grid points.
  */
  gAux = malloc(sizeof(*gAux)*par->ncell);
  for(gi=0;gi<par->ncell;gi++){
    gAux[gi].mol = malloc(sizeof(*(gAux[gi].mol))*par->nSpecies);
    for(molI=0;molI<par->nSpecies;molI++){
      gAux[gi].mol[molI].specNumDens = NULL;
      gAux[gi].mol[molI].dust        = malloc(sizeof(*(gAux[gi].mol[molI].dust))       *md[molI].nline);
      gAux[gi].mol[molI].knu         = malloc(sizeof(*(gAux[gi].mol[molI].knu))        *md[molI].nline);

      /* This next is repetition. I do it in order to be able to use the same sourcefunc.c functions for the interpolated grid values as for the 'standard' algorithm. With a sensible arrangement of memory for the grid values, this would be unnecessary.
      */
      for(li=0;li<md[molI].nline;li++){
        gAux[gi].mol[molI].dust[li] = gp[gi].mol[molI].dust[li];
        gAux[gi].mol[molI].knu[ li] = gp[gi].mol[molI].knu[ li];
      }
    }
  }

  if(img[im].doline){
    for(gi=0;gi<par->ncell;gi++){
      for(molI=0;molI<par->nSpecies;molI++){
        gAux[gi].mol[molI].specNumDens = malloc(sizeof(*(gAux[gi].mol[molI].specNumDens))*md[molI].nlev);

        for(ei=0;ei<md[molI].nlev;ei++)
          gAux[gi].mol[molI].specNumDens[ei] = gp[gi].mol[molI].binv\
            *gp[gi].mol[molI].nmol*gp[gi].mol[molI].pops[ei];
      }
    }
  }

  /* This is the start of loop 2/3, which loops over the rays. We trace each ray, then load into the image cube those for which the number of rays per pixel exceeds a minimum. The remaining image pixels we handle via an interpolation algorithm in loop 3.
  */
  defaultErrorHandler = gsl_set_error_handler_off();
  /*
The GSL documentation does not recommend leaving the error handler at the default within multi-threaded code.

While this is off however, gsl_* calls will not exit if they encounter a problem. We may need to pay some attention to trapping their errors.
  */

  omp_set_dynamic(0);
  #pragma omp parallel num_threads(par->nThreads)
  {
    /* Declaration of thread-private pointers.
    */
    int threadI = omp_get_thread_num();
    int ii, si, ri;
    gridInterp gips[numInterpPoints];

    if(par->traceRayAlgorithm==1){
      /* Allocate memory for the interpolation points:
      */
      for(ii=0;ii<numInterpPoints;ii++){
        gips[ii].mol = malloc(sizeof(*(gips[ii].mol))*par->nSpecies);
        for(si=0;si<par->nSpecies;si++){
          gips[ii].mol[si].specNumDens = malloc(sizeof(*(gips[ii].mol[si].specNumDens))*md[si].nlev);
          gips[ii].mol[si].dust        = malloc(sizeof(*(gips[ii].mol[si].dust))       *md[si].nline);
          gips[ii].mol[si].knu         = malloc(sizeof(*(gips[ii].mol[si].knu))        *md[si].nline);
        }
      }
    }

    #pragma omp for schedule(dynamic)
    for(ri=0;ri<numActiveRays;ri++){
      if(par->traceRayAlgorithm==0)
        traceray(rays[ri], cmbMolI, cmbLineI, im, par, gp, md, img, gAux\
          , nlinetot, allLineMolIs, allLineLineIs, cutoff, nStepsThruCell, oneOnNSteps);
      else if(par->traceRayAlgorithm==1)
        traceray_smooth(rays[ri], cmbMolI, cmbLineI, im, par, gp, md, img, gAux\
          , nlinetot, allLineMolIs, allLineLineIs, dc, numCells, epsilon, gips\
          , numSegments, oneOnNumSegments, nStepsThruCell, oneOnNSteps);

      if (threadI == 0){ /* i.e., is master thread */
        if(!silent) progressbar((double)(ri)*oneOnNumActiveRaysMinus1, 13);
      }
    }

    if(par->traceRayAlgorithm==1){
      for(ii=0;ii<numInterpPoints;ii++)
        freePop2(par->nSpecies, gips[ii].mol);
    }
  } /* End of parallel block. */

  gsl_set_error_handler(defaultErrorHandler);

  /* For pixels with more than a cutoff number of rays, just average those rays into the pixel:
  */
  for(ri=0;ri<numActiveRays;ri++){
    if(img[im].pixel[rays[ri].ppi].numRays >= minNumRaysForAverage){
      for(ichan=0;ichan<img[im].nchan;ichan++){
        img[im].pixel[rays[ri].ppi].intense[ichan] += rays[ri].intensity[ichan];
        img[im].pixel[rays[ri].ppi].tau[    ichan] += rays[ri].tau[      ichan];
      }
    }
  }
  numPixelsForInterp = 0;
  for(ppi=0;ppi<totalNumImagePixels;ppi++){
    if(img[im].pixel[ppi].numRays >= minNumRaysForAverage){
      oneOnNumRays = 1.0/(double)img[im].pixel[ppi].numRays;
      for(ichan=0;ichan<img[im].nchan;ichan++){
        img[im].pixel[ppi].intense[ichan] *= oneOnNumRays;
        img[im].pixel[ppi].tau[    ichan] *= oneOnNumRays;
      }
    }else
      numPixelsForInterp++;
  }

  if(numPixelsForInterp>0){
    /* Now we enter main loop 3/3, in which we loop over image pixels, and for any we need to interpolate, we do so. But first we need to invoke qhull to get a Delaunay triangulation of the projected points.
    */
    pt_array=malloc(sizeof(coordT)*2*numActiveRays);

    for(ri=0;ri<numActiveRays;ri++) {
      pt_array[ri*2+0] = rays[ri].x;
      pt_array[ri*2+1] = rays[ri].y;
    }

    sprintf(flags,"qhull d Qbb");
    if(qh_new_qhull(2, numActiveRays, pt_array, ismalloc, flags, NULL, NULL)) {
      if(!silent) bail_out("Qhull failed to triangulate");
      exit(1);
    }

    point[2]=0.;
    for(ppi=0;ppi<totalNumImagePixels;ppi++){
      if(img[im].pixel[ppi].numRays < minNumRaysForAverage){
        xi = (int)(ppi%(unsigned int)img[im].pxls);
        yi = floor(ppi/(double)img[im].pxls);
        point[0] = size*(0.5 + xi - imgCentreXPixels);
        point[1] = size*(0.5 + yi - imgCentreYPixels);

        qh_setdelaunay (3, 1, point);
        facet = qh_findbestfacet (point, qh_ALL, &bestdist, &isoutside);
        if(isoutside){
          c=0;
          FOREACHvertex_( facet->vertices ) {
            id = qh_pointid(vertex->point);
            triangle[c][0] = rays[id].x;
            triangle[c][1] = rays[id].y;
            ids[c]=id;
            c++;
          }

          calcTriangleBaryCoords(triangle, (double)point[0], (double)point[1], barys);

          /* Interpolate: */
          for(ichan=0;ichan<img[im].nchan;ichan++){
            img[im].pixel[ppi].intense[ichan] += barys[0]*rays[ids[0]].intensity[ichan]\
                                               + barys[1]*rays[ids[1]].intensity[ichan]\
                                               + barys[2]*rays[ids[2]].intensity[ichan];
            img[im].pixel[ppi].tau[    ichan] += barys[0]*rays[ids[0]].tau[ichan]\
                                               + barys[1]*rays[ids[1]].tau[ichan]\
                                               + barys[2]*rays[ids[2]].tau[ichan];
          }
        } /* end if !isoutside */
      } /* end if(img[im].pixel[ppi].numRays < minNumRaysForAverage) */
    } /* end loop over image pixels */

    qh_freeqhull(!qh_ALL);
    qh_memfreeshort (&curlong, &totlong);
    free(pt_array);
  } /* end if(numPixelsForInterp>0) */

  img[im].trans=tmptrans;

  freeGAux((unsigned long)par->ncell, par->nSpecies, gAux);
  if(par->traceRayAlgorithm==1)
    free(dc);
  for(ri=0;ri<numActiveRays;ri++){
    free(rays[ri].tau);
    free(rays[ri].intensity);
  }
  free(rays);
  free(allLineMolIs);
  free(allLineLineIs);
}

/*....................................................................*/
void calcTriangleBaryCoords(double vertices[3][2], double x, double y, double barys[3]){
  double mat[2][2], vec[2], det;
  /*
The barycentric coordinates (L0,L1,L2) of a point r[] in a triangle v[] are given by

	(v[0][0]-v[2][0]  v[1][0]-v[2][0]) (L0)   (r[0]-v[2][0])
	(                                )*(  ) = (            ),
	(v[0][1]-v[2][1]  v[1][1]-v[2][1]) (L1)   (r[1]-v[2][1])

with L2 = 1 - L0 - L1.
  */
  mat[0][0] = vertices[0][0] - vertices[2][0];
  mat[0][1] = vertices[1][0] - vertices[2][0];
  mat[1][0] = vertices[0][1] - vertices[2][1];
  mat[1][1] = vertices[1][1] - vertices[2][1];
  vec[0] = x - vertices[2][0];
  vec[1] = y - vertices[2][1];
  det = mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0];
  /*** We're assuming that the triangle is not pathological, i.e that det!=0. */
  barys[0] = ( mat[1][1]*vec[0] - mat[0][1]*vec[1])/det;
  barys[1] = (-mat[1][0]*vec[0] + mat[0][0]*vec[1])/det;
  barys[2] = 1.0 - barys[0] - barys[1];
}

