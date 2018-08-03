/*
 *  raytrace.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
TODO:
  - In raytrace(), look at rearranging the code to do the qhull step before choosing the rays. This would allow cells with all vertices outside the image boundaries to be excluded. If the image is much smaller than the model, this could lead to significant savings in time. The only downside might be memory useage...
  - In raytrace(), for par->traceRayAlgorithm==1, in theory we could deduce the cell geometry via the grid-point neighbour linkage, without needing to call delaunay() again.
 */

#include "lime.h"
#include "raythrucells.h"

typedef struct {
  double x,y, *intensity, *tau;
  unsigned int ppi;
  _Bool isInsideImage;
} rayData;

struct baryVelBuffType {
  int numVertices,numEdges,(*edgeVertexIndices)[2];
  double **vertexVels,**edgeVels,*entryCellBary,*midCellBary,*exitCellBary,*shapeFns;
};

typedef struct{
  double x[DIM], xCmpntRay, B[3];
  struct populations *mol;
  struct continuumLine cont;
} gridInterp;

/*....................................................................*/
void calcGridContDustOpacity(configInfo *par, const double freq\
  , double *lamtab, double *kaptab, const int nEntries, struct grid *gp){

  int id;
  double gtd;
  gsl_spline *spline = NULL;
  gsl_interp_accel *acc = NULL;
  double *kappatab = NULL;
  double *knus=NULL, *dusts=NULL;
  double *freqs=NULL;

  kappatab = malloc(sizeof(*kappatab)*1);
  knus     = malloc(sizeof(*knus)    *1);
  dusts    = malloc(sizeof(*dusts)   *1);
  freqs    = malloc(sizeof(*freqs)   *1);

  freqs[0] = freq;

  if(par->dust == NULL)
    kappatab[0] = 0.;
  else{
    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_cspline,nEntries);
    gsl_spline_init(spline,lamtab,kaptab,nEntries);
    kappatab[0] = interpolateKappa(freq, lamtab, kaptab, nEntries, spline, acc);
  }

  for(id=0;id<par->ncell;id++){
    gasIIdust(gp[id].x[0],gp[id].x[1],gp[id].x[2],&gtd);
    calcDustData(par, gp[id].dens, freqs, gtd, kappatab, 1, gp[id].t, knus, dusts); /* in aux.c. */
    gp[id].cont.knu = knus[0];
    gp[id].cont.dust = dusts[0];
  }

  if(par->dust != NULL){
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
  }
  free(knus);
  free(dusts);
  free(freqs);
  free(kappatab);
}

/*....................................................................*/
void
calcLineAmpSample(const double x[3], const double dx[3], const double ds\
  , const double binv, double *projVels, const int nSteps\
  , const double oneOnNSteps, const double deltav, double *vfac){
  /*
The bulk velocity of the model material can vary significantly with position, thus so can the value of the line-shape function at a given frequency and direction. The present function calculates 'vfac', an approximate average of the line-shape function along a path of length ds in the direction of the line of sight.

Note that this is called from within the multi-threaded block.
  */
  int i;
  double v,val;

  *vfac=0.;
  for(i=0;i<nSteps;i++){
    v = deltav - projVels[i]; /* projVels contains the component of the local bulk velocity in the direction of dx, whereas deltav is the recession velocity of the channel we are interested in (corrected for bulk source velocity and line displacement from the nominal frequency). Remember also that, since dx points away from the observer, positive values of the projected velocity also represent recessions. Line centre occurs when v==0, i.e. when deltav==dotProduct3D(dx,vel). That is the reason for the subtraction here. */
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
void
calcLineAmpInterp(const double projVelRay, const double binv\
  , const double deltav, double *vfac){
  /*
The bulk velocity of the model material can vary significantly with position, thus so can the value of the line-shape function at a given frequency and direction. The present function calculates 'vfac', an approximate average of the line-shape function along a path of length ds in the direction of the line of sight.

Note that this is called from within the multi-threaded block.
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
calcLineAmpErf(const double projVelOld, const double projVelNew, const double binv\
  , const double deltav, const double segmentLen, double *vfac){
  /*
The bulk velocity of the model material can vary significantly with position, thus so can the value of the line-shape function at a given frequency and direction. The present function calculates 'vfac', an approximate average of the line-shape function along a path of length ds in the direction of the line of sight.

Note that this is called from within the multi-threaded block.
  */
  double vbOld,vbNew;

  vbOld = binv*(deltav - projVelOld); /* projVelOld etc are components of the local bulk velocity in the direction of the ray, whereas deltav is the recession velocity of the channel we are interested in (corrected for bulk source velocity and line displacement from the nominal frequency). Remember also that, since the ray points away from the observer, positive values of the projected velocity also represent recessions. Line centre occurs when v==0, i.e. when deltav==projVelRay. That is the reason for the subtraction here. */
  vbNew = binv*(deltav - projVelNew);

  if (fabs(vbNew-vbOld)>(2.0*BIN_WIDTH))
    *vfac = geterf(vbOld,vbNew);
  else
    *vfac = gaussline(0.5*(vbOld+vbNew),1.0);
}

/*....................................................................*/
void
line_plane_intersect(struct grid *gp, double *ds, int posn, int *nposn\
  , double *dx, double *x, double cutoff){
  /*
This function returns ds as the (always positive-valued) distance between the present value of x and the next Voronoi face in the direction of vector dx, and nposn as the id of the grid cell that abuts that face. 

Note that this is called from within the multi-threaded block.
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
traceray(rayData ray, const int im\
  , configInfo *par, struct grid *gp, molData *md, imageInfo *img\
  , const double cutoff, const int nSteps, const double oneOnNSteps){
  /*
For a given image pixel position, this function evaluates the intensity of the total light emitted/absorbed along that line of sight through the (possibly rotated) model. The calculation is performed for several frequencies, one per channel of the output image.

Note that the algorithm employed here is similar to that employed in the function calculateJBar() which calculates the average radiant flux impinging on a grid cell: namely the notional photon is started at the side of the model near the observer and 'propagated' in the receding direction until it 'reaches' the far side. This is rather non-physical in conception but it makes the calculation easier.

Note that this is called from within the multi-threaded block.

Reads gp attributes x, dir, numNeigh, neigh, cont
if(par->polarization): B
if(img[im].doline): mol[molI].binv, mol[molI].specNumDens
if(!if(par->useVelFuncInRaytrace)): vel
  */
  int ichan,stokesId,di,i,posn,nposn,molI,lineI;
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
    line_plane_intersect(gp,&ds,posn,&nposn,dx,x,cutoff); /* Reads gp attributes numNeigh, x, dir, neigh. Returns a new ds equal to the distance to the next Voronoi face, and nposn, the ID of the grid cell that abuts that face. */

    if(par->polarization){ /* Should also imply img[im].doline==0. */
      sourceFunc_pol(gp[posn].B, gp[posn].cont, img[im].rotMat, snu_pol, &alpha);
      dtau=alpha*ds;
      calcSourceFn(dtau, par, &remnantSnu, &expDTau);
      remnantSnu *= ds;

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
      if(img[im].doline && par->useVelFuncInRaytrace){
        for(i=0;i<nSteps;i++){
          d = i*ds*oneOnNSteps;
          velocity(x[0]+(dx[0]*d),x[1]+(dx[1]*d),x[2]+(dx[2]*d),vel);
          projVels[i] = dotProduct3D(dx,vel);
        }
      }

      /* Calculate first the continuum stuff because it is the same for all channels:
      */
      contJnu = 0.0;
      contAlpha = 0.0;
      sourceFunc_cont(gp[posn].cont, &contJnu, &contAlpha);

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
                if(par->useVelFuncInRaytrace) /* because only in this case do we have projVels. */
                  calcLineAmpSample(x,dx,ds,gp[posn].mol[molI].binv,projVels,nSteps,oneOnNSteps,deltav,&vfac);
                else
                  vfac = gaussline(deltav-dotProduct3D(dx,gp[posn].vel),gp[posn].mol[molI].binv);

                /* Increment jnu and alpha for this Voronoi cell by the amounts appropriate to the spectral line. */
                sourceFunc_line(&md[molI],vfac,&(gp[posn].mol[molI]),lineI,&jnu,&alpha);
              }
            }
          }
        } /* end if(img[im].doline) */

        dtau=alpha*ds;
        /* Should we check for overly strong masers as in calculateJBar()?
        if(dtau < -30) dtau = -30;  
        */
        calcSourceFn(dtau, par, &remnantSnu, &expDTau);
        remnantSnu *= jnu*ds;
#ifdef FASTEXP
        brightnessIncrement = FastExp(ray.tau[ichan])*remnantSnu;
#else
        brightnessIncrement =    exp(-ray.tau[ichan])*remnantSnu;
#endif
        ray.intensity[ichan] += brightnessIncrement;
        ray.tau[ichan]+=dtau;
      } /* end loop over channels */
    } /* end if(par->polarization) */

    /* Move the working point to the edge of the next Voronoi cell. */
    for(di=0;di<DIM;di++) x[di]+=ds*dx[di];
    col+=ds;
    posn=nposn;
  } while(col < 2.0*fabs(zp));
}

/*....................................................................*/
void
doBaryInterp(const intersectType intcpt, struct grid *gp\
  , double xCmpntsRay[3], unsigned long gis[3]\
  , molData *md, const int numMols, gridInterp *gip){
  /*
The present routine takes (i) N values V_i at the vertices of a simplex, and (ii) the barycentric coordinates of a point, and returns a linear interpolation of the vertex values for the point location. This is essentially the same technique as the use of linear shape functions in Finite Element analysis. The idea is that if you define N shape functions Q_i, the ith shape function having the property that it is zero-valued at each of the vertices except the ith, and has value unity there, then the interpolation value at point r_ is given by

	       _N
	       \
	f(r_) =  >    V_i*Q_i(r_).
	       /_i=1

For a linear interpolation, each shape function Q_i(r_) is in fact just equal to the ith barycentric coordinate B_i of r_.

In the present case, N==3, the simplex is a triangular face of a Delaunay cell, and the point at which we desire the interpolated value is the intersection of a ray with that face. Several grid quantities of interest are interpolated.

For a readable definition of barycentric coordinates, see the wikipedia article of that name.

Note that this is called from within the multi-threaded block.
  */

  int di,molI,levelI;

  (*gip).xCmpntRay = intcpt.bary[0]*xCmpntsRay[0]\
                   + intcpt.bary[1]*xCmpntsRay[1]\
                   + intcpt.bary[2]*xCmpntsRay[2];
  for(di=0;di<DIM;di++){
    (*gip).x[di] = intcpt.bary[0]*gp[gis[0]].x[di]\
                 + intcpt.bary[1]*gp[gis[1]].x[di]\
                 + intcpt.bary[2]*gp[gis[2]].x[di];
  }

  for(di=0;di<3;di++){ /* 3 not DIM, because a B field only makes sense in 3 dimensions. */
/* ****** Maybe some test of DIM==3?? Would need to be systematic across the entire code..? */
    (*gip).B[di] = intcpt.bary[0]*gp[gis[0]].B[di]\
                 + intcpt.bary[1]*gp[gis[1]].B[di]\
                 + intcpt.bary[2]*gp[gis[2]].B[di];
  }

  for(molI=0;molI<numMols;molI++){
    (*gip).mol[molI].binv = intcpt.bary[0]*gp[gis[0]].mol[molI].binv\
                          + intcpt.bary[1]*gp[gis[1]].mol[molI].binv\
                          + intcpt.bary[2]*gp[gis[2]].mol[molI].binv;

    for(levelI=0;levelI<md[molI].nlev;levelI++){
      (*gip).mol[molI].specNumDens[levelI]\
        = intcpt.bary[0]*gp[gis[0]].mol[molI].specNumDens[levelI]\
        + intcpt.bary[1]*gp[gis[1]].mol[molI].specNumDens[levelI]\
        + intcpt.bary[2]*gp[gis[2]].mol[molI].specNumDens[levelI];
    }
  }

  (*gip).cont.dust = intcpt.bary[0]*gp[gis[0]].cont.dust\
                   + intcpt.bary[1]*gp[gis[1]].cont.dust\
                   + intcpt.bary[2]*gp[gis[2]].cont.dust;

  (*gip).cont.knu  = intcpt.bary[0]*gp[gis[0]].cont.knu\
                   + intcpt.bary[1]*gp[gis[1]].cont.knu\
                   + intcpt.bary[2]*gp[gis[2]].cont.knu;
}

/*....................................................................*/
void
doSegmentInterp(gridInterp gips[3], const int iA, molData *md\
  , const int numMols, const double oneOnNumSegments, const int si){
  /*
Note that this is called from within the multi-threaded block.
  */

  const double fracA = (si + 0.5)*oneOnNumSegments, fracB = 1.0 - fracA;
  const int iB = 1 - iA;
  int di,molI,levelI;

  gips[2].xCmpntRay = fracA*gips[iB].xCmpntRay + fracB*gips[iA].xCmpntRay; /* This does not seem to be used. */

  for(di=0;di<DIM;di++){
    gips[2].x[di] = fracA*gips[iB].x[di] + fracB*gips[iA].x[di];
  }
  for(di=0;di<3;di++){ /* 3 not DIM, because a B field only makes sense in 3 dimensions. */
    gips[2].B[di] = fracA*gips[iB].B[di] + fracB*gips[iA].B[di];
  }

  for(molI=0;molI<numMols;molI++){
    gips[2].mol[molI].binv = fracA*gips[iB].mol[molI].binv + fracB*gips[iA].mol[molI].binv;

    for(levelI=0;levelI<md[molI].nlev;levelI++){
      gips[2].mol[molI].specNumDens[levelI] = fracA*gips[iB].mol[molI].specNumDens[levelI]\
                                            + fracB*gips[iA].mol[molI].specNumDens[levelI];
    }
  }

  gips[2].cont.dust = fracA*gips[iB].cont.dust + fracB*gips[iA].cont.dust;
  gips[2].cont.knu  = fracA*gips[iB].cont.knu  + fracB*gips[iA].cont.knu;
}

/*....................................................................*/
void
calc2ndOrderShapeFunctions(struct baryVelBuffType *ptrToBuff, const int rayNum){
  /*
Note that this is called from within the multi-threaded block.
  */
  int vi,ei;
  char message[STR_LEN_0];
  double *barys=NULL;

  if(rayNum==0)
    barys = (*ptrToBuff).entryCellBary;
  else if(rayNum==1)
    barys = (*ptrToBuff).exitCellBary;
  else if(rayNum==2)
    barys = (*ptrToBuff).midCellBary;
  else{
    if(!silent){
      sprintf(message, "Bad ray number %d", rayNum);
      bail_out(message);
    }
    exit(1);
  }

  for(vi=0;vi<(*ptrToBuff).numVertices;vi++)
    (*ptrToBuff).shapeFns[vi] = barys[vi]*(2.0*barys[vi] - 1.0);

  for(ei=0;ei<(*ptrToBuff).numEdges;ei++){
    (*ptrToBuff).shapeFns[vi] = 4.0*barys[(*ptrToBuff).edgeVertexIndices[ei][0]]\
                                   *barys[(*ptrToBuff).edgeVertexIndices[ei][1]];

    vi++;
  }
}

/*....................................................................*/
void
doBaryInterpVel(const int numDims, struct baryVelBuffType *ptrToBuff\
, double vels[numDims]){
  /*
Note that this is called from within the multi-threaded block.
  */
  int di,vi,ei;

  for(di=0;di<numDims;di++)
    vels[di] = 0.0;

  for(vi=0;vi<(*ptrToBuff).numVertices;vi++)
    for(di=0;di<numDims;di++)
      vels[di] += (*ptrToBuff).vertexVels[vi][di]*(*ptrToBuff).shapeFns[vi];

  for(ei=0;ei<(*ptrToBuff).numEdges;ei++){
    for(di=0;di<numDims;di++)
      vels[di] += (*ptrToBuff).edgeVels[ei][di]*(*ptrToBuff).shapeFns[vi];
    vi++;
  }
}

/*....................................................................*/
void
doBaryInterpsVel(const int numDims, struct baryVelBuffType *ptrToBuff\
  , _Bool doRay[3], double rayVels[3][numDims]){
  /*
In this function we take (vector) velocities at 4 locations describing a tetrahedral Delaunay cell, plus 6 additional velocity values at the half-way points along each of the 6 edges of the cell, plus entry and exit (barycentric) locations of a ray passing through the cell, and return an interpolated value of the velocities at the entry and exit points, as well as at a point half-way between them.

This is similar to doBaryInterp() except that a quadratic interpolation is done rather than a linear one. For this we need values not just at the N vertices but also half-way along the N(N-1)/2 edges. The shape function Q_i which has value unity at the ith vertex and zero at all the other vertices, and also at the half-edge points, is given by

	Q_i(r_) = B_i(2*B_i - 1).

For the point ij half-way between vertices i and j the appropriate shape function is

	Q_ij(r_) = 4*B_i*B_j.

The interpolated value at r_ is, as before, the sum of the values at the vertices and the half-edge points, each multiplied by its associated shape function.

In the present case, N==4, so there are 6 half-edge points.

The remaining difference between the present function and doBaryInterp() is that here we only interpolate velocity.

Note that this is called from within the multi-threaded block.
  */

  int i;

  for(i=0;i<3;i++){
    if(doRay[i]){
      calc2ndOrderShapeFunctions(ptrToBuff, i);
      doBaryInterpVel(numDims, ptrToBuff, rayVels[i]);
    }
  }
}

/*....................................................................*/
void
doSegmentInterpVector(const int numDims, const double rayVels[3][numDims]\
  , const double x, double vel[numDims]){
  /*
This function is supplied with values of velocity at the beginning, end and midpoint of the line segment which represents the passage of a ray through a cell (stored in that order in the input array 'rayVels'); what it does is perform a 2nd-order (i.e. quadratic) interpolation to obtain an estimate of the velocity at the given displacement along the line segment (returned in the array 'vel').

Given y0, y1 and y2 at x0, x1 and x2, the Lagrange interpolating polynomial is

	         x-x1    x-x2         x-x0    x-x2         x-x0    x-x1
	y ~ y0*-------*------- + y1*-------*------- + y2*-------*-------.
	        x0-x1   x0-x2        x1-x0   x1-x2        x2-x0   x2-x1

If we calculate x as a fractional value along the line segment then, for the y values we have in hand, this reduces to

	y ~ y0*(x-1)*(2x-1) + y1*x*(2x-1) + y2*4x*(1-x).

***** Appears to be unused.
  */

  double shapeFns[3];
  int di;

  shapeFns[0] = (x - 1.0)*(2.0*x - 1.0);
  shapeFns[1] = x*(2.0*x - 1.0);
  shapeFns[2] = 4.0*x*(1.0 - x);

  for(di=0;di<numDims;di++)
    vel[di] = rayVels[0][di]*shapeFns[0] + rayVels[1][di]*shapeFns[1] + rayVels[2][di]*shapeFns[2];
}

/*....................................................................*/
double doSegmentInterpScalar(const double ys[3], const double x){
  /*
This function is supplied with values of y at the beginning, end and midpoint of the line segment; what it does is perform a 2nd-order (i.e. quadratic) interpolation to obtain an estimate of y at the given fractional displacement along the line segment.

Given y0, y1 and y2 at x0, x1 and x2, the Lagrange interpolating polynomial is

	         x-x1    x-x2         x-x0    x-x2         x-x0    x-x1
	y ~ y0*-------*------- + y1*-------*------- + y2*-------*-------.
	        x0-x1   x0-x2        x1-x0   x1-x2        x2-x0   x2-x1

If x is the fractional distance along the line segment, then for the y values we have in hand, this reduces to

	y ~ y0*(x-1)*(2x-1) + y1*x*(2x-1) + y2*4x*(1-x).

Note that this is called from within the multi-threaded block.
  */

  double shapeFns[3],y;

  shapeFns[0] = (x - 1.0)*(2.0*x - 1.0);
  shapeFns[1] = x*(2.0*x - 1.0);
  shapeFns[2] = 4.0*x*(1.0 - x);

  y = ys[0]*shapeFns[0] + ys[1]*shapeFns[1] + ys[2]*shapeFns[2];

  return y;
}

/*....................................................................*/
void
traceray_smooth(rayData ray, const int im\
  , configInfo *par, struct grid *gp, double *vertexCoords, molData *md\
  , imageInfo *img, struct simplex *dc, const unsigned long numCells\
  , const double epsilon, gridInterp gips[3], struct baryVelBuffType *ptrToBuff\
  , const int numSegments, const double oneOnNumSegments){
  /*
For a given image pixel position, this function evaluates the intensity of the total light emitted/absorbed along that line of sight through the (possibly rotated) model. The calculation is performed for several frequencies, one per channel of the output image.

Note that the algorithm employed here to solve the RTE is similar to that employed in the function calculateJBar() which calculates the average radiant flux impinging on a grid cell: namely the notional photon is started at the side of the model near the observer and 'propagated' in the receding direction until it 'reaches' the far side. This is rather non-physical in conception but it makes the calculation easier.

This version of traceray implements a new algorithm in which the population values are interpolated linearly from those at the vertices of the Delaunay cell which the working point falls within.

A note about the object 'gips': this is an array with 3 elements, each one a struct of type 'gridInterp'. This struct is meant to store as many of the grid-point quantities (interpolated from the appropriate values at actual grid locations) as are necessary for solving the radiative transfer equations along the ray path. The first 2 entries give the values for the entry and exit points to a Delaunay cell, but which is which can change, and is indicated via the variables entryI and exitI (this is a convenience to avoid copying the values, since the values for the exit point of one cell are obviously just those for entry point of the next). The third entry stores values interpolated along the ray path within a cell.

Note that this is called from within the multi-threaded block.
  */
  const int numFaces = DIM+1,nVertPerFace=3,numRayInterpSamp=3;
  int ichan,stokesId,di,status,lenChainPtrs=0,entryI,exitI,vi,vvi,ci,ei,fi;
  int si,molI,lineI,k,i;
  double xp,yp,zp,x[DIM],dir[DIM],projVelRay=0.0,vel[DIM],projVelOffset=0.0,projVel2ndDeriv;
  double xCmpntsRay[nVertPerFace],ds,snu_pol[3],dtau,contJnu,contAlpha;
  double jnu,alpha,lineRedShift,vThisChan,deltav,vfac,remnantSnu,expDTau;
  double brightnessIncrement,projVelOld=0.0,projVelNew=0.0;
  intersectType entryIntcptFirstCell, *cellExitIntcpts=NULL;
  unsigned long *chainOfCellIds=NULL,dci,dci0,dci1;
  unsigned long gis[2][nVertPerFace],gi,gi0,gi1,trialGi;
  double rayVels[numRayInterpSamp][DIM],projRayVels[numRayInterpSamp];
  _Bool doRay[numRayInterpSamp],matchFound,neighNotFound;
  struct interCellKeyType{
    int exitedFaceIs[nVertPerFace],fiEnteredCell;
  } *interCellKey=NULL;

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
  status = followRayThroughCells(DIM, x, dir, vertexCoords, dc, numCells, epsilon\
    , NULL, &entryIntcptFirstCell, &chainOfCellIds, &cellExitIntcpts, &lenChainPtrs);

  if(status!=0){
    free(chainOfCellIds);
    free(cellExitIntcpts);
    return;
  }

  if(img[im].doline && img[im].doInterpolateVels){
    /*
There is a problem when we want to copy cell-centric barycentric coords (BCs) for the entry face of all but the first cell in the chain when all we have to work with is the exit intercept (which includes face-centric BCs). This occurs because the order of the BCs in the exit cell corresponds to the order of the vertices in the cell it exits from, but we need the order for the cell entered. Thus we construct here an array of length lenChainPtrs-1 which gives a key to the exited-face vertices from the entered-face ones, and also stores the opposite-face vertex of the entered cell.
    */
    interCellKey = malloc(sizeof(*interCellKey)*(lenChainPtrs-1));
    dci1 = chainOfCellIds[0];
    for(ci=0;ci<lenChainPtrs-1;ci++){
      dci0 = dci1;
      dci1 = chainOfCellIds[ci+1];

      /* Obtain the indices of the grid points on the vertices of the exited face.
      */
      vvi = 0;
      for(fi=0;fi<numFaces;fi++){
        if(fi!=cellExitIntcpts[ci].fi){
          gis[0][vvi++] = dc[dci0].vertx[fi];
        }
      }

      /* Now we run through all vertices of the entered cell and attempt to match the gis.
      */
      interCellKey[ci].fiEnteredCell = -1; /* This flags that no opposite vertex was found - indicates a bug somewhere. */
      vvi = 0;
      for(vi=0;vi<(*ptrToBuff).numVertices;vi++){
        trialGi = dc[dci1].vertx[vi];
        matchFound = 0; /* default. */
        for(fi=0;fi<nVertPerFace;fi++){
          if(trialGi==gis[0][fi]){
            matchFound = 1;
        break;
          }
        }

        if(matchFound){
          interCellKey[ci].exitedFaceIs[vvi] = fi;
          vvi++;
        }else{
          if(interCellKey[ci].fiEnteredCell!=-1){
            if(!silent) bail_out("More than one opposite vertex found! This is some sort of bug.");
            exit(1);
          }
          interCellKey[ci].fiEnteredCell = vi;
        }
      }
    }
  } /* end if(img[im].doline) */

  entryI = 0;
  exitI  = 1;
  dci = chainOfCellIds[0];

  /* Obtain the indices of the grid points on the vertices of the entry face.
  */
  vvi = 0;
  for(fi=0;fi<numFaces;fi++){
    if(fi!=entryIntcptFirstCell.fi){
      gis[entryI][vvi++] = dc[dci].vertx[fi];
    }
  }

  /* Calculate, for each of the 3 vertices of the entry face, the displacement components in the direction of 'dir'. *** NOTE *** that if all the rays are parallel, we could precalculate these for all the vertices.
  */
  for(vi=0;vi<nVertPerFace;vi++)
    xCmpntsRay[vi] = dotProduct3D(dir, gp[gis[entryI][vi]].x);

  /* Calculate the values (via linear interpolation) of necessary grid quantities for the entry point to the first cell:
  */
  doBaryInterp(entryIntcptFirstCell, gp, xCmpntsRay, gis[entryI]\
    , md, par->nSpecies, &gips[entryI]);

  for(ci=0;ci<lenChainPtrs;ci++){
    /* For each cell we have 2 data structures which give information about respectively the entry and exit points of the ray, including the barycentric coordinates of the intersection point between the ray and the appropriate face of the cell. (If we follow rays in 3D space then the cells will be tetrahedra and the faces triangles.) If we know the value of a quantity Q for each of the vertices, then the linear interpolation of the Q values for any face is (for a 3D space) bary[0]*Q[0] + bary[1]*Q[1] + bary[2]*Q[2], where the indices are taken to run over the vertices of that face. Thus we can calculate the interpolated values Q_entry and Q_exit. Further linear interpolation along the path between entry and exit is straightforward.
    */

    dci = chainOfCellIds[ci];

    /* Obtain the indices of the grid points on the vertices of the exit face. */
    vvi = 0;
    for(fi=0;fi<numFaces;fi++){
      if(fi!=cellExitIntcpts[ci].fi){
        gis[exitI][vvi++] = dc[dci].vertx[fi];
      }
    }

    /* Calculate, for each of the 3 vertices of the exit face, the displacement components in the direction of 'dir'. *** NOTE *** that if all the rays are parallel, we could precalculate these for all the vertices.
    */
    for(vi=0;vi<nVertPerFace;vi++)
      xCmpntsRay[vi] = dotProduct3D(dir, gp[gis[exitI][vi]].x);

    /* Calculate the values (via linear interpolation) of necessary grid quantities for the exit point to the cell:
    */
    doBaryInterp(cellExitIntcpts[ci], gp, xCmpntsRay, gis[exitI]\
      , md, par->nSpecies, &gips[exitI]);

    if(img[im].doline && img[im].doInterpolateVels){
      /*
Calculate the values (via 2nd-order interpolation) of the systemic velocity for three points: the entry point to the cell, the exit point, and the point half-way between. (There is a fair bit of setting up to do before the actual call to doBaryInterpsVel).

First lot of setups: copying over the 4 vertex and the 6 mid-edge velocities.
      */
      for(vi=0;vi<(*ptrToBuff).numVertices;vi++){
        gi = dc[dci].vertx[vi];
        for(di=0;di<DIM;di++)
          (*ptrToBuff).vertexVels[vi][di] = gp[gi].vel[di];
      }
      for(ei=0;ei<(*ptrToBuff).numEdges;ei++){
        gi0 = dc[dci].vertx[(*ptrToBuff).edgeVertexIndices[ei][0]];
        gi1 = dc[dci].vertx[(*ptrToBuff).edgeVertexIndices[ei][1]];
        /* Find index of gi0 neighbours which points to gi1: */
        neighNotFound = 1; /* default */
        for(k=0;k<gp[gi0].numNeigh;k++){
          if(gp[gi0].neigh[k]->id==(int)gi1){ /* **** eventually change id to unsigned long! *** */
            neighNotFound = 0;
            break;
          }
        }
        if(neighNotFound){
          if(!silent) bail_out("Neighbour not found.");
          exit(1);
        }
        for(di=0;di<DIM;di++)
          (*ptrToBuff).edgeVels[ei][di] = gp[gi0].v2[3*k+di]; /* v2 is currently the mid-edge velocity. This is not very robust: should change it to allow a variable (but odd) number of velocity samples per edge, then pick the N_VEL_SEG_PER_HALFth. */
      }

      /*
More setups. Now we want to populate the barycentric coords of the entry and exit points to the cell. We need 4 values for each, because what we want is cell-specific (CS) BC. What we have calculated already are face-specific (FS) BC which omit the zero-valued BC for the vertex opposite the intersected face. To copy the 3 FS values to the array of 4 CS, we keep the order (since the order of the FS ones is the same as the order of vertices in dc[dci].vertx, which is what we want, just with the vertex opposite the intersected face omitted), just inserting 0 for the vertex opposite the intersected face.

However! While this works fine all the time for the exit intercepts, and for cell ci==0 for the entry intercept, for cells ci>0 we need the FS entry BC for the cell, which to be sure are the same values as the FS exit BC of the previous cell; but the vertex indices are all screwed up, because what we have are indices appropriate to cell ci-1 but we need them for cell ci. To remedy all this we have pre-calculated a key to the old cell indices from the new.
      */
      if(ci==0){
        vvi = 0;
        for(vi=0;vi<(*ptrToBuff).numVertices;vi++){
          if(vi==entryIntcptFirstCell.fi){
            (*ptrToBuff).entryCellBary[vi] = 0.0;
          }else{
            (*ptrToBuff).entryCellBary[vi] = entryIntcptFirstCell.bary[vvi];
            vvi++;
          }
        }

        doRay[0] = 1;
        doRay[1] = 1;
        doRay[2] = 1;

      }else{
        vvi = 0;
        for(vi=0;vi<(*ptrToBuff).numVertices;vi++){
          if(vi==interCellKey[ci-1].fiEnteredCell){
            (*ptrToBuff).entryCellBary[vi] = 0.0;
          }else{
            (*ptrToBuff).entryCellBary[vi] = cellExitIntcpts[ci-1].bary[interCellKey[ci-1].exitedFaceIs[vvi]];
            vvi++;
          }
        }

        doRay[0] = 0;
        doRay[1] = 1;
        doRay[2] = 1;

        for(di=0;di<DIM;di++)
          rayVels[0][di] = rayVels[1][di];
      }

      vvi = 0;
      for(vi=0;vi<(*ptrToBuff).numVertices;vi++){
        if(vi==cellExitIntcpts[ci].fi){
          (*ptrToBuff).exitCellBary[vi] = 0.0;
        }else{
          (*ptrToBuff).exitCellBary[vi] = cellExitIntcpts[ci].bary[vvi];
          vvi++;
        }

        (*ptrToBuff).midCellBary[vi] = 0.5*((*ptrToBuff).entryCellBary[vi] + (*ptrToBuff).exitCellBary[vi]);
      }

      doBaryInterpsVel(DIM, ptrToBuff, doRay, rayVels);

      for(i=0;i<numRayInterpSamp;i++)
        projRayVels[i] = dotProduct3D(dir, rayVels[i]);

      /*
We're going to calculate the line amplitude increment per segment for the 2nd-order interpolation scheme as follows. We have 3 values of projRayVels across the path the ray takes through the cell: one at cell entry, one at cell exit and the third at the half-way point. This quantity is the scalar projected value of velocity in the ray direction. Since these values were obtained via a 2nd-order interpolation within the cell, we may use them to derive 2nd-order interpolated values along the ray path, and be certain that these values are the same as we would have obtained if we had done the full 2nd-order interpolation within the cell, using the vertex and edge values etc. A quadratic is a quadratic, right? So these 3 projRayVels values are kind of shorthand for the full cell data. Anyway, we are going to break this ray path up into equal-width segments. Within each segment, we will approximate variation of the projected velocity by a linear function of distance, which allows us to express the line amplitude integral within this segment by an error function. But which linear function? The slope is easy, it will be the same as the 1st derivative of the vel quadratic across the whole in-cell path. In fact for the erf lookup we just need the start and end values of proj vel. We can (and do) get projected velocities at the beginning and end of each segment. The best linear function though is offset vertically, such that the integrated square of the difference between the 'true' 2nd-order function and this 1st-order function that we assume so we can use erfs is as small as possible. This occurs when the start and end proj vel values at the start and end of the segment are offset by -y"*x^2/6, where y" is the 2nd derivative and x is the segment length fraction.
      */

      projVelOld = doSegmentInterpScalar(projRayVels, 0);
      projVel2ndDeriv = (projRayVels[0] + projRayVels[1] - 2.0*projRayVels[2])*4.0; /* The times 4 is actually a divide by deltaX^2, because in this case deltaX is nominally 0.5, i.e. half-way across the path through the cell. */
      projVelOffset = -projVel2ndDeriv*oneOnNumSegments*oneOnNumSegments/6.0;
    } /* end if(img[im].doline) && img[im].doInterpolateVels */

    /* At this point we have interpolated all the values of interest to both the entry and exit points of the cell. Now we break the path between entry and exit into several segments and calculate all these values at the midpoint of each segment.

At the moment I will fix the number of segments, but it might possibly be faster to rather have a fixed segment length (in barycentric coordinates) and vary the number depending on how many of these lengths fit in the path between entry and exit.
    */
    ds = (gips[exitI].xCmpntRay - gips[entryI].xCmpntRay)*oneOnNumSegments;

    for(si=0;si<numSegments;si++){
      doSegmentInterp(gips, entryI, md, par->nSpecies, oneOnNumSegments, si);

      if(par->polarization){ /* Should also imply img[im].doline==0. */
        sourceFunc_pol(gips[2].B, gips[2].cont, img[im].rotMat, snu_pol, &alpha);
        dtau=alpha*ds;
        calcSourceFn(dtau, par, &remnantSnu, &expDTau);
        remnantSnu *= ds;

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
        if(img[im].doline){
          if(img[im].doInterpolateVels){
            projVelNew = doSegmentInterpScalar(projRayVels, (si + 1.0)*oneOnNumSegments);
          }else{
            velocity(gips[2].x[0], gips[2].x[1], gips[2].x[2], vel);
            projVelRay = dotProduct3D(dir, vel);
          }
        }

        /* Calculate first the continuum stuff because it is the same for all channels:
        */
        contJnu = 0.0;
        contAlpha = 0.0;
        sourceFunc_cont(gips[2].cont, &contJnu, &contAlpha);

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

                  if(img[im].doInterpolateVels)
                    calcLineAmpErf(projVelOld, projVelNew, gips[2].mol[molI].binv, deltav-projVelOffset, oneOnNumSegments, &vfac);
                  else
                    calcLineAmpInterp(projVelRay, gips[2].mol[molI].binv, deltav, &vfac);

                  /* Increment jnu and alpha for this Voronoi cell by the amounts appropriate to the spectral line.
                  */
                  sourceFunc_line(&md[molI], vfac, &(gips[2].mol[molI]), lineI, &jnu, &alpha);
                } /* end if within freq range. */
              } /* end loop over lines this mol. */
            } /* end loop over all mols. */
          } /* end if doLine. */

          dtau = alpha*ds;
//???          if(dtau < -30) dtau = -30; // as in calculateJBar()?
          calcSourceFn(dtau, par, &remnantSnu, &expDTau);
          remnantSnu *= jnu*ds;
#ifdef FASTEXP
          brightnessIncrement = FastExp(ray.tau[ichan])*remnantSnu;
#else
          brightnessIncrement =    exp(-ray.tau[ichan])*remnantSnu;
#endif
          ray.intensity[ichan] += brightnessIncrement;
          ray.tau[ichan] += dtau;
        } /* End loop over channels. */

        if(img[im].doInterpolateVels) projVelOld = projVelNew;
      } /* End if(par->polarization). */
    } /* End loop over segments within cell. */

    entryI = exitI;
    exitI = 1 - exitI;
  } /* End loop over cells in the chain traversed by the ray. */

  if(img[im].doline && img[im].doInterpolateVels)
    free(interCellKey);

  free(chainOfCellIds);
  free(cellExitIntcpts);
}

/*....................................................................*/
_Bool
locateRayOnImage(double x[2], const double size, const double imgCentreXPixels\
  , const double imgCentreYPixels, imageInfo *img, const int im\
  , int *xi, int *yi, unsigned int *ppi){

  _Bool isInsideImage;

  /* Calculate which pixel the projected position (x[0],x[1]) falls within.
  */
  *xi = floor(x[0]/size + imgCentreXPixels);
  *yi = floor(x[1]/size + imgCentreYPixels);
  if(*xi<0 || *xi>=img[im].pxls || *yi<0 || *yi>=img[im].pxls){
    isInsideImage = 0;
    *ppi = 0; /* Under these circumstances it ought never to be accessed, but it is not good to leave it without a value at all. */
  }else{
    isInsideImage = 1;
    *ppi = (unsigned int)(*yi)*(unsigned int)img[im].pxls + (unsigned int)(*xi);
  }

  return isInsideImage;
}

/*....................................................................*/
void
assignRayOnImage(double x[2], const double size, const double imgCentreXPixels\
  , const double imgCentreYPixels, imageInfo *img, const int im\
  , const int maxNumRaysPerPixel, rayData *rays, int *numActiveRays){
  /*
The present function does several things, as follows:
	- Calculates the image position in pixel coordinates of the proposed ray position specified by x[].
	- If the proposed ray is inside the image bounds, and the count of rays for that image pixel does not exceed the maximum allowed, the function:
	  * adds 1 to the count of rays for that pixel of the image.
	  * increments *numActiveRays;
	  * stores information for the new ray in rays[*numActiveRays].

Returned information is thus:
	- An updated array img[im].pixel[ppi].numRays.
	- An updated list of accepted rays.
	- An updated value of *numActiveRays.
  */

  int xi,yi,ichan;
  _Bool isInsideImage;
  unsigned int ppi;

  isInsideImage = locateRayOnImage(x, size, imgCentreXPixels\
    , imgCentreYPixels, img, im, &xi, &yi, &ppi);

  /* See if we want to keep the ray. For the time being we will include those outside the image bounds, because we need to interpolate between outside-image rays and inside-image ones, but a cleverer algorithm would exclude all but those which are affected by this (i.e. that immediately abut the image bounds). Note that maxNumRaysPerPixel<1 is used to flag that there is no upper limit to the number of rays per pixel.
  */
  if(!isInsideImage || maxNumRaysPerPixel<1 || img[im].pixel[ppi].numRays<maxNumRaysPerPixel){
    if(isInsideImage)
      img[im].pixel[ppi].numRays++;

    rays[*numActiveRays].isInsideImage = isInsideImage;
    rays[*numActiveRays].ppi = ppi;
    rays[*numActiveRays].x = x[0];
    rays[*numActiveRays].y = x[1];
    rays[*numActiveRays].tau       = malloc(sizeof(double)*img[im].nchan);
    rays[*numActiveRays].intensity = malloc(sizeof(double)*img[im].nchan);
    for(ichan=0;ichan<img[im].nchan;ichan++) {
      rays[*numActiveRays].tau[ichan] = 0.0;
      rays[*numActiveRays].intensity[ichan] = 0.0;
    }

    (*numActiveRays)++;
  }
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

/*....................................................................*/
double *
extractGridXs(const unsigned short numDims, const unsigned long numPoints\
  , struct grid *gp){

  double *xValues=NULL;
  unsigned long i_ul;
  unsigned short i_us;

  xValues = malloc(sizeof(*xValues)*numDims*numPoints);
  for(i_ul=0;i_ul<numPoints;i_ul++){
    for(i_us=0;i_us<numDims;i_us++)
      xValues[numDims*i_ul+i_us] = gp[i_ul].x[i_us];
  }

  return xValues;
}

/*....................................................................*/
struct simplex *
convertCellType(const unsigned short numDims, const unsigned long numCells\
  , struct cell *dc, struct grid *gp){

  struct simplex *cells=NULL;
  unsigned long dci,gi;
  unsigned short vi,di;

  cells = malloc(sizeof(*cells)*numCells);
  for(dci=0;dci<numCells;dci++){
    cells[dci].id = dci;

    for(vi=0;vi<numDims+1;vi++)
      cells[dci].vertx[vi] = dc[dci].vertx[vi]->id;

    for(di=0;di<numDims;di++)
      cells[dci].centre[di] = 0.0;
    for(vi=0;vi<numDims+1;vi++){
      gi = cells[dci].vertx[vi];
      for(di=0;di<numDims;di++)
        cells[dci].centre[di] += gp[gi].x[di];
    }
    for(di=0;di<numDims;di++)
      cells[dci].centre[di] *= (1.0/(double)(numDims+1));
  }
  for(dci=0;dci<numCells;dci++){
    for(vi=0;vi<numDims+1;vi++){
      if(dc[dci].neigh[vi]==NULL)
        cells[dci].neigh[vi] = NULL;
      else
        cells[dci].neigh[vi] = &cells[dc[dci].neigh[vi]->id];
    }
  }

  return cells;
}

/*....................................................................*/
void
get2DCells(rayData *rays, const int numActiveRays, struct simplex **cells2D\
  , unsigned long *numCells){

  const int numDims=2,numFaces=numDims+1;
  const double oneOnNFaces=1.0/(double)numFaces;
  coordT *pt_array;
  int ri,i,di,vi;
  char flags[255];
  boolT ismalloc = False;
  int curlong, totlong;
  facetT *facet,*neighbor,**neighborp;
  vertexT *vertex,**vertexp;
  unsigned long fi,ffi,id,dci;
  _Bool neighbourNotFound;
  char message[STR_LEN_0];
  double sum;

  pt_array = malloc(sizeof(*pt_array)*numDims*numActiveRays);

  for(ri=0;ri<numActiveRays;ri++) {
    pt_array[ri*numDims+0] = rays[ri].x;
    pt_array[ri*numDims+1] = rays[ri].y;
  }

  sprintf(flags,"qhull d Qbb Qt");
  if(qh_new_qhull(numDims, numActiveRays, pt_array, ismalloc, flags, NULL, NULL)) {
    if(!silent) bail_out("Qhull failed to triangulate");
    exit(1);
  }

  (*numCells) = 0;
  FORALLfacets {
    if(!facet->upperdelaunay)
      (*numCells)++;
  }

  (*cells2D) = malloc(sizeof(**cells2D)*(*numCells));

  fi = 0;
  FORALLfacets {
    if (!facet->upperdelaunay) {
      (*cells2D)[fi].id = (unsigned long)facet->id; /* Do NOT expect this to be equal to fi. */
      fi++;
    }
  }

  fi = 0;
  FORALLfacets {
    if (!facet->upperdelaunay) {
      i = 0;
      FOREACHneighbor_(facet) {
        if(neighbor->upperdelaunay){
          (*cells2D)[fi].neigh[i] = NULL;
        }else{
          /* Have to find the member of *cells2D with the same id as neighbour.*/
          ffi = 0;
          neighbourNotFound=1;
          while(ffi<(*numCells) && neighbourNotFound){
            if((*cells2D)[ffi].id==(unsigned long)neighbor->id){
              (*cells2D)[fi].neigh[i] = &(*cells2D)[ffi];
              neighbourNotFound = 0;
            }
            ffi++;
          }

          if(ffi>=(*numCells) && neighbourNotFound){
            if(!silent){
              sprintf(message, "Something weird going on. Cannot find a cell with ID %lu", (unsigned long)(neighbor->id));
              bail_out(message);
            }
            exit(1);
          }
        }
        i++;
      }

      i = 0;
      FOREACHvertex_( facet->vertices ) {
        id = (unsigned long)qh_pointid(vertex->point);
        (*cells2D)[fi].vertx[i] = id;
        i++;
      }

      fi++;
    }
  }

  /* We need to process the list of cells a bit further - calculate their centres, and reset the id values to be the same as the index of the cell in the list. (This last because we are going to construct other lists to indicate which cells have been visited etc.)
  */
  for(dci=0;dci<(*numCells);dci++){
    for(di=0;di<numDims;di++){
      sum = 0.0;
      for(vi=0;vi<numFaces;vi++){
        id = (*cells2D)[dci].vertx[vi];
        sum += pt_array[id*numDims+di];
      }
      (*cells2D)[dci].centre[di] = sum*oneOnNFaces;
    }

    (*cells2D)[dci].id = dci;
  }

  qh_freeqhull(!qh_ALL);
  qh_memfreeshort (&curlong, &totlong);
  free(pt_array);
}

/*....................................................................*/
void
raytrace(int im, configInfo *par, struct grid *gp, molData *md\
  , imageInfo *img, double *lamtab, double *kaptab, const int nEntries){
  /*
This function constructs an image cube by following sets of rays (at least 1 per image pixel) through the model, solving the radiative transfer equations as appropriate for each ray. The ray locations within each pixel are chosen randomly within the pixel, but the number of rays per pixel is set equal to the number of projected model grid points falling within that pixel, down to a minimum equal to par->alias.

Note that the argument 'md', and the grid element '.mol', are only accessed for line images.
  */
  const int maxNumRaysPerPixel=20; /**** Arbitrary - could make this a global, or an argument. Set it to zero to indicate there is no maximum. */
  const double cutoff = par->minScale*1.0e-7;
  const int numFaces=1+DIM,numInterpPoints=3,numSegments=5,minNumRaysForAverage=2;
  const double oneOnNFaces=1.0/(double)numFaces, oneOnNumSegments = 1.0/(double)numSegments;
  const double epsilon = 1.0e-6; // Needs thinking about. Double precision is much smaller than this.
  const int nStepsThruCell=10;
  const double oneOnNSteps=1.0/(double)nStepsThruCell;

  double pixelSize,imgCentreXPixels,imgCentreYPixels,minfreq,absDeltaFreq,x,xs[2],sum,oneOnNumRays;
  unsigned int totalNumImagePixels,ppi,numPixelsForInterp;
  int ichan,numCircleRays,numActiveRaysInternal,numActiveRays,lastChan;
  int gi,molI,lineI,i,di,xi,yi,ri,vi,ei,i0,i1;
  int cmbMolI,cmbLineI;
  rayData *rays;
  struct cell *dc=NULL;
  struct simplex *cells=NULL;
  unsigned long numCells,dci,numPointsInAnnulus;
  double local_cmb,cmbFreq,circleSpacing,scale,angle,rSqu;
  double *vertexCoords=NULL;
  gsl_error_handler_t *defaultErrorHandler=NULL;
  struct baryVelBuffType velBuff,*ptrToBuff=NULL;
#ifndef NO_PROGBARS
  double progFraction,oneOnNumActiveRaysMinus1;
#endif

  pixelSize = img[im].distance*img[im].imgres;
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
We need to calculate or choose a single value of 'local' CMB flux, also single values (i.e. one of each per grid point) of dust and knu, all corresponding the the nominal image frequency. The sensible thing would seem to be to calculate them afresh for each new image; and for continuum images, this is what in fact has always been done. For line images however local_cmb and the dust/knu values were calculated for the frequency of each spectral line and stored respectively in the molData struct and the struct populations element of struct grid. These multiple values (of dust/knu at least) are required during the main solution kernel of LIME, so for line images at least they were kept until the present point, just so one from their number could be chosen. :-/

At the present point in the code, for line images, instead of calculating the 'continuum' values of local_cmb/dust/knu, the algorithm chose the nearest 'line' frequency and calculates the required numbers from that. The intent is to preserve (for the present at least) the former numerical behaviour, while changing the way the information is parcelled out among the variables and structs. I.e. a dedicated 'continuum' pair of dust/knu values is now available for each grid point in addition to the array of spectral line values. This decoupling allows better management of memory and avoids the deceptive use of spectral-line variables for continuum use.
  */
  if(img[im].doline){
    if (img[im].trans>=0){
      cmbMolI  = img[im].molI;
      cmbLineI = img[im].trans;

    }else{ /* User didn't set trans. Find the nearest line to the image frequency. */
      minfreq = fabs(img[im].freq - md[0].freq[0]);;
      cmbMolI = 0;
      cmbLineI = 0;
      for(molI=0;molI<par->nSpecies;molI++){
        for(lineI=0;lineI<md[molI].nline;lineI++){
          if((molI==0 && lineI==0)) continue;

          absDeltaFreq = fabs(img[im].freq - md[molI].freq[lineI]);
          if(absDeltaFreq < minfreq){
            minfreq = absDeltaFreq;
            cmbMolI = molI;
            cmbLineI = lineI;
          }
        }
      }
    }
    cmbFreq = md[cmbMolI].freq[cmbLineI];

  }else{ /* continuum image */
    cmbFreq = img[im].freq;
  }

  local_cmb = planckfunc(cmbFreq,LOCAL_CMB_TEMP);
  calcGridContDustOpacity(par, cmbFreq, lamtab, kaptab, nEntries, gp); /* Reads gp attributes x, dens, and t and writes attributes cont.dust and cont.knu. */

  for(ppi=0;ppi<totalNumImagePixels;ppi++){
    for(ichan=0;ichan<img[im].nchan;ichan++){
      img[im].pixel[ppi].intense[ichan] = 0.0;
      img[im].pixel[ppi].tau[    ichan] = 0.0;
    }
  }

  for(ppi=0;ppi<totalNumImagePixels;ppi++)
    img[im].pixel[ppi].numRays = 0;

  /*
The set of rays which we plan to follow have starting points which are defined by the set of (non-sink) grid points, projected onto a plane parallel to the observer's X and Y axes. In addition to these points we will add another tranche located on a circle in this plane with its centre at (X,Y) == (0,0) and radius equal to the model radius. Addition of these circle points seems to be necessary to make qhull behave properly. We choose evenly-spaced points on this circle and choose the spacing such that it is the same as the average nearest-neighbour spacing of the grid points in the outer 1/3 annulus of the circular projected model, assuming the grid points were evenly distributed within the annulus.

How to calculate this distance? Well if we have N points randomly but evenly distributed inside an annulus of radius (2/3 to 1)*R it is not hard to show that the mean NN spacing is G(3/2)*R/sqrt(9*N/5), where G() is the gamma function. In fact G(3/2)=sqrt(pi)/2. This will add about 12*sqrt(pi*N/5) points.
  */
  numPointsInAnnulus = 0;
  for(gi=0;gi<par->pIntensity;gi++){
    /* Note that we are *NOT* rotating the model before doing this projection. The reason is that we just want a rough estimate which avoids the centre, which usually has a concentration of points. */
    rSqu = gp[gi].x[0]*gp[gi].x[0] + gp[gi].x[1]*gp[gi].x[1];
    if(rSqu > (4.0/9.0)*par->radiusSqu) numPointsInAnnulus += 1;
  }
  if(numPointsInAnnulus>0){
    circleSpacing = (1.0/6.0)*par->radius*sqrt(5.0*M_PI/(double)numPointsInAnnulus);
    numCircleRays = (int)(2.0*M_PI*par->radius/circleSpacing);
  }else{
    numCircleRays = 0;
  }

  /* The following is the first of the 3 main loops in raytrace. Here we loop over the (internal or non-sink) grid points. We're doing 2 things: loading the rotated, projected coordinates into the rays list, and counting the rays per image pixel.
  */
  rays = malloc(sizeof(rayData)*(par->pIntensity+numCircleRays)); /* We may need to reallocate this later. */
  numActiveRaysInternal = 0;
  for(gi=0;gi<par->pIntensity;gi++){
    /* Apply the inverse (i.e. transpose) rotation matrix. (We use the inverse matrix here because here we want to rotate grid coordinates to the observer frame, whereas inside traceray() we rotate observer coordinates to the grid frame.)
    */
    for(i=0;i<2;i++){
      xs[i]=0.0;
      for(di=0;di<DIM;di++){
        xs[i] += gp[gi].x[di]*img[im].rotMat[di][i];
      }
    }

    assignRayOnImage(xs, pixelSize, imgCentreXPixels, imgCentreYPixels, img, im, maxNumRaysPerPixel, rays, &numActiveRaysInternal);
  } /* End loop 1, over grid points. */

  /* Add the circle rays:
  */
  numActiveRays = numActiveRaysInternal;
  if(numCircleRays>0){
    scale = 2.0*M_PI/(double)numCircleRays;
    for(i=0;i<numCircleRays;i++){
      angle = i*scale;
      xs[0] = par->radius*cos(angle);
      xs[1] = par->radius*sin(angle);
      assignRayOnImage(xs, pixelSize, imgCentreXPixels, imgCentreYPixels, img, im, maxNumRaysPerPixel, rays, &numActiveRays);
    }
  }
#ifndef NO_PROGBARS
  oneOnNumActiveRaysMinus1 = 1.0/(double)(numActiveRaysInternal-1);
#endif

  if(numActiveRays<par->pIntensity+numCircleRays)
    rays = realloc(rays, sizeof(rayData)*numActiveRays);

  if(par->traceRayAlgorithm==1){
    delaunay(DIM, gp, (unsigned long)par->ncell, 1, 0, &dc, &numCells); /* mallocs dc if getCells==T */
    /*
Required elements of gp:
		.id
		.x

Sets elements of gp:
		.sink
		.numNeigh
		.neigh
    */

//**** Actually we can figure out the cell geometry from the grid neighbours.

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

    vertexCoords = extractGridXs(DIM, (unsigned long)par->ncell, gp); /* Reads gp[*].x */
    cells = convertCellType(DIM, numCells, dc, gp); /* Reads gp[*].x */
    free(dc);

    if(img[im].doline && img[im].doInterpolateVels){
      /* Set up the buffer we will use to store various quantities when doing barycentric interpolation of velocities:
      */
      velBuff.numVertices = DIM+1;
      velBuff.numEdges = velBuff.numVertices*(velBuff.numVertices-1)/2;

      velBuff.entryCellBary = malloc(sizeof(*velBuff.entryCellBary)*velBuff.numVertices);
      velBuff.midCellBary   = malloc(sizeof(*velBuff.midCellBary)  *velBuff.numVertices);
      velBuff.exitCellBary  = malloc(sizeof(*velBuff.exitCellBary) *velBuff.numVertices);
      velBuff.vertexVels    = malloc(sizeof(*velBuff.vertexVels)   *velBuff.numVertices);
      for(i=0;i<velBuff.numVertices;i++)
        velBuff.vertexVels[i] = malloc(sizeof(**velBuff.vertexVels)*DIM);

      velBuff.edgeVertexIndices = malloc(sizeof(*velBuff.edgeVertexIndices)*velBuff.numEdges);
      velBuff.edgeVels          = malloc(sizeof(*velBuff.edgeVels)         *velBuff.numEdges);
      for(i=0;i<velBuff.numEdges;i++)
        velBuff.edgeVels[i] = malloc(sizeof(**velBuff.edgeVels)*DIM);

      velBuff.shapeFns = malloc(sizeof(*velBuff.shapeFns)*(velBuff.numVertices+velBuff.numEdges));

      /* Populate the edge indices key for the cells:
      */
      ei = 0;
      for(i0=0;i0<velBuff.numVertices-1;i0++){
        for(i1=i0+1;i1<velBuff.numVertices;i1++){
          velBuff.edgeVertexIndices[ei][0] = i0;
          velBuff.edgeVertexIndices[ei][1] = i1;
          ei++;
        }
      }

      ptrToBuff = &velBuff;
    }

  }else if(par->traceRayAlgorithm!=0){
    if(!silent) bail_out("Unrecognized value of par.traceRayAlgorithm");
    exit(1);
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
#ifndef NO_PROGBARS
    int threadI = omp_get_thread_num();
#endif
    int ii, si, ri;
    gridInterp gips[numInterpPoints];

    if(par->traceRayAlgorithm==1){
      /* Allocate memory for the interpolation points:
      */
      if(img[im].doline){
        for(ii=0;ii<numInterpPoints;ii++){
          gips[ii].mol = malloc(sizeof(*(gips[ii].mol))*par->nSpecies);
          for(si=0;si<par->nSpecies;si++){
            gips[ii].mol[si].specNumDens\
              = malloc(sizeof(*(gips[ii].mol[si].specNumDens))*md[si].nlev);
            gips[ii].mol[si].cont    = NULL;
            gips[ii].mol[si].pops    = NULL;
            gips[ii].mol[si].partner = NULL;
          }
        }
      }else{ /* continuum image */
        for(ii=0;ii<numInterpPoints;ii++)
          gips[ii].mol = NULL;
      }
    }

    #pragma omp for schedule(dynamic)
    for(ri=0;ri<numActiveRaysInternal;ri++){
      if(par->traceRayAlgorithm==0)
        traceray(rays[ri], im, par, gp, md, img\
          , cutoff, nStepsThruCell, oneOnNSteps);

      else if(par->traceRayAlgorithm==1)
        traceray_smooth(rays[ri], im, par, gp, vertexCoords, md, img\
          , cells, numCells, epsilon, gips, ptrToBuff\
          , numSegments, oneOnNumSegments);

#ifndef NO_PROGBARS
      if (threadI == 0){ /* i.e., is master thread */
        progFraction = (double)(ri)*oneOnNumActiveRaysMinus1;
        if(!silent) progressbar(progFraction, 13);
      }
#endif
    }

    if(par->traceRayAlgorithm==1){
      for(ii=0;ii<numInterpPoints;ii++)
        freePopulation(par->nSpecies, gips[ii].mol);
    }
  } /* End of parallel block. */

  gsl_set_error_handler(defaultErrorHandler);
  if(!silent) printDone(13);

  if(par->traceRayAlgorithm==1){
    free(cells);
    free(vertexCoords);
    if(img[im].doline && img[im].doInterpolateVels){
      free(velBuff.shapeFns);
      for(i=0;i<velBuff.numEdges;i++)
        free(velBuff.edgeVels[i]);
      free(velBuff.edgeVels);
      free(velBuff.edgeVertexIndices);
      for(i=0;i<velBuff.numVertices;i++)
        free(velBuff.vertexVels[i]);
      free(velBuff.vertexVels);
      free(velBuff.exitCellBary);
      free(velBuff.midCellBary);
      free(velBuff.entryCellBary);
    }
  }

  /* We take, in formal terms, all the 'active' or accepted rays on the model-radius circle to be outside the model; thus we set their intensity and tau to zero.
  */
  for(ri=numActiveRaysInternal;ri<numActiveRays;ri++){
    for(ichan=0;ichan<img[im].nchan;ichan++){
      rays[ri].intensity[ichan] = 0.0;
      rays[ri].tau[      ichan] = 0.0;
    }
  }

  /* For pixels with more than a cutoff number of rays, just average those rays into the pixel:
  */
  for(ri=0;ri<numActiveRays;ri++){
    if(rays[ri].isInsideImage && img[im].pixel[rays[ri].ppi].numRays >= minNumRaysForAverage){
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
    double *grid2DCoords=NULL,rasterStarts[2],rasterDirs[2]={0.0,1.0};
    struct simplex *cells2D=NULL;
    unsigned long num2DCells,gis[3];
    intersectType entryIntcptFirstCell,*cellExitIntcpts=NULL;
    unsigned long *chainOfCellIds=NULL,*rasterCellIDs=NULL;
    int lenChainPtrs=0,status=0,startYi,si;
    double triangle[3][2],barys[3],y,deltaY;
    _Bool *rasterPixelIsInCells=NULL;

    rasterCellIDs        = malloc(sizeof(*rasterCellIDs)*img[im].pxls);
    rasterPixelIsInCells = malloc(sizeof(*rasterPixelIsInCells)*img[im].pxls);

    grid2DCoords = malloc(sizeof(double)*2*numActiveRays);
    for(ri=0;ri<numActiveRays;ri++) {
      grid2DCoords[ri*2+0] = rays[ri].x;
      grid2DCoords[ri*2+1] = rays[ri].y;
    }

    get2DCells(rays, numActiveRays, &cells2D, &num2DCells);

    rasterStarts[1] = pixelSize*(0.5 - imgCentreYPixels);
    for(xi=0;xi<img[im].pxls;xi++){
      x = pixelSize*(0.5 + xi - imgCentreXPixels);
      rasterStarts[0] = x;

      for(yi=0;yi<img[im].pxls;yi++){
        rasterPixelIsInCells[yi] = 0; /* default - signals that the pixel is outside the cell mesh. */
        rasterCellIDs[yi] = 0;
      }

      status = followRayThroughCells(2, rasterStarts, rasterDirs, grid2DCoords\
        , cells2D, num2DCells, epsilon, NULL, &entryIntcptFirstCell, &chainOfCellIds\
        , &cellExitIntcpts, &lenChainPtrs);

      if(status!=0)
    continue;

      startYi = img[im].pxls; /* default */
      for(yi=0;yi<img[im].pxls;yi++){
        deltaY = pixelSize*yi;
        if(deltaY>=entryIntcptFirstCell.dist){
          startYi = yi;
      break;
        }
      }

      /* Obtain the cell ID for each raster pixel:
      */
      si = 0;
      for(yi=startYi;yi<img[im].pxls;yi++){
        deltaY = pixelSize*yi;

        while(si<lenChainPtrs && deltaY>=cellExitIntcpts[si].dist)
          si++;

        if(si>=lenChainPtrs)
      break;

        rasterCellIDs[yi] = chainOfCellIds[si];
        rasterPixelIsInCells[yi] = 1;
      }

      /* Now interpolate for each pixel of the raster:
      */
      for(yi=0;yi<img[im].pxls;yi++){
        ppi = yi*img[im].pxls + xi;
        if(img[im].pixel[ppi].numRays >= minNumRaysForAverage)
      continue;

        y = pixelSize*(0.5 + yi - imgCentreYPixels);

        if(rasterPixelIsInCells[yi]){
          dci = rasterCellIDs[yi]; /* Just for short. */
          for(vi=0;vi<3;vi++){
            gis[vi] = cells2D[dci].vertx[vi];
            triangle[vi][0] = rays[gis[vi]].x;
            triangle[vi][1] = rays[gis[vi]].y;
          }

          calcTriangleBaryCoords(triangle, x, y, barys);

          for(ichan=0;ichan<img[im].nchan;ichan++){
            img[im].pixel[ppi].intense[ichan] += barys[0]*rays[gis[0]].intensity[ichan]\
                                               + barys[1]*rays[gis[1]].intensity[ichan]\
                                               + barys[2]*rays[gis[2]].intensity[ichan];
            img[im].pixel[ppi].tau[    ichan] += barys[0]*rays[gis[0]].tau[ichan]\
                                               + barys[1]*rays[gis[1]].tau[ichan]\
                                               + barys[2]*rays[gis[2]].tau[ichan];
          } /* End loop over ichan */
        } /* End if rasterPixelIsInCells */
      } /* End loop over yi */

      free(chainOfCellIds);
      free(cellExitIntcpts);
    } /* End loop over xi */

    free(cells2D);
    free(grid2DCoords);
    free(rasterPixelIsInCells);
    free(rasterCellIDs);
  } /* end if(numPixelsForInterp>0) */

  for(ri=0;ri<numActiveRays;ri++){
    free(rays[ri].tau);
    free(rays[ri].intensity);
  }
  free(rays);

  /*
Add and subtract appropriate amounts of cmb.

Some explanation is probably helpful here to explain what is going on. If we think of a ray at a given frequency passing through the model, the starting value of its intensity I(0) will be the cosmic background value, which in the bands of interest to LIME can be assumed to be the familiar ~2.7K black-body value. This is the value encoded in local_cmb. According to the backwards-propagation algorithm for solving the RTE described in Hogerheijde & van der Tak, Astron. Astrophys. 362, 697 (2000), the final term in the sum giving the intensity is I(0)*exp(-tau), where tau is the accumulated opacity of the model. For rays passing through areas of zero molecular column density (e.g. outside the model radius), the final or total radiation intensity can be assumed to equal I(0). However, LIME users prefer images not have scalar offsets, no matter how faithful to reality the offset is, thus we also subtract away a constant I(0) from the whole image; thus image areas outside the model radius are (in the present function) left at, or returned to, zero, and I(0) (aka local_cmb) is subtracted from all the in-radius pixels after the final RTE addition.

Note further that users also do not like the resulting zero-valued pixels (!), hence the addition of IMG_MIN_ALLOWED done to all such pixels in functions write4Dfits(). Such is life.
  */
  if(par->polarization){ /* just add cmb to Stokes I, which is the first 'channel' */
    lastChan = 0;
  }else{
    lastChan = img[im].nchan;
  }

#ifdef FASTEXP
  for(ppi=0;ppi<totalNumImagePixels;ppi++){
    for(ichan=0;ichan<lastChan;ichan++){
      img[im].pixel[ppi].intense[ichan] += (FastExp(img[im].pixel[ppi].tau[ichan])-1.0)*local_cmb;
    }
  }
#else
  for(ppi=0;ppi<totalNumImagePixels;ppi++){
    for(ichan=0;ichan<lastChan;ichan++){
      img[im].pixel[ppi].intense[ichan] += (exp(   -img[im].pixel[ppi].tau[ichan])-1.0)*local_cmb;
    }
  }
#endif
}

