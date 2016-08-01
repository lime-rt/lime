/*
 *  grid.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2016 The LIME development team
 *
TODO:
  - There is no need to malloc nmol if all the images are non-line.
 */

#include "lime.h"


void gridLineInit(configInfo *par, molData *md, struct grid *gp){
  int i,id, ilev;

  for(i=0;i<par->nSpecies;i++){
    /* Calculate Doppler and thermal line broadening */
    for(id=0;id<par->ncell;id++) {
      gp[id].mol[i].dopb = sqrt(gp[id].dopb*gp[id].dopb+2.*KBOLTZ/md[i].amass*gp[id].t[0]);
      gp[id].mol[i].binv = 1./gp[id].mol[i].dopb;
    }

    /* Allocate space for populations etc */
    for(id=0;id<par->ncell; id++){
      gp[id].mol[i].pops = malloc(sizeof(double)*md[i].nlev);
      gp[id].mol[i].dust = malloc(sizeof(double)*md[i].nline);
      gp[id].mol[i].knu  = malloc(sizeof(double)*md[i].nline);
      for(ilev=0;ilev<md[i].nlev;ilev++) gp[id].mol[i].pops[ilev]=0.0;
    }
  }
}

void calcGridMolDensities(configInfo *par, struct grid *g){
  int id,ispec,i;

  for(id=0;id<par->ncell; id++){
    for(ispec=0;ispec<par->nSpecies;ispec++){
      g[id].nmol[ispec] = 0.0;
      for(i=0;i<par->numDensities;i++)
        g[id].nmol[ispec] += g[id].abun[ispec]*g[id].dens[i]*par->nMolWeights[i];
    }
  }
}

void calcGridDustOpacity(configInfo *par, molData *md, struct grid *gp){
  FILE *fp;
  char string[80];
  int i=0,k,j,iline,id,si,di;
  double loglam, *lamtab, *kaptab, *kappatab, gtd, densityForDust;
  gsl_spline *spline;

  for(si=0;si<par->nSpecies;si++){
    kappatab = malloc(sizeof(*kappatab)*md[si].nline);

    if(par->dust == NULL){
      for(i=0;i<md[si].nline;i++) kappatab[i]=0.;
    } else {
      gsl_interp_accel *acc=gsl_interp_accel_alloc();
      if((fp=fopen(par->dust, "r"))==NULL){
        if(!silent) bail_out("Error opening dust opacity data file!");
        exit(1);
      }
      while(fgetc(fp) != EOF){
        fgets(string,80,fp);
        i++;
      }
      rewind(fp);
      if(i>0){
        lamtab=malloc(sizeof(*lamtab)*i);
        kaptab=malloc(sizeof(*kaptab)*i);
      } else {
        if(!silent) bail_out("No opacities read");
        exit(1);
      }
      for(k=0;k<i;k++){
        fscanf(fp,"%lf %lf\n", &lamtab[k], &kaptab[k]);
        lamtab[k]=log10(lamtab[k]/1e6);
        kaptab[k]=log10(kaptab[k]);
      }
      fclose(fp);
      spline=gsl_spline_alloc(gsl_interp_cspline,i);
      gsl_spline_init(spline,lamtab,kaptab,i);
      for(j=0;j<md[si].nline;j++) {
        loglam=log10(CLIGHT/md[si].freq[j]);
        if(loglam < lamtab[0]){
          kappatab[j]=0.1*pow(10.,kaptab[0] + (loglam-lamtab[0]) * (kaptab[1]-kaptab[0])/(lamtab[1]-lamtab[0]));
        } else if(loglam > lamtab[i-1]){
          kappatab[j]=0.1*pow(10.,kaptab[i-2] + (loglam-lamtab[i-2]) * (kaptab[i-1]-kaptab[i-2])/(lamtab[i-1]-lamtab[i-2]));
        } else kappatab[j]=0.1*pow(10.,gsl_spline_eval(spline,loglam,acc));
      }
      gsl_spline_free(spline);
      gsl_interp_accel_free(acc);
      free(kaptab);
      free(lamtab);
    }

    for(id=0;id<par->ncell;id++){
      densityForDust = 0.0;
      for(di=0;di<par->numDensities;di++)
        densityForDust += gp[id].dens[di]*par->dustWeights[di];

      for(iline=0;iline<md[si].nline;iline++){
        gasIIdust(gp[id].x[0],gp[id].x[1],gp[id].x[2],&gtd);
        gp[id].mol[si].knu[iline]=kappatab[iline]*2.4*AMU*densityForDust/gtd;
        //Check if input model supplies a dust temperature. Otherwise use the kinetic temperature
        if(gp[id].t[1]==-1) {
          gp[id].mol[si].dust[iline]=planckfunc(iline,gp[id].t[0],md,si);
        } else {
          gp[id].mol[si].dust[iline]=planckfunc(iline,gp[id].t[1],md,si);
        }
      }
    }

    free(kappatab);
  }

  return;
}

void calcGridCollRates(configInfo *par, molData *md, struct grid *g){
  int i,id,ipart,itrans,itemp,tnint=-1;
  struct cpData part;
  double fac;

  for(i=0;i<par->nSpecies;i++){
    for(id=0;id<par->ncell;id++){
      g[id].mol[i].partner = malloc(sizeof(struct rates)*md[i].npart);
    }

    for(ipart=0;ipart<md[i].npart;ipart++){
      part = md[i].part[ipart];
      for(id=0;id<par->ncell;id++){
        for(itrans=0;itrans<part.ntrans;itrans++){
          if((g[id].t[0]>part.temp[0])&&(g[id].t[0]<part.temp[part.ntemp-1])){
            for(itemp=0;itemp<part.ntemp-1;itemp++){
              if((g[id].t[0]>part.temp[itemp])&&(g[id].t[0]<=part.temp[itemp+1])){
                tnint=itemp;
              }
            }
            fac=(g[id].t[0]-part.temp[tnint])/(part.temp[tnint+1]-part.temp[tnint]);
            g[id].mol[i].partner[ipart].t_binlow = tnint;
            g[id].mol[i].partner[ipart].interp_coeff = fac;

	  } else if(g[id].t[0]<=part.temp[0]) {
	    g[id].mol[i].partner[ipart].t_binlow = 0;
	    g[id].mol[i].partner[ipart].interp_coeff = 0.0;
	  } else {
	    g[id].mol[i].partner[ipart].t_binlow = part.ntemp-2;
	    g[id].mol[i].partner[ipart].interp_coeff = 1.0;
	  }
        } /* End loop over transitions. */
      } /* End loop over grid points. */
    } /* End loop over collision partners. */
  } /* End loop over radiating molecules. */
}

void
qhull(configInfo *par, struct grid *gp){
  int i,j,k,id;
  char flags[255];
  boolT ismalloc = False;
  facetT *facet;
  vertexT *vertex,**vertexp;
  coordT *pt_array;
  int simplex[DIM+1];
  int curlong, totlong;

  pt_array=malloc(DIM*sizeof(coordT)*par->ncell);

  for(i=0;i<par->ncell;i++) {
    for(j=0;j<DIM;j++) {
      pt_array[i*DIM+j]=gp[i].x[j];
    }
  }

  sprintf(flags,"qhull d Qbb");
  if (!qh_new_qhull(DIM, par->ncell, pt_array, ismalloc, flags, NULL, NULL)) {
    /* Identify points */
    FORALLvertices {
      id=qh_pointid(vertex->point);
      gp[id].numNeigh=qh_setsize(vertex->neighbors);
      if(  gp[id].neigh != NULL )
        {
          free( gp[id].neigh );
        }
      gp[id].neigh=malloc(sizeof(struct grid *)*gp[id].numNeigh);
      for(k=0;k<gp[id].numNeigh;k++) {
        gp[id].neigh[k]=NULL;
      }
    }
    
    /* Identify neighbors */
    FORALLfacets {
      if (!facet->upperdelaunay) {
        j=0;
        FOREACHvertex_ (facet->vertices) simplex[j++]=qh_pointid(vertex->point);
        for(i=0;i<DIM+1;i++){
          for(j=0;j<DIM+1;j++){
            k=0;
            if(i!=j){
              while(gp[simplex[i]].neigh[k] != NULL && gp[simplex[i]].neigh[k]->id != gp[simplex[j]].id) {
                k++;
              }
              gp[simplex[i]].neigh[k]=&gp[simplex[j]];
            }
          }
        }
      }
    }
  } else {
    if(!silent) bail_out("Qhull failed to triangulate");
    exit(1);
  }

  for(i=0;i<par->ncell;i++){
    j=0;
    for(k=0;k<gp[i].numNeigh;k++){
      if(gp[i].neigh[k] != NULL)
        {
          j++;
        }
    }
    gp[i].numNeigh=j;
  }
  qh_freeqhull(!qh_ALL);
  qh_memfreeshort (&curlong, &totlong);
  free(pt_array);
}

void
distCalc(configInfo *par, struct grid *gp){
  int i,k,l;

  for(i=0;i<par->ncell;i++){
    if( gp[i].dir != NULL )
      {
        free( gp[i].dir );
      }
    if( gp[i].ds != NULL )
      {
        free( gp[i].ds );
      }
    gp[i].dir=malloc(sizeof(point)*gp[i].numNeigh);
    gp[i].ds =malloc(sizeof(double)*gp[i].numNeigh);
    memset(gp[i].dir, 0., sizeof(point) * gp[i].numNeigh);
    memset(gp[i].ds, 0., sizeof(double) * gp[i].numNeigh);
    for(k=0;k<gp[i].numNeigh;k++){
      for(l=0;l<3;l++) gp[i].dir[k].x[l] = gp[i].neigh[k]->x[l] - gp[i].x[l];
      gp[i].ds[k]=sqrt(gp[i].dir[k].x[0]*gp[i].dir[k].x[0]+gp[i].dir[k].x[1]*gp[i].dir[k].x[1]+gp[i].dir[k].x[2]*gp[i].dir[k].x[2]);
      for(l=0;l<3;l++) gp[i].dir[k].xn[l] = gp[i].dir[k].x[l]/gp[i].ds[k];
    }
    gp[i].nphot=ininphot*gp[i].numNeigh;
  }
}


void
write_VTK_unstructured_Points(configInfo *par, struct grid *g){
  FILE *fp;
  double length;
  int i,j,l=0;
  char flags[255];
  boolT ismalloc = False;
  facetT *facet;
  vertexT *vertex,**vertexp;
  coordT *pt_array;
  int curlong, totlong;

  pt_array=malloc(sizeof(coordT)*DIM*par->ncell);

  for(i=0;i<par->ncell;i++) {
    for(j=0;j<DIM;j++) {
      pt_array[i*DIM+j]=g[i].x[j];
    }
  }

  if((fp=fopen(par->gridfile, "w"))==NULL){
    if(!silent) bail_out("Error writing grid file!");
    exit(1);
  }
  fprintf(fp,"# vtk DataFile Version 3.0\n");
  fprintf(fp,"Lime grid\n");
  fprintf(fp,"ASCII\n");
  fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");
  fprintf(fp,"POINTS %d float\n",par->ncell);
  for(i=0; i<par->ncell; i++) {
    fprintf(fp,"%e %e %e\n", g[i].x[0], g[i].x[1], g[i].x[2]);
  }
  fprintf(fp, "\n");

  sprintf(flags,"qhull d Qbb T0");

  if (!qh_new_qhull(DIM, par->ncell, pt_array, ismalloc, flags, NULL, NULL)) {
    FORALLfacets {
      if (!facet->upperdelaunay) l++;
    }
    fprintf(fp,"CELLS %d %d\n",l, 5*l);
    FORALLfacets {
      if (!facet->upperdelaunay) {
        fprintf(fp,"4 ");
        FOREACHvertex_ (facet->vertices) {
          fprintf(fp, "%d ", qh_pointid(vertex->point));
        }
        fprintf(fp, "\n");
      }
    }
  }
  qh_freeqhull(!qh_ALL);
  qh_memfreeshort (&curlong, &totlong);
  fprintf(fp,"\nCELL_TYPES %d\n",l);
  for(i=0;i<l;i++){
    fprintf(fp, "10\n");
  }
  fprintf(fp,"POINT_DATA %d\n",par->ncell);
  fprintf(fp,"SCALARS H2_density float 1\n");
  fprintf(fp,"LOOKUP_TABLE default\n");
  for(i=0;i<par->ncell;i++){
    fprintf(fp, "%e\n", g[i].dens[0]);
  }
  fprintf(fp,"SCALARS Mol_density float 1\n");
  fprintf(fp,"LOOKUP_TABLE default\n");
  for(i=0;i<par->ncell;i++){
    fprintf(fp, "%e\n", g[i].abun[0]*g[i].dens[0]);
  }
  fprintf(fp,"SCALARS Gas_temperature float 1\n");
  fprintf(fp,"LOOKUP_TABLE default\n");
  for(i=0;i<par->ncell;i++){
    fprintf(fp, "%e\n", g[i].t[0]);
  }
  fprintf(fp,"VECTORS velocity float\n");
  for(i=0;i<par->ncell;i++){
    length=sqrt(g[i].vel[0]*g[i].vel[0]+g[i].vel[1]*g[i].vel[1]+g[i].vel[2]*g[i].vel[2]);
    if(length > 0.){
      fprintf(fp, "%e %e %e\n", g[i].vel[0]/length,g[i].vel[1]/length,g[i].vel[2]/length);
    } else {
      fprintf(fp, "%e %e %e\n", g[i].vel[0],g[i].vel[1],g[i].vel[2]);
    }
  }

  fclose(fp);
  free(pt_array);
}

void
dumpGrid(configInfo *par, struct grid *g){
  if(par->gridfile) write_VTK_unstructured_Points(par, g);
}

void
getArea(configInfo *par, struct grid *g, const gsl_rng *ran){
  int i,j,k,b;//=-1;
  double *angle,best;
  /*	double wsum; */
  double x,y,z,pt_phi,sinPtPhi,pt_theta;
  /* Lots of circles approach -- badly broken, needs to be fixed  */
  /*
     for(i=0;i<par->pIntensity;i++){
     angle=malloc(sizeof(*angle)*g[i].numNeigh);
     g[i].w=malloc(sizeof(double)*g[i].numNeigh);
     for(k=0;k<g[i].numNeigh;k++){
     best=0;
     for(j=0;j<g[i].numNeigh;j++){
     angle[j]=( g[i].dir[k].xn[0]*g[i].dir[j].xn[0]
     +g[i].dir[k].xn[1]*g[i].dir[j].xn[1]
     +g[i].dir[k].xn[2]*g[i].dir[j].xn[2]);
     if(angle[j]>best && angle[j]<1){
     best=angle[j];
     b=j;
     }
     }
     g[i].w[k]=PI*pow(tan(acos(best)),2);
     wsum+=g[i].w[k];
     }
     for(k=0;k<g[i].numNeigh;k++) g[i].w[k]=g[i].w[k]/wsum;
     free(angle);
     }
     */


  for(i=0;i<par->pIntensity;i++){
    angle=malloc(sizeof(*angle)*g[i].numNeigh);
    g[i].w=malloc(sizeof(double)*g[i].numNeigh);
    memset(g[i].w, 0, sizeof(double) * g[i].numNeigh);
    for(k=0;k<1000;k++){
      pt_theta=gsl_rng_uniform(ran)*2*PI;
      pt_phi=gsl_rng_uniform(ran)*PI;
      sinPtPhi=sin(pt_phi);
      x=cos(pt_theta)*sinPtPhi;
      y=sin(pt_theta)*sinPtPhi;
      z=cos(pt_phi);
      j=0;
      best=( x*g[i].dir[j].xn[0]
            +y*g[i].dir[j].xn[1]
            +z*g[i].dir[j].xn[2]);
      b=j;
      for(j=1;j<g[i].numNeigh;j++){
        angle[j]=( x*g[i].dir[j].xn[0]
                  +y*g[i].dir[j].xn[1]
                  +z*g[i].dir[j].xn[2]);
        if(angle[j]>best){
          best=angle[j];
          b=j;
        }
      }
      g[i].w[b]+=0.001;
    }
    free(angle);
  }
}


void
getMass(configInfo *par, struct grid *g, const gsl_rng *ran){
  double mass=0.,dist;
  double vol=0.,dp,dpbest,*farea,suma;
  int i,k,j,best=-1;
  typedef struct {coordT *pt_array;int vps;int *flag;} S;
  S *pts;
  char flags[255];
  boolT ismalloc = False;
  facetT *facet, *neighbor, **neighborp;
  vertexT *vertex;
  coordT *pt_array;
  int curlong, totlong;

  pts=malloc(sizeof(S)*par->pIntensity);
  for(i=0;i<par->pIntensity;i++){
    pts[i].vps=0;
  }
  pt_array=malloc(DIM*sizeof(coordT)*par->ncell);
  for(i=0;i<par->ncell;i++) {
    for(j=0;j<DIM;j++) {
      pt_array[i*DIM+j]=g[i].x[j];
    }
  }

  sprintf(flags,"qhull v Qbb");
  if (!qh_new_qhull(DIM, par->ncell, pt_array, ismalloc, flags, NULL, NULL)) {
    qh_setvoronoi_all();
    FORALLvertices {
      i=qh_pointid(vertex->point);
      if(i<par->pIntensity){
        pts[i].vps=0;
        FOREACHneighbor_(vertex) {
          if (!neighbor->upperdelaunay) pts[i].vps++;
        }
        if(pts[i].vps > 0){
          pts[i].pt_array=malloc(DIM*sizeof(coordT)*pts[i].vps);
          pts[i].flag=malloc(DIM*sizeof(int)*pts[i].vps);
        } else {
          if(!silent) bail_out("Qhull error");
          exit(0);
        }
        k=0;
        FOREACHneighbor_(vertex) {
          if (!neighbor->upperdelaunay) {
            for(j=0;j<DIM;j++) pts[i].pt_array[k*DIM+j]=neighbor->center[j]+(gsl_rng_uniform(ran)*2-1)*par->radius/1e6;
            k+=1;
          }
        }
      }
    }
  } else {
    if(!silent) bail_out("Qhull error");
    exit(0);
  }
  qh_freeqhull(!qh_ALL);
  qh_memfreeshort (&curlong, &totlong);

  /* Calculate Voronoi volumes */
  sprintf(flags,"qhull ");
  for(i=0;i<par->pIntensity;i++){
    g[i].w=malloc(sizeof(double)*g[i].numNeigh);
    vol=0.;
    suma=0.;
    // farea=malloc(sizeof(double)*pts[i].vps);
    if(pts[i].vps>0){
      farea=malloc(sizeof(*farea)*pts[i].vps);
    } else {
      if(!silent) bail_out("Qhull error");
      exit(0);
    }
    if (!qh_new_qhull(DIM, pts[i].vps, pts[i].pt_array, ismalloc, flags, NULL, NULL)) {
      FORALLfacets {
        dpbest=0.;
        for(j=0;j<g[i].numNeigh;j++){
          dp=facet->normal[0]*g[i].dir[j].xn[0]
            +facet->normal[1]*g[i].dir[j].xn[1]
            +facet->normal[2]*g[i].dir[j].xn[2];
          if(fabs(dp)>dpbest){
            dpbest=fabs(dp);
            best=j;
          }
        }
        if (!facet->normal)
          continue;
        if (facet->upperdelaunay && qh ATinfinity)
          continue;
        farea[best]=qh_facetarea (facet);
        suma+=farea[best];
        if (!qh DELAUNAY) {
          qh_distplane (qh interior_point, facet, &dist);
          vol += -dist * farea[best]/ qh hull_dim;
        }
      }
      for(j=0;j<g[i].numNeigh;j++) g[i].w[j]=farea[j]/suma; //if(g[i].w[j]<1e-2) g[i].w[j]=1e-2;
      free(pts[i].flag);
      free(pts[i].pt_array);
      mass+=vol*g[i].dens[0];
    }
    free(farea);
  }
  qh_freeqhull(!qh_ALL);
  qh_memfreeshort (&curlong, &totlong);
  free(pts);
  free(pt_array);
  if(!silent) quotemass(mass*2.37*1.67e-27/1.989e30);
}


void mallocAndSetDefaultGrid(struct grid **gp, const unsigned int numPoints){
  unsigned int i;

  *gp = malloc(sizeof(struct grid)*numPoints);
  for(i=0;i<numPoints; i++){
    (*gp)[i].a0 = NULL;
    (*gp)[i].a1 = NULL;
    (*gp)[i].a2 = NULL;
    (*gp)[i].a3 = NULL;
    (*gp)[i].a4 = NULL;
    (*gp)[i].mol = NULL;
    (*gp)[i].dir = NULL;
    (*gp)[i].neigh = NULL;
    (*gp)[i].w = NULL;
    (*gp)[i].ds = NULL;
    (*gp)[i].dens=NULL;
    (*gp)[i].abun=NULL;
    (*gp)[i].nmol=NULL;
    (*gp)[i].t[0]=-1;
    (*gp)[i].t[1]=-1;
    (*gp)[i].conv=0;
  }
}

void readOrBuildGrid(configInfo *par, struct grid **gp){
  double lograd;		/* The logarithm of the model radius		*/
  double logmin;	    /* Logarithm of par->minScale				*/
  double r,theta,phi,sinPhi,x,y,z,semiradius;	/* Coordinates								*/
  double temp;
  int k=0,i,j;            /* counters									*/
  int flag;
  struct gridInfoType gridInfoRead;
  int status;
  char **collPartNames;
  int numCollPartRead;
  char message[80];
  const int maxNumAttempts=1000;
  _Bool numRandomsThisPoint;
  int numSecondRandoms=0;
  char errStr[80];

  par->dataFlags = 0;
  if(par->gridInFile!=NULL){
    status = readGrid(par->gridInFile, lime_FITS, &gridInfoRead, gp\
      , &collPartNames, &numCollPartRead, &(par->dataFlags));

    if(status){
      if(!silent){
        sprintf(message, "Read of grid file failed with status return %d", status);
        bail_out(message);
      }
      exit(1);
    }

    /* Test that dataFlags obeys the rules. */
    /* No other bit may be set if DS_bit_x is not: */
    if(anyBitSet(par->dataFlags, (DS_mask_all & ~(1 << DS_bit_x))) && !bitIsSet(par->dataFlags, DS_bit_x)){
      if(!silent) bail_out("You may not read a grid file without X, ID or IS_SINK data.");
      exit(1);
    }

    /* DS_bit_ACOEFF may not be set if either DS_bit_neighbours or DS_bit_velocity is not: */
    if(bitIsSet(par->dataFlags, DS_bit_ACOEFF)\
    && !(bitIsSet(par->dataFlags, DS_bit_neighbours) && bitIsSet(par->dataFlags, DS_bit_velocity))){
      if(!silent) bail_out("You may not read a grid file with ACOEFF but no VEL or neighbour data.");
      exit(1);
    }

    /* DS_bit_populations may not be set unless all the others are set as well: */
    if(bitIsSet(par->dataFlags, DS_bit_populations)\
    && !allBitsSet(par->dataFlags, (DS_mask_all & ~(1 << DS_bit_populations)))){
      if(!silent) bail_out("You may not read a grid file with pop data unless all other data is present.");
      exit(1);
    }

    /* Test gridInfoRead values against par values and overwrite the latter, with a warning, if necessary.
    */
    if(gridInfoRead.nSinkPoints>0 && par->sinkPoints>0){
      if((int)gridInfoRead.nSinkPoints!=par->sinkPoints){
        if(!silent) warning("par->sinkPoints will be overwritten");
        par->sinkPoints = (int)gridInfoRead.nSinkPoints;
      }
      if((int)gridInfoRead.nInternalPoints!=par->pIntensity){
        if(!silent) warning("par->pIntensity will be overwritten");
        par->pIntensity = (int)gridInfoRead.nInternalPoints;
      }
      par->ncell = par->sinkPoints + par->pIntensity;
    }
    if(gridInfoRead.nDims!=DIM){ /* At present this situation is already detected and handled inside readGridExtFromFits(), but it may not be in future. The test here has no present functionality but saves trouble later if we change grid.x from an array to a pointer. */
      if(!silent){
        sprintf(message, "Grid file had %d dimensions but there should be %d.", (int)gridInfoRead.nDims, DIM);
        bail_out(message);
      }
      exit(1);
    }
    if(gridInfoRead.nSpecies>0 && par->nSpecies>0 && (int)gridInfoRead.nSpecies!=par->nSpecies){
      if(!silent){
        sprintf(message, "Grid file had %d species but you have provided moldata files for %d."\
          , (int)gridInfoRead.nSpecies, par->nSpecies);
        bail_out(message);
      }
      exit(1);
/**** should compare name to name - at some later time after we have read these from the moldata files? */
    }
    if(gridInfoRead.nDensities>0 && par->numDensities>0 && (int)gridInfoRead.nDensities!=par->numDensities){
      if(!silent){
        sprintf(message, "Grid file had %d densities but you have provided %d."\
          , (int)gridInfoRead.nDensities, par->numDensities);
        bail_out(message);
      }
      exit(1);
    }

/*
**** Ideally we should also have a test on nACoeffs.

**** Ideally we should also have a test on the mols entries - at some later time after we have read the corresponding values from the moldata files?
*/
    if(allBitsSet(par->dataFlags, DS_mask_4) && !silent)
      warning("Sorry, LIME is not yet able to further process level populations read from file.");
  } /* End of read grid file. Whether and what we subsequently calculate will depend on the value of par->dataStageI returned. */

  if(!anyBitSet(par->dataFlags, DS_mask_x)){ /* This should only happen if we did not read a file. Generate the grid point locations. */
    mallocAndSetDefaultGrid(gp, (unsigned int)par->ncell);

    gsl_rng *ran = gsl_rng_alloc(gsl_rng_ranlxs2);	/* Random number generator */
#ifdef TEST
    gsl_rng_set(ran,342971);
#else
    gsl_rng_set(ran,time(0));
#endif  
  
    lograd=log10(par->radius);
    logmin=log10(par->minScale);

    /* Sample pIntensity number of points */
    for(k=0;k<par->pIntensity;k++){
      flag=0;
      numRandomsThisPoint=0;
      do{
        temp=gsl_rng_uniform(ran);

        if(numRandomsThisPoint==1)
          numSecondRandoms++;
        numRandomsThisPoint++;

        /* Pick a point and check if we like it or not */
        j=0;
        while(!flag && j<maxNumAttempts){
          if(par->sampling==0){
            r=pow(10,logmin+gsl_rng_uniform(ran)*(lograd-logmin));
            theta=2.*PI*gsl_rng_uniform(ran);
            phi=PI*gsl_rng_uniform(ran);
            sinPhi=sin(phi);
            x=r*cos(theta)*sinPhi;
            y=r*sin(theta)*sinPhi;
            if(DIM==3) z=r*cos(phi);
            else z=0.;
          } else if(par->sampling==1){
            x=(2*gsl_rng_uniform(ran)-1)*par->radius;
            y=(2*gsl_rng_uniform(ran)-1)*par->radius;
            if(DIM==3) z=(2*gsl_rng_uniform(ran)-1)*par->radius;
            else z=0;
          } else if(par->sampling==2){
            r=pow(10,logmin+gsl_rng_uniform(ran)*(lograd-logmin));
            theta=2.*PI*gsl_rng_uniform(ran);
            if(DIM==3) {
              z=2*gsl_rng_uniform(ran)-1.;
              semiradius=r*sqrt(1.-z*z);
              z*=r;
            } else {
              z=0.;
              semiradius=r;
            }
            x=semiradius*cos(theta);
            y=semiradius*sin(theta);
          } else {
            if(!silent) bail_out("Don't know how to sample model");
            exit(1);
          }
          if((x*x+y*y+z*z)<par->radiusSqu) flag=pointEvaluation(par,temp,x,y,z);
          j++;
        }
      } while(!flag);
      /* Now pointEvaluation has decided that we like the point */

      /* Assign values to the k'th grid point */
      /* I don't think we actually need to do this here... */
      (*gp)[k].id=k;
      (*gp)[k].x[0]=x;
      (*gp)[k].x[1]=y;
      if(DIM==3) (*gp)[k].x[2]=z;
      else (*gp)[k].x[2]=0.;

      (*gp)[k].sink=0;

      if(!silent) progressbar((double) k/((double)par->pIntensity-1), 4);
    }
    /* end model grid point assignment */
    if(!silent) printDone(4);

    if(!silent && numSecondRandoms>0){
      sprintf(errStr, ">1 random point needed for %d grid points out of %d.", numSecondRandoms, par->pIntensity);
      warning(errStr);
    }

    /* Add surface sink particles */
    for(k=par->pIntensity;k<par->ncell;k++){
      theta=gsl_rng_uniform(ran)*2*PI;
      if(DIM==3) z=2*gsl_rng_uniform(ran)-1.;
      else z=0.;
      semiradius=sqrt(1.-z*z);
      x=semiradius*cos(theta);
      y=semiradius*sin(theta);;
      (*gp)[k].id=k;
      (*gp)[k].x[0]=par->radius*x;
      (*gp)[k].x[1]=par->radius*y;
      (*gp)[k].x[2]=par->radius*z;
      (*gp)[k].sink=1;
    }
    /* end grid allocation */

    gsl_rng_free(ran);

    smooth(par, *gp);
    if(!silent) printDone(5);

    par->dataFlags |= DS_mask_1;
  }

  if(onlyBitsSet(par->dataFlags, DS_mask_1)) /* Only happens if (i) we read no file and have constructed this data within LIME, or (ii) we read a file at dataStageI==1. */
    writeGridIfRequired(par, *gp, NULL, lime_FITS);

  if(!allBitsSet(par->dataFlags, DS_mask_neighbours)){
    qhull(par, *gp); /* Mallocs and sets .neigh, sets .numNeigh */

    par->dataFlags |= DS_mask_2;
  }
  distCalc(par, *gp); /* Mallocs and sets .dir & .ds, sets .nphot. We don't store these values so we have to calculate them whether we read a file or not. */

  if(onlyBitsSet(par->dataFlags, DS_mask_2)) /* Only happens if (i) we read no file and have constructed this data within LIME, or (ii) we read a file at dataStageI==2. */
    writeGridIfRequired(par, *gp, NULL, lime_FITS);

  if(!allBitsSet(par->dataFlags, DS_mask_velocity)){
    for(i=0;i<par->pIntensity;i++)
      velocity((*gp)[i].x[0],(*gp)[i].x[1],(*gp)[i].x[2], (*gp)[i].vel);
    for(i=par->pIntensity;i<par->ncell;i++){
      for(j=0;j<DIM;j++)
        (*gp)[i].vel[j]=0.;
    }

    par->dataFlags |= DS_mask_velocity;
  }

  if(!allBitsSet(par->dataFlags, DS_mask_density)){
    for(i=0;i<par->ncell; i++)
      (*gp)[i].dens = malloc(sizeof(double)*par->numDensities);
    for(i=0;i<par->pIntensity;i++)
      density((*gp)[i].x[0],(*gp)[i].x[1],(*gp)[i].x[2], (*gp)[i].dens);
    for(i=par->pIntensity;i<par->ncell;i++)
      (*gp)[i].dens[0]=1e-30;

    par->dataFlags |= DS_mask_density;
  }

  checkGridDensities(par, *gp);

  if(!allBitsSet(par->dataFlags, DS_mask_abundance)){
    for(i=0;i<par->ncell; i++){
      (*gp)[i].abun = malloc(sizeof(double)*par->nSpecies);
    }
    for(i=0;i<par->pIntensity;i++)
      abundance((*gp)[i].x[0],(*gp)[i].x[1],(*gp)[i].x[2], (*gp)[i].abun);
    for(i=par->pIntensity;i<par->ncell;i++)
      (*gp)[i].abun[0]=0;

    par->dataFlags |= DS_mask_abundance;
  }

  for(i=0;i<par->ncell; i++) /* We don't store the nmol values so we have to do this malloc whether we read a file or not. */
    (*gp)[i].nmol = malloc(sizeof(double)*par->nSpecies); //**** mind you, it would be better to malloc them just before calculating them in molinit.

  if(!allBitsSet(par->dataFlags, DS_mask_turb_doppler)){
    for(i=0;i<par->pIntensity;i++)
      doppler((*gp)[i].x[0],(*gp)[i].x[1],(*gp)[i].x[2],&(*gp)[i].dopb);	
    for(i=par->pIntensity;i<par->ncell;i++)
      (*gp)[i].dopb=0.;

    par->dataFlags |= DS_mask_turb_doppler;
  }

  if(!allBitsSet(par->dataFlags, DS_mask_temperatures)){
    for(i=0;i<par->pIntensity;i++)
      temperature((*gp)[i].x[0],(*gp)[i].x[1],(*gp)[i].x[2], (*gp)[i].t);
    for(i=par->pIntensity;i<par->ncell;i++){
      (*gp)[i].t[0]=par->tcmb;
      (*gp)[i].t[1]=par->tcmb;
    }

    par->dataFlags |= DS_mask_temperatures;
  }

  if(!allBitsSet(par->dataFlags, DS_mask_ACOEFF)){
    calcInterpCoeffs(par,*gp); /* Mallocs and sets .a0, .a1 etc. */

    par->dataFlags |= DS_mask_3;
  }

  if(onlyBitsSet(par->dataFlags, DS_mask_3)) /* Only happens if (i) we read no file and have constructed this data within LIME, or (ii) we read a file at dataStageI==3. */
    writeGridIfRequired(par, *gp, NULL, lime_FITS);

  dumpGrid(par,*gp);
}


