/*
 *  grid.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2016 The LIME development team
 *
 */

#include "lime.h"


void
gridAlloc(configInfo *par, struct grid **g){
  int i;

  *g=malloc(sizeof(struct grid)*(par->pIntensity+par->sinkPoints));
  memset(*g, 0., sizeof(struct grid) * (par->pIntensity+par->sinkPoints));

  for(i=0;i<(par->pIntensity+par->sinkPoints); i++){
    (*g)[i].v1 = NULL;
    (*g)[i].v2 = NULL;
    (*g)[i].v3 = NULL;
    (*g)[i].mol = NULL;
    (*g)[i].dir = NULL;
    (*g)[i].neigh = NULL;
    (*g)[i].w = NULL;
    (*g)[i].ds = NULL;
    (*g)[i].dens=malloc(sizeof(double)*par->numDensities);
    (*g)[i].abun=malloc(sizeof(double)*par->nSpecies);
    (*g)[i].t[0]=-1;
    (*g)[i].t[1]=-1;
  }
}

void gridLineInit(configInfo *par, molData *md, struct grid *gp){
  int i,id, ilev;

  for(i=0;i<par->nSpecies;i++){
    /* Calculate Doppler and thermal line broadening */
    for(id=0;id<par->ncell;id++) {
      gp[id].mol[i].dopb = sqrt(gp[id].dopb_turb*gp[id].dopb_turb+2.*KBOLTZ/md[i].amass*gp[id].t[0]);
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
      g[id].mol[ispec].nmol = 0.0;
      for(i=0;i<par->numDensities;i++)
        g[id].mol[ispec].nmol += g[id].abun[ispec]*g[id].dens[i]*par->nMolWeights[i];
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

      gasIIdust(gp[id].x[0],gp[id].x[1],gp[id].x[2],&gtd);
      for(iline=0;iline<md[si].nline;iline++){
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
  double fac, uprate, downrate=0.0;

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

/*....................................................................*/
void delaunay(const int numDims, struct grid *gp, const unsigned long numPoints\
  , const _Bool getCells, struct cell **dc, unsigned long *numCells){
  /*
The principal purpose of this function is to perform a Delaunay triangulation for the set of points defined by the input argument 'g'. This is achieved via routines in the 3rd-party package qhull.

A note about qhull nomenclature: a vertex is what you think it is - i.e., a point; but a facet means in this context a Delaunay triangle (in 2D) or tetrahedron (in 3D).

Required elements of structs:
	struct grid *gp:
		.id
		.x

Elements of structs are set as follows:
	struct grid *gp:
		.numNeigh
		.neigh (this is malloc'd too large and at present not realloc'd.)

	cellType *dc (if getCells>0):
		.id
		.neigh
		.vertx
  */

  coordT *pt_array=NULL;
  unsigned long ppi,id,pointIdsThisFacet[numDims+1],idI,idJ,fi,ffi;
  int i,j,k;
  char flags[255];
  boolT ismalloc = False;
  vertexT *vertex,**vertexp;
  facetT *facet, *neighbor, **neighborp;
  int curlong, totlong;
  _Bool neighbourNotFound;

  pt_array=malloc(sizeof(coordT)*numDims*numPoints);
  for(ppi=0;ppi<numPoints;ppi++) {
    for(j=0;j<numDims;j++) {
      pt_array[ppi*numDims+j]=gp[ppi].x[j];
    }
  }

  sprintf(flags,"qhull d Qbb");
  if (qh_new_qhull(numDims, (int)numPoints, pt_array, ismalloc, flags, NULL, NULL)) {
    if(!silent) bail_out("Qhull failed to triangulate");
    exit(1);
  }

  /* Identify points */
  FORALLvertices {
    id=(unsigned long)qh_pointid(vertex->point);
    /* Note this is NOT the same value as vertex->id. Only the id gained via the call to qh_pointid() is the same as the index of the point in the input list. */

    gp[id].numNeigh=qh_setsize(vertex->neighbors);
    /* Note that vertex->neighbors refers to facets abutting the vertex, not other vertices. In general there seem to be more facets surrounding a point than vertices (in fact there seem to be exactly 2x as many). In any case, mallocing to N_facets gives extra room. */

    if(gp[id].neigh!=NULL)
      free( gp[id].neigh );
    gp[id].neigh=malloc(sizeof(struct grid *)*gp[id].numNeigh);
    for(k=0;k<gp[id].numNeigh;k++) {
      gp[id].neigh[k]=NULL;
    }
  }

  /* Identify the Delaunay neighbors of each point. This is a little involved, because the only direct information we have about which vertices are linked to which others is stored in qhull's facetT objects.
  */
  *numCells = 0;
  FORALLfacets {
    if (!facet->upperdelaunay) {
      /* Store the point IDs in a list for convenience. These ID values are conveniently ordered such that qh_pointid() returns ppi for gp[ppi]. 
      */
      j=0;
      FOREACHvertex_ (facet->vertices) pointIdsThisFacet[j++]=(unsigned long)qh_pointid(vertex->point);

      for(i=0;i<numDims+1;i++){
        idI = pointIdsThisFacet[i];
        for(j=0;j<numDims+1;j++){
          idJ = pointIdsThisFacet[j];
          if(i!=j){
            /* Cycle through all the non-NULL links of gp[idI], storing the link if it is new.
            */
            k=0;
            while(gp[idI].neigh[k] != NULL && gp[idI].neigh[k]->id != gp[idJ].id)
              k++;
            gp[idI].neigh[k]=&gp[idJ];
          }
        }
      }
      (*numCells)++;
    }
  }

  for(ppi=0;ppi<numPoints;ppi++){
    j=0;
    for(k=0;k<gp[ppi].numNeigh;k++){
      if(gp[ppi].neigh[k] != NULL)
        j++;
    }
    gp[ppi].numNeigh=j;
  }

  if(getCells){
    (*dc) = malloc(sizeof(**dc)*(*numCells));
    fi = 0;
    FORALLfacets {
      if (!facet->upperdelaunay) {
        (*dc)[fi].id = (unsigned long)facet->id; /* Do NOT expect this to be equal to fi. */
        fi++;
      }
    }

    fi = 0;
    FORALLfacets {
      if (!facet->upperdelaunay) {
        i = 0;
        FOREACHneighbor_(facet) {
          if(neighbor->upperdelaunay){
            (*dc)[fi].neigh[i] = NULL;
          }else{
            /* Have to find the member of *dc with the same id as neighbour.*/
            ffi = 0;
            neighbourNotFound=1;
            while(ffi<(*numCells) && neighbourNotFound){
              if((*dc)[ffi].id==(unsigned long)neighbor->id){
                (*dc)[fi].neigh[i] = &(*dc)[ffi];
                neighbourNotFound = 0;
              }
              ffi++;
            }
            if(ffi>=(*numCells) && neighbourNotFound){
              if(!silent) bail_out("Something weird going on.");
              exit(1);
            }
          }
          i++;
        }

        i = 0;
        FOREACHvertex_( facet->vertices ) {
          id = (unsigned long)qh_pointid(vertex->point);
          (*dc)[fi].vertx[i] = &gp[id];
          i++;
        }

        fi++;
      }
    }
  }

  qh_freeqhull(!qh_ALL);
  qh_memfreeshort (&curlong, &totlong);
  free(pt_array);
}

void
distCalc(configInfo *par, struct grid *g){
  int i,k,l;

  for(i=0;i<par->ncell;i++){
    if( g[i].dir != NULL )
      {
        free( g[i].dir );
      }
    if( g[i].ds != NULL )
      {
        free( g[i].ds );
      }
    g[i].dir=malloc(sizeof(point)*g[i].numNeigh);
    g[i].ds =malloc(sizeof(double)*g[i].numNeigh);
    memset(g[i].dir, 0., sizeof(point) * g[i].numNeigh);
    memset(g[i].ds, 0., sizeof(double) * g[i].numNeigh);
    for(k=0;k<g[i].numNeigh;k++){
      for(l=0;l<3;l++) g[i].dir[k].x[l] = g[i].neigh[k]->x[l] - g[i].x[l];
      g[i].ds[k]=sqrt(g[i].dir[k].x[0]*g[i].dir[k].x[0]+g[i].dir[k].x[1]*g[i].dir[k].x[1]+g[i].dir[k].x[2]*g[i].dir[k].x[2]);
      for(l=0;l<3;l++) g[i].dir[k].xn[l] = g[i].dir[k].x[l]/g[i].ds[k];
    }
    g[i].nphot=ininphot*g[i].numNeigh;
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


void
buildGrid(configInfo *par, struct grid *g){
  double lograd;		/* The logarithm of the model radius		*/
  double logmin;	    /* Logarithm of par->minScale				*/
  double r,theta,phi,sinPhi,x,y,z,semiradius;	/* Coordinates								*/
  double temp;
  int k=0,i,j;            /* counters									*/
  int flag;
  struct cell *dc=NULL; /* Not used at present. */
  unsigned long numCells;
  const int maxNumAttempts=1000;
  _Bool numRandomsThisPoint;
  int numSecondRandoms=0;
  char errStr[80];

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
    g[k].id=k;
    g[k].x[0]=x;
    g[k].x[1]=y;
    if(DIM==3) g[k].x[2]=z;
    else g[k].x[2]=0.;

    g[k].sink=0;
    /* This next step needs to be done, even though it looks stupid */
    g[k].dir=malloc(sizeof(point)*1);
    g[k].ds =malloc(sizeof(double)*1);
    g[k].neigh =malloc(sizeof(struct grid *)*1);
    if(!silent) progressbar((double) k/((double)par->pIntensity-1), 4);
  }
  /* end model grid point assignment */
  if(!silent) done(4);

  if(!silent && numSecondRandoms>0){
    sprintf(errStr, ">1 random point needed for %d grid points out of %d.", numSecondRandoms, par->pIntensity);
    warning(errStr);
  }

  /* Add surface sink particles */
  for(i=0;i<par->sinkPoints;i++){
    theta=gsl_rng_uniform(ran)*2*PI;
    if(DIM==3) z=2*gsl_rng_uniform(ran)-1.;
    else z=0.;
    semiradius=sqrt(1.-z*z);
    x=semiradius*cos(theta);
    y=semiradius*sin(theta);;
    g[k].id=k;
    g[k].x[0]=par->radius*x;
    g[k].x[1]=par->radius*y;
    g[k].x[2]=par->radius*z;
    g[k].sink=1;
    g[k].abun[0]=0;
    g[k].dens[0]=1e-30;//************** what is the low but non zero value for?
    g[k].t[0]=par->tcmb;
    g[k].t[1]=par->tcmb;
    g[k++].dopb_turb=0.;
  }
  /* end grid allocation */

  /* Check that the user has supplied all necessary functions:
  */
  density(    0.0,0.0,0.0, g[0].dens);
  temperature(0.0,0.0,0.0, g[0].t);
  doppler(    0.0,0.0,0.0,&g[0].dopb_turb);
  abundance(  0.0,0.0,0.0, g[0].abun);
  /* Note that velocity() is the only one of the 5 mandatory functions which is still needed (in raytrace) unless par->doPregrid. Therefore we test it already in parseInput(). */

  delaunay(DIM, g, (unsigned long)par->ncell, 0, &dc, &numCells);
  distCalc(par, g);
  smooth(par,g);

  for(i=0;i<par->pIntensity;i++){
    density(    g[i].x[0],g[i].x[1],g[i].x[2], g[i].dens);
    temperature(g[i].x[0],g[i].x[1],g[i].x[2], g[i].t);
    doppler(    g[i].x[0],g[i].x[1],g[i].x[2],&g[i].dopb_turb);
    abundance(  g[i].x[0],g[i].x[1],g[i].x[2], g[i].abun);
  }

  checkGridDensities(par, g);

  if(par->polarization){
    for(i=0;i<par->pIntensity;i++)
      magfield(g[i].x[0],g[i].x[1],g[i].x[2], g[i].B);
  }else{
    for(i=0;i<par->pIntensity;i++){
      g[i].B[0]=0.0;
      g[i].B[1]=0.0;
      g[i].B[2]=0.0;
    }
  }

  //	getArea(par,g, ran);
  //	getMass(par,g, ran);
  getVelocities(par,g);
  dumpGrid(par,g);
  free(dc);

  gsl_rng_free(ran);
  if(!silent) done(5);
}


