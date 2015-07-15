/*
 *  grid.c
 *  LIME, The versatile 3D line modeling environment
 *
 *  Created by Christian Brinch on 11/16/06.
 *  Copyright 2006-2014, Christian Brinch,
 *  <brinch@nbi.dk>
 *  Niels Bohr institutet
 *  University of Copenhagen
 *	All rights reserved.
 *
 */

#include "lime.h"


void
gridAlloc(inputPars *par, struct grid **g){
  int i;
  double temp[99];

  *g=malloc(sizeof(struct grid)*(par->pIntensity+par->sinkPoints));
  memset(*g, 0., sizeof(struct grid) * (par->pIntensity+par->sinkPoints));

  if(par->doPregrid || par->restart) par->collPart=1;
  else{
    for(i=0;i<99;i++) temp[i]=-1;
    density(AU,AU,AU,temp);
    i=0;
    par->collPart=0;
    while(temp[i++]>-1) par->collPart++;
  }

  for(i=0;i<(par->pIntensity+par->sinkPoints); i++){
    (*g)[i].a0 = NULL;
    (*g)[i].a1 = NULL;
    (*g)[i].a2 = NULL;
    (*g)[i].a3 = NULL;
    (*g)[i].a4 = NULL;
    (*g)[i].mol = NULL;
    (*g)[i].dir = NULL;
    (*g)[i].neigh = NULL;
    (*g)[i].w = NULL;
    (*g)[i].ds = NULL;
    (*g)[i].dens=malloc(sizeof(double)*par->collPart);
    (*g)[i].abun=malloc(sizeof(double)*par->nSpecies);
    (*g)[i].nmol=malloc(sizeof(double)*par->nSpecies);
    (*g)[i].t[0]=-1;
    (*g)[i].t[1]=-1;
  }
}

void
freePopulation(const inputPars *par, const molData* m, struct populations* pop ) {
  if( pop !=NULL )
    {
      int j,k;
      for( j=0; j<par->nSpecies; j++ )
        {
          if( pop[j].pops != NULL )
            {
              free( pop[j].pops );
            }
          if( pop[j].knu != NULL )
            {
              free( pop[j].knu );
            }
          if( pop[j].dust != NULL )
            {
              free( pop[j].dust );
            }
          if( pop[j].partner != NULL )
            {
              if( m != NULL )
                {
                  for(k=0; k<m[j].npart; k++)
                    {
                      if( pop[j].partner[k].up != NULL )
                        {
                          free(pop[j].partner[k].up);
                        }
                      if( pop[j].partner[k].down != NULL )
                        {
                          free(pop[j].partner[k].down);
                        }
                    }
                }
              free( pop[j].partner );
            }
        }
      free(pop);
    }
}
void
freeGrid(const inputPars *par, const molData* m ,struct grid* g){
  int i;
  if( g != NULL )
    {
      for(i=0;i<(par->pIntensity+par->sinkPoints); i++){
        if(g[i].a0 != NULL)
          {
            free(g[i].a0);
          }
        if(g[i].a1 != NULL)
          {
            free(g[i].a1);
          }
        if(g[i].a2 != NULL)
          {
            free(g[i].a2);
          }
        if(g[i].a3 != NULL)
          {
            free(g[i].a3);
          }
        if(g[i].a4 != NULL)
          {
            free(g[i].a4);
          }
        if(g[i].dir != NULL)
          {
            free(g[i].dir);
          }
        if(g[i].neigh != NULL)
          {
            free(g[i].neigh);
          }
        if(g[i].w != NULL)
          {
            free(g[i].w);
          }
        if(g[i].dens != NULL)
          {
            free(g[i].dens);
          }
        if(g[i].nmol != NULL)
          {
            free(g[i].nmol);
          }
        if(g[i].abun != NULL)
          {
            free(g[i].abun);
          }
        if(g[i].ds != NULL)
          {
            free(g[i].ds);
          }
        if(g[i].mol != NULL)
          {
            freePopulation( par, m, g[i].mol );
          }
      }
      free(g);
    }
}

void
qhull(inputPars *par, struct grid *g){
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
      pt_array[i*DIM+j]=g[i].x[j];
    }
  }

  sprintf(flags,"qhull d Qbb");
  if (!qh_new_qhull(DIM, par->ncell, pt_array, ismalloc, flags, NULL, NULL)) {
    /* Identify points */
    FORALLvertices {
      id=qh_pointid(vertex->point);
      g[id].numNeigh=qh_setsize(vertex->neighbors);
      if(  g[id].neigh != NULL )
        {
          free( g[id].neigh );
        }
      g[id].neigh=malloc(sizeof(struct grid *)*g[id].numNeigh);
      for(k=0;k<g[id].numNeigh;k++) {
        g[id].neigh[k]=NULL;
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
              while(g[simplex[i]].neigh[k] != NULL && g[simplex[i]].neigh[k]->id != g[simplex[j]].id) {
                k++;
              }
              g[simplex[i]].neigh[k]=&g[simplex[j]];
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
    for(k=0;k<g[i].numNeigh;k++){
      if(g[i].neigh[k] != NULL)
        {
          j++;
        }
    }
    g[i].numNeigh=j;
  }
  qh_freeqhull(!qh_ALL);
  qh_memfreeshort (&curlong, &totlong);
  free(pt_array);
}

void
distCalc(inputPars *par, struct grid *g){
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
write_VTK_unstructured_Points(inputPars *par, struct grid *g){
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
dumpGrid(inputPars *par, struct grid *g){
  if(par->gridfile) write_VTK_unstructured_Points(par, g);
}

void
getArea(inputPars *par, struct grid *g, const gsl_rng *ran){
  int i,j,k,b;//=-1;
  double *angle,best;
  /*	double wsum; */
  double x,y,z,pt_phi,sinPtPhi,pt_theta;
  /* Lots of circles approach -- badly broken, needs to be fixed  */
  /*
     for(i=0;i<par->pIntensity;i++){
     angle=malloc(sizeof(double)*g[i].numNeigh);
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
    angle=malloc(sizeof(double)*g[i].numNeigh);
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
getMass(inputPars *par, struct grid *g, const gsl_rng *ran){
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
      farea=malloc(sizeof(double)*pts[i].vps);
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
buildGrid(inputPars *par, struct grid *g){
  double lograd;		/* The logarithm of the model radius		*/
  double logmin;	    /* Logarithm of par->minScale				*/
  double r,theta,phi,sinPhi,x,y,z,semiradius;	/* Coordinates								*/
  double temp;
  int k=0,i;            /* counters									*/
  int flag;

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
    temp=gsl_rng_uniform(ran);
    flag=0;
    /* Pick a point and check if we like it or not */
    do{
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
    g[k].dens[0]=1e-30;
    g[k].t[0]=par->tcmb;
    g[k].t[1]=par->tcmb;
    g[k++].dopb=0.;
  }
  /* end grid allocation */

  qhull(par, g);
  distCalc(par, g);
  smooth(par,g);

  for(i=0;i<par->pIntensity;i++){
    density(    g[i].x[0],g[i].x[1],g[i].x[2], g[i].dens);
    temperature(g[i].x[0],g[i].x[1],g[i].x[2], g[i].t);
    doppler(    g[i].x[0],g[i].x[1],g[i].x[2],&g[i].dopb);	
    abundance(  g[i].x[0],g[i].x[1],g[i].x[2], g[i].abun);
  }

  //	getArea(par,g, ran);
  //	getMass(par,g, ran);
  getVelosplines(par,g);
  dumpGrid(par,g);

  gsl_rng_free(ran);
  if(!silent) done(5);
}


