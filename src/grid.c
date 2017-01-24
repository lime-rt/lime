/*
 *  grid.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2016 The LIME development team
 *
TODO:
	- In readOrBuildGrid(), test for the presence of the 5 mandatory functions (actually 4, since velocity() is already tested in aux.c:parseInput() ) before doing smoothing.
 */

#include "lime.h"
#include "tree_random.h"
#include "gridio.h"

/*....................................................................*/
void mallocAndSetDefaultGrid(struct grid **gp, const unsigned int numPoints){
  unsigned int i;

  *gp = malloc(sizeof(struct grid)*numPoints);
  for(i=0;i<numPoints; i++){
    (*gp)[i].v1 = NULL;
    (*gp)[i].v2 = NULL;
    (*gp)[i].v3 = NULL;
    (*gp)[i].mol = NULL;
    (*gp)[i].dir = NULL;
    (*gp)[i].neigh = NULL;
    (*gp)[i].w = NULL;
    (*gp)[i].ds = NULL;
    (*gp)[i].dens=NULL;
    (*gp)[i].abun=NULL;
    (*gp)[i].t[0]=-1;
    (*gp)[i].t[1]=-1;
    (*gp)[i].B[0]=-1;
    (*gp)[i].B[1]=-1;
    (*gp)[i].B[2]=-1;
    (*gp)[i].conv=0;
  }
}

/*....................................................................*/
void gridPopsInit(configInfo *par, molData *md, struct grid *gp){
  int i,id,ilev;

  for(i=0;i<par->nSpecies;i++){
    /* Allocate space for populations etc */
    for(id=0;id<par->ncell; id++){
      gp[id].mol[i].pops = malloc(sizeof(double)*md[i].nlev);
      for(ilev=0;ilev<md[i].nlev;ilev++)
        gp[id].mol[i].pops[ilev] = 0.0;
    }
  }
}

/*....................................................................*/
void calcGridMolDoppler(configInfo *par, molData *md, struct grid *gp){
  int i,id;

  for(i=0;i<par->nSpecies;i++){
    /* Calculate Doppler and thermal line broadening */
    for(id=0;id<par->ncell;id++) {
      gp[id].mol[i].dopb = sqrt(gp[id].dopb_turb*gp[id].dopb_turb\
                                + 2.*KBOLTZ/md[i].amass*gp[id].t[0]);
      gp[id].mol[i].binv = 1./gp[id].mol[i].dopb;
    }
  }
}

/*....................................................................*/
void calcGridMolDensities(configInfo *par, struct grid *gp){
  int id,ispec,i;

  for(id=0;id<par->ncell;id++){
    for(ispec=0;ispec<par->nSpecies;ispec++){
      gp[id].mol[ispec].nmol = 0.0;
      for(i=0;i<par->numDensities;i++)
        gp[id].mol[ispec].nmol += gp[id].abun[ispec]*gp[id].dens[i]\
                                  *par->nMolWeights[i];
    }
  }
}

/*....................................................................*/
void calcGridMolSpecNumDens(configInfo *par, molData *md, struct grid *gp){
  int gi,ispec,ei;

  for(gi=0;gi<par->ncell;gi++){
    for(ispec=0;ispec<par->nSpecies;ispec++){
      for(ei=0;ei<md[ispec].nlev;ei++){
        gp[gi].mol[ispec].specNumDens[ei] = gp[gi].mol[ispec].binv\
          *gp[gi].mol[ispec].nmol*gp[gi].mol[ispec].pops[ei];
      }
    }
  }
}

/*....................................................................*/
void readDustFile(char *dustFileName, double **lamtab, double **kaptab\
  , int *nEntries){

  /* NOTE! The calling routine must free lamtab and kaptab after use.
  */
  FILE *fp;
  int i=0,k;
  char string[80];

  /* Open the file and count the values it contains.
  */
  if((fp=fopen(dustFileName, "r"))==NULL){
    if(!silent) bail_out("Error opening dust opacity data file!");
    exit(1);
  }
  while(fgetc(fp) != EOF){
    fgets(string,80,fp);
    i++;
  }
  rewind(fp);

  /* Now read the values.
  */
  if(i>0){
    *lamtab=malloc(sizeof(**lamtab)*i);
    *kaptab=malloc(sizeof(**kaptab)*i);
  } else {
    if(!silent) bail_out("No opacities read");
    exit(1);
  }
  for(k=0;k<i;k++){
    fscanf(fp,"%lf %lf\n", &(*lamtab)[k], &(*kaptab)[k]);
    (*lamtab)[k]=log10((*lamtab)[k]/1e6);
    (*kaptab)[k]=log10((*kaptab)[k]);
  }
  fclose(fp);

  *nEntries = i;
}

/*....................................................................*/
double interpolateKappa(const double freq, double *lamtab, double *kaptab\
  , const int nEntries, gsl_spline *spline, gsl_interp_accel *acc){
  /* Note that the multiplications by 0.1 below are to convert cm^2/g to m^2/kg. */

  double loglam, kappa;

  loglam=log10(CLIGHT/freq);
  if(loglam < lamtab[0])
    kappa = 0.1*pow(10.,kaptab[0] + (loglam-lamtab[0])\
          *(kaptab[1]-kaptab[0])/(lamtab[1]-lamtab[0]));
  else if(loglam > lamtab[nEntries-1])
    kappa = 0.1*pow(10.,kaptab[nEntries-2] + (loglam-lamtab[nEntries-2])\
          *(kaptab[nEntries-1]-kaptab[nEntries-2])\
          /(lamtab[nEntries-1]-lamtab[nEntries-2]));
  else
    kappa = 0.1*pow(10.,gsl_spline_eval(spline,loglam,acc));

  return kappa;
}

/*....................................................................*/
void calcGridContDustOpacity(configInfo *par, const double freq\
  , double *lamtab, double *kaptab, const int nEntries, struct grid *gp){

  double kappa,densityForDust,gtd;
  int id,di;
  gsl_spline *spline = NULL;
  gsl_interp_accel *acc = NULL;

  if(par->dust == NULL)
    kappa = 0.;
  else{
    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_cspline,nEntries);
    gsl_spline_init(spline,lamtab,kaptab,nEntries);
    kappa = interpolateKappa(freq, lamtab, kaptab, nEntries, spline, acc);
  }

  for(id=0;id<par->ncell;id++){
    densityForDust = 0.0;
    for(di=0;di<par->numDensities;di++)
      densityForDust += gp[id].dens[di]*par->dustWeights[di];

    gasIIdust(gp[id].x[0],gp[id].x[1],gp[id].x[2],&gtd);
    gp[id].cont.knu = kappa*2.4*AMU*densityForDust/gtd;
    /* Check if input model supplies a dust temperature. Otherwise use the kinetic temperature. */
    if(gp[id].t[1]<=0.0) { /* Flags that the user has not set it. */
      gp[id].cont.dust = planckfunc(freq,gp[id].t[0]);
    } else {
      gp[id].cont.dust = planckfunc(freq,gp[id].t[1]);
    }
  }

  if(par->dust != NULL){
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
  }
}

/*....................................................................*/
void calcGridLinesDustOpacity(configInfo *par, molData *md, double *lamtab\
  , double *kaptab, const int nEntries, struct grid *gp){

  int iline,id,si,di;
  double *kappatab,gtd,densityForDust,dustToGas,t;
  gsl_spline *spline = NULL;
  gsl_interp_accel *acc = NULL;

  if(par->dust != NULL){
    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_cspline,nEntries);
    gsl_spline_init(spline,lamtab,kaptab,nEntries);
  }

  for(id=0;id<par->ncell;id++){
    for(si=0;si<par->nSpecies;si++){
      free(gp[id].mol[si].cont);
      gp[id].mol[si].cont = malloc(sizeof(*(gp[id].mol[si].cont))*md[si].nline);
    }
  }

  for(si=0;si<par->nSpecies;si++){
    kappatab = malloc(sizeof(*kappatab)*md[si].nline);

    if(par->dust == NULL){
      for(iline=0;iline<md[si].nline;iline++)
        kappatab[iline] = 0.;
    }else{
      for(iline=0;iline<md[si].nline;iline++)
        kappatab[iline] = interpolateKappa(md[si].freq[iline]\
                        , lamtab, kaptab, nEntries, spline, acc);
    }

    for(id=0;id<par->ncell;id++){
      densityForDust = 0.0;
      for(di=0;di<par->numDensities;di++)
        densityForDust += gp[id].dens[di]*par->dustWeights[di];

      gasIIdust(gp[id].x[0],gp[id].x[1],gp[id].x[2],&gtd);
      dustToGas = 2.4*AMU*densityForDust/gtd;
      for(iline=0;iline<md[si].nline;iline++){
        gp[id].mol[si].cont[iline].knu = kappatab[iline]*dustToGas;
        /* Check if input model supplies a dust temperature. Otherwise use the kinetic temperature. */
        if(gp[id].t[1]<=0.0){ /* Flags that the user has not set it. */
          t = gp[id].t[0];
        }else
          t = gp[id].t[1];
        gp[id].mol[si].cont[iline].dust = planckfunc(md[si].freq[iline],t);
      }
    }

    free(kappatab);
  }

  if(par->dust != NULL){
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
  }
}

/*....................................................................*/
void calcGridCollRates(configInfo *par, molData *md, struct grid *gp){
  int i,id,ipart,itrans,itemp,tnint=-1;
  struct cpData part;
  double fac;

  for(i=0;i<par->nSpecies;i++){
    for(id=0;id<par->ncell;id++){
      gp[id].mol[i].partner = malloc(sizeof(struct rates)*md[i].npart);
    }

    for(ipart=0;ipart<md[i].npart;ipart++){
      part = md[i].part[ipart];
      for(id=0;id<par->ncell;id++){
        for(itrans=0;itrans<part.ntrans;itrans++){
          if((gp[id].t[0]>part.temp[0])&&(gp[id].t[0]<part.temp[part.ntemp-1])){
            for(itemp=0;itemp<part.ntemp-1;itemp++){
              if((gp[id].t[0]>part.temp[itemp])&&(gp[id].t[0]<=part.temp[itemp+1])){
                tnint=itemp;
              }
            }
            fac=(gp[id].t[0]-part.temp[tnint])/(part.temp[tnint+1]-part.temp[tnint]);
            gp[id].mol[i].partner[ipart].t_binlow = tnint;
            gp[id].mol[i].partner[ipart].interp_coeff = fac;

	  } else if(gp[id].t[0]<=part.temp[0]) {
	    gp[id].mol[i].partner[ipart].t_binlow = 0;
	    gp[id].mol[i].partner[ipart].interp_coeff = 0.0;
	  } else {
	    gp[id].mol[i].partner[ipart].t_binlow = part.ntemp-2;
	    gp[id].mol[i].partner[ipart].interp_coeff = 1.0;
	  }
        } /* End loop over transitions. */
      } /* End loop over grid points. */
    } /* End loop over collision partners. */
  } /* End loop over radiating molecules. */
}

/*....................................................................*/
void
delaunay(const int numDims, struct grid *gp, const unsigned long numPoints\
  , const _Bool getCells, _Bool checkSink, struct cell **dc, unsigned long *numCells){
  /*
The principal purpose of this function is to perform a Delaunay triangulation for the set of points defined by the input argument 'g'. This is achieved via routines in the 3rd-party package qhull.

A note about qhull nomenclature: a vertex is what you think it is - i.e., a point; but a facet means in this context a Delaunay triangle (in 2D) or tetrahedron (in 3D). This nomenclature arises because the Delaunay cells are indeed facets (or rather projections of facets) of the convex hull constructed in 1 higher dimension.

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
  char message[80];

  /* pt_array contains the grid point locations in the format required by qhull.
  */
  pt_array=malloc(sizeof(coordT)*numDims*numPoints);
  for(ppi=0;ppi<numPoints;ppi++) {
    for(j=0;j<numDims;j++) {
      pt_array[ppi*numDims+j]=gp[ppi].x[j];
    }
  }

  /* Run qhull to generate the Delaunay mesh. (After this, all the information of importance is stored in variables defined in the qhull header.)
  */
  sprintf(flags,"qhull d Qbb Qt");
  if (qh_new_qhull(numDims, (int)numPoints, pt_array, ismalloc, flags, NULL, NULL)) {
    if(!silent) bail_out("Qhull failed to triangulate");
    exit(1);
  }

  if(checkSink){
    FORALLfacets {
      if(!facet->upperdelaunay){
        FOREACHneighbor_(facet) {
          if(neighbor->upperdelaunay){ /* This should indicate that facet lies on the edge of the model. */
            FOREACHvertex_(neighbor->vertices){
              ppi = (unsigned long)qh_pointid(vertex->point);
              if(ppi<numPoints)
                gp[ppi].sink = 1;
            }
          }
        }
      }
    }
  }

  /* Malloc .neigh for each grid point. At present it is not known how many neighbours a point will have, so all the mallocs are larger than needed.
  */
  FORALLvertices {
    id=(unsigned long)qh_pointid(vertex->point);
    /* Note this is NOT the same value as vertex->id. Only the id gained via the call to qh_pointid() is the same as the index of the point in the input list. */

    gp[id].numNeigh=qh_setsize(vertex->neighbors);
    /* Note that vertex->neighbors refers to facets abutting the vertex, not other vertices. In general there seem to be more facets surrounding a point than vertices (in fact there seem to be exactly 2x as many). In any case, mallocing to N_facets gives extra room. */

    free(gp[id].neigh);
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

  /* Count the actual number of neighbours per point.
  */
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

/*....................................................................*/
unsigned long
reorderGrid(const unsigned long numPoints, struct grid *gp){
  /*
The algorithm works its way up the list of points with one index and down with another. The 'up' travel stops at the 1st sink point it finds, the 'down' at the 1st non-sink point. If at that point the 'up' index is lower in value than the 'down', the points are swapped. This is just a tiny bit tricky because we also need to make sure all the neigh pointers are swapped. That's why we make an ordered list of indices and perform the swaps on that as well.
  */

  unsigned long indices[numPoints],upI,dnI,nExtraSinks=0,i,ngi;
  int j;
  struct grid tempGp;

  for(upI=0;upI<numPoints;upI++) indices[upI] = upI;

  upI = 0;
  dnI = numPoints-1;
  while(1){
    while(upI<numPoints && !gp[upI].sink) upI++;
    while(dnI>=0        &&  gp[dnI].sink) dnI--;

  if(upI>=dnI) break;

    nExtraSinks++;

    i = indices[dnI];
    indices[dnI] = indices[upI];
    indices[upI] = i;

    tempGp = gp[dnI];
    gp[dnI] = gp[upI];
    gp[upI] = tempGp;

    /* However we want to retain the .id values as sequential.
    */
    gp[dnI].id = dnI;
    gp[upI].id = upI;
  }

  /*
Now we sort out the .neigh values. An example of how this should work is as follows. Suppose we swapped points 30 and 41. We have fixed up the .id values, but the swap is still shown in the 'indices' array. Thus we will have

	gp[30].id == 30 (but all the other data is from 41)
	gp[41].id == 41 (but all the other data is from 30)
	indices[30] == 41
	indices[41] == 30

Suppose further that the old value of gp[i].neigh[j] is &gp[30]. We detect that we need to fix it (change it to &gp[41]) because ngi=&gp[30].id=30 != indices[ngi=30]=41. gp[i].neigh[j] is then reset to &gp[indices[30]] = &gp[41], i.e. to point to the same data as used to be in location 30.
  */
  for(i=0;i<numPoints;i++){
    for(j=0;j<gp[i].numNeigh;j++){
      ngi = (unsigned long)gp[i].neigh[j]->id;
      if(ngi != indices[ngi])
        gp[i].neigh[j] = &gp[indices[ngi]];
    }
  }

  return nExtraSinks;
}

/*....................................................................*/
void distCalc(configInfo *par, struct grid *gp){
  int i,k,l;

  for(i=0;i<par->ncell;i++){
    free(gp[i].dir);
    free(gp[i].ds);
    gp[i].dir=malloc(sizeof(point) *gp[i].numNeigh);
    gp[i].ds =malloc(sizeof(double)*gp[i].numNeigh);
    memset(gp[i].dir, 0., sizeof(point) * gp[i].numNeigh);
    memset(gp[i].ds, 0., sizeof(double) * gp[i].numNeigh);
    for(k=0;k<gp[i].numNeigh;k++){
      for(l=0;l<3;l++)
        gp[i].dir[k].x[l] = gp[i].neigh[k]->x[l] - gp[i].x[l];
        gp[i].ds[k] = sqrt(  gp[i].dir[k].x[0]*gp[i].dir[k].x[0]\
                         + gp[i].dir[k].x[1]*gp[i].dir[k].x[1]\
                         + gp[i].dir[k].x[2]*gp[i].dir[k].x[2]);
      for(l=0;l<3;l++)
        gp[i].dir[k].xn[l] = gp[i].dir[k].x[l]/gp[i].ds[k];
    }
    gp[i].nphot=RAYS_PER_POINT;
  }
}


/*....................................................................*/
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

/*....................................................................*/
void
dumpGrid(configInfo *par, struct grid *g){
  if(par->gridfile) write_VTK_unstructured_Points(par, g);
}

/*....................................................................*/
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


/*....................................................................*/
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

/*....................................................................*/
int pointEvaluation(configInfo *par, const double uniformRandom, double *r){
  double fracDensity;

  fracDensity = gridDensity(par, r);

  if(uniformRandom < fracDensity) return 1;
  else return 0;
}

/*....................................................................*/
void randomsViaRejection(configInfo *par, const unsigned int desiredNumPoints, gsl_rng *randGen\
  , double (*outRandLocations)[DIM]){

  double lograd; /* The logarithm of the model radius. */
  double logmin; /* Logarithm of par->minScale. */
  double r,theta,phi,sinPhi,z,semiradius;
  double uniformRandom;
  int j,di;
  unsigned int i_u;
  int pointIsAccepted;
  double x[DIM];
  const int maxNumAttempts=1000;
  int numRandomsThisPoint,numSecondRandoms=0;
  char errStr[80];

  lograd=log10(par->radius);
  logmin=log10(par->minScale);

  /* Sample pIntensity number of points */
  for(i_u=0;i_u<desiredNumPoints;i_u++){
    pointIsAccepted=0;
    numRandomsThisPoint=0;
    do{
      uniformRandom=gsl_rng_uniform(randGen);

      if(numRandomsThisPoint==1)
        numSecondRandoms++;
      numRandomsThisPoint++;

      /* Pick a point and check if we like it or not */
      j=0;
      while(!pointIsAccepted && j<maxNumAttempts){
        if(par->sampling==0){
          r=pow(10,logmin+gsl_rng_uniform(randGen)*(lograd-logmin));
          theta=2.*PI*gsl_rng_uniform(randGen);
          phi=PI*gsl_rng_uniform(randGen);
          sinPhi=sin(phi);
          x[0]=r*cos(theta)*sinPhi;
          x[1]=r*sin(theta)*sinPhi;
          if(DIM==3) x[2]=r*cos(phi);
        } else if(par->sampling==1){
          x[0]=(2*gsl_rng_uniform(randGen)-1)*par->radius;
          x[1]=(2*gsl_rng_uniform(randGen)-1)*par->radius;
          if(DIM==3) x[2]=(2*gsl_rng_uniform(randGen)-1)*par->radius;
        } else if(par->sampling==2){
          r=pow(10,logmin+gsl_rng_uniform(randGen)*(lograd-logmin));
          theta=2.*PI*gsl_rng_uniform(randGen);
          if(DIM==3) {
            z=2*gsl_rng_uniform(randGen)-1.;
            semiradius=r*sqrt(1.-z*z);
            z*=r;
            x[2]=z;
          } else {
            semiradius=r;
          }
          x[0]=semiradius*cos(theta);
          x[1]=semiradius*sin(theta);
        } else {
          if(!silent) bail_out("Don't know how to sample model");
          exit(1);
        }
        pointIsAccepted = pointEvaluation(par, uniformRandom, x);
        j++;
      }
    } while(!pointIsAccepted);
    /* Now pointEvaluation has decided that we like the point */

    for(di=0;di<DIM;di++)
      outRandLocations[i_u][di]=x[di];

    if(!silent) progressbar((double)i_u/((double)desiredNumPoints-1), 4);
  }

  if(!silent && numSecondRandoms>0){
    sprintf(errStr, ">1 random point needed for %d grid points out of %u.", numSecondRandoms, desiredNumPoints);
    warning(errStr);
  }
}

/*....................................................................*/
void writeGridIfRequired(configInfo *par, struct grid *gp, molData *md, const int fileFormatI){
  int status = 0;
  char **collPartNames=NULL; /*** this is a placeholder until we start reading these. */
  char message[80];
  int dataStageI=0;

  /* Work out the data stage:
  */
  if(!allBitsSet(par->dataFlags, DS_mask_1)){
    if(!silent) warning("Trying to write at data stage 0.");
    return;
  }

  if(      allBitsSet(par->dataFlags, DS_mask_4)){
    dataStageI = 4;
  }else if(allBitsSet(par->dataFlags, DS_mask_3)){
    dataStageI = 3;
  }else if(allBitsSet(par->dataFlags, DS_mask_2)){
    dataStageI = 2;
  }else{
    dataStageI = 1;
  }

  if(par->writeGridAtStage[dataStageI-1]){
    struct gridInfoType gridInfo;
    unsigned short i_us;
    struct keywordType *primaryKwds=malloc(sizeof(struct keywordType)*1);

    gridInfo.nInternalPoints = par->pIntensity;
    gridInfo.nSinkPoints     = par->sinkPoints;
    gridInfo.nLinks          = 0; /* This quantity is calculated when writing to file. */
    gridInfo.nNNIndices      = 0; /* This quantity is calculated when writing to file. */
    gridInfo.nDims           = DIM;
    gridInfo.nSpecies        = par->nSpecies;
    gridInfo.nDensities      = par->numDensities;
    gridInfo.nLinkVels       = NUM_VEL_COEFFS;
    if(md==NULL)
      gridInfo.mols = NULL;
    else{
      gridInfo.mols = malloc(sizeof(*(gridInfo.mols))*gridInfo.nSpecies);
      for(i_us=0;i_us<gridInfo.nSpecies;i_us++){
        gridInfo.mols[i_us].molName = md[i_us].molName; /* NOTE*** this just copies the pointer. Should be safe enough if we just want to read this value. */
        gridInfo.mols[i_us].nLevels = md[i_us].nlev;
        gridInfo.mols[i_us].nLines  = md[i_us].nline;
      }
    }

    initializeKeyword(&primaryKwds[0]);
    primaryKwds[0].datatype = TDOUBLE;
    primaryKwds[0].keyname = "RADIUS  ";
    primaryKwds[0].doubleValue = par->radius;
    primaryKwds[0].comment = "[m] Model radius.";

    status = writeGrid(par->gridOutFiles[dataStageI-1], fileFormatI\
      , gridInfo, primaryKwds, 1, gp, collPartNames, par->dataFlags);

    free(primaryKwds);
    free(gridInfo.mols);

    if(status){
      sprintf(message, "writeGrid at data stage %d returned with status %d", dataStageI, status);
      if(!silent) bail_out(message);
      exit(1);
    }
  }

  free(collPartNames);
}

/*....................................................................*/
void
readOrBuildGrid(configInfo *par, struct grid **gp){
  const gsl_rng_type *ranNumGenType = gsl_rng_ranlxs2;
  int i,j,k,di,si,levelI=0,status=0,numCollPartRead;
  double theta,semiradius,z;
  double *outRandDensities=NULL;
  double (*outRandLocations)[DIM]=NULL;
  treeRandConstType rinc;
  treeRandVarType rinv;
  struct cell *dc=NULL; /* Not used at present. */
  unsigned long numCells;
  struct gridInfoType gridInfoRead;
  char **collPartNames;
  char message[80];
  double x[DIM];
  treeType tree;

  par->dataFlags = 0;
  if(par->gridInFile!=NULL){
    const int numDesiredKwds=1;
    struct keywordType *desiredKwds=malloc(sizeof(struct keywordType)*numDesiredKwds);

    initializeKeyword(&desiredKwds[0]);
    desiredKwds[0].datatype = TDOUBLE;
    desiredKwds[0].keyname = "RADIUS  ";
    /* Currently not doing anything with the read keyword. */

    status = readGrid(par->gridInFile, lime_FITS, &gridInfoRead, desiredKwds\
      , numDesiredKwds, gp, &collPartNames, &numCollPartRead, &(par->dataFlags));

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

    /* DS_bit_populations may not be set unless all the others (except DS_bit_magfield) are set as well: */
    if(bitIsSet(par->dataFlags, DS_bit_populations)\
    && !allBitsSet(par->dataFlags & DS_mask_all_but_mag, DS_mask_populations)){
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
  } /* End of read grid file. Whether and what we subsequently calculate will depend on the value of par->dataStageI returned. */

  if(!anyBitSet(par->dataFlags, DS_mask_x)){ /* This should only happen if we did not read a file. Generate the grid point locations. */
    mallocAndSetDefaultGrid(gp, (unsigned int)par->ncell);

    rinc.randGen = gsl_rng_alloc(ranNumGenType);	/* Random number generator */
#ifdef TEST
    gsl_rng_set(rinc.randGen,342971);
#else
    gsl_rng_set(rinc.randGen,time(0));
#endif  

    outRandDensities = malloc(sizeof(double   )*par->pIntensity); /* Not used at present; and in fact they are not useful outside this routine, because they are not the values of the physical density at that point, just what densityFunc3D() returns, which is not necessarily the same thing. */
    outRandLocations = malloc(sizeof(*outRandLocations)*par->pIntensity);

    if(par->samplingAlgorithm==0){
      randomsViaRejection(par, (unsigned int)par->pIntensity, rinc.randGen, outRandLocations);

    } else if(par->samplingAlgorithm==1){
      rinc.par = *par;
      rinc.verbosity = 0;
      rinc.numInRandoms      = TREE_N_RANDOMS;
      rinc.maxRecursion      = TREE_MAX_RECURSION;
      rinc.maxNumTrials      = TREE_MAX_N_TRIALS;
      rinc.dither            = TREE_DITHER;
      rinc.maxNumTrialsDbl = (double)rinc.maxNumTrials;
      rinc.doShuffle = 1;
      rinc.doQuasiRandom = 1;

      for(di=0;di<DIM;di++){
        rinv.fieldOrigin[di] = -par->radius;
        rinv.fieldWidth[di] = 2.0*par->radius;
      }
      rinv.expectedDesNumPoints = (double)par->pIntensity;
      rinc.desiredNumPoints = (unsigned int)par->pIntensity;

      rinv.numHighPoints = par->numGridDensMaxima;
      if(par->numGridDensMaxima>0){
        rinv.highPointLocations = malloc(sizeof(*(rinv.highPointLocations))*par->numGridDensMaxima);
        rinv.highPointDensities = malloc(sizeof(double   )*par->numGridDensMaxima);
        for(i=0;i<par->numGridDensMaxima;i++){
          for(di=0;di<DIM;di++){
            rinv.highPointLocations[i][di] = par->gridDensMaxLoc[i][di];
          }
          rinv.highPointDensities[i] = par->gridDensMaxValues[i];
        }
      }else{
        rinv.highPointLocations = NULL;
        rinv.highPointDensities = NULL;
      }

      initializeTree(&rinc, &rinv, gridDensity, &tree);
      constructTheTree(&rinc, &rinv, levelI, gridDensity, &tree);
      fillTheTree(&rinc, &tree, gridDensity, outRandLocations, outRandDensities);

      free(tree.leaves);
      freeRinv(rinv);
      free(rinc.inRandLocations);

    } else {
      if(!silent) bail_out("Unrecognized sampling algorithm.");
      exit(1);
    }

    for(k=0;k<par->pIntensity;k++){
      /* Assign values to the k'th grid point */
      (*gp)[k].id=k;
      (*gp)[k].x[0]=outRandLocations[k][0];
      (*gp)[k].x[1]=outRandLocations[k][1];
      if(DIM==3) (*gp)[k].x[2]=outRandLocations[k][2];
      (*gp)[k].sink=0;
    }

    /* end model grid point assignment */
    if(!silent) printDone(4);

    /* Add surface sink particles */
    for(k=par->pIntensity;k<par->ncell;k++){
      theta=gsl_rng_uniform(rinc.randGen)*2*PI;

      if(DIM==3) {
        z=2*gsl_rng_uniform(rinc.randGen)-1.;
        semiradius=sqrt(1.-z*z);
        x[2]=z;
      } else {
        semiradius=1.0;
      }

      x[0]=semiradius*cos(theta);
      x[1]=semiradius*sin(theta);;
      (*gp)[k].id=k;
      (*gp)[k].x[0]=par->radius*x[0];
      (*gp)[k].x[1]=par->radius*x[1];
      if(DIM==3) (*gp)[k].x[2]=par->radius*x[2];
      (*gp)[k].sink=1;
    }
    /* end grid allocation */

    free(outRandLocations);
    free(outRandDensities);
    gsl_rng_free(rinc.randGen);

    if(par->samplingAlgorithm==0){
      smooth(par,*gp);
      if(!silent) printDone(5);
    }

    par->dataFlags |= DS_mask_1;
  }

  if(onlyBitsSet(par->dataFlags, DS_mask_1)) /* Only happens if (i) we read no file and have constructed this data within LIME, or (ii) we read a file at dataStageI==1. */
    writeGridIfRequired(par, *gp, NULL, lime_FITS);

  if(!allBitsSet(par->dataFlags, DS_mask_neighbours)){
    unsigned long nExtraSinks;

    delaunay(DIM, *gp, (unsigned long)par->ncell, 0, 1, &dc, &numCells);

    /* We just asked delaunay() to flag any grid points with IDs lower than par->pIntensity (which means their distances from model centre are less than the model radius) but which are nevertheless found to be sink points by virtue of the geometry of the mesh of Delaunay cells. Now we need to reshuffle the list of grid points, then reset par->pIntensity, such that all the non-sink points still have IDs lower than par->pIntensity.
    */ 
    nExtraSinks = reorderGrid((unsigned long)par->ncell, *gp);
    par->pIntensity -= nExtraSinks;
    par->sinkPoints += nExtraSinks;

    par->dataFlags |= DS_mask_neighbours;
  }
  distCalc(par, *gp); /* Mallocs and sets .dir & .ds, sets .nphot. We don't store these values so we have to calculate them whether we read a file or not. */

  if(onlyBitsSet(par->dataFlags, DS_mask_2)) /* Only happens if (i) we read no file and have constructed this data within LIME, or (ii) we read a file at dataStageI==2. */
    writeGridIfRequired(par, *gp, NULL, lime_FITS);

  if(!allBitsSet(par->dataFlags, DS_mask_velocity)){
    for(i=0;i<par->pIntensity;i++)
      velocity((*gp)[i].x[0],(*gp)[i].x[1],(*gp)[i].x[2],(*gp)[i].vel);

    /* Set velocity values also for sink points (otherwise Delaunay ray-tracing has problems) */
    for(i=par->pIntensity;i<par->ncell;i++)
      velocity((*gp)[i].x[0],(*gp)[i].x[1],(*gp)[i].x[2],(*gp)[i].vel);

    par->dataFlags |= DS_mask_velocity;
  }

  if(!allBitsSet(par->dataFlags, DS_mask_density)){
    for(i=0;i<par->ncell; i++)
      (*gp)[i].dens = malloc(sizeof(double)*par->numDensities);
    for(i=0;i<par->pIntensity;i++)
      density((*gp)[i].x[0],(*gp)[i].x[1],(*gp)[i].x[2],(*gp)[i].dens);
    for(i=par->pIntensity;i<par->ncell;i++){
      for(j=0;j<par->numDensities;j++)
        (*gp)[i].dens[j]=1e-30;//************** what is the low but non zero value for?
    }

    par->dataFlags |= DS_mask_density;
  }

  checkGridDensities(par, *gp);

  if(!allBitsSet(par->dataFlags, DS_mask_abundance)){
    for(i=0;i<par->ncell; i++)
      (*gp)[i].abun = malloc(sizeof(double)*par->nSpecies);
    for(i=0;i<par->pIntensity;i++)
      abundance((*gp)[i].x[0],(*gp)[i].x[1],(*gp)[i].x[2],(*gp)[i].abun);
    for(i=par->pIntensity;i<par->ncell;i++){
      for(si=0;si<par->nSpecies;si++)
        (*gp)[i].abun[si]=0;
    }

    par->dataFlags |= DS_mask_abundance;
  }

  if(!allBitsSet(par->dataFlags, DS_mask_turb_doppler)){
    for(i=0;i<par->pIntensity;i++)
      doppler((*gp)[i].x[0],(*gp)[i].x[1],(*gp)[i].x[2],&(*gp)[i].dopb_turb);	
    for(i=par->pIntensity;i<par->ncell;i++)
      (*gp)[i].dopb_turb=0.;

    par->dataFlags |= DS_mask_turb_doppler;
  }

  if(!allBitsSet(par->dataFlags, DS_mask_temperatures)){
    for(i=0;i<par->pIntensity;i++)
      temperature((*gp)[i].x[0],(*gp)[i].x[1],(*gp)[i].x[2],(*gp)[i].t);
    for(i=par->pIntensity;i<par->ncell;i++){
      (*gp)[i].t[0]=par->tcmb;
      (*gp)[i].t[1]=par->tcmb;
    }

    par->dataFlags |= DS_mask_temperatures;
  }

  if(!allBitsSet(par->dataFlags, DS_mask_magfield)){
    if(par->polarization){
      for(i=0;i<par->pIntensity;i++)
        magfield((*gp)[i].x[0],(*gp)[i].x[1],(*gp)[i].x[2],(*gp)[i].B);

      par->dataFlags |= DS_mask_magfield;

    }else{
      for(i=0;i<par->pIntensity;i++){
        (*gp)[i].B[0]=0.0;
        (*gp)[i].B[1]=0.0;
        (*gp)[i].B[2]=0.0;
      }
    }

    for(i=par->pIntensity;i<par->ncell;i++){
      (*gp)[i].B[0]=0.0;
      (*gp)[i].B[1]=0.0;
      (*gp)[i].B[2]=0.0;
    }
  }

  if(!allBitsSet(par->dataFlags, DS_mask_ACOEFF)){
    getVelocities(par,*gp); /* Mallocs and sets .v1, .v2, .v3 */

    par->dataFlags |= DS_mask_ACOEFF;
  }

  if(onlyBitsSet(par->dataFlags & DS_mask_all_but_mag, DS_mask_3)) /* Only happens if (i) we read no file and have constructed this data within LIME, or (ii) we read a file at dataStageI==3. */
    writeGridIfRequired(par, *gp, NULL, lime_FITS);

  dumpGrid(par,*gp);
  free(dc);
}


