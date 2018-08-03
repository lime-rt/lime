/*
 *  grid_aux.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#include "lime.h"
#include "gridio.h" /* For struct gridInfoType, struct keywordType etc */

/*....................................................................*/
void mallocAndSetDefaultGrid(struct grid **gp, const size_t numPoints, const size_t numSpecies){
  size_t i,j;

  *gp = malloc(sizeof(**gp)*numPoints);
  for(i=0;i<numPoints; i++){
    (*gp)[i].v1 = NULL;
    (*gp)[i].v2 = NULL;
    (*gp)[i].v3 = NULL;
    (*gp)[i].dir = NULL;
    (*gp)[i].neigh = NULL;
    (*gp)[i].w = NULL;
    (*gp)[i].ds = NULL;
    (*gp)[i].dens=NULL;
    (*gp)[i].t[0]=-1.0;
    (*gp)[i].t[1]=-1.0;
    (*gp)[i].B[0]=0.0;
    (*gp)[i].B[1]=0.0;
    (*gp)[i].B[2]=0.0;
    (*gp)[i].conv=0;

    if(numSpecies > 0){
      (*gp)[i].mol = malloc(sizeof(*(*gp)[i].mol)*numSpecies);
      for(j=0;j<numSpecies;j++){
        (*gp)[i].mol[j].pops        = NULL;
        (*gp)[i].mol[j].specNumDens = NULL;
        (*gp)[i].mol[j].partner     = NULL;
        (*gp)[i].mol[j].cont        = NULL;
        (*gp)[i].mol[j].dopb = 0.0;
        (*gp)[i].mol[j].binv = 0.0;
        (*gp)[i].mol[j].nmol = 0.0;
        (*gp)[i].mol[j].abun = 0.0;
      }
    }else
      (*gp)[i].mol = NULL;
  }
}

/*....................................................................*/
void checkGridDensities(configInfo *par, struct grid *gp){
  /* This checks that none of the density samples is too small. */
  int i;
  static _Bool warningAlreadyIssued=0;
  char errStr[STR_LEN_1];

  if(!silent){ /* Warn if any densities too low. */
    i = 0;
    while(i<par->pIntensity && !warningAlreadyIssued){
      if(gp[i].dens[0]<TYPICAL_ISM_DENS){
        warningAlreadyIssued = 1;
        snprintf(errStr, STR_LEN_1, "gp[%d].dens[0] at %.1e is below typical values for the ISM (~%.1e).", i, gp[i].dens[0], TYPICAL_ISM_DENS);
        warning(errStr);
        warning("This could give you convergence problems. NOTE: no further warnings will be issued.");
      }
      i++;
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
void calcGridMolDensities(configInfo *par, struct grid **gp){
  int id,ispec,i;

  for(id=0;id<par->ncell;id++){
    for(ispec=0;ispec<par->nSpecies;ispec++){
      (*gp)[id].mol[ispec].nmol = 0.0;
      for(i=0;i<par->numDensities;i++)
        (*gp)[id].mol[ispec].nmol += (*gp)[id].mol[ispec].abun*(*gp)[id].dens[i]\
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
		.sink
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
  char message[STR_LEN_1];//****[80];

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

    if(gp[id].numNeigh<=0){
      if(!silent){
        snprintf(message, STR_LEN_1, "qhull failed silently, grid point %lu has 0 neighbours. Smoother gridDensity() might help.", id);
        bail_out(message);
      }
exit(1);
    }

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
                snprintf(message, STR_LEN_1, "Something weird going on. Cannot find a cell with ID %lu", (unsigned long)(neighbor->id));
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
    gp[i].dir=malloc(sizeof(*(gp[i].dir)) *gp[i].numNeigh);
    gp[i].ds =malloc(sizeof(double)*gp[i].numNeigh);
    memset(gp[i].dir, 0., sizeof(*(gp[i].dir)) * gp[i].numNeigh);
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
  if(par->nSpecies>0){
    for(i=0;i<par->ncell;i++){
      fprintf(fp, "%e\n", g[i].mol[0].abun*g[i].dens[0]);
    }
  }else{
    for(i=0;i<par->ncell;i++){
      fprintf(fp, "%e\n", 0.0);
    }
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
getEdgeVelocities(configInfo *par, struct grid *gp){
  int i,k,j,l;
  double vel[3], x[3];
  
  for(i=0;i<par->ncell;i++){
    gp[i].v1=malloc(3*gp[i].numNeigh*sizeof(double));
    gp[i].v2=malloc(3*gp[i].numNeigh*sizeof(double));
    gp[i].v3=malloc(3*gp[i].numNeigh*sizeof(double));

    for(k=0;k<gp[i].numNeigh;k++){
      for(j=0;j<3;j++) x[j]=gp[i].x[j];		
      for(l=0;l<5;l++){
        velocity(x[0],x[1],x[2],vel);	

        if (l==1) {
	  gp[i].v1[3*k]=vel[0]; gp[i].v1[3*k+1]=vel[1]; gp[i].v1[3*k+2]=vel[2];
        }
        if (l==2) {
          gp[i].v2[3*k]=vel[0]; gp[i].v2[3*k+1]=vel[1]; gp[i].v2[3*k+2]=vel[2];
        }
        if (l==3) {
          gp[i].v3[3*k]=vel[0]; gp[i].v3[3*k+1]=vel[1]; gp[i].v3[3*k+2]=vel[2];
        }
		
        for(j=0;j<3;j++) x[j]=x[j]+(gp[i].dir[k].xn[j]*gp[i].ds[k])/4.;
      }
    }
  }

  par->edgeVelsAvailable = 1;
}

/*....................................................................*/
int setupAndWriteGrid(configInfo *par, struct grid *gp, molData *md, char *outFileName){
  const int numKwds=3;
  int i,status = 0;
  struct gridInfoType gridInfo;
  unsigned short i_us;
  struct keywordType *primaryKwds=malloc(sizeof(struct keywordType)*numKwds);

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
    for(i_us=0;i_us<gridInfo.nSpecies;i_us++){ /* md should only ==NULL if gridInfo.nSpecies==0. */
      copyInparStr(md[i_us].molName, &gridInfo.mols[i_us].molName);
      gridInfo.mols[i_us].nLevels = md[i_us].nlev;
      gridInfo.mols[i_us].nLines  = md[i_us].nline;
    }
  }

  i = 0;
  initializeKeyword(&primaryKwds[i]);
  primaryKwds[i].datatype = lime_DOUBLE;
  sprintf(primaryKwds[i].keyname, "RADIUS  ");
  primaryKwds[i].doubleValue = par->radius;
  sprintf(primaryKwds[i].comment, "[m] Model radius.");

  i++;
  initializeKeyword(&primaryKwds[i]);
  primaryKwds[i].datatype = lime_DOUBLE;
  sprintf(primaryKwds[i].keyname, "MINSCALE");
  primaryKwds[i].doubleValue = par->minScale;
  sprintf(primaryKwds[i].comment, "[m] Minimum model scale.");

  i++;
  initializeKeyword(&primaryKwds[i]);
  primaryKwds[i].datatype = lime_INT;
  sprintf(primaryKwds[i].keyname, "NSOLITER");
  primaryKwds[i].intValue = par->nSolveItersDone;
  sprintf(primaryKwds[i].comment, "Number of RTE iterations performed.");

  status = writeGrid(outFileName\
    , gridInfo, primaryKwds, numKwds, gp, par->collPartNames, par->dataFlags);

  freeKeywords(primaryKwds, numKwds);
  freeGridInfo(&gridInfo);

  return status;
}

/*....................................................................*/
void writeGridIfRequired(configInfo *par, struct grid *gp, molData *md, const int dataStageI){
  if(par->writeGridAtStage[dataStageI-1]){
    int status=0;
    char message[80];

    status = setupAndWriteGrid(par, gp, md, par->gridOutFiles[dataStageI-1]);

    if(status){
      sprintf(message, "writeGrid at data stage %d returned with status %d", dataStageI, status);
      if(!silent) bail_out(message);
exit(1);
    }
  }
}

