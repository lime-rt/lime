/*
 *  raythrucells.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2016 The LIME development team
 *
TODO:
 */

#include "raythrucells.h"

void __attribute__((weak))
error(int errCode, char *message){
  printf("Error: %s\n", message);
  exit(1);
}

/*....................................................................*/
faceType *
extractFace(const int numDims, double *vertexCoords, struct simplex *dc\
  , const unsigned long dci, const int fi){
  /* Given a simplex dc[dci] and the face index (in the range {0...numDims}) fi, this returns the desired information about that face. Note that the ordering of the elements of face.r[] is the same as the ordering of the vertices of the simplex, dc[dci].vertx[]; just the vertex fi is omitted.

Note that the element 'centre' of the faceType struct is mean to contain the spatial coordinates of the centre of the simplex, not of the face. This is designed to facilitate orientation of the face and thus to help determine whether rays which cross it are entering or exiting the simplex.
 */

  const int numFaces=numDims+1;
  int vi, vvi, di;
  static faceType face;
  unsigned long gi;

  vvi = 0;
  for(vi=0;vi<numFaces;vi++){
    if(vi!=fi){
      gi = dc[dci].vertx[vi];
      for(di=0;di<numDims;di++){
        face.r[vvi][di] = vertexCoords[numDims*gi+di];
      }
      vvi++;
    }
  }

  for(di=0;di<numDims;di++)
    face.simplexCentre[di] = dc[dci].centre[di];

  return &face;
}

/*....................................................................*/
int
getNewEntryFaceI(const int numDims, const unsigned long dci, const struct simplex newCell){
  /* Finds the index of the old cell in the face list of the new cell. */

  const int numFaces=numDims+1;
  _Bool matchFound = 0;
  int ffi = 0, newEntryFaceI;

  while(ffi<numFaces && matchFound==0){
    if(newCell.neigh[ffi]!=NULL && newCell.neigh[ffi]->id==dci){
      matchFound = 1;
      newEntryFaceI = ffi;
    }
    ffi++;
  }

  /* Sanity check: */
  if(!matchFound)
    error(RTC_ERR_OLD_NOT_FOUND, "Cannot find old cell ID in new cell data.");

  return newEntryFaceI;
}

/*....................................................................*/
facePlusBasisType
calcFaceInNMinus1(const int numDims, const int numVertices, faceType *face){
  /*
Each of the faces of a polytope in N spatial dimensions is itself a polytope in N-1 dimensions. Each face can therefore be represented via the coordinates of its vertices expressed in an (N-1)-dimensional frame oriented so as to be parallel to the face. The function of the present routine is to perform this decomposition and to return both the (N-1)-dimensional frame (specified via the coordinates in N dimensions of its N-1 basis vectors) and the coordinates in it of each of the M face vertices. (If the face is simplicial, M==N.)

The calling routine should call freeFacePlusBasis() after it is finished with the returned object.

Note that numVertices (a.k.a. M) is expected to be >= numDims (a.k.a. N). The thing would not make much sense otherwise.
  */

  facePlusBasisType facePlusBasis;
  double vs[numVertices-1][numDims],dotValue,norm;
  int di,vi,i,j,ddi;

  /* The first part of the routine is finding the components in N-space of the N-1 orthonormal axes parallel to the face. It is essentially a modified (i.e numerically stabler) Gram-schmidt orthonormalization of the first N-1 vectors from face vertex 0 to each of its other vertices. In pseudo-code, suppose we have N-1 vectors *v and want to orthonormalize them to *u:

	for(i=0;i<N-1;i++){
	  u[i] = v[i]
	  for(j=0;j<i;j++){
	    u[i] -= (u[i].u[j])*u[j]
	  }
	  u[i] /= |u[i]|
	}

  */
  for(di=0;di<numDims;di++)
    facePlusBasis.origin[di] = (*face).r[0][di];

  for(vi=0;vi<numVertices-1;vi++){
    for(di=0;di<numDims;di++)
      vs[vi][di] = (*face).r[vi+1][di] - (*face).r[0][di];
  }

  for(i=0;i<numDims-1;i++){
    for(di=0;di<numDims;di++)
      facePlusBasis.axes[i][di] = vs[i][di];

    for(j=0;j<i;j++){
      dotValue = 0.0;
      for(di=0;di<numDims;di++)
        dotValue += facePlusBasis.axes[i][di]*facePlusBasis.axes[j][di];
      for(di=0;di<numDims;di++)
        facePlusBasis.axes[i][di] -= dotValue*facePlusBasis.axes[j][di];
    }

    norm = 0.0;
    for(di=0;di<numDims;di++){
      norm += facePlusBasis.axes[i][di]*facePlusBasis.axes[i][di];
    }
    norm = 1.0/sqrt(norm);
    for(di=0;di<numDims;di++)
      facePlusBasis.axes[i][di] *= norm;
  }

  /* Now we calculate the coords of each vertex in the N-1 system:
  */
  vi = 0;
  for(ddi=0;ddi<numDims-1;ddi++)
    facePlusBasis.r[vi][ddi] = 0.0;

  for(vi=1;vi<numVertices;vi++){
    for(ddi=0;ddi<numDims-1;ddi++){

      dotValue = 0.0;
      for(di=0;di<numDims;di++)
        dotValue += vs[vi-1][di]*facePlusBasis.axes[ddi][di];

      facePlusBasis.r[vi][ddi] = dotValue;
    }
  }

  return facePlusBasis;
}

/*....................................................................*/
intersectType
intersectLineWithFace(const int numDims, double *x, double *dir, faceType *face, const double epsilon){
  /*
This function calculates the intersection between a line and the face of a simplex oriented in that space. Obviously the number of dimensions of the space must be >=2. The intersection point may be expressed as

	px_ = x_ + a*dir_

where px_, x_ and dir_ are vectors and 'a' is a scalar. The routine returns the following information:
	- The value of 'a'.
	- The so-called barycentric coordinates (BC) of px_ in the face.

The scalar 'a' is found as follows. We need the additional point y_ which can be any of the face vertices. We state the vector identity

	a*dir_ = (y_ - x_) + (px_ - y_).

Clearly (px_ - y_) is parallel to the face. If we take the scalar product of both sides with the vector n_ which is normal to the face we arrive at

	a*dir_.n_ = (y_ - x_).n_ + (px_ - y_).n_.

	          = (y_ - x_).n_ + 0.
Thus
	     (y_ - x_).n_
	a = --------------.
	       dir_.n_

Notes:
	* This routine works best when the sides of the face are not too disparate in size.

	* There is, of course, no guarantee that the line actually intersects the face, even if the line and the face are non-parallel. There are also borderline cases, i.e. where the line passes close to a vertex, in which an exact calculation would show that the intersection occurs (or doesn't occur), but the imprecise computed value claims that it doesn't (or does). Intersection may be judged via the values of the barycentric coordinates (BC): if the BC all fall in the interval [0,1], the line intersects the face; if one of the BC is negative, it doesn't.
  */
  const double oneOnEpsilon=1.0/epsilon;
  double vs[numDims-1][numDims],norm[numDims],normDotDx, numerator, pxInFace[numDims-1], tMat[numDims-1][numDims-1], bVec[numDims-1], det;
  int i,j,k,di,ci,ri,vi,ciOfMax,ciOfMin;
  double testSumForCW=0.0,maxSingularValue,singularValue;
  facePlusBasisType facePlusBasis;
  char errStr[80];
  intersectType intcpt;

  for(vi=0;vi<numDims-1;vi++){
    for(di=0;di<numDims;di++)
      vs[vi][di] = (*face).r[vi+1][di]-(*face).r[0][di];
  }

  /* First calculate a normal vector to the face (note that it doesn't need to be of length==1).
  */
  if(numDims==2){
    norm[0] = -vs[0][1];
    norm[1] =  vs[0][0];

  }else if(numDims==3){
    /* Calculate norm via cross product. */
    for(di=0;di<numDims;di++){
      j = (di+1)%numDims;
      k = (di+2)%numDims;
      norm[di] = vs[0][j]*vs[1][k] - vs[0][k]*vs[1][j];
    }

  }else{ /* Assume numDims>3 */
    /* Calculate norm via SVD. */
    int status=0;
    gsl_matrix *matrix = gsl_matrix_alloc(numDims-1, numDims);
    gsl_matrix *svv    = gsl_matrix_alloc(numDims,   numDims);
    gsl_vector *svs    = gsl_vector_alloc(numDims);
    gsl_vector *work   = gsl_vector_alloc(numDims);

    for(ci=0;ci<numDims;ci++){
      for(ri=0;ri<numDims-1;ri++)
        gsl_matrix_set(matrix, ri, ci, vs[ri][ci]);
    }

    status = gsl_linalg_SV_decomp(matrix, svv, svs, work);
    if(status){
      sprintf(errStr, "SVD decomposition failed (GSL error %d).", status);
      error(RTC_ERR_SVD_FAIL, errStr);
    }

    /*
Since we have N-1 equations in N unknowns, we would expect at least one of the N elements of svs to be zero (within rounding error). The column of svv which corresponds to this value should then be the normal vector which we require. We'll just check however that not more than 1 value of svs is zero.

The GSL doco says that SV_decomp returns sorted svs values, but I prefer not to rely on that.
    */

    ci = 0;
    ciOfMax = ci;
    maxSingularValue = gsl_vector_get(svs,ci);
    for(ci=1;ci<numDims;ci++){
      singularValue = gsl_vector_get(svs,ci);
      if(singularValue>maxSingularValue){
        ciOfMax = ci;
        maxSingularValue = singularValue;
      }
    }

    ciOfMin = -1; /* Impossible default. */
    for(ci=0;ci<numDims;ci++){
      if(ci==ciOfMax) continue;

      singularValue = gsl_vector_get(svs,ci);
      if(singularValue*oneOnEpsilon<maxSingularValue){
        if(ciOfMin>=0){
          /* This is an error because it indicates that >1 singular values are 'small'. */
          sprintf(errStr, "Simplex face does not span an N-1 subspace.");
          error(RTC_ERR_NON_SPAN, errStr);
        }

        ciOfMin = ci;
      }
    }

    for(di=0;di<numDims;di++)
      norm[di] = gsl_matrix_get(svv,di,ciOfMin);

    gsl_vector_free(work);
    gsl_vector_free(svs);
    gsl_matrix_free(svv);
    gsl_matrix_free(matrix);
  }

  /* Since we don't know a priori whether the vertices of the face are listed CW or ACW seen from inside the simplex, we will work this out by dotting the normal vector with a vector from the centre of the simplex to vertex 0 of the face. (A simplex is always convex so this ought to work.)
  */
  testSumForCW = 0.0;
  for(di=0;di<numDims;di++)
    testSumForCW += norm[di]*((*face).r[0][di] - (*face).simplexCentre[di]);

  if(testSumForCW<0.0){
    for(di=0;di<numDims;di++)
      norm[di] *= -1.0;
  }

  /* Calculate the scalar (or dot) product between norm and dir.
  */
  normDotDx = 0.0;
  for(di=0;di<numDims;di++)
    normDotDx += norm[di]*dir[di];

  if(normDotDx>0.0){ /* it is an exit face. */
    intcpt.orientation = 1;
  }else if(normDotDx<0.0){ /* it is an entry face. */
    intcpt.orientation = -1;
  }else{ /* normDotDx==0.0, i.e. line and face are parallel. */
    intcpt.orientation = 0;
    for(di=0;di<numDims;di++)
      intcpt.bary[di] = 0.0;
    intcpt.dist    = 0.0;
    intcpt.collPar = 0.0;

    return intcpt;
  }

  /* If we've got to here, we can be sure the line and the face are not parallel, and that we therefore expect meaningful results for the calculation of 'a'.
  */
  numerator = 0.0;
  for(di=0;di<numDims;di++)
    numerator += norm[di]*((*face).r[0][di] - x[di]); /* n_.(y_ - x_) */

  intcpt.dist = numerator/normDotDx;

  /* In order to calculate the barycentric coordinates, we need to set up a N-1 coordinate basis in the plane of the face.
  */
  facePlusBasis = calcFaceInNMinus1(numDims, numDims, face);

  /* Now we want to express the intersection point in these coordinates:
  */
  for(i=0;i<numDims-1;i++){
    pxInFace[i] = 0.0;
    for(di=0;di<numDims;di++)
      pxInFace[i] += (x[di] + intcpt.dist*dir[di] - facePlusBasis.origin[di])*facePlusBasis.axes[i][di];
  }

  /*
The barycentric coordinates x_ = {L_1,L_2,...,L_{N-1}} are given by

	T x_ = b_

where T is an (N-1)*(N-1) matrix with entries

	T_{i,j} = facePlusBasis.r[j+1][i] - facePlusBasis.r[0][i]

and
	b_i = pxInFace[i] - facePlusBasis.r[0][i].

The final BC L_0 is given by

	         _N-1
	         \
	L_0 = 1 - >    L_i.
	         /_i=1
  */

  if(numDims==2 || numDims==3){
    for(i=0;i<numDims-1;i++){
      for(j=0;j<numDims-1;j++)
        tMat[i][j] = facePlusBasis.r[j+1][i] - facePlusBasis.r[0][i];
      bVec[i] = pxInFace[i] - facePlusBasis.r[0][i];
    }

    if(numDims==2)
      intcpt.bary[1] = bVec[0]/tMat[0][0];

    else{ /* numDims==3 */
      det = tMat[0][0]*tMat[1][1] - tMat[0][1]*tMat[1][0];
      /*** We're assuming that the triangle is not pathological, i.e that det!=0. */
      intcpt.bary[1] = ( tMat[1][1]*bVec[0] - tMat[0][1]*bVec[1])/det;
      intcpt.bary[2] = (-tMat[1][0]*bVec[0] + tMat[0][0]*bVec[1])/det;
    }

  }else{ /* Assume numDims>3 */
    int dummySignum,status=0;
    gsl_matrix *gslT = gsl_matrix_alloc(numDims-1, numDims-1);
    gsl_vector *gsl_x = gsl_vector_alloc(numDims-1);
    gsl_vector *gsl_b = gsl_vector_alloc(numDims-1);
    gsl_permutation *p = gsl_permutation_alloc(numDims-1);

    for(i=0;i<numDims-1;i++){
      for(j=0;j<numDims-1;j++)
        gsl_matrix_set(gslT, i, j, facePlusBasis.r[j+1][i] - facePlusBasis.r[0][i]);
      gsl_vector_set(gsl_b, i, pxInFace[i] - facePlusBasis.r[0][i]);
    }

    status = gsl_linalg_LU_decomp(gslT,p,&dummySignum);
    if(status){
      sprintf(errStr, "LU decomposition failed (GSL error %d).", status);
      error(RTC_ERR_LU_DECOMP_FAIL, errStr);
    }

    status = gsl_linalg_LU_solve(gslT,p,gsl_b,gsl_x);
    if(status){
      sprintf(errStr, "LU solver failed (GSL error %d).", status);
      error(RTC_ERR_LU_SOLVE_FAIL, errStr);
    }

    for(i=0;i<numDims-1;i++)
      intcpt.bary[i+1] = gsl_vector_get(gsl_x,i);

    gsl_permutation_free(p);
    gsl_vector_free(gsl_b);
    gsl_vector_free(gsl_x);
    gsl_matrix_free(gslT);
  }

  intcpt.bary[0] = 1.0;
  for(i=1;i<numDims;i++)
    intcpt.bary[0] -= intcpt.bary[i];

  /* Finally, calculate the 'collision parameter':
  */
  di = 0;
  if(intcpt.bary[di] < 0.5)
    intcpt.collPar = intcpt.bary[di];
  else
    intcpt.collPar = 1.0 - intcpt.bary[di];

  for(di=1;di<numDims;di++){
    if(intcpt.bary[di] < 0.5){
      if(intcpt.bary[di] < intcpt.collPar)
        intcpt.collPar = intcpt.bary[di];
    }else{ /* intcpt.bary[di]>=0.5 */
      if(1.0 - intcpt.bary[di] < intcpt.collPar)
        intcpt.collPar = 1.0 - intcpt.bary[di];
    }
  }

  return intcpt;
}

/*....................................................................*/
int
buildRayCellChain(const int numDims, double *x, double *dir, double *vertexCoords\
  , struct simplex *dc, _Bool **cellVisited, unsigned long dci, int entryFaceI\
  , int levelI, int nCellsInChain, const double epsilon\
  , faceType **facePtrs[N_DIMS+1], unsigned long **chainOfCellIds, intersectType **cellExitIntcpts\
  , int *lenChainPtrs){
  /*
This function is designed to follow a ray (defined by a starting locus 'x' and a direction vector 'dir') through a convex connected set of cells (assumed simplicial). The function returns an integer status value directly, and two lists (plus their common length) via the argument interface: chainOfCellIds and cellExitIntcpts. Taken together, these lists define a chain of cells traversed by the ray.

The task of determining which cells are traversed by the ray is simple in principle, but complications arise in computational practice due to the finite precision of floating-point calculations. Where the ray crosses a cell face near to one of its edges, numerical calculation of the 'impact parameter' may return an answer which is erroneous either way: i.e., a ray which just misses a face may be reported as hitting it, and vice versa. To deal with this, a distinction is made between impacts which are (i) 'good', that is far from any face edge; (ii) 'bad', i.e. clearly missing the face; and (iii) 'marginal', i.e. closer to the edge than some preset cutoff which is represented in the argument list by the number 'epsilon'. Note that a tally is kept of those cells which have already been tested for the current ray, and any face which abuts a neighbouring cell which has been visited already will be flagged as 'bad'.

The function therefore looks at all the exit faces of the cell and sorts them into these three categories. How it proceeds then depends on the relative numbers of each, as described below.

	- If there is more than 1 good exit face, an exception is generated. This is a sign that epsilon has been chosen with too small a value.

	- If there is just 1 good exit face, the marginal ones are ignored. The current cell data are appended to the chain and the cell abutting the exit face becomes the new working cell.

	- If there are no good exit faces, we look at the marginal ones. If there is only 1 of these, we loop as before. If there are more than 1, the function is called recursively for each cell on the far side of an exit face.

Thus there are two alternate modes of operation within the function: a straightforward loop along the cells in a single chain, which will continue so long as there is only a single exit face of type either 'good' or 'marginal'; and a recursive launch of the function at a fork in the chain into each of its several possible branches.

The function terminates under the following conditions:
	- It detects that the edge of the model is reached (returns success).
	- There are no exit faces (returns failure).
	- There are no good exit faces and either
		* all of the recursive calls to marginal faces have been unsuccessful (returns failure), or
		* one of these has been successful (returns success).

At a successful termination, therefore, details of all the cells to the edge of the model are correctly stored in chainOfCellIds and cellExitIntcpts, and the number of these cells is returned in lenChainPtrs.

***** Note that it is assumed here that a mis-indentification of the actual cell traversed by a ray will not ultimately matter to the external calling routine. This is only reasonable if whatever function or property is being sampled by the ray does not vary in a stepwise manner at any cell boundary. *****
  */

  const int numFaces=numDims+1;
  _Bool followingSingleChain;
  const int bufferSize=1024;
  int numGoodExits, numMarginalExits, fi, goodExitFis[numFaces], marginalExitFis[numFaces], exitFi, i, status, newEntryFaceI;
  faceType *pFace;
  intersectType intcpt[numFaces];

  followingSingleChain = 1; /* default */
  do{ /* Follow the chain through 'good' cells, i.e. ones for which entry and exit are nicely distant from face edges. (Here we also follow marginal exits if there are no good ones, and only 1 marginal one.) */
    (*cellVisited)[dci] = 1;

    /* Store the current cell ID (we leave storing the exit face for later, when we know what it is). (If there is not enough room in chainOfCellIds and cellExitIntcpts, realloc them to new value of lenChainPtrs.) */
    if(nCellsInChain>=(*lenChainPtrs)){
      *lenChainPtrs += bufferSize;
      *chainOfCellIds  = realloc(*chainOfCellIds,  sizeof(**chainOfCellIds) *(*lenChainPtrs));
      *cellExitIntcpts = realloc(*cellExitIntcpts, sizeof(**cellExitIntcpts)*(*lenChainPtrs));
    }
    (*chainOfCellIds)[nCellsInChain] = dci;

    /* calculate num good and bad exits */
    numGoodExits = 0;
    numMarginalExits = 0;
    for(fi=0;fi<numFaces;fi++){
      if(fi!=entryFaceI && (dc[dci].neigh[fi]==NULL || !(*cellVisited)[dc[dci].neigh[fi]->id])){
        /* Store points for this face: */
        if(facePtrs==NULL){
          pFace = extractFace(numDims, vertexCoords, dc, dci, fi);
        }else{
          pFace = &(*facePtrs)[dci][fi];
        }

        /* Now calculate the intercept: */
        intcpt[fi] = intersectLineWithFace(numDims, x, dir, pFace, epsilon);
        intcpt[fi].fi = fi; /* Ultimately we need this so we can relate the bary coords for the face back to the Delaunay cell. */

        if(intcpt[fi].orientation>0){ /* it is an exit face. */
          if(intcpt[fi].collPar-epsilon>0.0){
            goodExitFis[numGoodExits] = fi;
            numGoodExits++;
          }else if (intcpt[fi].collPar+epsilon>0.0){
            marginalExitFis[numMarginalExits] = fi;
            numMarginalExits++;
          }
        }
      }
    }

    if(numGoodExits>1){
      error(RTC_ERR_BUG, "Some sort of bug: more than 1 firm candidate found for ray exit from cell.");

    }else if(numGoodExits==1 || numMarginalExits==1){
      if(numGoodExits==1)
        exitFi = goodExitFis[0];
      else
        exitFi = marginalExitFis[0];

      /* Store the exit face details: */
      (*cellExitIntcpts)[nCellsInChain] = intcpt[exitFi];

      nCellsInChain++;

      if(dc[dci].neigh[exitFi]==NULL){ /* Signals that we have reached the edge of the model. */
        /* Realloc the ptrs to their final sizes: */
        *chainOfCellIds  = realloc(*chainOfCellIds,  sizeof(**chainOfCellIds) *nCellsInChain);
        *cellExitIntcpts = realloc(*cellExitIntcpts, sizeof(**cellExitIntcpts)*nCellsInChain);
        *lenChainPtrs = nCellsInChain;

        return 0;
      }

      entryFaceI = getNewEntryFaceI(numDims, dci, *(dc[dci].neigh[exitFi]));
      dci = dc[dci].neigh[exitFi]->id;

    }else{
      followingSingleChain = 0;
    }
  }while(followingSingleChain);

  /* Now we have run out of good (or at least single) exit-face options, let's try the marginal ones. */

  if(numMarginalExits<1)
    return 3; /* Unsuccessful end of this chain. */

  /* If we have got to this point, we must have numMarginalExits>1; thus we have a fork in the chain, and must explore each branch. We recurse here because a recursive scheme is the best way to do that. */

  i = 0;
  status = 4; /* default */
  while(i<numMarginalExits && status>0){
    exitFi = marginalExitFis[i];
    (*cellExitIntcpts)[nCellsInChain] = intcpt[exitFi];

    if(dc[dci].neigh[exitFi]==NULL){ /* Signals that we have reached the edge of the model. */
      /* Realloc the ptrs to their final sizes: */
      nCellsInChain++;
      *chainOfCellIds  = realloc(*chainOfCellIds,  sizeof(**chainOfCellIds) *nCellsInChain);
      *cellExitIntcpts = realloc(*cellExitIntcpts, sizeof(**cellExitIntcpts)*nCellsInChain);
      *lenChainPtrs = nCellsInChain;

      status = 0;

    }else{
      newEntryFaceI = getNewEntryFaceI(numDims, dci, *(dc[dci].neigh[exitFi]));

      /* Now we dive into the branch: */
      status = buildRayCellChain(numDims, x, dir, vertexCoords, dc, cellVisited\
        , dc[dci].neigh[exitFi]->id, newEntryFaceI, levelI+1, nCellsInChain+1\
        , epsilon, facePtrs, chainOfCellIds, cellExitIntcpts, lenChainPtrs);
    }
    i++;
  }

  return status;
}

/*....................................................................*/
int
followRayThroughCells(const int numDims, double *x, double *dir, double *vertexCoords\
  , struct simplex *dc, const unsigned long numCells, const double epsilon\
  , faceType **facePtrs[N_DIMS+1], intersectType *entryIntcpt, unsigned long **chainOfCellIds\
  , intersectType **cellExitIntcpts, int *lenChainPtrs){
  /*
The present function follows a ray through a connected, convex set of cells (assumed to be simplices) and returns information about the chain of cells it passes through. If the ray is found to pass through 1 or more cells, the function returns 0, indicating success; if not, it returns a non-zero value. The chain description consists of three pieces of information: (i) intercept information for the entry face of the first cell encountered; (ii) the IDs of the cells in the chain; (iii) intercept information for the exit face of the ith cell.

The calling routine should free chainOfCellIds, cellExitIntcpts & cellVisited after use.

The argument facePtrs may be set to NULL, in which case the function will construct each face from the list of cells etc as it needs it. This saves on memory but takes more time. If the calling routine supplies these values it needs to do something like as follows:

	faceType *pFace,*facePtrs[N_DIMS+1]=malloc(sizeof(*(*facePtrs[N_DIMS+1]))*numFaces); // numFaces must of course be calculated beforehand.
	for(i=0;i<numFaces;i++){
	  for(j=0;j<numDims+1;j++){
	    pFace = extractFace(numDims, vertexCoords, dc, i, j);
	    facePtrs[i][j] = *pFace;
	  }
	}
	status = followRayThroughCells(... &facePtrs, ...);

Note finally that if facePtrs is supplied non-NULL, vertexCoords may be left at NULL. If  is filled, it should be malloc'd as

	vertexCoords = malloc(sizeof(double)*numDims*numPoints);

and filled as

	for(i=0;i<numPoints;i++)
	  for(j=0;j<numDims;j++)
	    vertexCoords[numDims*i+j] = // grid point i, coordinate j


  */

  const int numFaces=numDims+1, maxNumEntryFaces=100;
  int numEntryFaces, fi, entryFis[maxNumEntryFaces], i, status;
  faceType *pFace;
  unsigned long dci, entryDcis[maxNumEntryFaces];
  intersectType intcpt, entryIntcpts[maxNumEntryFaces];
  _Bool *cellVisited=NULL;

  /* Choose a set of starting faces by testing all the 'external' faces of cells which have some. */
  numEntryFaces = 0;
  for(dci=0;dci<numCells;dci++){
    for(fi=0;fi<numFaces;fi++){
      if(dc[dci].neigh[fi]==NULL){ /* means that this face lies on the outside of the model. */
        /* Store points for this face: */
        if(facePtrs==NULL){
          pFace = extractFace(numDims, vertexCoords, dc, dci, fi);
        }else{
          pFace = &(*facePtrs)[dci][fi];
        }

        /* Now calculate the intercept: */
        intcpt = intersectLineWithFace(numDims, x, dir, pFace, epsilon);
        intcpt.fi = fi; /* Ultimately we need this so we can relate the bary coords for the face back to the Delaunay cell. */

        if(intcpt.orientation<0){ /* it is an entry face. */
          if(intcpt.collPar+epsilon>0.0){
            if(numEntryFaces>maxNumEntryFaces)
              error(RTC_ERR_TOO_MANY_ENTRY, "Too many entry faces.");

            entryDcis[   numEntryFaces] = dci;
            entryFis[    numEntryFaces] = fi;
            entryIntcpts[numEntryFaces] = intcpt;
            entryIntcpts[numEntryFaces] = intcpt;
            numEntryFaces++;
          }
        }
      }
    }
  }

  if(numEntryFaces<=0)
    return 2;

  *lenChainPtrs = 1024; /* This can be increased within followCellChain(). */
  *chainOfCellIds  = malloc(sizeof(**chainOfCellIds) *(*lenChainPtrs));
  *cellExitIntcpts = malloc(sizeof(**cellExitIntcpts)*(*lenChainPtrs));
  cellVisited = malloc(sizeof(*cellVisited)*numCells);
  for(dci=0;dci<numCells;dci++)
    cellVisited[dci] = 0;

  i = 0;
  status = 1; /* default */
  while(i<numEntryFaces && status>0){
    status = buildRayCellChain(numDims, x, dir, vertexCoords, dc, &cellVisited\
      , entryDcis[i], entryFis[i], 0, 0, epsilon, facePtrs, chainOfCellIds, cellExitIntcpts, lenChainPtrs);
    i++;
  }

  if(status==0)
    *entryIntcpt = entryIntcpts[i-1];
  /* Note that the order of the bary coords, and the value of fi, are with reference to the vertx list of the _entered_ cell. This can't of course be any other way, because this ray enters this first cell from the exterior of the model, where there are no cells. For all the intersectType objects in the list cellExitIntcpts, the bary coords etc are with reference to the exited cell. */

//*** this is not too good because *entryIntcpt, *chainOfCellIds, *cellExitIntcpts and lenChainPtrs are left at unsuitable values if status!=0.

  free(cellVisited);

  return status;
}

