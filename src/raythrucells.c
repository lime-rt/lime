/*
 *  raythrucells.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2016 The LIME development team
 *
 */

#include "lime.h"

/*....................................................................*/
int followRayThroughDelCells(double x[DIM], double dir[DIM], struct grid *gp\
  , struct cell *dc, const unsigned long numCells, const double epsilon\
  , intersectType *entryIntcpt, unsigned long **chainOfCellIds\
  , intersectType **cellExitIntcpts, int *lenChainPtrs){
  /*
The present function follows a ray through a connected, convex set of Delaunay cells and returns information about the chain of cells it passes through. If the ray is found to pass through 1 or more cells, the function returns 0, indicating success; if not, it returns a non-zero value. The chain description consists of three pieces of information: (i) intercept information for the entry face of the first cell encountered; (ii) the IDs of the cells in the chain; (iii) intercept information for the exit face of the ith cell.
  */

  const int numFaces=DIM+1, maxNumEntryFaces=100;
  int numEntryFaces, fi, entryFis[maxNumEntryFaces], i, status;
  faceType face;
  unsigned long dci, entryDcis[maxNumEntryFaces];//, trialDci;
  intersectType intcpt, entryIntcpts[maxNumEntryFaces];
  _Bool *cellVisited=NULL;

  /* Choose a set of starting faces by testing all the 'external' faces of cells which have some. */
  numEntryFaces = 0;
  for(dci=0;dci<numCells;dci++){
    for(fi=0;fi<numFaces;fi++){
      if(dc[dci].neigh[fi]==NULL){ /* means that this face lies on the outside of the model. */
        /* Store points for this face: */
        face = extractFace(gp, dc, dci, fi);

        /* Now calculate the intercept: */
        intersectLineTriangle(x, dir, face, &intcpt);
        intcpt.fi = fi; /* Ultimately we need this so we can relate the bary coords for the face back to the Delaunay cell. */

        if(intcpt.orientation<0){ /* it is an entry face. */
          if(intcpt.collPar+epsilon>0.0){
            if(numEntryFaces>maxNumEntryFaces){
              if(!silent) bail_out("Too many entry faces.");
              exit(1);
            }

            entryDcis[   numEntryFaces] = dci;
            entryFis[    numEntryFaces] = fi;
            entryIntcpts[numEntryFaces] = intcpt;
            numEntryFaces++;
          }
        }
      }
    }
  }

  if(numEntryFaces<=0)
    return 1;

  *lenChainPtrs = 1024; /* This can be increased within followCellChain(). */
  *chainOfCellIds  = malloc(sizeof(**chainOfCellIds) *(*lenChainPtrs));
  *cellExitIntcpts = malloc(sizeof(**cellExitIntcpts)*(*lenChainPtrs));
  cellVisited = malloc(sizeof(*cellVisited)*numCells);
  for(dci=0;dci<numCells;dci++)
    cellVisited[dci] = 0;

  i = 0;
  status = 1; /* default */
  while(i<numEntryFaces && status>0){
    status = buildRayCellChain(x, dir, gp, dc, &cellVisited, entryDcis[i], entryFis[i], 0, 0, epsilon, chainOfCellIds, cellExitIntcpts, lenChainPtrs);
    i++;
  }

  if(status==0)
    *entryIntcpt = entryIntcpts[i-1];
  /* Note that the order of the bary coords, and the value of fi, are with reference to the vertx list of the _entered_ cell. This can't of course be any other way, because this ray enters this first cell from the exterior of the model, where there are no cells. For all the intersectType objects in the list cellExitIntcpts, the bary coords etc are with reference to the exited cell. */

  free(cellVisited);

  return status;
}

/*....................................................................*/
int buildRayCellChain(double x[DIM], double dir[DIM], struct grid *gp\
  , struct cell *dc, _Bool **cellVisited, unsigned long dci, int entryFaceI\
  , int levelI, int nCellsInChain, const double epsilon\
  , unsigned long **chainOfCellIds, intersectType **cellExitIntcpts, int *lenChainPtrs){
  /*
This function is designed to follow a ray (defined by a starting locus 'x' and a direction vector 'dir') through a convex connected set of Delaunay cells. The function returns an integer status value directly, and two lists (plus their common length) via the argument interface: chainOfCellIds and cellExitIntcpts. Taken together, these lists define a chain of Delaunay cells traversed by the ray.

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

***** Note that it is assumed here that a mis-indentification of the actual cell traversed by a ray will not ultimately matter to the external calling routine. This is only reasonable if whatever function or property is being sampled by the ray does not vary in a step function across cell boundaries. *****
  */

  const int numFaces=DIM+1;
  _Bool followingSingleChain;
  const int bufferSize=1024;
  int numGoodExits, numMarginalExits, fi, goodExitFis[numFaces], marginalExitFis[numFaces], exitFi, i, status, newEntryFaceI;
  faceType face;
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
        face = extractFace(gp, dc, dci, fi);

        /* Now calculate the intercept: */
        intersectLineTriangle(x, dir, face, &intcpt[fi]);
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
      if(!silent) bail_out("Some sort of bug: more than 1 firm candidate found for ray exit from cell.");
      exit(1);

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

      entryFaceI = getNewEntryFaceI(dci, *(dc[dci].neigh[exitFi]));
      dci = dc[dci].neigh[exitFi]->id;

    } else
      followingSingleChain = 0;
  } while(followingSingleChain);

  /* Now we have run out of good (or at least single) exit-face options, let's try the marginal ones. */

  if(numMarginalExits<1)
    return 1; /* Unsuccessful end of this chain. */

  /* If we have got to this point, we must have numMarginalExits>1; thus we have a fork in the chain, and must explore each branch. We recurse here because a recursive scheme is the best way to do that. */

  i = 0;
  status = 1; /* default */
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
      newEntryFaceI = getNewEntryFaceI(dci, *(dc[dci].neigh[exitFi]));

      /* Now we dive into the branch: */
      status = buildRayCellChain(x, dir, gp, dc, cellVisited, dc[dci].neigh[exitFi]->id\
        , newEntryFaceI, levelI+1, nCellsInChain+1, epsilon\
        , chainOfCellIds, cellExitIntcpts, lenChainPtrs);
    }
    i++;
  }

  return status;
}

/*....................................................................*/
int getNewEntryFaceI(const unsigned long dci, const struct cell newCell){
  /* Finds the index of the old cell in the face list of the new cell. */

  const int numFaces=DIM+1;
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
  if(!matchFound){
    bail_out("Cannot find old cell ID in new cell data.");
    exit(1);
  }

  return newEntryFaceI;
}

/*....................................................................*/
faceType extractFace(struct grid *gp, struct cell *dc, const unsigned long dci, const int fi){
  /* Given a simplex dc[dci] and the face index (in the range {0...DIM}) fi, this returns the desired information about that face. Note that the ordering of the elements of face.r[] is the same as the ordering of the vertices of the simplex, dc[dci].vertx[]; just the vertex fi is omitted.

Note that the element 'centre' of the faceType struct is mean to contain the spatial coordinates of the centre of the simplex, not of the face. This is designed to facilitate orientation of the face and thus to help determine whether rays which cross it are entering or exiting the simplex.
 */

  const int numFaces=DIM+1;
  int vi, vvi, di;
  struct grid point;
  faceType face;

  unsigned long gi;

  vvi = 0;
  for(vi=0;vi<numFaces;vi++){
    if(vi!=fi){
      gi = dc[dci].vertx[vi]->id;
      point = gp[gi];
      for(di=0;di<DIM;di++){
        face.r[vvi][di] = point.x[di];
      }
      vvi++;
    }
  }

  for(di=0;di<DIM;di++){
    face.centre[di] = dc[dci].centre[di];
  }

  return face;
}

/*....................................................................*/
void intersectLineTriangle(double x[3], double dir[3], faceType face\
  , intersectType *intcpt){
  /*
This function calculates the intersection between a line in 3 space and a triangle oriented in that space. The intersection point may be expressed as

	px = x + a*dir

where px, x and dir are vectors and 'a' is a scalar. The routine returns the following information:
	- The value of 'a'.
	- The so-called barycentric coordinates (BC) of px in the triangle.

This routine works best when the sides of the triangle are not too disparate in size.

Notes:
	* The line and the triangle may be parallel. In this case, the function returns 1; otherwise 0.

	* There is, of course, no guarantee that the line actually intersects the triangle, even if the line and the triangular plane intersect. There are also borderline cases, i.e. where the line passes close to a vertex, in which an exact calculation would show that the intersection occurs (or doesn't occur), but the imprecise computed value claims that it doesn't (or does). Intersection may be judged via the values of the barycentric coordinates: if these all fall in the interval [0,1], the line intersects the triangle; if one of the BC is negative, it doesn't.

***** NOTE ***** that the intersectType now includes a field for the index of the face. This is not known to or supplied by the present function, so it sets it to -1 just to avoid an uninitialized value.
  */

  const int numDims=3;
  double norm[numDims], normDotDx, numerator, px2D[numDims-1], mat[numDims-1][numDims-1], vec[numDims-1], det;
  int j, k, di;
  double testSumForCW=0.0;
  triangle2D face2D;

  intcpt->fi = -1;

  /* First calculate a normal vector to the triangle from the cross product. Since we don't know a priori whether the vertices of the face are listed CW or ACW seen from inside the cell, we will work this out by dotting the norm vector with a vector from the centre of the cell to vertex 0.
  */
  for(di=0;di<numDims;di++){
    j = (di+1)%numDims;
    k = (di+2)%numDims;
    norm[di] = (face.r[1][j]-face.r[0][j])*(face.r[2][k]-face.r[0][k])\
             - (face.r[1][k]-face.r[0][k])*(face.r[2][j]-face.r[0][j]);

    testSumForCW += norm[di]*(face.r[0][di] - face.centre[di]);
  }

  if(testSumForCW<0.0){
    for(di=0;di<numDims;di++)
      norm[di] *= -1.0;
  }

  /* Calculate the scalar (or dot) product between norm and dir.
  */
  normDotDx = norm[0]*dir[0] + norm[1]*dir[1] + norm[2]*dir[2];

  if(normDotDx>0.0){ /* it is an exit face. */
    intcpt->orientation = 1;
  }else if(normDotDx<0.0){ /* it is an entry face. */
    intcpt->orientation = -1;
  }else{ /* normDotDx==0.0, i.e. line and plane are parallel. */
    intcpt->orientation = 0;
    intcpt->bary[0] = 0.0;
    intcpt->bary[1] = 0.0;
    intcpt->bary[2] = 0.0;
    intcpt->dist    = 0.0;
    intcpt->collPar = 0.0;

    return;
  }

  /* If we've got to here, we can be sure the line and plane are not parallel, and that we therefore expect meaningful results. */

  numerator = norm[0]*(face.r[0][0] - x[0])\
            + norm[1]*(face.r[0][1] - x[1])\
            + norm[2]*(face.r[0][2] - x[2]);

  intcpt->dist = numerator/normDotDx;

  /* In order to calculate the barycentric coordinates, we need to set up a 2D coordinate basis in the plane of the triangle. I'll define the X axis as 10^, where for short I use the generic notation ji_ to mean the vector (tri[j]-tri[i]), and ji^ to mean the unit vector in the same direction. The Y axis is the unit vector in the direction (20_ - (20_.10^)*10^).
  */
  face2D = calcTriangle2D(face);

  /* Now we want to express the intersection point in these coordinates:
  */
  px2D[0] = (x[0] + intcpt->dist*dir[0] - face.r[0][0])*face2D.xAxis[0]\
          + (x[1] + intcpt->dist*dir[1] - face.r[0][1])*face2D.xAxis[1]\
          + (x[2] + intcpt->dist*dir[2] - face.r[0][2])*face2D.xAxis[2];
  px2D[1] = (x[0] + intcpt->dist*dir[0] - face.r[0][0])*face2D.yAxis[0]\
          + (x[1] + intcpt->dist*dir[1] - face.r[0][1])*face2D.yAxis[1]\
          + (x[2] + intcpt->dist*dir[2] - face.r[0][2])*face2D.yAxis[2];

  /*
The barycentric coordinates (L0,L1,L2) are given by

	(tri2D[0][0]-tri2D[2][0]  tri2D[1][0]-tri2D[2][0]) (L0)   (px2D[0]-tri2D[2][0])
	(                                                )*(  ) = (                   ),
	(tri2D[0][1]-tri2D[2][1]  tri2D[1][1]-tri2D[2][1]) (L1)   (px2D[1]-tri2D[2][1])

with L2 = 1 - L0 - L1.
  */
  mat[0][0] = face2D.r[0][0]-face2D.r[2][0];
  mat[0][1] = face2D.r[1][0]-face2D.r[2][0];
  mat[1][0] = face2D.r[0][1]-face2D.r[2][1];
  mat[1][1] = face2D.r[1][1]-face2D.r[2][1];
  vec[0] = px2D[0]-face2D.r[2][0];
  vec[1] = px2D[1]-face2D.r[2][1];
  det = mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0];
  /*** We're assuming that the triangle is not pathological, i.e that det!=0. */
  intcpt->bary[0] = ( mat[1][1]*vec[0] - mat[0][1]*vec[1])/det;
  intcpt->bary[1] = (-mat[1][0]*vec[0] + mat[0][0]*vec[1])/det;
  intcpt->bary[2] = 1.0 - intcpt->bary[0] - intcpt->bary[1];

  di = 0;
  if(intcpt->bary[di] < 0.5)
    intcpt->collPar = intcpt->bary[di];
  else
    intcpt->collPar = 1.0 - intcpt->bary[di];

  for(di=1;di<numDims;di++){
    if(intcpt->bary[di] < 0.5){
      if(intcpt->bary[di] < intcpt->collPar)
        intcpt->collPar = intcpt->bary[di];
    }else{ /* intcpt->bary[di]>=0.5 */
      if(1.0 - intcpt->bary[di] < intcpt->collPar)
        intcpt->collPar = 1.0 - intcpt->bary[di];
    }
  }

  return;
}

/*....................................................................*/
triangle2D calcTriangle2D(faceType face){
/**** all this could be precalculated (norm and mat too). */

  triangle2D face2D;
  double lengthX, lengthY, dotResult;

  face2D.xAxis[0] = face.r[1][0]-face.r[0][0];
  face2D.xAxis[1] = face.r[1][1]-face.r[0][1];
  face2D.xAxis[2] = face.r[1][2]-face.r[0][2];
  lengthX = sqrt(face2D.xAxis[0]*face2D.xAxis[0]\
               + face2D.xAxis[1]*face2D.xAxis[1]\
               + face2D.xAxis[2]*face2D.xAxis[2]);
  face2D.xAxis[0] /= lengthX;
  face2D.xAxis[1] /= lengthX;
  face2D.xAxis[2] /= lengthX;

  dotResult = face2D.xAxis[0]*(face.r[2][0]-face.r[0][0])\
            + face2D.xAxis[1]*(face.r[2][1]-face.r[0][1])\
            + face2D.xAxis[2]*(face.r[2][2]-face.r[0][2]);

  face2D.yAxis[0] = (face.r[2][0]-face.r[0][0]) - dotResult*face2D.xAxis[0];
  face2D.yAxis[1] = (face.r[2][1]-face.r[0][1]) - dotResult*face2D.xAxis[1];
  face2D.yAxis[2] = (face.r[2][2]-face.r[0][2]) - dotResult*face2D.xAxis[2];
  lengthY = sqrt(face2D.yAxis[0]*face2D.yAxis[0]\
               + face2D.yAxis[1]*face2D.yAxis[1]\
               + face2D.yAxis[2]*face2D.yAxis[2]);
  face2D.yAxis[0] /= lengthY;
  face2D.yAxis[1] /= lengthY;
  face2D.yAxis[2] /= lengthY;

  /* The expressions for the triangle vertices in this basis are simple:
  */
  face2D.r[0][0] = 0.0;
  face2D.r[0][1] = 0.0;
  face2D.r[1][0] = lengthX;
  face2D.r[1][1] = 0.0;
  face2D.r[2][0] = dotResult;
  face2D.r[2][1] = lengthY;

  return face2D;
}

/*....................................................................*/
void doBaryInterp(const intersectType intcpt, struct grid *gp\
  , struct gAuxType *gAux, double xCmpntsRay[3], unsigned long gis[3]\
  , molData *md, const int numMols, gridInterp *gip){

  int di, molI, levelI, lineI;

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
        = intcpt.bary[0]*gAux[gis[0]].mol[molI].specNumDens[levelI]\
        + intcpt.bary[1]*gAux[gis[1]].mol[molI].specNumDens[levelI]\
        + intcpt.bary[2]*gAux[gis[2]].mol[molI].specNumDens[levelI];
    }

    for(lineI=0;lineI<md[molI].nline;lineI++){
      (*gip).mol[molI].dust[lineI] = intcpt.bary[0]*gp[gis[0]].mol[molI].dust[lineI]\
                                   + intcpt.bary[1]*gp[gis[1]].mol[molI].dust[lineI]\
                                   + intcpt.bary[2]*gp[gis[2]].mol[molI].dust[lineI];

      (*gip).mol[molI].knu[lineI]  = intcpt.bary[0]*gp[gis[0]].mol[molI].knu[lineI]\
                                   + intcpt.bary[1]*gp[gis[1]].mol[molI].knu[lineI]\
                                   + intcpt.bary[2]*gp[gis[2]].mol[molI].knu[lineI];
    }
  }
}

/*....................................................................*/
void doSegmentInterp(gridInterp gips[3], const int iA, molData *md\
  , const int numMols, const double oneOnNumSegments, const int si){

  const double fracA = (si + 0.5)*oneOnNumSegments, fracB = 1.0 - fracA;
  const int iB = 1 - iA;
  int di, molI, levelI, lineI;

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

    for(lineI=0;lineI<md[molI].nline;lineI++){
      gips[2].mol[molI].dust[lineI] = fracA*gips[iB].mol[molI].dust[lineI]\
                                    + fracB*gips[iA].mol[molI].dust[lineI];

      gips[2].mol[molI].knu[ lineI] = fracA*gips[iB].mol[molI].knu[ lineI]\
                                    + fracB*gips[iA].mol[molI].knu[ lineI];
    }
  }
}

