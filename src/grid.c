/*
 *  grid.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
TODO:
	- In readOrBuildGrid(), test for the presence of the 5 mandatory functions (actually 4, since velocity() is already tested in aux.c:parseInput() ) before doing smoothing.
 */

#include "lime.h"
#include "tree_random.h"
#include "gridio.h"

/*....................................................................*/
void
dumpGrid(configInfo *par, struct grid *g){
  if(par->gridfile) write_VTK_unstructured_Points(par, g);
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
  char errStr[STR_LEN_0];
#ifndef NO_PROGBARS
  double progFraction;
#endif

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
          theta=2.*M_PI*gsl_rng_uniform(randGen);
          phi=M_PI*gsl_rng_uniform(randGen);
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
          theta=2.*M_PI*gsl_rng_uniform(randGen);
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

#ifndef NO_PROGBARS
    progFraction = (double)i_u/((double)desiredNumPoints-1);
    if(!silent) progressbar(progFraction, 4);
#endif
  }

  if(!silent && numSecondRandoms>0){
    snprintf(errStr, STR_LEN_0, ">1 random point needed for %d grid points out of %u.", numSecondRandoms, desiredNumPoints);
    warning(errStr);
  }
}

/*....................................................................*/
void
treePrintMessage(const int status, const char message[TREE_STRLEN]){
  char errStr[STR_LEN_0];

  if(silent)
return;

  if(     status==TREE_MSG_MESSAGE)
    printMessage((char*)message);
  else if(status==TREE_MSG_WARN)
    warning((char*)message);
  else if(status==TREE_MSG_ERROR)
    bail_out((char*)message);
  else{
    snprintf(errStr, STR_LEN_0, "Message status %d not understood.\n", status);
    bail_out(errStr);
exit(1);
  }
}

/*....................................................................*/
void
buildGrid(configInfo *par, struct grid **gp){
  const gsl_rng_type *ranNumGenType = gsl_rng_ranlxs2;
  int i,j,k,di,si;
  double theta,semiradius,z,dummyT[2],dummyScalar;
  double *outRandDensities=NULL,*dummyPointer=NULL,x[DIM];
  double (*outRandLocations)[DIM]=NULL;
  treeRandConstType rinc;
  gsl_rng *randGen;
  struct cell *dc=NULL; /* Not used at present. */
  unsigned long numCells;

  /* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
Check for the existence of any mandatory functions we have not supplied grid values for.

Note that we need density and temperature values whether par->doMolCalcs or not.
  */
  if(!allBitsSet(par->dataFlags, DS_mask_density)){
    if(bitIsSet(defaultFuncFlags, USERFUNC_density)){
      if(!silent) bail_out("You need to supply a density() function.");
exit(1);
    }
  }

  if(!allBitsSet(par->dataFlags, DS_mask_temperatures)){
    if(bitIsSet(defaultFuncFlags, USERFUNC_temperature)){
      if(!silent) bail_out("You need to supply a temperature() function.");
exit(1);
    }
  }

  if(par->doMolCalcs){
    if(!allBitsSet(par->dataFlags, DS_mask_abundance)){
      if(bitIsSet(defaultFuncFlags, USERFUNC_abundance)){
        if(bitIsSet(defaultFuncFlags, USERFUNC_molNumDensity)){
          if(!silent) bail_out("You must provide either an abundance() or a molNumDensity() function.");
exit(1);
        }

        par->useAbun = 0;

      }else{
        if(!bitIsSet(defaultFuncFlags, USERFUNC_molNumDensity)){
          if(!silent) warning("abundance() function takes precendence, molNumDensity() ignored.");
        }

        par->useAbun = 1;
      }
    }

    if(!allBitsSet(par->dataFlags, DS_mask_turb_doppler)){
      if(bitIsSet(defaultFuncFlags, USERFUNC_doppler)){
        if(!silent) bail_out("You need to supply a doppler() function.");
exit(1);
      }
    }

    if(!allBitsSet(par->dataFlags, DS_mask_velocity)){
      if(bitIsSet(defaultFuncFlags, USERFUNC_velocity)){
        if(!silent) bail_out("You need to supply a velocity() function.");
exit(1);
      }
    }

    if(!allBitsSet(par->dataFlags, DS_mask_ACOEFF)){
      if(bitIsSet(defaultFuncFlags, USERFUNC_velocity)){
        if(!silent) warning("There were no edge velocities in the file, and you haven't supplied a velocity() function.");
      }
    }
  } /* End if par->doMolCalcs */

  /* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
Generate the grid point locations.
  */
  if(!anyBitSet(par->dataFlags, DS_mask_x)){ /* This should only happen if we did not read a file. Generate the grid point locations. */
    mallocAndSetDefaultGrid(gp, (size_t)par->ncell, (size_t)par->nSpecies);

    outRandDensities = malloc(sizeof(double   )*par->pIntensity); /* Not used at present; and in fact they are not useful outside this routine, because they are not the values of the physical density at that point, just what densityFunc3D() returns, which is not necessarily the same thing. */
    outRandLocations = malloc(sizeof(*outRandLocations)*par->pIntensity);

    randGen = gsl_rng_alloc(ranNumGenType);	/* Random number generator */
    if(fixRandomSeeds)
      gsl_rng_set(randGen,342971);
    else
      gsl_rng_set(randGen,time(0));

    if(par->samplingAlgorithm==0){
      randomsViaRejection(par, (unsigned int)par->pIntensity, randGen, outRandLocations);

    } else if(par->samplingAlgorithm==1){
      setConstDefaults(&rinc);

      if(fixRandomSeeds)
        rinc.randSeed = 342971;
      else
        rinc.randSeed = time(0);

      rinc.numDims = DIM;
      rinc.par = *par;
      rinc.desiredNumPoints = (unsigned int)par->pIntensity;
      for(di=0;di<DIM;di++){
        rinc.wholeFieldOrigin[di] = -par->radius;
        rinc.wholeFieldWidth[di] = 2.0*par->radius;
      }
      rinc.verbosity = 0;
      rinc.monitorFunc = NULL;

      rinc.totalNumHighPoints = par->numGridDensMaxima;

      if(rinc.totalNumHighPoints>0){
        rinc.allHighPointLoc   = malloc(sizeof(*(rinc.allHighPointLoc  ))*rinc.totalNumHighPoints);
        rinc.allHighPointDensy = malloc(sizeof(*(rinc.allHighPointDensy))*rinc.totalNumHighPoints);
        for(i=0;i<rinc.totalNumHighPoints;i++){
          for(di=0;di<rinc.numDims;di++){
            rinc.allHighPointLoc[i][di] = par->gridDensMaxLoc[i][di];
          }
          rinc.allHighPointDensy[i] = par->gridDensMaxValues[i];
        }
      }else{
        rinc.allHighPointLoc = NULL;
        rinc.allHighPointDensy = NULL;
      }

      treeGenerateRandoms(&rinc, gridDensity, outRandLocations, outRandDensities);

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
      theta=gsl_rng_uniform(randGen)*2*M_PI;

      if(DIM==3) {
        z=2*gsl_rng_uniform(randGen)-1.;
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
    gsl_rng_free(randGen);

    if(par->samplingAlgorithm==0){
      smooth(par,*gp);
      if(!silent) printDone(5);
    }

    par->dataFlags |= DS_mask_1;
  }

  if(onlyBitsSet(par->dataFlags, DS_mask_1)) /* Only happens if (i) we read no file and have constructed this data within LIME, or (ii) we read a file at dataStageI==1. */
    writeGridIfRequired(par, *gp, NULL, 1);

  /* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
Generate the remaining values if needed. **Note** that we check a few of them to make sure the user has set the appropriate values.
  */
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
    writeGridIfRequired(par, *gp, NULL, 2);

  if(!allBitsSet(par->dataFlags, DS_mask_density)){
    /* Note that we have checked in parseInput() that the user has defined sufficient values. */
    for(i=0;i<par->ncell; i++)
      (*gp)[i].dens = malloc(sizeof(double)*par->numDensities);
    for(i=0;i<par->pIntensity;i++)
      density((*gp)[i].x[0],(*gp)[i].x[1],(*gp)[i].x[2],(*gp)[i].dens);
    for(i=par->pIntensity;i<par->ncell;i++){
      for(j=0;j<par->numDensities;j++)
        (*gp)[i].dens[j]=EPS; //************** what is the low but non zero value for? Probably to make sure no ills happen in case something gets divided by this?
    }

    par->dataFlags |= DS_mask_density;
  }

  if(par->doMolCalcs)
    checkGridDensities(par, *gp); /* Check that none of the density samples is too small. */

  if(!allBitsSet(par->dataFlags, DS_mask_temperatures)){
    if(!bitIsSet(defaultFuncFlags, USERFUNC_temperature)){
      /* Check that the user has defined gas temperatures at least (if the dust temp was not defined, it is taken to be the same as the gas temp).
      */
      dummyT[0] = -1.0; /* a non-physical temperature. */
      temperature(0.0,0.0,0.0, dummyT);
      if(dummyT[0]<0.0){
        if(!silent) bail_out("You need to set gas temperatures in your model.");
exit(1);
      }
    }

    for(i=0;i<par->pIntensity;i++)
      temperature((*gp)[i].x[0],(*gp)[i].x[1],(*gp)[i].x[2],(*gp)[i].t);
    for(i=par->pIntensity;i<par->ncell;i++){
      (*gp)[i].t[0]=par->tcmb;
      (*gp)[i].t[1]=par->tcmb;
    }

    par->dataFlags |= DS_mask_temperatures;
  }

  if(onlyBitsSet(par->dataFlags, DS_mask_3)) /* Only happens if (i) we read no file and have constructed this data within LIME, or (ii) we read a file at dataStageI==3. */
    writeGridIfRequired(par, *gp, NULL, 3); /* Sufficient information for a continuum image. */

  if(par->doMolCalcs){
    if(!allBitsSet(par->dataFlags, DS_mask_abundance)){
      /* Means we didn't read abun values from file, we have to calculate them via the user-supplied fuction. */
      dummyPointer = malloc(sizeof(*dummyPointer)*par->nSpecies);
      if(par->useAbun){
        if(!bitIsSet(defaultFuncFlags, USERFUNC_abundance)){
          /* Check that the user set reasonable values for all species.
          */
          for(si=0;si<par->nSpecies;si++)
            dummyPointer[si] = -1.0; /* non-physical values. */
          abundance(0.0,0.0,0.0, dummyPointer);
          for(si=0;si<par->nSpecies;si++){
            if(dummyPointer[si]<0.0){
              if(!silent) bail_out("You need to set abundances for all species in your model.");
exit(1);
            }
          }
        }

        for(i=0;i<par->pIntensity;i++){
          abundance((*gp)[i].x[0],(*gp)[i].x[1],(*gp)[i].x[2],dummyPointer);
          for(si=0;si<par->nSpecies;si++)
            (*gp)[i].mol[si].abun = dummyPointer[si];
        }
        for(i=par->pIntensity;i<par->ncell;i++){
          for(si=0;si<par->nSpecies;si++)
            (*gp)[i].mol[si].abun = 0.0;
        }
      }else{
        if(!bitIsSet(defaultFuncFlags, USERFUNC_molNumDensity)){
          /* Check that the user set reasonable values for all species.
          */
          for(si=0;si<par->nSpecies;si++)
            dummyPointer[si] = -1.0; /* non-physical values. */
          molNumDensity(0.0,0.0,0.0, dummyPointer);
          for(si=0;si<par->nSpecies;si++){
            if(dummyPointer[si]<0.0){
              if(!silent) bail_out("You need to set molNumDensity for all species in your model.");
exit(1);
            }
          }
        }

        for(i=0;i<par->pIntensity;i++){
          molNumDensity((*gp)[i].x[0],(*gp)[i].x[1],(*gp)[i].x[2],dummyPointer);
          for(si=0;si<par->nSpecies;si++)
            (*gp)[i].mol[si].nmol = dummyPointer[si];
        }
        for(i=par->pIntensity;i<par->ncell;i++){
          for(si=0;si<par->nSpecies;si++)
            (*gp)[i].mol[si].nmol = 0.0;
        }
      }
      free(dummyPointer);

      par->dataFlags |= DS_mask_abundance;
    }

    if(!allBitsSet(par->dataFlags, DS_mask_turb_doppler)){
      if(!bitIsSet(defaultFuncFlags, USERFUNC_doppler)){
        /* Check that the user set reasonable values.
        */
        dummyScalar = -1.0; /* a non-physical value. */
        doppler(0.0,0.0,0.0, &dummyScalar);
        if(dummyScalar<0.0){
          if(!silent) bail_out("You need to set gas turbulence doppler values in your model.");
exit(1);
        }
      }

      for(i=0;i<par->pIntensity;i++)
        doppler((*gp)[i].x[0],(*gp)[i].x[1],(*gp)[i].x[2],&(*gp)[i].dopb_turb);	
      for(i=par->pIntensity;i<par->ncell;i++)
        (*gp)[i].dopb_turb=0.;

      par->dataFlags |= DS_mask_turb_doppler;
    }

    if(!allBitsSet(par->dataFlags, DS_mask_velocity)){
      /* There seems to be no way we can test if the user has set velocities properly because -ve component values are of course possible. */
      for(i=0;i<par->pIntensity;i++)
        velocity((*gp)[i].x[0],(*gp)[i].x[1],(*gp)[i].x[2],(*gp)[i].vel);

      /* Set velocity values also for sink points (otherwise Delaunay ray-tracing has problems) */
      for(i=par->pIntensity;i<par->ncell;i++)
        velocity((*gp)[i].x[0],(*gp)[i].x[1],(*gp)[i].x[2],(*gp)[i].vel);

      par->dataFlags |= DS_mask_velocity;
    }

    if(!allBitsSet(par->dataFlags, DS_mask_ACOEFF)){
      if(!bitIsSet(defaultFuncFlags, USERFUNC_velocity)){
        getEdgeVelocities(par,*gp); /* Mallocs and sets .v1, .v2, .v3, which are only used within calculateJBar(), which is only called if par->doMolCalcs. This also sets par->edgeVelsAvailable. */

        par->dataFlags |= DS_mask_ACOEFF;
      }
    }
  } /* End if(par->doMolCalcs) */

  if(!allBitsSet(par->dataFlags, DS_mask_magfield)){
    if(par->polarization){
      /* There seems to be no way we can test if the user has set B field values properly because -ve component values are of course possible. */
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

  if(onlyBitsSet(par->dataFlags & DS_mask_all_but_mag, DS_mask_4)) /* Only happens if (i) we read no file and have constructed this data within LIME, or (ii) we read a file at dataStageI==4. */
    writeGridIfRequired(par, *gp, NULL, 4);

  dumpGrid(par,*gp);
  free(dc);
}

