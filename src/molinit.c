/*
 *  molinit.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015 The LIME development team
 *
TODO:
	- Get rid of, or regularize somehow, the printf statements - and clean up all the other new messages which are going to dick with the stdout when curses are selected? (Sigh.)
 */

#include "lime.h"

char *collpartnames[] = {"H2","p-H2","o-H2","electrons","H","He","H+"}; /* definition from LAMDA */

void calcMolCMBs(configInfo *par, molData *md){
  int si, iline;

  for(si=0;si<par->nSpecies;si++){
    md[si].cmb	     = malloc(sizeof(double)*md[si].nline);
    md[si].local_cmb = malloc(sizeof(double)*md[si].nline);

    /* fix the normalization at 230GHz */
    md[si].norm=planckfunc(0,par->tcmb,md,0);
    md[si].norminv=1./md[si].norm;
    for(iline=0;iline<md[si].nline;iline++){
      if(par->tcmb>0.)
        md[si].cmb[iline] = planckfunc(iline,par->tcmb,md,si)/md[si].norm;
      else
        md[si].cmb[iline]=0.;

      md[si].local_cmb[iline] = planckfunc(iline,2.728,md,si)/md[si].norm;
    }
  }
}

double
planckfunc(int iline, double temp, molData *md, int s){
  double bb=10.,wn;
  if(temp<eps) bb = 0.0;
  else {
    wn=md[s].freq[iline]/CLIGHT;
    if (HPLANCK*md[s].freq[iline]>100.*KBOLTZ*temp) 
      bb=2.*HPLANCK*wn*wn*md[s].freq[iline]*exp(-HPLANCK*md[s].freq[iline]/KBOLTZ/temp);
    else 
      bb=2.*HPLANCK*wn*wn*md[s].freq[iline]/(exp(HPLANCK*md[s].freq[iline]/KBOLTZ/temp)-1);
  }
  return bb;
}

void readMolData(configInfo *par, molData *md, int **allUniqueCollPartIds, int *numCollPartsFound){
  /* NOTE! allUniqueCollPartIds is malloc'd in the present function, but not freed. The calling program must free it elsewhere.
  */
  int i,j,k,ilev,idummy,iline,numPartsAcceptedThisMol,ipart,collPartId,itemp,itrans;
  double dummy;
  _Bool cpFound,previousCpFound;
  const int sizeI=200;
  char string[sizeI], specref[90], partstr[90];
  FILE *fp;

  *allUniqueCollPartIds = malloc(sizeof(**allUniqueCollPartIds)*MAX_N_COLL_PART);
  *numCollPartsFound = 0;

  for(i=0;i<par->nSpecies;i++){
    if((fp=fopen(par->moldatfile[i], "r"))==NULL) {
      if(!silent) bail_out("Error opening molecular data file");
      exit(1);
    }

    /* Read the header of the data file */
    fgets(string, sizeI, fp);
    fgets(specref, 90, fp);
    fgets(string, sizeI, fp);
    fscanf(fp, "%lf\n", &md[i].amass);
    fgets(string, sizeI, fp);
    fscanf(fp, "%d\n", &md[i].nlev);
    fgets(string, sizeI, fp);

    md[i].amass *= AMU;

    md[i].eterm=malloc(sizeof(double)*md[i].nlev);
    md[i].gstat=malloc(sizeof(double)*md[i].nlev);

    /* Read the level energies and statistical weights */
    for(ilev=0;ilev<md[i].nlev;ilev++){
      fscanf(fp, "%d %lf %lf", &idummy, &md[i].eterm[ilev], &md[i].gstat[ilev]);
      fgets(string, sizeI, fp);
    }

    /* Read the number of transitions and allocate array space */
    fgets(string, sizeI, fp);
    fscanf(fp, "%d\n", &md[i].nline);
    fgets(string, sizeI, fp);

    md[i].lal     = malloc(sizeof(int)   *md[i].nline);
    md[i].lau     = malloc(sizeof(int)   *md[i].nline);
    md[i].aeinst  = malloc(sizeof(double)*md[i].nline);
    md[i].freq    = malloc(sizeof(double)*md[i].nline);
    md[i].beinstu = malloc(sizeof(double)*md[i].nline);
    md[i].beinstl = malloc(sizeof(double)*md[i].nline);

    /* Read transitions, Einstein A, and frequencies */
    for(iline=0;iline<md[i].nline;iline++){
      fscanf(fp, "%d %d %d %lf %lf %lf\n", &idummy, &md[i].lau[iline]\
        , &md[i].lal[iline], &md[i].aeinst[iline], &md[i].freq[iline], &dummy);
      md[i].freq[iline]*=1e9;
      md[i].lau[iline]-=1;
      md[i].lal[iline]-=1;
    }

    /* Calculate Einsten B's */
    for(iline=0;iline<md[i].nline;iline++){
      /*		md[i].freq[iline]=(md[i].eterm[md[i].lau[iline]]-md[i].eterm[md[i].lal[iline]])*100*CLIGHT; */
      md[i].beinstu[iline]=md[i].aeinst[iline]*(CLIGHT/md[i].freq[iline])*(CLIGHT/md[i].freq[iline])/(HPLANCK*md[i].freq[iline])/2.;
      md[i].beinstl[iline]=md[i].gstat[md[i].lau[iline]]/md[i].gstat[md[i].lal[iline]]*md[i].beinstu[iline];
    }

    /* Collision rates below here */
    fgets(string, sizeI, fp);
    fscanf(fp,"%d\n", &md[i].npart);

    md[i].part = malloc(sizeof(*(md[i].part))*md[i].npart);

    /* Not all the collision partners listed in the moldata file may have associated density functions. Those which don't can play no role and should therefore be ignored. We will try not to store them in md, although due to the demands of backward-compatibility, this will sometimes not be possible, e.g. if the user has not set values for par->collPartIds and if at the same time there are fewer density functions than the total number of collision partners specified in the moldata files. To cover these cases we introduce a new struct cpData attribute: densityIndex, with default value -1 signalling that there is no density function for the associated CP.
    */
    k = 0; /* Index to only those CPs which are found to be associated with a density function. */
    for(ipart=0;ipart<md[i].npart;ipart++){
      fgets(string, sizeI, fp);
      fscanf(fp,"%d\n", &collPartId);

      /* We want to test if the comment after the coll partner ID number is longer than the buffer size. To do this, we write a character - any character, as long as it is not \0 - to the last element of the buffer before reading into it:
      */
      string[sizeof(string)-1] = 'x';
      if(fgets(string, sizeI, fp)==NULL){
        if(!silent) bail_out("Read of collision-partner comment line failed.");
        exit(1);
      } else{
        if(string[sizeof(string)-1]=='\0' && string[sizeof(string)-2]!='\n'){
          /* The presence now of a final \0 means the comment string was either just long enough for the buffer, or too long; the absence of \n in the 2nd-last place means it was too long.
          */
          if(!silent) bail_out("Collision-partner comment line is too long.");
          exit(1);
        }
      }

      /* Look for this CP in par->collPartIds
      */
      cpFound = 0;
      if(par->collPartIds!=NULL){
        for(j=0;j<par->numDensities;j++)
          if(collPartId==par->collPartIds[j]) cpFound = 1;
      }

      if(par->collPartIds==NULL || cpFound){
        /* Check to see if we have encountered collPartId already in a previous moldata file. If not, add it to the list of unique coll parts.
        */
        j = 0;
        previousCpFound = 0;
        while(j<(*numCollPartsFound) && !previousCpFound){
          if(collPartId==(*allUniqueCollPartIds)[j])
            previousCpFound = 1;
          j++;
        }

        if(!previousCpFound){
          if((*numCollPartsFound)>=MAX_N_COLL_PART){
            if(!silent) bail_out("Unrecognized collision partner ID found.");
            exit(1);
          }

          (*allUniqueCollPartIds)[*numCollPartsFound] = collPartId;
          (*numCollPartsFound)++;
        }

        md[i].part[k].collPartId = collPartId;
        md[i].part[k].densityIndex = -1; /* Default, signals that there is no density function for this CP. */

        if(par->lte_only){
          readDummyCollPart(fp, sizeI);
          md[i].part[k].ntemp  = -1;
          md[i].part[k].ntrans = -1;
          md[i].part[k].down  = NULL;
          md[i].part[k].temp  = NULL;
          md[i].part[k].lcl   = NULL;
          md[i].part[k].lcu   = NULL;

        }else{ /* Add the CP data to md[i].part, we will need it to solve the population levels. */
          fgets(string, sizeI, fp);
          fscanf(fp,"%d\n", &md[i].part[k].ntrans);
          fgets(string, sizeI, fp);
          fscanf(fp,"%d\n", &md[i].part[k].ntemp);
          fgets(string, sizeI, fp);

          md[i].part[k].temp = malloc(sizeof(double)*md[i].part[k].ntemp);
          md[i].part[k].lcl  = malloc(sizeof(int)   *md[i].part[k].ntrans);
          md[i].part[k].lcu  = malloc(sizeof(int)   *md[i].part[k].ntrans);

          for(itemp=0;itemp<md[i].part[k].ntemp;itemp++){
            fscanf(fp, "%lf", &md[i].part[k].temp[itemp]);
          }

          fscanf(fp,"\n");
          fgets(string, sizeI, fp);

          md[i].part[k].down = malloc(sizeof(double)\
            *md[i].part[k].ntrans*md[i].part[k].ntemp);

          for(itrans=0;itrans<md[i].part[k].ntrans;itrans++){
            fscanf(fp, "%d %d %d", &idummy, &md[i].part[k].lcu[itrans], &md[i].part[k].lcl[itrans]);
            md[i].part[k].lcu[itrans]-=1;
            md[i].part[k].lcl[itrans]-=1;
            for(itemp=0;itemp<md[i].part[k].ntemp;itemp++){
              j = itrans*md[i].part[k].ntemp+itemp;
              fscanf(fp, "%lf", &md[i].part[k].down[j]);
              md[i].part[k].down[j] /= 1.0e6;
            }
            fscanf(fp,"\n");
          }
        } /* End if(par->lte_only) */

        k++;
      } /* End if CP found in par->collPartIds. */
    } /* End loop over collision partners this molecule. */
    numPartsAcceptedThisMol = k;

    if(numPartsAcceptedThisMol!=md[i].npart){
      md[i].part = realloc(md[i].part, sizeof(*(md[i].part))*numPartsAcceptedThisMol);
      md[i].npart = numPartsAcceptedThisMol;
    }

    /* Print out collision partner information.
    */
    strcpy(partstr, collpartnames[md[i].part[0].collPartId-1]);
    for(ipart=1;ipart<md[i].npart;ipart++){
      strcat( partstr, ", ");
      strcat( partstr, collpartnames[md[i].part[ipart].collPartId-1]);
    }
    if(!silent) {
      collpartmesg(specref, md[i].npart);
      collpartmesg2(partstr, ipart);
      collpartmesg3(par->numDensities, 0);
    }

    fclose(fp);
  } /* end loop over molecule index i */

  if((*numCollPartsFound)<=0){
    if(!silent) bail_out("No recognized collision partners read from file.");
    exit(1);
  }
}

void setUpDensityAux(configInfo *par, int *allUniqueCollPartIds, const int numUniqueCollParts){
  /*
The present function, which needs to be called only if we have to calculate the energy level populations at the grid points, deals with the user-settable vectors par->collPartIds and par->nMolWeights. The former of these is used to associate density values with collision-partner species, and the latter is used in converting, for each radiating species, its abundance to a number density, stored respectively in the grid struct attributes abun and nmol. The function deals specifically with the case in which the user has either not set par->collPartIds or par->nMolWeights at all (which they may choose to do), or has set the incorrectly. In either case the respective parameter will have been freed and set to NULL in checkUserDensWeights(). The function tries its best to guess likely values for the parameters, in line with the algorithm used in the code before par->collPartIds and par->nMolWeights were introduced.
  */
  int i;

  if(par->collPartIds==NULL){
    /*
To preserve backward compatibility I am going to try to make the same guesses as were made before par->collPartIds was introduced, but this is made tricky by the fact that the switch block in the previous code did not cover all possibilities. I'm going to add some warnings too. We want users to be able to run their old model.c files for as long as possible, but at the same time urge them to make use of the new facility for specifying par->collPartIds.
    */
    if(par->numDensities > numUniqueCollParts){
      if(!silent) bail_out("Too many density profiles defined.");
      exit(1);
    }

    if(numUniqueCollParts==1){
      if(allUniqueCollPartIds[0]==CP_H2 || allUniqueCollPartIds[0]==CP_p_H2 || allUniqueCollPartIds[0]==CP_o_H2){
        par->collPartIds = malloc(sizeof(int)*par->numDensities); /* par->numDensities must ==1 at this point. */
        par->collPartIds[0] = allUniqueCollPartIds[0];
      }else{
        if(!silent) bail_out("No H2 collision partner, and user didn't set par.collPartIds.");
        exit(1);
      }

    }else if(numUniqueCollParts==2){
      if((allUniqueCollPartIds[0]==CP_p_H2 && allUniqueCollPartIds[1]==CP_o_H2)\
      || (allUniqueCollPartIds[1]==CP_p_H2 && allUniqueCollPartIds[0]==CP_o_H2)){
        par->collPartIds = malloc(sizeof(int)*par->numDensities);
        for(i=0;i<par->numDensities;i++) /* At this point par->numDensities can be only ==1 (previously signalled via 'flag') or ==2. */
          par->collPartIds[i] = allUniqueCollPartIds[i];

        if(par->numDensities==1 && !silent) warning("Calculating molecular density with respect to first collision partner only.");

      }else{
        if(!silent) bail_out("No H2 collision partners, and user didn't set par.collPartIds.");
        exit(1);
      }

    }else if(numUniqueCollParts==par->numDensities){ /* At this point, numUniqueCollParts must be >2. */
      par->collPartIds = malloc(sizeof(int)*par->numDensities);
      for(i=0;i<par->numDensities;i++)
        par->collPartIds[i] = allUniqueCollPartIds[i];

      if(!silent) warning("Calculating molecular density with respect to first two collision partners only.");

    }else{ /* numUniqueCollParts>2 && par->numDensities<numUniqueCollParts */
      if(!silent) bail_out("More than 2 collision partners, but number of density returns doesn't match.");
      exit(1);
    }

    if(!silent) {
      warning("User didn't set par.collPartIds, I'm having to guess them. Guessed:");
#ifdef NO_NCURSES
      for(i=0;i<par->numDensities;i++){
        printf("  Collision partner %d assigned code %d (=%s)\n", i, par->collPartIds[i], collpartnames[par->collPartIds[i]-1]);
      }
      printf("\n");
#endif
    }

    if(par->nMolWeights==NULL){
      /*
The same backward-compatible guesses are made here as for par->collPartIds in the foregoing section of code. I've omitted warnings and errors because they have already been issued during the treatment of par->collPartIds.
      */
      if(numUniqueCollParts==1){
        if(allUniqueCollPartIds[0]==CP_H2 || allUniqueCollPartIds[0]==CP_p_H2 || allUniqueCollPartIds[0]==CP_o_H2){
          par->nMolWeights = malloc(sizeof(double)*par->numDensities);
          par->nMolWeights[0] = 1.0;
        }

      }else if(numUniqueCollParts==2){
        if((allUniqueCollPartIds[0]==CP_p_H2 && allUniqueCollPartIds[1]==CP_o_H2)\
        || (allUniqueCollPartIds[1]==CP_p_H2 && allUniqueCollPartIds[0]==CP_o_H2)){
          par->nMolWeights = malloc(sizeof(double)*par->numDensities);
          for(i=0;i<par->numDensities;i++) /* At this point par->numDensities can be only ==1 (previously signalled via 'flag') or ==2. */
            par->nMolWeights[i] = 1.0;
        }

      }else if(numUniqueCollParts==par->numDensities){ /* At this point, numUniqueCollParts must be >2. */
        par->nMolWeights = malloc(sizeof(double)*par->numDensities);
        for(i=0;i<par->numDensities;i++)
          par->nMolWeights[i] = 0.0;
        par->nMolWeights[0] = 1.0;
        par->nMolWeights[1] = 1.0;
      }

    }else{
      if(!silent){
        warning("Your choices for par.nMolWeights have been let stand, but it");
        warning("is risky to set them without also setting par.collPartIds.");
      }
    }

  }else if(par->nMolWeights==NULL){ /* We get here only if the user has not supplied these values (or not supplied the right number of them) in their model.c. */
    if(par->numDensities==1){
      if(par->collPartIds[0]==CP_H2 || par->collPartIds[0]==CP_p_H2 || par->collPartIds[0]==CP_o_H2){
        par->nMolWeights = malloc(sizeof(double)*par->numDensities);
        par->nMolWeights[0] = 1.0;
      }else{
        if(!silent) bail_out("No H2 collision partner, and user didn't set par.nMolWeights.");
        exit(1);
      }

    }else if(par->numDensities==2){
      if((par->collPartIds[0]==CP_p_H2 && par->collPartIds[1]==CP_o_H2)\
      || (par->collPartIds[1]==CP_p_H2 && par->collPartIds[0]==CP_o_H2)){
        par->nMolWeights = malloc(sizeof(double)*par->numDensities);
        for(i=0;i<par->numDensities;i++){
          par->nMolWeights[i] = 1.0;
        }

      }else{
        if(!silent) bail_out("No H2 collision partners, and user didn't set par.nMolWeights.");
        exit(1);
      }

    }else{ /* par->numDensities>2 */
      par->nMolWeights = malloc(sizeof(double)*par->numDensities);
      for(i=0;i<par->numDensities;i++){
        par->nMolWeights[i] = 0.0;
      }
      par->nMolWeights[0] = 1.0;
      par->nMolWeights[1] = 1.0;
    }

    if(!silent) warning("User didn't set par.nMolWeights, having to guess them.");
  }
}

void assignMolCollPartsToDensities(configInfo *par, molData *md){
  /*
If we have reached this point, par->collPartIds (and par->nMolWeights) should have been malloc'd and filled with sensible values. Here we set up indices which allow us to associate a density function with each collision partner of each radiating molecule. This information is made use of in stateq.c.
  */
  int i,j,ipart;

  for(i=0;i<par->nSpecies;i++){
    for(ipart=0;ipart<md[i].npart;ipart++){
      md[i].part[ipart].densityIndex = -1; /* Default, signals that there is no density function for this CP. */
      for(j=0;j<par->numDensities;j++){
        if(md[i].part[ipart].collPartId==par->collPartIds[j]){
          md[i].part[ipart].densityIndex = j;
        }
      }
      if(md[i].part[ipart].densityIndex==-1){
        if(!silent) bail_out("No density function has been found for molecule/coll. part. combination.");
        exit(1);
      }
    }
  }
}

void readDummyCollPart(FILE *fp, const int strLen){
  char string[strLen];
  int ntrans, ntemp, itemp, itrans, idummy, dummyLcu, dummyLcl;
  double dummyTemp, dummyDown;

  fgets(string, strLen, fp);
  fscanf(fp,"%d\n", &ntrans);
  fgets(string, strLen, fp);
  fscanf(fp,"%d\n", &ntemp);
  fgets(string, strLen, fp);

  for(itemp=0;itemp<ntemp;itemp++)
    fscanf(fp, "%lf", &dummyTemp);

  fscanf(fp,"\n");
  fgets(string, strLen, fp);

  for(itrans=0;itrans<ntrans;itrans++){
    fscanf(fp, "%d %d %d", &idummy, &dummyLcu, &dummyLcl);
    for(itemp=0;itemp<ntemp;itemp++){
      fscanf(fp, "%lf", &dummyDown);
    }
    fscanf(fp,"\n");
  }
}

