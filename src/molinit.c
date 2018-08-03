/*
 *  molinit.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
TODO:
	- Clean up all the new messages which are going to dick with the stdout when curses are selected? (Sigh.)
 */

#include "lime.h"

/*....................................................................*/
void mallocAndSetDefaultMolData(const int nSpecies, molData **md){
  int i;

  (*md)=malloc(sizeof(molData)*nSpecies);
  for(i=0;i<nSpecies;i++){
    (*md)[i].nlev  = -1;
    (*md)[i].nline = -1;
    (*md)[i].npart = -1;
    (*md)[i].amass = -1.0;
    (*md)[i].part    = NULL;
    (*md)[i].lal     = NULL;
    (*md)[i].lau     = NULL;
    (*md)[i].aeinst  = NULL;
    (*md)[i].gir     = NULL;
    (*md)[i].freq    = NULL;
    (*md)[i].beinstu = NULL;
    (*md)[i].beinstl = NULL;
    (*md)[i].eterm   = NULL;
    (*md)[i].gstat   = NULL;
    (*md)[i].cmb     = NULL;
  }
}

/*....................................................................*/
void readDummyCollPart(FILE *fp, const int strLen){
  char string[strLen];
  int ntrans=0,ntemp=0,itemp,itrans,idummy,dummyLcu,dummyLcl;
  double dummyTemp,dummyDown;

  /* Stops compiler warnings when -Wunused-variable */
  (void)dummyTemp;
  (void)dummyDown;
  (void)idummy;
  (void)dummyLcu;
  (void)dummyLcl;
  (void)string;

  checkFgets(fgets(string, strLen, fp), "dcp 0");
  checkFscanf(fscanf(fp,"%d\n", &ntrans), 1, "dcp 1");
  checkFgets(fgets(string, strLen, fp), "dcp 2");
  checkFscanf(fscanf(fp,"%d\n", &ntemp), 1, "dcp 3");
  checkFgets(fgets(string, strLen, fp), "dcp 4");

  for(itemp=0;itemp<ntemp;itemp++)
    checkFscanf(fscanf(fp, "%lf", &dummyTemp), 1, "dcp 5");

  checkFscanf(fscanf(fp,"\n"), 0, "dcp 6");
  checkFgets(fgets(string, strLen, fp), "dcp 7");

  for(itrans=0;itrans<ntrans;itrans++){
    checkFscanf(fscanf(fp, "%d %d %d", &idummy, &dummyLcu, &dummyLcl), 3, "dcp 8");
    for(itemp=0;itemp<ntemp;itemp++){
      checkFscanf(fscanf(fp, "%lf", &dummyDown), 1, "dcp 9");
    }
    checkFscanf(fscanf(fp,"\n"), 0, "dcp 10");
  }
}

/*....................................................................*/
void
checkFirstLineMolDat(FILE *fp, char *moldatfile){
  const int sizeI=200;
  char string[sizeI],message[STR_LEN_0+1];
  char *expectedLine="!MOLECULE";

  if(fgets(string, sizeI, fp)==NULL){
    if(!silent){
      snprintf(message, STR_LEN_0, "Moldat file %s seems to be empty.", moldatfile);
      bail_out(message);
    }
exit(1);
  }

  if(strncmp(string, expectedLine, strlen(expectedLine))!=0){
    if(!silent){
      snprintf(message, STR_LEN_0, "Bad format first line of moldat file %s.", moldatfile);
      bail_out(message);
    }
exit(1);
  }
}

/*....................................................................*/
void readMolData(configInfo *par, molData *md, int **allUniqueCollPartIds, int *numUniqueCollPartsFound){
  /* NOTE! allUniqueCollPartIds is malloc'd in the present function, but not freed. The calling program must free it elsewhere.
  */
  int i,j,k,ilev,idummy,iline,numPartsAcceptedThisMol,ipart,collPartId=-1,itemp,itrans;
  double dummy;
  _Bool cpWasFoundInUserList,previousCpFound;
  const int sizeI=200;
  char string[sizeI],message[STR_LEN_0+1];
  FILE *fp;

  /* Stops compiler warnings when -Wunused-variable */
  (void)idummy;
  (void)dummy;

  *allUniqueCollPartIds = malloc(sizeof(**allUniqueCollPartIds)*MAX_N_COLL_PART);
  *numUniqueCollPartsFound = 0;

  for(i=0;i<par->nSpecies;i++){
    if((fp=fopen(par->moldatfile[i], "r"))==NULL) {
      if(!silent) bail_out("Error opening molecular data file");
exit(1);
    }

    checkFirstLineMolDat(fp, par->moldatfile[i]);

    checkFgets(fgets(md[i].molName, 80, fp), "molName");
    md[i].molName[strcspn(md[i].molName, "\r\n")] = 0;
    checkFgets(fgets(string, sizeI, fp), "blank line");
    checkFscanf(fscanf(fp, "%lf\n", &md[i].amass), 1, "amass");
    checkFgets(fgets(string, sizeI, fp), "blank line");
    checkFscanf(fscanf(fp, "%d\n", &md[i].nlev), 1, "nlev");
    checkFgets(fgets(string, sizeI, fp), "blank line");

    md[i].amass *= AMU;

    md[i].eterm=malloc(sizeof(double)*md[i].nlev);
    md[i].gstat=malloc(sizeof(double)*md[i].nlev);

    /* Read the level energies and statistical weights */
    for(ilev=0;ilev<md[i].nlev;ilev++){
      checkFscanf(fscanf(fp, "%d %lf %lf", &idummy, &md[i].eterm[ilev], &md[i].gstat[ilev]), 3, "eterm, gstat");
      checkFgets(fgets(string, sizeI, fp), "blank line");
    }

    /* Read the number of transitions and allocate array space */
    checkFgets(fgets(string, sizeI, fp), "blank line");
    checkFscanf(fscanf(fp, "%d\n", &md[i].nline), 1, "nline");
    checkFgets(fgets(string, sizeI, fp), "blank line");

    md[i].lal     = malloc(sizeof(int)   *md[i].nline);
    md[i].lau     = malloc(sizeof(int)   *md[i].nline);
    md[i].aeinst  = malloc(sizeof(double)*md[i].nline);
    md[i].freq    = malloc(sizeof(double)*md[i].nline);
    md[i].beinstu = malloc(sizeof(double)*md[i].nline);
    md[i].beinstl = malloc(sizeof(double)*md[i].nline);

    /* Read transitions, Einstein A, and frequencies */
    for(iline=0;iline<md[i].nline;iline++){
      checkFscanf(fscanf(fp, "%d %d %d %lf %lf %lf\n", &idummy, &md[i].lau[iline]\
        , &md[i].lal[iline], &md[i].aeinst[iline], &md[i].freq[iline], &dummy), 6, "lau, lal etc");
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
    checkFgets(fgets(string, sizeI, fp), "blank line");
    checkFscanf(fscanf(fp,"%d\n", &md[i].npart), 1, "npart");

    md[i].part = malloc(sizeof(*(md[i].part))*md[i].npart);

    /* Not all the collision partners listed in the moldata file may have associated density functions. Those which don't can play no role and should therefore be ignored. We will try not to store them in md, although due to the demands of backward-compatibility, this will sometimes not be possible, e.g. if the user has not set values for par->collPartIds and if at the same time there are fewer density functions than the total number of collision partners specified in the moldata files. To cover these cases we introduce a new struct cpData attribute: densityIndex, with default value -1 signalling that there is no density function for the associated CP.
    */
    k = 0; /* Index to only those CPs which are found to be associated with a density function. */
    for(ipart=0;ipart<md[i].npart;ipart++){
      checkFgets(fgets(string, sizeI, fp), "blank line");
      checkFscanf(fscanf(fp,"%d\n", &collPartId), 1, "collPartId integer");

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
          if(!silent){
            snprintf(message, STR_LEN_0, "Collision-partner comment line must be shorter than %d characters.", sizeI-1);
            bail_out(message);
          }
exit(1);
        }
      }

      /* Look for this CP in par->collPartIds. Note that since we have already called collparts.checkUserDensWeights() by this stage, therefore either we have 1 par->collPartIds entry per density function return, or par->collPartIds==NULL.
      */
      cpWasFoundInUserList = 0;
      if(par->collPartIds!=NULL){
        for(j=0;j<par->numDensities;j++)
          if(collPartId==par->collPartIds[j]) cpWasFoundInUserList = 1;

//**** Won't this go wrong when the user has supplied par->collPartNames but not par->collPartIds? See collparts.checkUserDensWeights().
      }

      if(par->collPartIds==NULL || cpWasFoundInUserList){
        /* Check to see if we have encountered collPartId already in a previous moldata file. If not, add it to the list of unique coll parts.
        */
        j = 0;
        previousCpFound = 0;
        while(j<(*numUniqueCollPartsFound) && !previousCpFound){
          if(collPartId==(*allUniqueCollPartIds)[j])
            previousCpFound = 1;
          j++;
        }

        if(!previousCpFound){
          if((*numUniqueCollPartsFound)>=MAX_N_COLL_PART){
            if(!silent){
              snprintf(message, STR_LEN_0, "More than %d unique collision partners found in the moldata files.", MAX_N_COLL_PART);
              bail_out(message);
            }
exit(1);
          }

          (*allUniqueCollPartIds)[*numUniqueCollPartsFound] = collPartId;
          (*numUniqueCollPartsFound)++;
        }

        setCollPartsDefaults(&(md[i].part[k]));
        md[i].part[k].collPartId = collPartId;

        if(par->lte_only){
          readDummyCollPart(fp, sizeI);

        }else{ /* Add the CP data to md[i].part, we will need it to solve the population levels. */
          checkFgets(fgets(string, sizeI, fp), "blank line");
          checkFscanf(fscanf(fp,"%d\n", &md[i].part[k].ntrans), 1, "ntrans");
          checkFgets(fgets(string, sizeI, fp), "blank line");
          checkFscanf(fscanf(fp,"%d\n", &md[i].part[k].ntemp), 1, "ntemp");
          checkFgets(fgets(string, sizeI, fp), "blank line");

          md[i].part[k].temp = malloc(sizeof(double)*md[i].part[k].ntemp);
          md[i].part[k].lcl  = malloc(sizeof(int)   *md[i].part[k].ntrans);
          md[i].part[k].lcu  = malloc(sizeof(int)   *md[i].part[k].ntrans);

          for(itemp=0;itemp<md[i].part[k].ntemp;itemp++){
            checkFscanf(fscanf(fp, "%lf", &md[i].part[k].temp[itemp]), 1, "temp");
          }

          checkFscanf(fscanf(fp,"\n"), 0, "LF");
          checkFgets(fgets(string, sizeI, fp), "blank line");

          md[i].part[k].down = malloc(sizeof(double)\
            *md[i].part[k].ntrans*md[i].part[k].ntemp);

          for(itrans=0;itrans<md[i].part[k].ntrans;itrans++){
            checkFscanf(fscanf(fp, "%d %d %d", &idummy, &md[i].part[k].lcu[itrans], &md[i].part[k].lcl[itrans]), 3, "lcu, lcl");
            md[i].part[k].lcu[itrans]-=1;
            md[i].part[k].lcl[itrans]-=1;
            for(itemp=0;itemp<md[i].part[k].ntemp;itemp++){
              j = itrans*md[i].part[k].ntemp+itemp;
              checkFscanf(fscanf(fp, "%lf", &md[i].part[k].down[j]), 1, "down");
              md[i].part[k].down[j] /= 1.0e6;  /* collision data is assuming that coll.part density is in [cm^-3], but LIME uses [m^-3] */
            }
            checkFscanf(fscanf(fp,"\n"), 0, "LF");
          }
        } /* End if(par->lte_only) */
        k++;
      }else{ /* read and discard to keep the file reading in sync */
        readDummyCollPart(fp, sizeI);
      } /* End if CP found in par->collPartIds. */
    } /* End loop over collision partners this molecule. */
    numPartsAcceptedThisMol = k;

    if(numPartsAcceptedThisMol!=md[i].npart){
      md[i].part = realloc(md[i].part, sizeof(*(md[i].part))*numPartsAcceptedThisMol);
      md[i].npart = numPartsAcceptedThisMol;
    }

    fclose(fp);
  } /* end loop over molecule index i */
}

/*....................................................................*/
void assignMolCollPartsToDensities(configInfo *par, molData *md){
  /*
If we have reached this point, par->collPartIds (and par->nMolWeights) should have been malloc'd and filled with sensible values. Here we set up indices which allow us to associate a density function with each collision partner of each radiating molecule. This information is made use of in solver.c.
  */
  int i,j,ipart;
  char message[STR_LEN_0+1];

  for(i=0;i<par->nSpecies;i++){
    for(ipart=0;ipart<md[i].npart;ipart++){
      md[i].part[ipart].densityIndex = -1; /* Default, signals that there is no density function for this CP. */
      for(j=0;j<par->numDensities;j++){
        if(md[i].part[ipart].collPartId==par->collPartIds[j]){
          md[i].part[ipart].densityIndex = j;
        }
      }
      if(md[i].part[ipart].densityIndex==-1){
        if(!silent){
          snprintf(message, STR_LEN_0, "No density function found for molecule %d coll. part. %d.", i, ipart);
          warning(message);
        }
      }
    }
  }
}

/*....................................................................*/
void _calcMolCMBs(configInfo *par, molData *md){
  int si, iline;

  for(si=0;si<par->nSpecies;si++){
    md[si].cmb	     = malloc(sizeof(double)*md[si].nline);
    for(iline=0;iline<md[si].nline;iline++){
      if(par->tcmb>0.)
        md[si].cmb[iline] = planckfunc(md[si].freq[iline],par->tcmb);
      else
        md[si].cmb[iline]=0.;
    }
  }
}

/*....................................................................*/
void setUpGir(configInfo *par, molData *md){
  int i,ilev,jlev;
  double dummy;
  FILE *fp;
  char message[STR_LEN_0+1];

  for(i=0;i<par->nSpecies;i++){
    md[i].gir = malloc(sizeof(double)*md[i].nlev*md[i].nlev);
    /* Read the pumping rate coefficients onto gir array */
    for(ilev=0;ilev<md[i].nlev;ilev++){
      for(jlev=0;jlev<md[i].nlev;jlev++){
        md[i].gir[ilev*md[i].nlev+jlev] = 0.;
      }
    }
    if((fp=fopen(par->girdatfile[i], "r")) != NULL){
      while (fscanf(fp, "%d %d %lf", &ilev, &jlev, &dummy) != EOF){
        if(ilev<=0 || jlev<=0){
          if(!silent){
            snprintf(message, STR_LEN_0, "Girdat file %d: level <1 found. (Values found on the bad line: %d, %d.)", i, ilev, jlev);
            bail_out(message);
exit(1);
          }
        }
        if(ilev>md[i].nlev || jlev>md[i].nlev){
          if(!silent){
            snprintf(message, STR_LEN_0, "Girdat file %d: level >%d found. (Values found on the bad line: %d, %d.)", i, md[i].nlev, ilev, jlev);
            bail_out(message);
exit(1);
          }
        }
        md[i].gir[(ilev-1)*md[i].nlev+jlev-1] = dummy;
      }
      fclose(fp);
    }
  }
}

/*....................................................................*/
void molInit(configInfo *par, molData *md){
  int i,j,jStart,numActiveCollParts;
  char partstr[90];
  int *allUniqueCollPartIds=NULL;
  int numUniqueCollPartsFound;

  readMolData(par, md, &allUniqueCollPartIds, &numUniqueCollPartsFound);
  setUpDensityAux(par, allUniqueCollPartIds, numUniqueCollPartsFound); /* In collparts.c */
  free(allUniqueCollPartIds);

  if(par->girdatfile!=NULL){
    setUpGir(par, md);
  }
  if(!par->lte_only){
    assignMolCollPartsToDensities(par, md);

    /* Print out collision partner information.
    */
    for(i=0;i<par->nSpecies;i++){
      for(j=0;j<md[i].npart;j++){
        if(md[i].part[j].densityIndex>=0)
          copyInparStr(par->collPartNames[md[i].part[j].densityIndex], &(md[i].part[j].name));
      }

      jStart = 0;
      while(md[i].part[jStart].densityIndex<0) jStart++;

      if(jStart<md[i].npart){
        numActiveCollParts = 1;
        strcpy(partstr, md[i].part[jStart].name);
        for(j=jStart+1;j<md[i].npart;j++){
          if(md[i].part[j].densityIndex<0)
        continue;

          strcat( partstr, ", ");
          strcat( partstr, md[i].part[j].name);
          numActiveCollParts++;
        }

        if(!silent) {
          collpartmesg(md[i].molName, numActiveCollParts);
          collpartmesg2(partstr);
          collpartmesg3(par->numDensities, 0);//**************** was the 2nd arg used in lime-1.5??
        }
      }else{
        if(!silent) {
          collpartmesg(md[i].molName, 0);
          collpartmesg3(par->numDensities, 0);//**************** was the 2nd arg used in lime-1.5??
        }
      }
    } /* end loop over molecule index i */
  }

  _calcMolCMBs(par, md);
}

