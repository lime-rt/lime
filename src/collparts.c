/*
 *  collparts.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *

Some explanation:
=================
We have three sources of information about the density and type of bulk gas (a.k.a. collision partner) species: (i) the number of elements returned by the density() function; (ii) the list par->collPartIds; (iii) and the transition rate values supplied in the moldatfiles. What we want to do is associate density functions with collision partners given in the moldatfiles. We would like to do this if possible without deviating from the original assumptions made in the earliest versions of LIME, and without loading the user down with too much crud.

Below is given a typical timeline of processing within LIME, showing for each stage the questions which the code asks about par->collPartIds and the moldatfiles, and the decisions taken in response.

main.initParImg():
------------------

  - inpars->collPartIds (i.e. a struct of type 'inputPars') is malloc'd to size MAX_N_COLL_PART and all the entries are set to 0. Zero does not correspond with any LAMDA species and can therefore be used to flag that the user has not reset the value to something sensible.

User's model file:
------------------

  - (optionally) some enties of inpars->collPartIds are set to >0 values.

run.parseInput():
-----------------

  par->collPartIds (a struct of type 'configInfo') is malloc'd to size MAX_N_COLL_PART and the inpars values are copied over.

collparts.checkUserDensWeights():
---------------------------------

  - The number of >0 (therefore user-set) values of par->collPartIds are counted. Counting starts from the 0th element and stops as soon as the first <=0 value is found.

  - The number of user-set values of par->collPartIds is compared to the number of density() function returns (by now stored in par->numDensities). If the numbers are unequal, par->collPartIds is freed and set to NULL (a warning is issued if the user set some values). If the numbers match, par->collPartIds is reallocated in size to that number (this is possibly unnecessary).

  - If the user has supplied par->collPartNames, with the number of entries correctly matching the number of density() returns (par->numDensities), in the case that par->collPartIds was set to NULL just above, then par->collPartIds is malloc'd to size par->numDensities; it is loaded then with an increasing sequence of integers.

  - par->collPartIds, if set, are checked for uniqueness.

molinit.readMolData():
----------------------

  - The collision partners listed in the set of moldatfiles are read from file. Here we first want to compile a list of unique (that is, without repetition) collision partners that occur in the files. We will only do this however if either (i) the user has supplied no values of par->collPartIds or (ii) the collision partner matches one supplied in par->collPartIds.

  - Note that ALL collision partners read from the ith file are stored in md[i].part.

collparts.setUpDensityAux():
----------------------------

  - The processing in this looks quite complicated, but the only thing that is done to par->collPartIds is if it ==NULL: it is then malloc'd to size par->numDensities and, under limited circumstances, its values are set to those read from the moldatfiles.

molinit.assignMolCollPartsToDensities():
----------------------------------------

  - par->collPartIds is compared here to the collision partner IDs in md[i].part[ipart].collPartId. If a CP is found in par->collPartIds but not in md[i].part, no action is taken. If any md[i].part[ipart].collPartId is not found in par->collPartIds, md[i].part[ipart].densityIndex is set to -1. This flags to solver.getFixedMatrix() to ignore the respective transition rates.

 */

#include "lime.h"

/*....................................................................*/
void setCollPartsDefaults(struct cpData *part){
  (*part).collPartId = -1;
  (*part).densityIndex = -1; /* Signals that there is no density function for this CP. */
  (*part).name = NULL; /* If it turns out to have a matching density function we will store the name. */
  (*part).ntemp  = -1;
  (*part).ntrans = -1;
  (*part).down  = NULL;
  (*part).temp  = NULL;
  (*part).lcl   = NULL;
  (*part).lcu   = NULL;
}

/*....................................................................*/
void checkUserDensWeights(configInfo *par){
  /*
This deals with four user-settable list parameters which relate to collision partners and their number densities: par->collPartIds, par->nMolWeights, par->collPartMolWeights and par->collPartNames. We have to see if these (optional) parameters were set, do some basic checks on them, and if they were set make sure they match the number of density values, which by this time should be stored in par->numDensities.

	* par->collPartIds: this list acts as a link between the N density function returns (I'm using here N as shorthand for par->numDensities) and the M collision partner ID integers found in the moldatfiles. This allows us to associate density functions with the collision partner transition rates provided in the moldatfiles.

	* par->collPartNames: essentially this has only cosmetic importance since it has no effect on the functioning of LIME, only on the names of the collision partners which are printed to stdout. Its main purpose is to reassure the user who has provided transition rates for a non-LAMDA collision species in their moldatfile that they are actually getting these values and not some mysterious reversion to LAMDA.

	The user can specify either, none, or both of these two parameters, with the following effects:

		Ids	Names	Effect
		----------------------
		0	0	LAMDA collision partners are assumed and the association between the density functions and the moldatfiles is essentially guessed.

		0	1	par->collPartIds is constructed to contain integers in a sequence from 1 to N. Naturally the user should write matching collision partner ID integers in their moldatfiles.

		1	0	LAMDA collision partners are assumed.

		1	1	User will get what they ask for.
		----------------------

	* par->collPartMolWeights: this MUST be present if par->collPartNames has been supplied, and it MUST then have the same number and order of elements as all the other collision-partner lists. If this parameter is supplied but par->collPartNames not, it will be ignored.

	* par->nMolWeights: this list gives the weights to be applied to the N density values when calculating molecular densities from abundances.
  */
  int i,j,numUserSetCPIds,numUserSetNMWs,numUserSetCPNames,numUserSetCPWeights;
  int *uniqueCPIds=NULL;
  double sum;

  par->collPartUserSetFlags = 0;

  /* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .*/
  /* Get the numbers of elements set by the user for each of the 4 parameters:
  */
  i = 0;
  while(i<MAX_N_COLL_PART && par->collPartIds[i]>0) i++;
  numUserSetCPIds = i;
  if(i>0) par->collPartUserSetFlags |= (1<<CPF_BIT_ids);

  i = 0;
  while(i<MAX_N_COLL_PART && par->nMolWeights[i]>=0.0) i++;
  numUserSetNMWs = i;
  if(i>0) par->collPartUserSetFlags |= (1<<CPF_BIT_weights);

  i = 0;
  while(i<MAX_N_COLL_PART && par->collPartNames[i]!=NULL) i++;
  numUserSetCPNames = i;
  if(i>0) par->collPartUserSetFlags |= (1<<CPF_BIT_names);

  i = 0;
  while(i<MAX_N_COLL_PART && par->collPartMolWeights[i]>=0) i++;
  numUserSetCPWeights = i;
  if(i>0) par->collPartUserSetFlags |= (1<<CPF_BIT_MolWeights);

  /* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .*/
  /* Perform checks on the numbers.
  */
  /* Check that we have either 0 par->collPartIds or the same number as the number of density values. If not, the par->collPartIds values the user set will be thrown away, the pointer will be reallocated, and new values will be written to it in setUpDensityAux(), taken from the values in the moldatfiles.
  */
  if(numUserSetCPIds != par->numDensities){
    free(par->collPartIds);
    par->collPartIds = NULL;
    /* Note that in the present case we will (for a line-emission image) look for the collision partners listed in the moldatfiles and set par->collPartIds from them. For that to happen, we require the number of collision partners found in the files to equal par->numDensities. */

    /* numUserSetCPIds==0 is ok, this just means the user has not set the parameter at all, but for other values we should issue some warnings, because if the user sets any at all, they should set the same number as there are returns from density():
    */
    if(numUserSetCPIds > 0){
      if(!silent) warning("par.collPartIds will be ignored - there should be 1 for each density.");
      numUserSetCPIds = 0;
    }

  }else{
    par->collPartIds = realloc(par->collPartIds, sizeof(*(par->collPartIds))*par->numDensities);
  }

  if(!silent && !par->useAbun && numUserSetNMWs>0)
    warning("You only need to set par->nMolWeights if you have provided an abundance function.");

  /* Check if we have either 0 par->nMolWeights or the same number as the number of density values.
  */
  if(par->useAbun && numUserSetNMWs != par->numDensities){
    free(par->nMolWeights);
    par->nMolWeights = NULL;
    /* Note that in the present case we will (for a line-emission image) look for the collision partners listed in the moldatfiles and set par->nMolWeights from them. */

    /* numUserSetNMWs==0 is ok, this just means the user has not set the parameter at all, but for other values we should issue some warnings, because if the user sets any at all, they should set the same number as there are returns from density():
    */
    if(numUserSetNMWs > 0){
      if(!silent) warning("par->nMolWeights will be ignored - there should be 1 for each density() return.");
      numUserSetNMWs = 0;
    }

  }else{
    par->nMolWeights = realloc(par->nMolWeights, sizeof(*(par->nMolWeights))*par->numDensities);
  }

  /* Re the interaction between par->collPartIds and par->collPartNames: the possible scenarios are given in the function header.
  */
  if(numUserSetCPNames != par->numDensities){
    for(i=0;i<MAX_N_COLL_PART;i++) free(par->collPartNames[i]);
    free(par->collPartNames);
    par->collPartNames = NULL;

    /* numUserSetCPNames==0 is ok, this just means the user has not set the parameter at all, but for other values we should issue some warnings, because if the user sets any at all, they should set the same number as there are returns from density():
    */
    if(numUserSetCPNames > 0){
      if(!silent) warning("par->collPartNames will be ignored - there should be 1 for each density() return.");
      numUserSetCPNames = 0;
    }

  }else{ /* If we get to here, then numUserSetCPNames==par->numDensities. */
    for(i=par->numDensities;i<MAX_N_COLL_PART;i++){
      free(par->collPartNames[i]);
    }
    par->collPartNames = realloc(par->collPartNames, sizeof(*(par->collPartNames))*par->numDensities);

    if(numUserSetCPIds<=0){ /* This implies that we must have set par->collPartIds=NULL a little above here. */
      par->collPartIds = malloc(sizeof(*(par->collPartIds))*par->numDensities);
      for(i=0;i<par->numDensities;i++)
        par->collPartIds[i] = i+1;
//**** Won't this go wrong when the collPartIds read from moldatfiles are compared to par->collPartIds in molinit.readMolData()??
      numUserSetCPIds = par->numDensities;
    }
  }

  /* The constraints on the list of CP molecular weights are similar, but NULL + warn that they will be ignored if there are no CP names.
  */
  if(numUserSetCPWeights != par->numDensities || numUserSetCPNames <= 0){
    free(par->collPartMolWeights);
    par->collPartMolWeights = NULL;

    /* numUserSetCPWeights==0 is ok, this just means the user has not set the parameter at all, but for other values we should issue some warnings, because if the user sets any at all, they should set the same number as there are returns from density():
    */
    if(numUserSetCPWeights > 0){
      if(numUserSetCPNames <= 0){
        if(!silent) warning("par->collPartMolWeights will be ignored - you must also set par->collPartNames.");
      }else{
        if(!silent) warning("par->collPartMolWeights will be ignored - there should be 1 for each density() return.");
      }
      numUserSetCPWeights = 0;
    }

  }else{
    par->collPartMolWeights = realloc(par->collPartMolWeights, sizeof(*(par->collPartMolWeights))*par->numDensities);
  }

  /* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .*/
  /* Now we do some sanity checks.
  */
  if(numUserSetCPIds>0){
    /* Check that they are unique.
    */
    uniqueCPIds = malloc(sizeof(int)*numUserSetCPIds);
    for(i=0;i<numUserSetCPIds;i++){
      for(j=0;j<i;j++){
        if(par->collPartIds[i]==uniqueCPIds[j]){
          if(!silent) bail_out("Your list of par.collPartIds is not unique.");
          exit(1);
        }
      }
      uniqueCPIds[i] = par->collPartIds[i];
    }
    free(uniqueCPIds);
  }

  if(par->useAbun && numUserSetNMWs>0){
    /* Check that they do not sum to zero.
    */
    sum = 0.0;
    for(i=0;i<numUserSetNMWs;i++){
      sum += par->nMolWeights[i];
    }
    if(sum<=0.0){
      if(!silent) bail_out("At least some of your par.nMolWeights must be non-zero!");
      exit(1);
    }
  }
}

/*....................................................................*/
void setUpDensityAux(configInfo *par, int *allUniqueCollPartIds, const int numUniqueCollParts){
  /*
The present function, which needs to be called only if we have to calculate the energy level populations at the grid points, deals with the user-settable vectors par->collPartIds and par->nMolWeights. The former of these is used to associate density values with collision-partner species, and the latter is used in converting, for each radiating species, its abundance to a number density, stored respectively in the grid struct attributes abun and nmol. The function deals specifically with the case in which the user has either not set par->collPartIds or par->nMolWeights at all (which they may choose to do), or has set them incorrectly. In either case the respective parameter will have been freed and set to NULL in checkUserDensWeights(). The function tries its best to guess likely values for the parameters, in line with the algorithm used in the code before par->collPartIds and par->nMolWeights were introduced.

allUniqueCollPartIds is the list of all unique CP ID integers detected in the supplied moldata files.
  */
  double lamdaMolWeights[] = {2.0159,2.0159,2.0159,5.486e-4,1.00794,4.0026,1.00739};
  char *lamdaNames[] = {"H2","p-H2","o-H2","electrons","H","He","H+"};
#ifndef NCURSES
  char message[STR_LEN_0+1];
#endif
  int i;

  if(par->collPartIds==NULL){ /* For this to happen means that the user set neither par->collPartIds nor par->collPartNames. */
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
//*** Huh? What if they both ==3? Don't think this warning here is valid.

    }else{ /* numUniqueCollParts>2 && par->numDensities<numUniqueCollParts */
      if(!silent) bail_out("More than 2 collision partners, but number of density returns doesn't match.");
      exit(1);
    }

    if(!silent) {
#ifdef NCURSES
      printMessage("User didn't set par.collPartIds, I'm having to guess them.");
#else
      printMessage("User didn't set par.collPartIds, I'm having to guess them. Guessed:");
      for(i=0;i<par->numDensities;i++){
        snprintf(message, STR_LEN_0, "  Collision partner %d assigned code %d (=%s)\n", i, par->collPartIds[i], lamdaNames[par->collPartIds[i]-1]);
        printMessage(message);
      }
      printMessage("\n");
#endif
    }

    if(par->useAbun){
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
    }

  }else if(par->useAbun && par->nMolWeights==NULL){ /* We get here only if the user has not supplied these values (or not supplied the right number of them) in their model.c. */
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

  /* If we get to here then par->collPartIds has definitely been malloc'd and its values set.
  */
//*** check that no par->collPartIds[i] are out of lamda range.
  if(par->collPartNames==NULL){ /* Then load it from the LAMDA names. */
    par->collPartNames=malloc(sizeof(*par->collPartNames)*par->numDensities);
    for(i=0;i<par->numDensities;i++)
      copyInparStr(lamdaNames[par->collPartIds[i]-1], &(par->collPartNames[i]));
  }

  if(par->collPartMolWeights==NULL){ /* Then load it from the LAMDA mol weights. */
    par->collPartMolWeights=malloc(sizeof(*par->collPartMolWeights)*par->numDensities);
    for(i=0;i<par->numDensities;i++)
      par->collPartMolWeights[i] = lamdaMolWeights[par->collPartIds[i]-1];
  }
}




