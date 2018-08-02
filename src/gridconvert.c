/*
 *  gridconvert.c
 *  This is an auxiliary program to convert LIME grid files to and from ASCII/FITS.
 *
 *  See ../COPYRIGHT
 *
Examples of valid patterns of options:

  gridconvert -fp -m 'hco+@xpol.dat'                               grid_5.ds test_pops.pop
  gridconvert -f                                                   grid_5.ds test_pregrid.asc
  gridconvert -p  -n 4000 -b 3000 -t 2.725                         test_pops.pop    test_fits_p.ds
  gridconvert     -n 4000 -b 3000 -t 2.725 -r 2.991957e+14 -c 'H2' test_pregrid.asc test_fits_a.ds
 */
#include <locale.h>
#include <argp.h>

#include "lime.h"
#include "gridio.h"

char *programName = "gridconvert";

int silent = 0;
int defaultFuncFlags = 0;
double defaultDensyPower = DENSITY_POWER;

#ifdef TEST
_Bool fixRandomSeeds = TRUE;
#else
_Bool fixRandomSeeds = FALSE;
#endif

const char *argp_program_version = VERSION;
const char *argp_program_bug_address = "https://github.com/lime-rt/lime";
/* Program documentation. */
static char doc[] = "gridconvert - a program to convert between LIME grid-file formats.";

/* A description of the arguments we accept. */
static char args_doc[] = "infilename outfilename";

/* The options we understand. */
static struct argp_option options[] = {
  {"fitsin",         'f',    0,    0, "Set if the input is FITS; left unset means the output is FITS." },
  {"popsnotpre",     'p',    0,    0, "Set if the non-FITS is pops format; left unset means it is pregrid format." },
  {"collpart",       'c', "str",   0, "Name of the bulk gas constituent." },
  {"ninternal",      'n', "int",   0, "Number of model points inside the spherical boundary." },
  {"nboundary",      'b', "int",   0, "Number of model points on the spherical boundary." },
  {"modelradius",    'r', "float", 0, "Radius (m) of the spherical model boundary." },
  {"cmbtemperature", 't', "float", 0, "Temperature (K) of the cosmic microwave background." },
  {"moldatfile",     'm', "str",   0, "Name of a LIME molecule data file." },
  { 0 }
};

/* Used by main to communicate with parse_opt. */
struct arguments
{
  _Bool fitsIsTheInput,popsNotPre;
  char *inFile,*outFile,*bulkSpeciesName,*moldatfile;
  int numInternalPoints,numBoundaryPoints;
  double modelRadius,tempCMB;
};

/*....................................................................*/
void myPrintUsage(){
  const int maxNumOptions=20,maxOptFieldLen=28;
  int i,numSpaces;
  char optionsStr[maxNumOptions+1],strBuffer[1024];

  for(i=0;i<maxNumOptions;i++){
    if(options[i].name==0)
  break;
    optionsStr[i] = (char)options[i].key;
  }
  if(i>=maxNumOptions && options[i].name!=0){
    warning("Bug! You need to increase maxNumOptions");
  }
  optionsStr[i] = '\0';

  printf("%s [-%s] <infile name> <outfile name>\n\n", programName, optionsStr);
  printf("OPTIONS SUMMARY:\n");

  for(i=0;i<maxNumOptions;i++){
    if(options[i].name==0)
  break;

    sprintf(strBuffer, " -%c, --%s", (char)options[i].key, options[i].name);
    numSpaces = maxOptFieldLen - strlen(strBuffer);
    if(numSpaces<1) numSpaces = 1;

    printf(" %s%*c%s\n", strBuffer, numSpaces, ' ', options[i].doc);
  }
}

/*....................................................................*/
/* Parse a single option. */
static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{
  /* Get the input argument from argp_parse, which we
     know is a pointer to our arguments structure. */
  struct arguments *arguments = state->input;

  switch (key)
    {
    case 'f':
      arguments->fitsIsTheInput = 1;
      break;
    case 'p':
      arguments->popsNotPre = 1;
      break;
    case 'c':
      arguments->bulkSpeciesName = arg;
      break;
    case 'n':
      arguments->numInternalPoints = atoi(arg);
      break;
    case 'b':
      arguments->numBoundaryPoints = atoi(arg);
      break;
    case 'r':
      arguments->modelRadius = atof(arg);
      break;
    case 't':
      arguments->tempCMB = atof(arg);
      break;
    case 'm':
      arguments->moldatfile = arg;
      break;

    case ARGP_KEY_ARG:
      if(state->arg_num == 0)
        arguments->inFile = arg;
      else if(state->arg_num == 1)
        arguments->outFile = arg;
      else
        argp_usage(state);

      break;

    case ARGP_KEY_END:
      break;

    default:
      return ARGP_ERR_UNKNOWN;
    }
  return 0;
}

/* Our argp parser. */
static struct argp argp = { options, parse_opt, args_doc, doc };

/*....................................................................*/
int main(int argc, char **argv) {
  /*
Options:
	-f --fitsin		# If set, the input file is expected to be a FITS grid file; otherwise, the output is expected to be.
	-p --popsnotpre		# Relates to the format of the non-FITS file. If set: it is expected to adhere to the (current) LIME 'predefgrid' format, otherwise it is assumed a 'popsin/popsout' file.
  */

  struct arguments arguments;
  char message[STR_LEN_1];
  FILE *fp;
  struct grid *gp=NULL;
  configInfo par;
  int status=0;
  molData *md=NULL;

  /* Default values. */
  arguments.inFile = NULL;
  arguments.outFile = NULL;
  arguments.fitsIsTheInput = 0;
  arguments.popsNotPre = 0;
  arguments.numInternalPoints = -1;
  arguments.numBoundaryPoints = -1;
  arguments.modelRadius = -1.0;
  arguments.tempCMB = -1.0;
  arguments.bulkSpeciesName = NULL;
  arguments.moldatfile = NULL;

  /* Parse our arguments; every option seen by parse_opt will be reflected in arguments.
  */
  if(argp_parse(&argp, argc, argv, 0, 0, &arguments)){
    bail_out("Argument parser returned with an error status.");
exit(1);
  }

  /* Do some checks of the arguments. I'm doing this here rather than using the argp machinery because the latter does not seem well adjusted to the situation we have here in which the expected number of arguments depends on the values of some of them.
  */
  if(arguments.inFile==NULL || arguments.outFile==NULL){
    snprintf(message, STR_LEN_1, "You must invoke this as '%s <options> <infile> <outfile>", programName);
    bail_out(message);
exit(1);
  }
  if(arguments.fitsIsTheInput){
    if(arguments.popsNotPre){
      if(arguments.moldatfile==NULL){
        sprintf(message, "If you want to write a 'pops' file you need to include option '-m <name of moldatfile>'");
        bail_out(message);
exit(1);
      }
    }
  }else{ /* Output to FITS */
    if(arguments.popsNotPre){
      if(arguments.numInternalPoints<0 || arguments.numBoundaryPoints<0 || arguments.tempCMB<0.0){
        if(arguments.numInternalPoints<0)
          sprintf(message, "If you want to read a 'pops' file you need to include option '-n <num internal pts>'");
        if(arguments.numBoundaryPoints<0)
          sprintf(message, "If you want to read a 'pops' file you need to include option '-b <num boundary pts>'");
        if(arguments.tempCMB<0.0)
          sprintf(message, "If you want to read a 'pops' file you need to include option '-t <CMB temp (K)>'");

        bail_out(message);
exit(1);
      }
    }else{
      if(arguments.numInternalPoints<0 || arguments.numBoundaryPoints<0 || arguments.modelRadius<0.0 || arguments.tempCMB<0.0 || arguments.bulkSpeciesName==NULL){
        if(arguments.numInternalPoints<0)
          sprintf(message, "If you want to read a pregrid file you need to include option '-n <num internal pts>'");
        if(arguments.numBoundaryPoints<0)
          sprintf(message, "If you want to read a pregrid file you need to include option '-b <num boundary pts>'");
        if(arguments.modelRadius<0.0)
          sprintf(message, "If you want to read a pregrid file you need to include option '-r <model radius (m)>'");
        if(arguments.tempCMB<0.0)
          sprintf(message, "If you want to read a pregrid file you need to include option '-t <CMB temp (K)>'");
        if(arguments.bulkSpeciesName==NULL)
          sprintf(message, "If you want to read a pregrid file you need to include option '-c <name of bulk species>'");

        bail_out(message);
exit(1);
      }
    }
  }

  /* Check that the infile exists:
  */
  if((fp=fopen(arguments.inFile, "r"))==NULL) {
    snprintf(message, STR_LEN_1, "Input file %s not found.", arguments.inFile);
    bail_out(message);
exit(1);
  }
  fclose(fp);

  if(arguments.fitsIsTheInput){
    int numCollPartRead=0,dataFlags=0,i;
    char **collPartNames;
    struct gridInfoType gridInfo;

    status = readGrid(arguments.inFile, &gridInfo, NULL\
      , 0, &gp, &collPartNames, &numCollPartRead, &dataFlags);
    /* We leave keywords==NULL because pregrid doesn't read any additional data, so there is no format defined for writing it. */

    if(status){
      if(!silent){
        snprintf(message, STR_LEN_1, "Read of grid file from FITS failed with status return %d", status);
        bail_out(message);
      }
exit(1);
    }

    par.ncell = (int)(gridInfo.nInternalPoints+gridInfo.nSinkPoints);
    par.nSpecies = 0;

    if(arguments.popsNotPre){
      int *allUniqueCollPartIds,numUniqueCollPartsFound;

      if(!allBitsSet(dataFlags, DS_mask_5)){
        bail_out("Can't write POPS file because there is not enough information in the input file.");
exit(1);
      }

      if(gridInfo.nSpecies>1){
        snprintf(message, STR_LEN_1, "There was data on %d species in your input file, but we can only handle 1.", gridInfo.nSpecies);
        warning(message);
      }
      if(gridInfo.nDensities>1){
        snprintf(message, STR_LEN_1, "There was data on %d coll. parts. in your input file, but we can only handle 1.", gridInfo.nDensities);
        warning(message);
      }

      /* Fill in elements of par that binpopsout() reads:
      */
      par.nSpecies = 1;
      par.moldatfile = malloc(sizeof(*par.moldatfile)*par.nSpecies);
      par.moldatfile[0] = arguments.moldatfile;
      par.collPartIds = NULL;
      par.numDensities = 1;
      par.lte_only = 0;
      par.binoutputfile = arguments.outFile;
      par.radius = arguments.modelRadius;
      par.nMolWeights = malloc(sizeof(*par.nMolWeights)*par.nSpecies);
      par.nMolWeights[0] = 1.0;

      mallocAndSetDefaultMolData(par.nSpecies, &md);

      readMolData(&par, md, &allUniqueCollPartIds, &numUniqueCollPartsFound);

      calcGridMolDoppler(&par, md, gp);
      calcGridMolDensities(&par, &gp);

      binpopsout(&par, gp, md);

      free(allUniqueCollPartIds);
      free(par.moldatfile);
      free(par.nMolWeights);

    }else{
      writeGridToAscii(arguments.outFile, gp, gridInfo.nInternalPoints, dataFlags);
    }

    if(collPartNames != NULL){
      for(i=0;i<numCollPartRead;i++)
        free(collPartNames[i]);
      free(collPartNames);
    }

    freeGrid(par.ncell, gridInfo.nSpecies, gp);
    freeGridInfo(&gridInfo);

  }else{ /* Read ASCII, write to FITS */
    if(arguments.popsNotPre){
      int popsdone=0;

      /* Fill in elements of par that popsin() reads:
      */
      par.restart    = arguments.inFile;
      par.pIntensity = arguments.numInternalPoints;
      par.sinkPoints = arguments.numBoundaryPoints;
      par.tcmb       = arguments.tempCMB;
      par.nSolveItersDone = 0; /* In reality there must have been some, but there is no facility in the pops file for storing this. */
      par.collPartNames = NULL;
      par.dataFlags = 0;

      mallocAndSetDefaultMolData(1, &md); /* popsin reallocs it. */

      popsin(&par, &gp, &md, &popsdone); /* Sets par.nSpecies among other elements. */

    }else{
      /* Fill in elements of par that predefinedGrid() reads:
      */
      par.gridfile = NULL;
      par.pregrid = arguments.inFile;
      par.pIntensity = arguments.numInternalPoints;
      par.sinkPoints = arguments.numBoundaryPoints;
      par.ncell = par.pIntensity + par.sinkPoints;
      par.radius = arguments.modelRadius;
      par.tcmb = arguments.tempCMB;
      par.nSpecies = 1;
      par.dataFlags = 0;
      par.nSolveItersDone = 0;
      par.collPartNames = malloc(sizeof(*par.collPartNames)*1);
      par.collPartNames[0] = malloc(sizeof(char)*STR_LEN_0);
      strcpy(par.collPartNames[0], arguments.bulkSpeciesName);
      par.doMolCalcs = TRUE;

      mallocAndSetDefaultGrid(&gp, (size_t)par.ncell, (size_t)par.nSpecies);
      predefinedGrid(&par, gp);
    }

    status = setupAndWriteGrid(&par, gp, md, arguments.outFile);

    if(status){
      if(!silent){
        snprintf(message, STR_LEN_1, "Write of grid file to FITS failed with status return %d", status);
        bail_out(message);
      }
exit(1);
    }

    freeGrid(par.ncell, par.nSpecies, gp);
  }

  freeMolData(par.nSpecies, md);

  return 0;
}

