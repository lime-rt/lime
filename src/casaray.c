/*
 *  casaray.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2017 The LIME development team
 *
TODO:
 */

#include "lime.h"
#include "casalime.h"
#include "py_utils.h"

int silent = 0;
statusType statusObj;

/*....................................................................*/
static PyObject* py_run(PyObject* self, PyObject* args){
  int status=0,nPars,nImgPars=0,nImages=0,numCollPartRead=0,nEntries=0,i,gi,si,ei;
  parTemplateType *parTemplates=NULL,*imgParTemplates=NULL;
  PyObject *pParClass,*pImgList,*pImgPars;
  inputPars inpars;
  image *inimg = NULL;
  char message[STR_LEN_0],**collPartNamesRead=NULL;
  molData *md=NULL;
  configInfo par;
  imageInfo *img=NULL;
  struct grid *gp=NULL;
  double *lamtab=NULL,*kaptab=NULL;

  /* Unpack the parameter object and check it against the parameter templates.
  */
  if(!PyArg_ParseTuple(args, "O", &pParClass))
return NULL;

  status = getParTemplates(pParClass, &parTemplates, &nPars);
  if(status){
    snprintf(message, STR_LEN_0, "getParTemplates() returned status value %d", status);
    PyErr_SetString(PyExc_ValueError, message);
    free(parTemplates);
return NULL;
  }

  /* Get the number of images and check that the image attributes, for all images, have the correct types.
  */
  pImgList = PyObject_GetAttrString(pParClass, "img");
  if(pImgList==NULL){
    free(parTemplates);
return NULL;
  }

  nImages = (int)PyList_Size(pImgList);
  if(nImages!=1){
    snprintf(message, STR_LEN_0, "You can only process 1 image, not %d", nImages);
    PyErr_SetString(PyExc_ValueError, message);
    free(parTemplates);
    Py_DECREF(pImgList);
return NULL;
  }

  pImgPars = PyList_GetItem(pImgList, (Py_ssize_t)0); /* Don't have to DECREF pImgPars, it is a borrowed reference. */

  status = getParTemplates(pImgPars, &imgParTemplates, &nImgPars);
  if(status){
    free(imgParTemplates);
    free(parTemplates);
    Py_DECREF(pImgList);
return NULL;
  }

  Py_DECREF(pImgList);

  status = mallocInputParStrs(&inpars); /* Note that the inimg equivalents are malloc'd in readParImg(). This separation seems a little messy. */
  if(status){
    PyErr_NoMemory();
return NULL;
  }

  status = readParImg(pParClass, parTemplates, nPars, imgParTemplates, nImgPars, &inpars, &inimg, &nImages, warning);
  /* This mallocs inimg. */

  free(imgParTemplates);
  free(parTemplates);

  if(status){
    snprintf(message, STR_LEN_0, "Function readParImg() returned with status %d.", status);
    PyErr_SetString(PyExc_AttributeError, message);
return NULL;
  }

  /* Calculate any non-user-set configuration parameters and do some sanity checks. This takes in inpars, inimg and converts them to the similar but much more comprehensive par, img. It also initializes the molecular data in md.
  */
  par.nSpecies = copyInpars(inpars, inimg, nImages, &par, &img); /* This mallocs img and some char* attributes of par. */

  pyFreeInputImgPars(inimg, nImages);
  freeInputPars(&inpars);

  setOtherEasyConfigValues(nImages, &par, &img);

  if(par.gridInFile==NULL){
    snprintf(message, STR_LEN_0, "You must define a grid input file.");
    PyErr_SetString(PyExc_ValueError, message);

    freeMolData(par.nSpecies, md);
    freeImgInfo(par.nImages, img);
    freeConfigInfo(&par);
return NULL;
  }

  /* Read the grid data into array-of-structs gp (1 struct per grid point).
  */
  readGridWrapper(&par, &gp, &collPartNamesRead, &numCollPartRead); /* This mallocs gp and collPartNamesRead. */

  if(!bitIsSet(par.dataFlags, DS_bit_populations)){
    snprintf(message, STR_LEN_0, "Grid input file doesn't have all needed data. Flags = 0x%x.", par.dataFlags);
    PyErr_SetString(PyExc_ValueError, message);

    freeGrid((unsigned int)par.ncell, (unsigned short)par.nSpecies, gp);

    freeMolData(par.nSpecies, md);
    freeImgInfo(par.nImages, img);
    freeConfigInfo(&par);
    freeArrayOfStrings(collPartNamesRead, numCollPartRead);
return NULL;
  }

  parseInput_new(&par, &img, TRUE); /* Sets par.numDensities */

  /* Allocate moldata array.
  */
  if(par.nSpecies>0){
    mallocAndSetDefaultMolData(par.nSpecies, &md);
  } /* otherwise leave it at NULL - we will not be using it. */

  checkUserDensWeights(&par); /* In collparts.c. Needs par.numDensities. */

  /* Calculate some necessary quantities.
  */
  if(par.dust != NULL)
    readDustFile(par.dust, &lamtab, &kaptab, &nEntries); /* This mallocs lamtab, kaptab. */

  i = 0;
  if(img[i].doline){
    molInit(&par, md);
    calcGridMolDoppler(&par, md, gp);

    if(par.useAbun)
      calcGridMolDensities(&par, &gp);

    for(gi=0;gi<par.ncell;gi++){
      for(si=0;si<par.nSpecies;si++){
        gp[gi].mol[si].specNumDens = malloc(sizeof(double)*md[si].nlev);
        for(ei=0;ei<md[si].nlev;ei++){
          gp[gi].mol[si].specNumDens[ei] = 0.0;
        }
      }
    }

    calcGridMolSpecNumDens(&par,md,gp);
  }

  /* raytracing
  */
  raytrace(i, &par, gp, md, img, lamtab, kaptab, nEntries);

  /* Fits cube output.
  */
  write4Dfits(i, 0, &par, img);

  /* Free everything.
  */
  if(par.dust != NULL){
    free(kaptab);
    free(lamtab);
  }
  freeGrid((unsigned int)par.ncell, (unsigned short)par.nSpecies, gp);
  freeMolData(par.nSpecies, md);
  freeImgInfo(par.nImages, img);
  freeConfigInfo(&par);

  freeArrayOfStrings(collPartNamesRead, numCollPartRead);

Py_RETURN_NONE;
}

/*....................................................................*/
/*
 * Bind Python function names to our C functions
 */
static PyMethodDef myModule_methods[] = {
  {"run",  py_run, METH_VARARGS},
  {NULL, NULL}
};

/*
 * Python calls this to let us initialize our module
 */
void initcasaray(void)
{
  (void) Py_InitModule("casaray", myModule_methods);
}

