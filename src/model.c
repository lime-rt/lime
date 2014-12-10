/*
 *  model.c
 *  LIME, The versatile 3D line modeling tool
 *
 *  Created by Christian Brinch on 11/05/07.
 *  Copyright 2006-2013, Christian Brinch,
 *  <brinch@nbi.dk>
 *  Niels Bohr institutet
 *  University of Copenhagen
 *	All rights reserved.
 *
 */

#include "Python.h"
#include "lime.h"

PyObject* density_py = NULL;
PyObject* velocity_py = NULL;
PyObject* temperature_py = NULL;
PyObject* doppler_py = NULL;
PyObject* abundance_py = NULL;
PyObject* py_module = NULL;
PyObject* math_module = NULL;

/******************************************************************************/

int
input(char* input_file, inputPars *par, image *img){

  FILE *f;
  char line[MAX_LINE];
  char keyword[MAX_LINE];
  char parameter[MAX_LINE];
  char value[MAX_LINE];
  int line_number = 0;
  errno = 0;

  /* Open the input file or exit if we can't open it. */
  f = fopen (input_file, "r");
  if (!f)
    {
      fprintf (stderr, "astrochem: error: Can't open %s: %s\n", input_file,
               strerror (errno));
      return EXIT_FAILURE;
    }

  /* Loop over the lines, and look for keywords (between brackets) and
     parameters/values (separated by "="). */

  while (fgets (line, MAX_LINE, f) != NULL)
    {
      line_number++;
      if (line[0] == '#')
        continue;               /* Skip comments */
      if (sscanf (line, "[ %512[a-zA-Z] ]", keyword) == 1)
        ;
      else if (sscanf (line, "%s = %s", parameter, value) == 2)
        {
          if (strcmp (keyword, "input") == 0)
            {
              if (strcmp (parameter, "radius") == 0)
                par->radius = strtof(value, NULL )*AU;
              else if (strcmp (parameter, "minScale") == 0)
                par->minScale = strtof(value, NULL)*AU;
              else if (strcmp (parameter, "pIntensity") == 0)
                par->pIntensity = strtol(value, NULL, 10 );
              else if (strcmp (parameter, "sinkPoints") == 0)
                par->sinkPoints = strtol(value, NULL, 10);
              else if (strcmp (parameter, "tcmb") == 0)
                par->tcmb = strtof(value, NULL);
              else if (strcmp (parameter, "dust") == 0)
                strcpy (par->dust, value);
              else if (strcmp (parameter, "moldatfile") == 0)
                strcpy (par->moldatfile[0], value);
              else if (strcmp (parameter, "pregrid") == 0)
                strcpy (par->pregrid, value);
              else if (strcmp (parameter, "restart") == 0)
                strcpy (par->restart, value);
              else if (strcmp (parameter, "lte_only") == 0)
                par->lte_only = strtol(value, NULL, 10);
              else if (strcmp (parameter, "blend") == 0)
                par->blend = strtol(value, NULL, 10);
              else if (strcmp (parameter, "antialias") == 0)
                par->antialias = strtol(value, NULL, 10);
              else if (strcmp (parameter, "polarization") == 0)
                par->polarization = strtol(value, NULL, 10);
              else if (strcmp (parameter, "sampling") == 0)
                par->sampling = strtol(value, NULL, 10);
              else
                {//TODO curses ?
                  fprintf (stderr, "lime: error: incorrect input in %s line %i : unknow parameter %s\n",
                           input_file, line_number, parameter );
                  return EXIT_FAILURE;
                }
            }
          else if (strcmp (keyword, "output") == 0)
            {
              if (strcmp (parameter, "outputfile") == 0)
                strcpy (par->outputfile, value);
              else if (strcmp (parameter, "binoutputfile") == 0)
                strcpy (par->binoutputfile, value);
              else if (strcmp (parameter, "gridfile") == 0)
                strcpy (par->gridfile, value);
              else
                { //TODO Curses
                  fprintf (stderr, "lime: error: incorrect input in %s line %i : unknow parameter %s\n",
                           input_file, line_number, parameter );
                  return EXIT_FAILURE;
                }
            }
          else if (strcmp (keyword, "image") == 0)
            {
              if (strcmp (parameter, "nchan") == 0)
                img[0].nchan = strtol( value , NULL, 10);
              else if (strcmp (parameter, "phi") == 0)
                img[0].phi = strtof( value, NULL );
              else if (strcmp (parameter, "source_vel") == 0)
                img[0].source_vel = strtof( value, NULL );
              else if (strcmp (parameter, "freq") == 0)
                img[0].freq = strtof( value, NULL );
              else if (strcmp (parameter, "bandwidth") == 0)
                img[0].bandwidth = strtof( value, NULL );
              else if (strcmp (parameter, "velres") == 0)
                img[0].velres = strtof( value, NULL );
              else if (strcmp (parameter, "trans") == 0)
                img[0].trans = strtol( value , NULL, 10);
              else if (strcmp (parameter, "pxls") == 0)
                img[0].pxls = strtol( value , NULL, 10);
              else if (strcmp (parameter, "imgres") == 0)
                img[0].imgres = strtof( value, NULL );
              else if (strcmp (parameter, "theta") == 0)
                img[0].theta = strtof( value, NULL );
              else if (strcmp (parameter, "distance") == 0)
                img[0].distance = strtof( value, NULL )*PC;
              else if (strcmp (parameter, "unit") == 0)
                img[0].unit = strtol( value , NULL, 10);
              else if (strcmp (parameter, "filename") == 0)
                strcpy (img[0].filename, value);
              else
                { //TODO Curses
                  fprintf (stderr, "lime: error: incorrect input in %s line %i : unknow parameter %s\n",
                           input_file, line_number, parameter );
                  return EXIT_FAILURE;
                }
            }
          else if (strcmp (keyword, "model") == 0)
            {
              if (strcmp (parameter, "module") == 0)
                strcpy( par->python_module_name, value );
              else if (strcmp (parameter, "path") == 0)
                strcpy( par->python_module_path, value );
              else if (strcmp (parameter, "path") == 0)
                strcpy( par->python_module_path, value );
              else if (strcmp (parameter, "density") == 0)
                strcpy( par->density_func_name, value );
              else if (strcmp (parameter, "velocity") == 0)
                strcpy( par->velocity_func_name, value );
              else if (strcmp (parameter, "temperature") == 0)
                strcpy( par->temperature_func_name, value );
              else if (strcmp (parameter, "doppler") == 0)
                strcpy( par->doppler_func_name, value );
              else if (strcmp (parameter, "abundance") == 0)
                strcpy( par->abundance_func_name, value );
              else
                { //TODO Curses
                  fprintf (stderr, "lime: error: incorrect input in %s line %i : unknow parameter %s\n",
                           input_file, line_number, parameter );
                  return EXIT_FAILURE;
                }
            }
          /* Unknown or unspecified keyword */
          else
            {
              fprintf (stderr, "lime: error: incorrect input in %s line %i : unknow keyword %s\n",
                       input_file, line_number, keyword);
              return EXIT_FAILURE;
            }
        }
      /* Error while reading a parameter/values */
      else
        {
          fprintf (stderr, "lime: error: incorrect input in %s line %i: Error while reading param = value \n",
                   input_file, line_number);
          return EXIT_FAILURE;
        }
    }
  /* Close the file */
  fclose (f);
  return EXIT_SUCCESS;
}


void
density(double x, double y, double z, double *density){
 if( density_py == NULL )
    {
      fprintf (stderr, "lime: error: density if not defined in model\n");
    }
  else
    {
      PyObject* pArgs = PyTuple_New(3);
      PyObject* pValue;
      pValue = PyFloat_FromDouble( x );
      if( !pValue )
        {
          Py_DECREF(pArgs);
          fprintf(stderr, "Cannot convert density argument\n");
          return;
        }
      PyTuple_SetItem(pArgs, 0, pValue );

      pValue = PyFloat_FromDouble( y );
      if( !pValue )
        {
          Py_DECREF(pArgs);
          fprintf(stderr, "Cannot convert density argument\n");
          return;
        }
      PyTuple_SetItem(pArgs, 1, pValue );

      pValue = PyFloat_FromDouble( z );
      if( !pValue )
        {
          Py_DECREF(pArgs);
          fprintf(stderr, "Cannot convert density argument\n");
          return;
        }
      PyTuple_SetItem(pArgs, 2, pValue );

      pValue = PyObject_CallObject( density_py, pArgs);
      Py_DECREF(pArgs);
      if (pValue != NULL)
        {
          density[0] = PyFloat_AsDouble(pValue);
          Py_DECREF(pValue);
        }
      else
        {
          PyErr_Print();
          fprintf(stderr, "Cannot compute density\n");
        }
    }
}

/******************************************************************************/

void
temperature(double x, double y, double z, double *temperature){
  int i,x0=0;
  double r;
  /*
   * Array containing temperatures as a function of radial
   * distance from origin (this is an example of a tabulated model)
   */
  double temp[2][10] = {
      {2.0e13, 5.0e13, 8.0e13, 1.1e14, 1.4e14, 1.7e14, 2.0e14, 2.3e14, 2.6e14, 2.9e14},
      {44.777, 31.037, 25.718, 22.642, 20.560, 19.023, 17.826, 16.857, 16.050, 15.364}
  };
  /*
   * Calculate coordinate distance from origin
   */
  r=sqrt(x*x+y*y+z*z);
  /*
   * Linear interpolation in temperature input
   */
  if(r > temp[0][0] && r<temp[0][9]){
    for(i=0;i<9;i++){
      if(r>temp[0][i] && r<temp[0][i+1]) x0=i;
    }
  }
  if(r<temp[0][0])
    temperature[0]=temp[1][0];
  else if (r>temp[0][9])
    temperature[0]=temp[1][9];
  else
    temperature[0]=temp[1][x0]+(r-temp[0][x0])*(temp[1][x0+1]-temp[1][x0])/(temp[0][x0+1]-temp[0][x0]);
}

/******************************************************************************/

void
abundance(double x, double y, double z, double *abundance){
  /*
   * Here we use a constant abundance. Could be a
   * function of (x,y,z).
   */
  abundance[0] = 1.e-9;
}

/******************************************************************************/

void
doppler(double x, double y, double z, double *doppler){
  /*
   * 200 m/s as the doppler b-parameter. This
   * can be a function of (x,y,z) as well.
   * Note that *doppler is a pointer, not an array.
   * Remember the * in front of doppler.
   */
  *doppler = 200.;
}

/******************************************************************************/

void
velocity(double x, double y, double z, double *vel){
  /*
   * Variables for spherical coordinates
   */
  double R, phi,r,theta;
  /*
   * Transform Cartesian coordinates into spherical coordinates
   */
  R=sqrt(x*x+y*y+z*z);
  theta=atan2(sqrt(x*x+y*y),z);
  phi=atan2(y,x);
  /*
   * Free-fall velocity in the radial direction onto a central
   * mass of 1.0 solar mass
   */
  r=-sqrt(2*6.67e-11*1.989e30/R);
  /*
   * Vector transformation back into Cartesian basis
   */
  vel[0]=r*sin(theta)*cos(phi);
  vel[1]=r*sin(theta)*sin(phi);
  vel[2]=r*cos(theta);
}

/******************************************************************************/


int python_call_initialize( const inputPars *par )
{
  Py_Initialize();
  PyObject* pName = PyString_FromString( par->python_module_name );
  PyObject* pNameMath = PyString_FromString( "math" );

  PyObject* sysPath = PySys_GetObject((char*)"path");
  PyObject* modulePath = PyString_FromString( par->python_module_path );
  PyList_Append(sysPath, modulePath );
  Py_DECREF( modulePath );

  py_module = PyImport_Import(pName);
  math_module = PyImport_Import(pNameMath);
  Py_DECREF(pName);
  Py_DECREF(pNameMath);

  if (py_module == NULL )
    {
      PyErr_Print();
      fprintf(stderr, "Failed to import \"%s\"\n", par->python_module_name );
      return EXIT_FAILURE;
    }
  if (math_module == NULL )
    {
      PyErr_Print();
      fprintf(stderr, "Failed to import math" );
      return EXIT_FAILURE;
    }

  density_py = PyObject_GetAttrString(py_module, par->density_func_name);
  if( density_py == NULL || !PyCallable_Check( density_py ) )
    {
      if (PyErr_Occurred())
        {
          PyErr_Print();
        }
      fprintf(stderr, "Cannot find density function \"%s\"\n", par->density_func_name );
      Py_XDECREF(density_py);
      Py_DECREF(py_module);
      Py_DECREF(math_module);
      return EXIT_FAILURE;
    }

  velocity_py = PyObject_GetAttrString(py_module, par->velocity_func_name);
  if( velocity_py == NULL || !PyCallable_Check( velocity_py ) )
    {
      if (PyErr_Occurred())
        {
          PyErr_Print();
        }
      fprintf(stderr, "Cannot find velocity function \"%s\"\n", par->velocity_func_name );
      Py_XDECREF(velocity_py);
      Py_DECREF( density_py );
      Py_DECREF(py_module);
      Py_DECREF(math_module);
      return EXIT_FAILURE;
    }

  temperature_py = PyObject_GetAttrString(py_module, par->temperature_func_name);
  if( temperature_py == NULL || !PyCallable_Check( temperature_py ) )
    {
      if (PyErr_Occurred())
        {
          PyErr_Print();
        }
      fprintf(stderr, "Cannot find temperature function \"%s\"\n", par->temperature_func_name );
      Py_XDECREF(temperature_py);
      Py_DECREF( velocity_py );
      Py_DECREF( density_py );
      Py_DECREF(py_module);
      Py_DECREF(math_module);
      return EXIT_FAILURE;
    }

  doppler_py = PyObject_GetAttrString(py_module, par->doppler_func_name);
  if( doppler_py == NULL || !PyCallable_Check( doppler_py ) )
    {
      if (PyErr_Occurred())
        {
          PyErr_Print();
        }
      fprintf(stderr, "Cannot find doppler function \"%s\"\n", par->doppler_func_name );
      Py_XDECREF(doppler_py);
      Py_DECREF( velocity_py );
      Py_DECREF( temperature_py );
      Py_DECREF( density_py );
      Py_DECREF(py_module);
      Py_DECREF(math_module);
      return EXIT_FAILURE;
    }

  abundance_py = PyObject_GetAttrString(py_module, par->abundance_func_name);
  if( abundance_py == NULL || !PyCallable_Check( abundance_py ) )
    {
      if (PyErr_Occurred())
        {
          PyErr_Print();
        }
      fprintf(stderr, "Cannot find abundance function \"%s\"\n", par->abundance_func_name );
      Py_XDECREF(abundance_py);
      Py_DECREF( doppler_py );
      Py_DECREF( velocity_py );
      Py_DECREF( temperature_py );
      Py_DECREF( density_py );
      Py_DECREF( py_module );
      Py_DECREF(math_module);
      return EXIT_FAILURE;
    }
  Py_DECREF( sysPath );
  return EXIT_SUCCESS;
}

void python_call_finalize()
{
  Py_CLEAR( abundance_py );
  Py_CLEAR( velocity_py );
  Py_CLEAR( temperature_py );
  Py_CLEAR( doppler_py );
  Py_CLEAR( density_py );
  Py_CLEAR( py_module );
  Py_CLEAR( math_module );
  Py_Finalize();
  abundance_py = NULL;
  velocity_py = NULL;
  temperature_py = NULL;
  doppler_py = NULL;
  density_py = NULL;
  py_module = NULL;
  math_module = NULL;
}
