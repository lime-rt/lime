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
PyObject* numpy_module = NULL;

/******************************************************************************/
void    python_call( PyObject* func_py, double x, double y, double z, double *output);

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
  python_call( density_py, x, y, z, density );
}

/******************************************************************************/

void
temperature(double x, double y, double z, double *temperature){
  python_call( temperature_py, x, y, z, temperature );
}

/******************************************************************************/

void
abundance(double x, double y, double z, double *abundance){
  python_call( abundance_py, x, y, z, abundance );
}

/******************************************************************************/

void
doppler(double x, double y, double z, double *doppler){
  python_call( doppler_py, x, y, z, doppler );
}

/******************************************************************************/

void
velocity(double x, double y, double z, double *vel){
  python_call( velocity_py, x, y, z, vel );
}

/******************************************************************************/


int python_call_initialize( const inputPars *par )
{
  Py_Initialize();
  PyObject* pName = PyString_FromString( par->python_module_name );
  PyObject* pNameMath = PyString_FromString( "math" );
  PyObject* pNameNumpy = PyString_FromString( "numpy" );

  PyObject* sysPath = PySys_GetObject((char*)"path");
  PyObject* modulePath = PyString_FromString( par->python_module_path );
  PyList_Append(sysPath, modulePath );
  Py_DECREF( modulePath );

  py_module = PyImport_Import(pName);
  math_module = PyImport_Import(pNameMath);
  numpy_module = PyImport_Import(pNameNumpy);
  Py_DECREF(pName);
  Py_DECREF(pNameMath);
  Py_DECREF(pNameNumpy);

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
  if (numpy_module == NULL )
    {
      PyErr_Print();
      fprintf(stderr, "Failed to import numpy, cannot use numpy in model" );
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
      Py_DECREF(numpy_module);
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
      Py_DECREF(numpy_module);
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
      Py_DECREF(numpy_module);
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
      Py_DECREF(numpy_module);
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
      Py_DECREF(numpy_module);
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
  Py_CLEAR( numpy_module );
  Py_Finalize();
  abundance_py = NULL;
  velocity_py = NULL;
  temperature_py = NULL;
  doppler_py = NULL;
  density_py = NULL;
  py_module = NULL;
  math_module = NULL;
  numpy_module = NULL;
}

void
python_call( PyObject* func_py, double x, double y, double z, double *output){
  if( func_py == NULL )
    {
      fprintf (stderr, "lime: error: func if not defined in model\n");
    }
  else
    {
      PyObject* pArgs = PyTuple_New(3);
      PyObject* pValue;
      pValue = PyFloat_FromDouble( x );
      if( !pValue )
        {
          Py_DECREF(pArgs);
          fprintf(stderr, "Cannot convert python func argument\n");
          return;
        }
      PyTuple_SetItem(pArgs, 0, pValue );

      pValue = PyFloat_FromDouble( y );
      if( !pValue )
        {
          Py_DECREF(pArgs);
          fprintf(stderr, "Cannot convert python func argument\n");
          return;
        }
      PyTuple_SetItem(pArgs, 1, pValue );

      pValue = PyFloat_FromDouble( z );
      if( !pValue )
        {
          Py_DECREF(pArgs);
          fprintf(stderr, "Cannot convert python func argument\n");
          return;
        }
      PyTuple_SetItem(pArgs, 2, pValue );

      pValue = PyObject_CallObject( func_py, pArgs);
      Py_DECREF(pArgs);
      if (pValue != NULL)
        {
          if( PyTuple_Check( pValue ) )
            {
              Py_ssize_t tuple_size = PyTuple_GET_SIZE( pValue );
              Py_ssize_t i;
              PyObject* tupleValue;
              for( i = 0; i<tuple_size; i++ )
                {
                  tupleValue = PyTuple_GET_ITEM( pValue, i );
                  if( PyFloat_Check( tupleValue ) )
                    {
                      output[i] = PyFloat_AsDouble(tupleValue);
                    }
                  else if( PyInt_Check( tupleValue ) )
                    {
                      output[i] = PyInt_AsLong( tupleValue );
                    }
                  else
                    {
                      fprintf(stderr, "Cannot convert returned tuple value\n");
                    }
                  Py_DECREF(tupleValue);
                }
            }
          else
            {
              if( PyFloat_Check( pValue ) )
                {
                  output[0] = PyFloat_AsDouble(pValue);
                }
              else if( PyInt_Check( pValue ) )
                {
                  output[0] = PyInt_AsLong( pValue );
                }
              else
                {
                  fprintf(stderr, "Cannot convert returned value\n");
                }
            }
          Py_DECREF(pValue);
        }
      else
        {
          PyErr_Print();
          fprintf(stderr, "Cannot compute python func\n");
        }
    }
}
