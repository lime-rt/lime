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

#include "lime.h"
#include "matheval.h"

/******************************************************************************/

void* density_evaluator = NULL;
void* velocity_evaluator = NULL;
void* doppler_evaluator = NULL;
void* abundance_evaluator = NULL;

int
input(char* input_file, inputPars *par, image *img){

  FILE *f;
  char line[MAX_LINE];
  char keyword[MAX_LINE];
  char parameter[MAX_LINE];
  char value[MAX_LINE];
  int line_number = 0;
  int i = 0, j = 0;
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
              if (strcmp (parameter, "density") == 0)
                {
                  if( density_evaluator != NULL )
                    {
                        evaluator_destroy ( density_evaluator );
                    }
                  density_evaluator = evaluator_create (value);
                }
              else if (strcmp (parameter, "doppler") == 0)
                {
                  if( doppler_evaluator != NULL )
                    {
                        evaluator_destroy ( doppler_evaluator );
                    }
                  doppler_evaluator = evaluator_create (value);
                }
              else if (strcmp (parameter, "abundance") == 0)
                {
                  if( abundance_evaluator != NULL )
                    {
                        evaluator_destroy ( abundance_evaluator );
                    }
                  abundance_evaluator = evaluator_create (value);
                }
              else if (strcmp (parameter, "velocity") == 0)
                {
                  if( velocity_evaluator != NULL )
                    {
                        evaluator_destroy ( velocity_evaluator );
                    }
                  velocity_evaluator = evaluator_create (value);
                }


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
}

/******************************************************************************/

void destroy_model_evaluator()
{
  if( density_evaluator != NULL )
    {
      evaluator_destroy ( density_evaluator );
    }
  if( abundance_evaluator != NULL )
    {
      evaluator_destroy ( abundance_evaluator );
    }
  if( doppler_evaluator != NULL )
    {
      evaluator_destroy ( doppler_evaluator );
    }
  if( velocity_evaluator != NULL )
    {
      evaluator_destroy ( velocity_evaluator );
    }
}

void
density(double x, double y, double z, double *density){
  if( density_evaluator == NULL )
    {
      if(!silent) bail_out("Density is not defined in model.c but is needed by LIME!");
    }
  else
    {
      density[0]= evaluator_evaluate_x_y_z( density_evaluator, x, y, z);
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
  if( abundance_evaluator == NULL )
    {
      if(!silent) bail_out("abundance is not defined in model.c but is needed by LIME!");
    }
  else
    {
      abundance[0]= evaluator_evaluate_x_y_z( abundance_evaluator, x, y, z);
    }
}

/******************************************************************************/

void
doppler(double x, double y, double z, double *doppler){
  if( doppler_evaluator == NULL )
    {
      if(!silent) bail_out("doppler is not defined in model.c but is needed by LIME!");
    }
  else
    {
      *doppler= evaluator_evaluate_x_y_z( doppler_evaluator, x, y, z);
    }
}

/******************************************************************************/

void
velocity(double x, double y, double z, double *vel){
  if( velocity_evaluator == NULL )
    {
      if(!silent) bail_out("velocity is not defined in model.c but is needed by LIME!");
    }

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
   r= evaluator_evaluate_x( velocity_evaluator, R);
  /*
   * Vector transformation back into Cartesian basis
   */
  vel[0]=r*sin(theta)*cos(phi);
  vel[1]=r*sin(theta)*sin(phi);
  vel[2]=r*cos(theta);
}

/******************************************************************************/


