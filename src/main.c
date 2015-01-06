/*
 *  main.c
 *  LIME, The versatile 3D line modeling environment
 *
 *  Created by Christian Brinch on 16/11/06.
 *  Copyright 2006-2014, Christian Brinch,
 *  <brinch@nbi.dk>
 *  Niels Bohr institutet
 *  University of Copenhagen
 *	All rights reserved.
 *
 *  LIME is derived from RATRAN by Michiel Hogerheijde and Floris van der Tak,
 *  Copyright 2000, Hogerheijde and van der Tak, A&A, 362, 697, 2000.
 *
 *  DISCLAIMER: LIME is provided as is and the author does not resume any
 *  responsibility for erroneous results due to bugs in the code.
 *
 *  Any publication with results obtain using LIME should refer to
 *  Brinch & Hogerheijde, A&A, 523, A25, 2010
 *
 */

#include <getopt.h>
#include "lime.h"

void usage (void);
void version (void);

int main ( int argc, char** argv) {
  int i;
  int initime=time(0);
  int popsdone=0;
  molData*     m = NULL;
  inputPars    par;
  struct grid* g = NULL;
  image*       img = NULL;

  char *input_file;

  silent = 0;
  /* Parse options and command line arguments. Diplay help
     message if no (or more than one) argument is given. */

  {
    int opt;

    static struct option longopts[] = {
      {"help", no_argument, NULL, 'h'},
      {"version", no_argument, NULL, 'V'},
      {"quiet", no_argument, NULL, 'q'},
      {0, 0, 0, 0}
    };

    while ((opt = getopt_long (argc, argv, "hVvq", longopts, NULL)) != -1)
      {
        switch (opt)
          {
          case 'h':
            usage ();
            exit (0);
            break;
          case 'V':
            version ();
            exit (0);
            break;
          case 'q':
            silent = 1;
            break;
          default:
            usage ();
            exit (1);
          }
      };
    argc -= optind;
    argv += optind;
    if (argc != 1)
      {
        usage ();
        exit (1);
      }
    input_file = argv[0];
  }

  if(!silent) greetings();
  if(!silent) screenInfo();

  parseInput( input_file ,&par,&img,&m);

  if( strlen(par.pregrid) > 0 )
    {
      gridAlloc(&par,&g);
      predefinedGrid(&par,g);
    }
  else if( strlen(par.restart) > 0 )
    {
      popsin(&par,&g,&m,&popsdone);
    }
  else
    {
      gridAlloc(&par,&g);
      buildGrid(&par,g);
    }

  for(i=0;i<par.nImages;i++){
    if(img[i].doline==1 && popsdone==0) {
      levelPops(m,&par,g,&popsdone);
    }
    if(img[i].doline==0) {
      continuumSetup(i,img,m,&par,g);
    }

    raytrace(i,&par,g,m,img);
    writefits(i,&par,m,img);
  }

  if(!silent) goodnight(initime,img[0].filename);

  freeGrid( &par, m, g);
  freeInput(&par, img, m);
  return 0;
}

/*
   Display help message.
 */

void
usage (void)
{
  fprintf (stdout, "Usage: lime [option...] [file]\n\n");
  fprintf (stdout, "Options:\n");
  fprintf (stdout, "   -h, --help         Display this help\n");
  fprintf (stdout, "   -V, --version      Print program version\n");
  fprintf (stdout, "   -q, --quiet        Suppress all messages\n");
  fprintf (stdout, "\n");
}

/*
   Display version.
 */

void
version (void)
{
  fprintf (stdout, "This is lime, version %s\n", VERSION );
  fprintf (stdout,
           "This is free software. You may redistribute copies of it under the terms\n");
  fprintf (stdout,
           "of the GNU General Public License. There is NO WARRANTY, to the extent\n");
  fprintf (stdout, "permitted by law.\n");
}
