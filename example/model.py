#!/usr/bin/env python
# Notes:
#	- The above pragma line is only required if you plan to run this module as a stand-alone script to run any test harnesses which may occur after "if __name__ == '__main__'".
#
#	- You should run lime with this script in the form
#		pylime model.py
#	  You will need the location of pylime in your PATH environment variable; also you should have the location of the par_classes.py module in your PYTHONPATH environment variable.

import math

# For definitions of the classes ModelParameters and ImageParameters:
from par_classes import *

# Note that the useful macros defined in lime.h are also provided here in the dictionary 'macros' provided as an argument to each funciton below. See the example code at the end for the full list of macro values provided.

#.......................................................................
def input(macros):
  par = ModelParameters()

  # Basic parameters. See cheat sheet for details.
  #
  par.radius            = 2000*macros["AU"]
  par.minScale          = 0.5*macros["AU"]
  par.pIntensity        = 4000
  par.sinkPoints        = 3000
  par.dust              = "jena_thin_e6.tab"
  par.moldatfile        = ["hco+@xpol.dat"] # must be a list, even when there is only 1 item.
  par.antialias         = 4
  par.sampling          = 2

  par.outputfile        = "populations.pop"
  par.binoutputfile     = "restart.pop"
  par.gridfile          = "grid.vtk"

  #    Setting elements of the following three arrays is optional. NOTE
  #    that, if you do set any of their values, you should set as many as
  #    the number of elements returned by your function density(). The
  #    ith element of the array in question will then be assumed to refer
  #    to the ith element in the density function return. The current
  #    maximum number of elements allowed is 7, which is the number of
  #    types of collision partner recognized in the LAMBDA database.
  #
  #    Note that there is no (longer) a hard connection between the
  #    number of density elements and the number of collision-partner
  #    species named in the moldata files. This means in practice that,
  #    if you set the values for par->collPartIds, you can, if you like,
  #    set some for which there are no transition rates supplied in the
  #    moldatfiles. This might happen for example if there is a molecule
  #    which contributes significantly to the total molecular density but
  #    for which there are no measured collision rates for the radiating
  #    species you are interested in.
  #
  #    You may also omit to mention in par->collPartIds a collision
  #    partner which is specified in the moldatfiles. In this case LIME
  #    will assume the density of the respective molecules is zero.
  #
  #    If you don't set any values for any or all of these arrays,
  #    i.e. if you omit any mention of them here (we preserve this
  #    possibility for purposes of backward compatibility), LIME will
  #    attempt to replicate the algorithms employed in version 1.5, which
  #    involve guessing which collision partner corresponds to which
  #    density element. Since this was not exactly a rigorous procedure,
  #    we recommend use of the arrays.
  #
  #    par->nMolWeights: this specifies how you want the number density
  #    of each radiating species to be calculated. At each grid point a
  #    sum (weighted by par->nMolWeights) of the density values is made,
  #    then this is multiplied by the abundance to return the number
  #    density.
  #
  #    par->dustWeights: this is similar, but the weighted sum of
  #    densities now feeds into the calculation of the dust opacity.
  #
  #    Note that there are convenient macros defined in ../src/lime.h for
  #    7 types of collision partner.
  #
  #    Below is an example of how you might use these parameters:

  par.collPartIds       = [macros["CP_H2"]] # must be a list, even when there is only 1 item.
  par.nMolWeights       = [1.0] # must be a list, even when there is only 1 item.
  par.dustWeights       = [1.0] # must be a list, even when there is only 1 item.

  # Definitions for image #0. Add further similar blocks for additional images.
  #
  par.img.append(ImageParameters()) # by default this list par.img has 0 entries.
  par.img[-1].nchan                  = 61             # Number of channels
  par.img[-1].velres                 = 500.0          # Channel resolution in m/s
  par.img[-1].trans                  = 3              # zero-indexed J quantum number
  par.img[-1].pxls                   = 100            # Pixels per dimension
  par.img[-1].imgres                 = 0.1            # Resolution in arc seconds
  par.img[-1].distance               = 140*macros["PC"] # source distance in m
  par.img[-1].source_vel             = 0.0            # source velocity in m/s
  par.img[-1].unit                   = 0              # 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau
  par.img[-1].filename               = "image0.fits"  # Output filename
#  par.img[-1].azimuth                = 0.0
#  par.img[-1].incl                   = 0.0
#  par.img[-1].posang                 = 0.0

  return par

#.......................................................................
#.......................................................................
# User-defined functions:

#.......................................................................
def density(macros, x, y, z):
  """
The value returned should be a list, each element of which is a density (in molecules per cubic metre) of a molecular species (or electrons). The molecule should be one of the 7 types currently recognized in the LAMDA database - see

	http://home.strw.leidenuniv.nl/~moldata/

Note that these species are expected to be the bulk constituent(s) of the physical system of interest rather than species which contribute significantly to spectral-line radiation. In LIME such species are often called 'collision partners'.

The identity of each collision partner is provided via the list parameter par.collPartIds. If you do provide this, obviously it must have the same number and ordering of elements as the density list you provide here; if you don't include it, LIME will try to guess the identities of the species you provide density values for.
  """

  rMin = 0.1*macros["AU"] # greater than zero to avoid a singularity at the origin.

  # Calculate radial distance from origin
  #
  r = math.sqrt(x*x+y*y+z*z)

  # Calculate a spherical power-law density profile
  # (Multiply with 1e6 to go to SI-units)
  #
  if r>rMin:
    rToUse = r
  else:
    rToUse = rMin # Just to prevent overflows at r==0!

  listOfDensities = [1.5e6*((rToUse/(300*macros["AU"]))**(-1.5))*1e6] # must be a list, even when there is only 1 item.

  return listOfDensities

#.......................................................................
def temperature(macros, x, y, z):
  """
This function should return a tuple of 2 temperatures (in kelvin). The 2nd is optional, i.e. you can return None for it, and LIME will do the rest.
  """

  # Array containing temperatures as a function of radial
  # distance from origin (this is an example of a tabulated model)
  #
  rToTemp = [
    [2.0e13, 5.0e13, 8.0e13, 1.1e14, 1.4e14, 1.7e14, 2.0e14, 2.3e14, 2.6e14, 2.9e14],
    [44.777, 31.037, 25.718, 22.642, 20.560, 19.023, 17.826, 16.857, 16.050, 15.364]
  ]

  # Calculate radial distance from origin
  #
  r = math.sqrt(x*x+y*y+z*z)

  # Linear interpolation in temperature input
  #
  xi = 0
  if r>rToTemp[0][0] and r<rToTemp[0][9]:
    for i in range(9):
      if r>rToTemp[0][i] and r<rToTemp[0][i+1]: xi=i

  if r<rToTemp[0][0]:
    temp0 = rToTemp[1][0]
  elif r>rToTemp[0][9]:
    temp0 = rToTemp[1][9]
  else:
    temp0 = rToTemp[1][xi]+(r-rToTemp[0][xi])*(rToTemp[1][xi+1]-rToTemp[1][xi])\
          / (rToTemp[0][xi+1]-rToTemp[0][xi])

  return (temp0, None)

#.......................................................................
def abundance(macros, x, y, z):
  """
This function should return a list of abundances (as fractions of the effective bulk density), 1 for each of the radiating species. Note that the number and identity of these species is set via the list of file names you provide in the par.moldatfile parameter, so make sure at least that the number of elements returned by abundance() is the same as the number in par.moldatfile!

Note that the 'effective bulk density' mentioned just above is calculated as a weighted sum of the values returned by the density() function, the weights being provided in the par.nMolWeights parameter.
  """

  # Here we use a constant abundance. Could be a
  # function of (x,y,z).
  #
  listOfAbundances = [1.0e-9] # must be a list, even when there is only 1 item.
  return listOfAbundances

#.......................................................................
def doppler(macros, x, y, z):
  """
This function returns the Doppler B parameter, defined in terms of a Doppler-broadened Gaussian linewidth as follows:

	             ( -[v-v0]^2 )
	flux(v) = exp(-----------).
	             (    B^2    )

Note that the present value refers only to the Doppler broadening due to bulk turbulence; LIME later adds in the temperature-dependent part (which also depends on molecular mass).
  """

  # 200 m/s as the doppler b-parameter. This
  # can be a function of (x,y,z) as well.
  # Note that *doppler is a pointer, not an array.
  # Remember the * in front of doppler.
  #
  dopplerBValue = 200.0
  return dopplerBValue

#.......................................................................
def velocity(macros, x, y, z):
  """
Gives the bulk gas velocity vector in m/s.
  """

  rMin = 0.1*macros["AU"] # greater than zero to avoid a singularity at the origin.

  # Calculate radial distance from origin
  #
  r = math.sqrt(x*x+y*y+z*z)

  if r>rMin:
    rToUse = r
  else:
    rToUse = rMin # Just to prevent overflows at r==0!

  # Free-fall velocity in the radial direction onto a central 
  # mass of 1.0 solar mass
  #
  ffSpeed = math.sqrt(2.0*macros["GRAV"]*1.989e30/rToUse)

  vel = [0,0,0] # just to initialize its size.
  vel[0] = -x*ffSpeed/rToUse
  vel[1] = -y*ffSpeed/rToUse
  vel[2] = -z*ffSpeed/rToUse

  return vel

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if __name__ == '__main__':
  # Put any private debugging tests here, which you can then run by calling the module directly from the unix command line.

  macros = {\
    "AMU"     :1.66053904e-27,\
    "CLIGHT"  :2.99792458e8,\
    "HPLANCK" :6.626070040e-34,\
    "KBOLTZ"  :1.38064852e-23,\
    "GRAV"    :6.67428e-11,\
    "AU"      :1.495978707e11,\
    "PC"      :3.08567758e16,\
    "PI"      :3.14159265358979323846,\
    "SPI"     :1.77245385091,\
    "CP_H2"   :1,\
    "CP_p_H2" :2,\
    "CP_o_H2" :3,\
    "CP_e"    :4,\
    "CP_H"    :5,\
    "CP_He"   :6,\
    "CP_Hplus":7\
  }

  par = input(macros)

  x = par.radius*0.1
  y = par.radius*0.07
  z = par.radius*0.12

  print density(    macros, x, y, z)[0]
  print temperature(macros, x, y, z)[0]
  print doppler(    macros, x, y, z)
  print velocity(   macros, x, y, z)

