LIME (Line Modeling Engine)
===========================

Created by [Christian Brinch](mailto:brinch@nbi.dk), Copyright 2006-2013

Niels Bohr institutet, University of Copenhagen
  
About LIME
----------

LIME is a 3D molecular excitation and radiation transfer code for
far-infrared and (sub-)millimeter wavelength. LIME will calculate
spectra of rotational transitions of atoms and molecules, given a
user-supplied physical model.  Details on the method can be found in
[C. Brinch and M. R. Hogerheijde, A&A 553, A25
(2010)](http://adsabs.harvard.edu/abs/2010A%26A...523A..25B)

Any scientific publication making use of the LIME code should also
reference this publication.

A comprehensive user manual is available online at the [LIME
website](http://www.nbi.dk/~brinch/lime.html).

Installation notes
------------------

The LIME code needs three library packages in order to compile: qhull,
gsl, and cfitsio.

### Mac OS X

The easiest way is to install
[MacPorts](http://www.macports.org). After MacPorts has been
installed, type in a terminal (as root or sudo)

```
$ port install qhull
$ port install cfitsio
$ port install gsl
```

Alternatively, the library packages may be installed with
[Fink](http://www.finkproject.org). After Fink has been installed,
type:

```
$ fink install qhull6.3.1-dev cfitsio gsl
```

### Linux / Unix / Mac OS X (alternative)

Download the sources from the following locations. Make sure to get
the latest versions.

- [qhull](http://www.qhull.org/download/)
- [cfitsio](http://heasarc.gsfc.nasa.gov/fitsio/)
- [gsl](http://www.gnu.org/software/gsl/)

The GNU scientific library (GSL), is present on most modern Unix and
Linux systems. Check for availability with `gsl-config --libs`. If
this command returns a library path, there is no need to download and
install it.

All three library packages are installed using the following

```
 $ configure --prefix=/path/to/LimePackage
 $ make
 $ make install
 ```

In some cases the qhull library will produce a segmentation fault unless it is 
compiled with the `-fno-strict-aliasing` flag. The example here does not 
require root. Some modifications to the Makefile may be required if another
location is set for the installation.

*Note for qhull2011.1 and later:* This version of qhull will sometimes
not compile unless the `-Wno-sign-conversion` flag is removed from the
qhull Makefile. The naming of the qhull library has changed between
version 2010.1 and 2011.1. Make sure to edit the qhull flag near the
top of the LIME Makefile accordingly. Also, the newest versions of
qhull does not include a configure script.

LIME is automatically compiled at runtime so there is no installation
required and also no configure script. Do not try to `make` or `make
install` LIME as this will produce an error.

Running the code
----------------

Source the file called source.me using the following commands.

For bash:

```
 $ . sourceme.bash
```

For csh:

```
$ source sourceme.csh
```

To find out whether you use bash or csh, do

```
$ echo $SHELL
```

The source command needs to be executed for each new session. The content of 
the file source.me can be placed in .bashrc or .cshrc for convenience in which 
case there is no need to source the source.me file before each session.

The code runs simply by typing

```
$ cd example/
$ lime model.c
```

The model file can have any name. The model file need not be written
in C.  Subroutines can also be written in Fortran (or other languages
that can be linked with C, e.g., python) as long as the names comply
with the standards of linking C and Fortran. See [this
page](http://tinyurl.com/y6sddr) for information on how to link C and
Fortran. If Fortran subroutines are used, the linking of LIME needs to
be done with the Fortran compiler. Modify the Makefile accordingly.


Changelog
---------

- 1.31: Bug was found and fixed in `getJbar()`

- 1.3: LIME can now run multiple instances in the same directory. Grid
  data structure has changed. A memory leak in qhull call has been
  closed.  model.c no longer has to contain functions that are not
  needed.  LIME can be restarted from previously calculated
  populations. gas- to-dust ratio has been implemented as a function
  rather than a constant.
      
- 1.23: Includes a tcp client so that LIME can automatically download
  LAMDA files. Makefile and lime script improved. Semantic updates.

- 1.22: Fixed a bug in the vtk output. Fixed a bug in the
  randomization of photon directions. vtk file now contains the
  normalized velocities.

- 1.21: A few memory leaks fixed. Output screen updated with more
  information.  LIME now always calculates molecular density with
  respect to total H2 A catch for negative optically depth has been
  implemented.
      
- 1.2: Minor bug fix in ray-tracer (Thanks to Ruud Visser). FITS
  header tweaked for CASA compatability. Photon transport slightly
  optimized.  The photons nitial directions are now
  randomized. Multiple lines no longer get mirrored. LIME should run
  faster now.
      
- 1.1: Raytracer has been rewritten and optimized and now includes an
  anti- aliasing filter. The FITS header has been corrected. LIME
  images can now be read directly by CASA. Continuum polarization
  implemented.

- 1.03: Bug fixes

- 1.02: Minor bug fixes and some clean-up of the code. The VTK file
  now contains density, temperature, molecular density and the
  velocity field. Faster evaluation of line-of-sight velocities

- 1.01: First full public release

[![Build Status](https://travis-ci.org/smaret/lime.svg?branch=master)](https://travis-ci.org/smaret/lime) [![Documentation Status](https://readthedocs.org/projects/lime/badge/?version=latest)](https://readthedocs.org/projects/lime/?badge=latest)
