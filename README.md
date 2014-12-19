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

#### Pre-requisites

The LIME code needs three libraries in order to compile: qhull, gsl,
cfitsio and ncurses. It also requires Python, Numpy and CMake, and a C
compiler.

On MacOSX these dependencies can be installed using
[MacPorts](http://www.macports.org). After MacPorts has been
installed, type in a terminal (as root or sudo)

```
$ sudo port install gsl gcc48 cfitsio atlas ncurses qhull py27-numpy
```

#### Compilation

LIME uses CMake as a build system. To compile it, type in a terminal:

```
mkdir build
cd build
cmake ../
make
```

On MacOSX, the default C compiler (Clang) causes LIME to segfault. On
this platform on needs to use the GNU C compiler instead. This can be
done by setting the ```CC`` variable *before* compiling the code:

```
export CC=gcc-mp-4.8
```

Running the code
----------------

The code runs simply by typing:

```
$ ./lime build/doc/example/input.ini
```

The ```input.ini``` file contains the code input parameters. The
physical parameters (density, temperature, velocity, abundances) are
read from a Python file. See the ```model.py``` file for more details.

[![Build Status](https://travis-ci.org/smaret/lime.svg)]
(https://travis-ci.org/smaret/lime)

