LIME (Line Modeling Engine)
===========================

Copyright (C) 2006-2014 Christian Brinch.

Copyright (C) 2015-2017 the LIME development team.

LIME was created by Christian Brinch, but is now maintained by several people. See the LIME repository on [GitHub](https://github.com/lime-rt/lime) for further details.
  
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

A comprehensive user manual is available online at the [ReadTheDocs](http://lime.readthedocs.org/) website.

Installation notes
------------------

The LIME code needs three library packages in order to compile: qhull,
gsl, and cfitsio. If these are not already present on your system, you will need to install them.

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

If one or more of these packages is not present on your system, you can download the sources from the following locations. Make sure to get
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
compiled with the `-fno-strict-aliasing` flag. Note that the example here does not 
require root privileges. Some modifications to the Makefile may be required if another
location is set for the installation.

*Note for qhull2011.1 and later:* This version of qhull will sometimes
not compile unless the `-Wno-sign-conversion` flag is removed from the
qhull Makefile. The naming of the qhull library has changed between
version 2010.1 and 2011.1. Make sure to edit the qhull flag near the
top of the LIME Makefile accordingly. Also, the newest versions of
qhull does not include a configure script.

LIME is automatically compiled at runtime so there is no installation
required. Do not try to `make` or `make install` LIME as this will produce an error.

Configuring LIME
----------------

We added a configure script with LIME version 1.9 to avoid the necessity to set extra environment variables or hack the Makefile etc. in order to deal with different names for cfitsio/qhull headers and libraries on different systems. You should run this script once after you install LIME on your machine, viz:

```
 $ cd <LIME directory>
 $ ./configure
```

Doing this generates a file `Makefile.defs`, without which LIME will not compile.


Running the code
----------------

The path to the `lime` script needs to be in your `PATH`
environment variable. If you are using bash, do:

```
 $ export PATH=/path/to/lime/:${PATH}
```

where `/path/to/lime/` is the directory where the LIME source code
is located. If you are using csh, do:

```
$ setenv PATH /path/to/lime/:${PATH}
```

To find out whether you use bash or csh, do

```
$ echo $SHELL
```

The code runs simply by typing (for example)

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

Note that with version 1.9, a model file can be written in python; however to use this you will need to compile LIME specially. Go to the LIME parent directory and type

```
 $ make pylime
```

This generates an executable named `pylime`. You run this with your python module as (for example)

```
 $ cd example
 $ ../pylime model.py
```

[![Build Status](https://travis-ci.org/lime-rt/lime.svg?branch=master)](https://travis-ci.org/lime-rt/lime) [![Documentation Status](https://readthedocs.org/projects/lime/badge/?version=latest)](https://readthedocs.org/projects/lime/?badge=latest)
