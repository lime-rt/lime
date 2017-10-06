The CASA interface to LIME.
+++++++++++++++++++++++++++

Note that wherever I have written below something having the form <some text between angle brackets>, the whole should be understood to be a placeholder for actual text which is not known a priori.


Purpose:
========

The purpose of this interface is to allow the user to run LIME from the CASA command line.

In order to keep the number of parameters under control, the functionality of LIME has been apportioned between two separate CASA tasks, limesolver and raytrace, described separately below.

limesolver:
-----------
This interface allows the user to choose a model from a library (currently containing 9 choices), to adjust its parameters, and to compute the resulting spatial distribution of excited states by means of LIME. The output from this task is a FITS or HDF5 file (currently the software is configured for HDF5) which contains all the necessary information to construct an image cube from the solution. Note that limesolver itself does not construct any images. Most of the (non-image) LIME parameters may also be set via this interface.

Note that there is currently no support to allow the user to supply their own bespoke model, as in 'traditional' LIME.

raytrace:
---------
This task reads the output file from limesolver and creates a single image via the standard set of LIME image parameters. Creation of multiple images in a single run is not supported.


Code relationships
==================

The functionality is currently distributed between 2 packages, LIME and ARTIST, which are under separate systems of development and version control. The relation between the packages may be envisaged in 3 tiers:

  CASA interface files (maintained as part of LIME)
             |
  Python modules to access respectively the LIME engine and the library of models (modellib) (part of ARTIST)
        |                             |
  The LIME engine (LIME)  <-  The modellib (ARTIST)


How to install, configure and run the code.
===========================================

The CASA interface to LIME is contained in the Artist package, which should be available as a tarball. In the description below I'll assume the user has unpacked this tarball into a directory ARTISTROOT/.

ARTIST:
-------
Instructions for compiling and installing ARTIST (with LIME included) are to be found in the file ARTISTROOT/README. Once these instructions have been followed successfully, the user ought to be able to import both 'lime' and 'modellib' in an interactive python session. Typing 'help(lime)' or 'help(modellib)' on that python command line should also confirm that these modules are correctly sourced from ARTISTROOT/lib/python<version>/site-packages/.

The CASA interface:
-------------------
There are two CASA interface files for each CASA task: a python file named 'task_<taskname>.py' and an XML file '<taskname>.xml'. These are contained in the directory ARTISTROOT/packages/lime/lime/casa/. Note that there may be other files there too.

To compile the CASA interfaces, the CASA script 'buildmytasks' should be run in the directory ARTISTROOT/packages/lime/lime/casa/. In the 5.1 release of CASA, this script is contained in the directory <casa root>/bin/. The 'buildmytasks' script generates 2 more files for each CASA task in the PWD (namely, <taskname>.py and <taskname>_cli.py) plus a file 'mytasks.py' (yes that actually is its name).

The user may find it more convenient to copy the 'task_<taskname>.py' and '<taskname>.xml' files from ARTISTROOT/packages/lime/lime/casa/ to some other directory before running buildmytasks on them.

The final thing you have to do to allow CASA to find the LIME tasks is to type, at the CASA prompt (i.e., after you have started CASA):

  execfile('<directory where it is>/mytasks.py')

If you then type

  tasklist

at the prompt, you should see limesolver and raytrace at the bottom of the resulting list.

Note that you can put the 'execfile' instruction in a script 'init.py' in your ~/.casa/ directory (just create one if you don't already have one) to save having to type it each time.


Running the tasks:
------------------
The tasks limesolver and raytrace are standard CASA tasks and are thus run in the same way as other CASA tasks (see e.g. https://casa.nrao.edu/docs/cookbook/casa_cookbook002.html#sec32).

It is of course possible to import the modules 'lime' and 'modellib' directly at the CASA prompt. These can then be run if desired without the interface wrappers, which allows the user more freedom, with the usual expense of some added complication. There is comprehensive information on these modules in the document ARTISTROOT/doc/PythonInterface.pdf.

Documentation describing the operations of LIME can be found in a variety of formats on either

  https://github.com/lime-rt/lime/blob/master/doc/usermanual.rst

or

  https://readthedocs.org/projects/lime/

but it is also expected to be available in HTML form under ARTISTROOT/packages/lime/lime/doc/_html/. The advantage of the latter is that it refers to the version of LIME which comes bundled with ARTIST, whereas the former describes the latest version of LIME, which may be in advance of the ARTIST one.

