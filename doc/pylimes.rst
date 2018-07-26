.. _pylimes:

Python flavours of LIME
=======================

1. Flavour 'pylime'
-------------------

This flavour was developed to allow the user to write their model file in python. For this, we have to compile and run in separate steps.

Compiling
~~~~~~~~~

Go to the package directory ``<LIME base dir>`` and type

::

    ./configure # if you haven't already done so.
    make pylime

If this completes ok without errors, the next thing to do is to add some directories to your PYTHONPATH environment variable. A script to do this is provided for your convenience. If your shell is bash, do

::

    . ./pylimerc.sh

Alternatively, if you use cshell, do

::

    source ./pylimerc.csh

The compiled executable is called ``pylime`` and will also be found in the ``<LIME base dir>`` directory.

.. note:: Running the ``pylimerc`` script also adds ``<LIME base dir>`` to your PATH environment variable - you don't have to do this as an extra step.

To run pylime, do

::

    pylime [options...] <model file>

Command-line options
~~~~~~~~~~~~~~~~~~~~

.. option:: -s

   Suppresses output messages.

.. option:: -t

   This runs LIME in a test mode, in which it is compiled with the debugging flag set; fixed random seeds are also employed in this mode, so the results of any two runs with the same model should be identical.

.. option:: -p nthreads

   Run in parallel mode with ``nthreads``. The default is a single thread, i.e. serial execution.

.. note::

   The number of threads may also be set with the :ref:`par->nThreads <par-nthreads>` parameter. This will override the value set via the -p option.

.. note::

   Curses-style output is not available for ``pylime``.

The model file
~~~~~~~~~~~~~~

This is written in python, but it follows a very similar format to the C-language file in original use with LIME. A template model file can be found at ``<LIME base dir>/example/model.py``.

Note that it is not necessary to recompile ``pylime`` every time you make a change to the model file.


2. Flavour 'pyshared'
---------------------

The purpose of this flavour is to emulate some features of the package ARTIST. The ``make`` target for this flavour is ``pyshared``. The result is two so-called shared objects, ``<LIME base dir>/liblime.so`` and ``<LIME base dir>/libmodellib.so``, which can be imported into python, either in a python script or in an interactive python session. (Wrapper modules to provide a richer interface are located in ``<LIME base dir>/python/``.) After doing

::

    cd <LIME base dir>
    ./configure # if you haven't already done so.
    make pyshared
    source ./pylimerc.csh # or the bash equivalent

you should be able to cd to anywhere else, run python, and do

::

    >>> import lime
    >>> import modellib

The available functions can be examined via python's help() function. A test script is available at ``<LIME base dir>/tests/pyshared_test.py``. This may also serve as a template for proper use of these modules from within python.

.. _modellib:

Modellib
~~~~~~~~

This repeats the functionality of the module of the same name in ARTIST. I.e., it offers a library of model templates, each of which can be tailored both by choice of parameters and by augmenting their grid-value functions.

.. warning:: The library of bespoke models in LIME has undergone very little testing. The ARTIST code was written by various people of varying ability; all that has been done here is to port that C++ code to C. Until someone gives it a thorough testing and bug-cleaning it should be regarded with suspicion.

As well as the bespoke models from the ARTIST version, you can also supply your own grid-value functions via a python file, in similar fashion to both traditional ``lime`` and ``pylime``. The call which directs ``modellib`` to use this is ``modellib.setUserModel("<name of the model file>")``. A template model file is available at ``<LIME base dir>/example/model_pyshared.py``.

.. note:: If you supply such a model file, you should only include grid-value functions, not parameters. You will see e.g. that the template file ``<LIME base dir>/example/model_pyshared.py`` has no ``input()`` function. For this flavour of LIME, parameters should be set from within python. See the test script ``<LIME base dir>/tests/pyshared_test.py`` for examples of how this is done.

The library of model templates:
*******************************

- **allen03a**: from 'Allen et al. 2003, ApJ, 599, 351'.

- **BonnorEbert56**: from 'Bonnor 1956, MNRAS, 116, 351 ; Ebert 1955, ZA (Zeitschrift fuer Astrophysik), 37, 217'.

- **CG97**: from 'Chiang & Goldreich 1997, ApJ, 490, 368'.

- **DDN01**: from 'Dullemond & Dominik 2001, ApJ, 560, 957'.

- **LiShu96**: from 'Li & Shu 1996, ApJ, 472, 211'.

- **Mamon88**: from 'Mamon et al. 1988, ApJ 328, 797'.

- **Mendoza09**: from 'Mendoza, Tejeda & Nagel, 2009, MNRAS, 393, 579'.

- **Shu77**: from 'Shu 1977, ApJ, 214, 488'.

- **Ulrich76**: from 'Ulrich 1976, ApJ, 210, 377'.


3. Flavour 'casalime'
---------------------

The final flavour of LIME offers similar functionality to ``pyshared``, but is designed to be used from the CASA command line. Originally the ``pyshared`` modules were used for this, but due to stupid clashes in threading and cfitsio, it was decided to redesign the CASA interface so that it launched LIME in a new process.

Compiling
~~~~~~~~~

::

    cd <LIME base dir>
    make casalime

This generates an executable called ``casalime``. As with ``lime`` and ``pylime`` flavours, you will want to make sure that ``<LIME base dir>`` is in your PATH environment variable, so CASA can find this executable. Also do

::

    source ./pylimerc.csh # or the bash equivalent

Testing
~~~~~~~

There is a test script ``<LIME base dir>/tests/casalime_test.py`` for checking that the tasks built ok.

CASA-specific compilation
~~~~~~~~~~~~~~~~~~~~~~~~~

The actual tasks which you run on the CASA command-line are called ``limesolver`` and ``raytrace``. More on how to use those :ref:`below <casa_tasks>`. For the moment we just want to get them running.

You will find the following four files under ``<LIME base dir>/casa``:

::

    limesolver.xml
    raytrace.xml
    task_limesolver.py
    task_raytrace.py

You can leave them there for the next step, but it is neater if you copy them somewhere else, to some convenient working directory. Suppose you have done that. CD to that working directory and invoke ``buildmytasks`` from the CASA distro you plan to use. That should generate the following new files:

::

    limesolver_cli.py
    limesolver.py
    mytasks.py
    raytrace_cli.py
    raytrace.py

The final step is to make sure that CASA can find these files when you start it up. If you don't already have a file ``~/.casa/init.py``, create one. Add the following line to it:

::

    execfile("<location of your task_* etc modules>/mytasks.py")

Once you've done that, you should be able to start CASA from anywhere and run the tasks ``limesolver`` and ``raytrace`` successfully.


.. _casa_tasks:

CASA tasks
~~~~~~~~~~

The CASA interface for setting task parameter values is not a very good tool for expressing the complicated and interrelated set of LIME parameters. Mostly for this reason, two simplifications have been made to flavour ``casalime``: the LIME functionality has been split between two tasks ``limesolver`` and ``raytrace``, and only 1 image at a time can be produced.

limesolver
**********

This generates the grid and solves the radiative transfer equations. It's not the job of ``limesolver`` to make images.

CASA tasks store parameter values via INP files. A template INP file is available at ``<LIME base dir>/casa/limesolver.template``. If you copy this, together with the files

::

    hco+@xpol.dat
    jena_thin_e6.tab
    model_pyshared.py

from ``<LIME base dir>/example`` to the directory you want to run CASA from, then you should be able from the CASA command line to do

::

    execfile('limesolver.template')
    go

for a nominal run of ``limesolver``. The output will be found in the same directory in the FITS file ``grid_5_mymodel.ds``. This conforms in format to the description in the header of the module ``<LIME base dir>/src/grid2fits.c``.

You will recognize most of the early parameters from LIME but those following ``modelID`` all pertain to :ref:`modellib <modellib>`.

raytrace
********

This task reads the grid file created by ``limesolver`` and makes a (single) image.

Once again there is a template INP file available: ``<LIME base dir>/casa/raytrace.template``. Perusal of this shows that the parameters are similar to the LIME ones, but two boolean parameters ``rotationStyle`` and ``doLine`` have been added.



