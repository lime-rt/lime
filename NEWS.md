NEWS
----

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


