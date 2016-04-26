# UV-heated outflow cavity walls
The Herschel Space Observatory detected widespread emission of highly rotationally excited CO (sometimes up to _J_=40-30 and higher) in embedded protostars from low to high mass. Two potential excitation mechanisms were suggested: stellar UV radiation and shocks. In Visser et al. (2012, hereafter V12), we presented LIME models for three low-mass protostars to analyse the contribution of UV excitation. In this scenario, the high-_J_ CO emission seen with Herschel arises in gas of a few 100 to a few 1000 K in the envelope-outflow interface, i.e., the walls of the outflow cavity.

This branch of the LIME GitHub repository contains the model setup for HH46 as used in the paper, along with modifications to the LIME source code, and updated to work with the current release of LIME. We make this available in the hope that all or part of it will be of use in future applications of LIME.


## Modifications to LIME
The general geometry is drawn in Fig. 2 of V12. The cavity wall is an abrupt change from envelope densities of 10<sup>4</sup>--10<sup>9</sup> cm<sup>-3</sup> to a cavity density of zero. This infinitely steep density gradient presented some difficulties in LIME. In addition, we wanted to use abundances from a chemical code, dust temperatures from RADMC, and gas temperatures from a PDR grid. All of this required a highly customized `model.c` file and modifications to the LIME source code. The following subsections provide some details.

#### Density in `model.c`
The envelope follows a standard power-law density profile. The cavity wall is an ellipsoid with a given semimajor and semiminor axis. The density inside the cavity is effectively zero. We passed the gas temperature as an extra parameter to the `density` function so that we could compute the thermal H<sub>2</sub> ortho/para ratio at each grid point.

#### Temperature in `model.c`
Dust temperatures and UV fluxes were pre-computed for a given stellar luminosity using Kees Dullemond's RADMC continuum radiative transfer code. The RADMC output was converted to `dusttemp.h` and `uvflux.h`, each containing 1D arrays for _r_ and _θ_ and a 2D array for _T_<sub>dust</sub> or _F_<sub>UV</sub>. The `temperature` function reads the arrays and uses the value for the point nearest the provided _x_, _y_ and _z_. At one point we tried a linear interpolation between the four nearest points, but it had no effect on the CO line fluxes.

The gas and dust temperatures are decoupled in the UV-irradiated cavity walls, similar to what happens in a PDR. Temperature grids for PDRs are plentiful in the literature. We used the one from Kaufman et al. (1999), which provides gas temperatures at the surface of a PDR as a function of UV flux and gas density. Based on Röllig et al. (2007), the gas temperature decreases roughly as exp(-0.6 _A_<sub>V</sub>) until it equals the dust temperature. See Section 5.1 of V12 for details. The Kaufman et al. grid is encoded in `tgasgrid.k99.h`.

#### Abundance in `model.c`
The CO abundances are obtained from a pre-computed grid of abundances for a range of densities, temperatures, UV fluxes and extinctions. Before building the LIME grid, the function `buildGrid` in `grid.c` initiates the input file `grid_in.dat` for the external abundance grid interpolation program. This step requires the species name to be set in a new parameter `par->species`. The `abundance` function in `model.c` appends to `grid_in.dat` the physical conditions for each LIME grid point and sets the abundance temporarily to -1 as a flag to `buildGrid`. This flag activates a call to the external program `grid_interpol.x`, whose output `grid_out.dat` is read by LIME and assigned to `g[i].abun[0]` for each grid point.

The CO abundance grid is included in this repo as an ASCII file in the directory `abundance_grid`. The ASCII file needs to be converted to binary format for use by `grid_interpol.x`. The conversion source code and interpolation source code are both included as well. Command-line instructions for compilation and conversion are provided near the bottom of this readme.

#### Modifications to LIME: general approach
The modifications to LIME for the outflow cavity wall setup are controlled by the compiler flag `CAVITY_WALLS`, which in turn is activated by running LIME with command-line option `-w`, i.e., `lime -w model.c`.

#### Grid point sampling and weighting
The steep density gradient at the cavity wall required careful sampling of the grid in LIME. After much trial and error, we settled on 30,000 grid points in three sets of 10,000. The corresponding `par->sampling` value is 23.

All points are sampled logarithmically in radius. The first 10,000 points follow uniform sampling in _φ_ and are weighted by the angular distance to the cavity wall. The second set is placed exactly along the cavity wall and no weighting is used (all points are kept). The third set is placed in a thin layer right behind the cavity wall, also without weighting. This sampling and weighting scheme uses functions `angletocavity` and `cavity_phi` from `model.c`.

#### Grid smoothing
The number of smoothing steps is reduced from 20 to 5. With 20 smoothing steps, the grid points placed at and near the cavity wall move away too much and the wall is no longer well sampled. The 5 smoothing steps push some points into the cavity, where the density is set to effectively zero.

#### Raytracing
Even with the careful sampling and weighting, and the reduced number of smoothing steps, it is unavoidable that some Voronoi cells extend across the cavity/envelope boundary. During raytracing, suppose a given step along the ray has spatial coordinates inside the cavity. If these coordinates fall in a Voronoi cell whose corresponding grid point lies on the _other_ side of the wall, then this raytracing step picks up a non-zero gas density and a non-zero contribution to the line emission. This might be okay if the points are sampled equally on either side of the wall, but there are more points in the envelope than in the cavity. The end result is an artificial increase in emission.

Our solution was to check at each step whether the coordinates are within the cavity, using the function `exclude` in `model.c`. If so, no emission is added for this step.


## Preliminary step: abundance grid
```
rm -rf hh46-co/grid_interpol.x
rm -rf hh46-co/mol_grid_CO.bin
cd abundance_grid
make
./convert.sh
cd ../hh46-co
ln -s ../abundance_grid/mol_grid_CO.bin .
ln -s ../abundance_grid/grid_interpol.x .
cd ..
```


## Running LIME
```
cd hh46-co
lime -w model.c
```
