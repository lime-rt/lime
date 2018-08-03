#!/usr/bin/python

# I've designed this as a quick-and-dirty test of the casalime executable.

import math
import task_limesolver as lsol
import task_raytrace as lray

_useUserModel = True

doSolve = True
doRay   = True

if _useUserModel:
  radius                  = 2000.0
  minScale                = 0.5
  distUnit                = 'au'
  tcmb                    = 2.728
  pIntensity              = 4000
  sinkPoints              = 3000
  samplingAlgorithm       = 0
  sampling                = 2
  lte_only                = False
  init_lte                = False
  nThreads                = 1
  nSolveIters             = 14
  moldatfile              = 'hco+@xpol.dat'
  dust                    = 'jena_thin_e6.tab'
  gridInFile              = ''
  gridOutFile             = 'grid_5_mymodel.ds'
  resetRNG                = False

  modelID                 = None
  T                       = None
  rhoc                    = None
  Mstar                   = None
  Rstar                   = None
  Tstar                   = None
  bgdens                  = None
  hph                     = None
  plsig1                  = None
  rin                     = None
  rout                    = None
  sig0                    = None
  mdisk                   = None
  soundspeed              = None
  h0                      = None
  ab0                     = None
  mdot                    = None
  tin                     = None
  ve                      = None
  mdota                   = None
  mu                      = None
  nu                      = None
  rc                      = None
  Tcloud                  = None
  age                     = None
  Rn                      = None

  userModelPath           = 'model_pyshared.py'

  abundance_func     = ''
  abundance_args     = []
  bmag_func          = ''
  bmag_args          = []
  density_func       = ''
  density_args       = []
  doppler_func       = ''
  doppler_args       = []
  tdust_func         = ''
  tdust_args         = []
  temperature_func   = ''
  temperature_args   = []
  velocity_func      = ''
  velocity_args      = []

else:
  radius                  = 800.0
  minScale                = 1.0
  distUnit                = 'au'
  tcmb                    = 2.728
  pIntensity              = 1000
  sinkPoints              = 1000
  samplingAlgorithm       = 0
  sampling                = 2
  lte_only                = False
  init_lte                = False
  nThreads                = 4
  nSolveIters             = 12
  moldatfile              = 'hco+@xpol.dat'
  dust                    = 'jena_thin_e6.tab'
  gridInFile              = ''
  gridOutFile             = 'grid_5_CG97.ds'
  resetRNG                = False

  modelID                 = 'CG97'
  T                       = None
  rhoc                    = None
  Mstar                   = 0.5
  Rstar                   = 2.0
  Tstar                   = 4000.0
  bgdens                  = 1.0e-4
  hph                     = 4.0
  plsig1                  = -1.0
  rin                     = 1.0
  rout                    = 100.0
  sig0                    = 0.01
  mdisk                   = None
  soundspeed              = None
  h0                      = None
  ab0                     = None
  mdot                    = None
  tin                     = None
  ve                      = None
  mdota                   = None
  mu                      = None
  nu                      = None
  rc                      = None
  Tcloud                  = None
  age                     = None
  Rn                      = None

  userModelPath           = ''

  abundance_func     = 'scalarConst'
  abundance_args     = [1.0e-9]
  bmag_func          = 'vectorConstR'
  bmag_args          = [0.0]
  density_func       = ''
  density_args       = []
  doppler_func       = 'scalarConst'
  doppler_args       = [100.0]
  tdust_func         = ''
  tdust_args         = []
  temperature_func   = ''
  temperature_args   = []
  velocity_func      = ''
  velocity_args      = []

if doSolve:
  lsol.limesolver(radius,minScale,distUnit,tcmb,sinkPoints,pIntensity,samplingAlgorithm,sampling
    ,lte_only,init_lte,nThreads,nSolveIters,moldatfile,dust,gridInFile,gridOutFile,resetRNG
    ,modelID,T,rhoc,Mstar,Rstar,Tstar,bgdens,hph,plsig1,rin,rout,sig0,mdisk,soundspeed,h0,ab0,mdot
    ,tin,ve,mdota,mu,nu,rc,Tcloud,age,Rn,userModelPath,abundance_func,abundance_args
    ,bmag_func,bmag_args,density_func,density_args,doppler_func,doppler_args,tdust_func,tdust_args
    ,temperature_func,temperature_args,velocity_func,velocity_args)

if _useUserModel:
  gridInFile              = 'grid_5_mymodel.ds'
  moldatfile              = 'hco+@xpol.dat'
  dust                    = 'jena_thin_e6.tab'

  filename                = 'image_casa_mymodel.fits'
  imgres                  = 0.1            # Resolution in arc seconds
  pxls                    = 100            # Pixels per dimension
  unit                    = 0
  freq                    = 0.0
  rotationStyle           = 0
  theta                   = 0.0
  phi                     = 0.0
  incl                    = 0.0
  posang                  = 0.0
  azimuth                 = 0.0
  distance                = 140
  distUnit                = 'pc'
  doLine                  = True
  nchan                   = 61
  velres                  = 500.0
  trans                   = 3
  molI                    = 0
  bandwidth               = 0.0
  source_vel              = 0.0

  nThreads                = 1
  traceRayAlgorithm       = 1
  doInterpolateVels       = True
  polarization            = False

else:
  gridInFile              = 'grid_5_CG97.ds'
  moldatfile              = 'hco+@xpol.dat'
  dust                    = 'jena_thin_e6.tab'

  filename                = 'image1_CG97_casa.fits'
  imgres                  = 0.02            # Resolution in arc seconds
  pxls                    = 101             # Pixels per dimension
  unit                    = 0
  freq                    = 0.0
  rotationStyle           = 0
  theta                   = 30.0
  phi                     = 0.0
  incl                    = 0.0
  posang                  = 0.0
  azimuth                 = 0.0
  distance                = 140
  distUnit                = 'pc'
  doLine                  = True
  nchan                   = 60
  velres                  = 100.0
  trans                   = 2
  molI                    = 0
  bandwidth               = 0.0
  source_vel              = 0.0

  nThreads                = 1
  traceRayAlgorithm       = 1
  doInterpolateVels       = True
  polarization            = False

if doRay:
  lray.raytrace(gridInFile,moldatfile,dust,filename,imgres,pxls,unit,freq,rotationStyle
    ,theta,phi,incl,posang,azimuth,distance,distUnit,nThreads,traceRayAlgorithm,doLine,nchan
    ,velres,trans,molI,bandwidth,source_vel,doInterpolateVels,polarization)

