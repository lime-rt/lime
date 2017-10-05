#!/usr/bin/python

# I've designed this as a quick-and-dirty test of the casalime modules.

import math
import task_limesolver as lsol
import task_raytrace as lray

AU = 1.49598e11    # AU to m
PC = 3.08568025e16 # PC to m

doSolve=True
doRay=True

radius                  = 800 * AU
minScale                = 1. * AU
tcmb                    = 2.7
sinkPoints              = 1000
pIntensity              = 1000
samplingAlgorithm       = 0
sampling                = 2
lte_only                = False
init_lte                = False
nThreads                = 1
nSolveIters             = 12
moldatfile              = ['hco+@xpol.dat'] # must be a list, even when there is only 1 item.
dust                    = 'jena_thin_e6.tab'
gridInFile              = ''
gridOutFiles            = ['','','','','grid_5.ds']
resetRNG                = False

modelID			= 'CG97'
T			= None
rhoc			= None
Mstar                   = 0.5
Rstar                   = 2.0
Tstar                   = 4000.0
bgdens                  = 1.0e-4
hph                     = 4.0
plsig1                  = -1.0
rin                     = 1.0
rout                    = 100.0
sig0                    = 0.01
mdisk			= None
cs			= None
h0			= None
ab0			= None
mdot			= None
tin			= None
ve			= None
mdota			= None
mu			= None
nu			= None
rc			= None
Tcloud			= None
age			= None
Rn			= None


if doSolve:
  lsol.limesolver(radius,minScale,tcmb,sinkPoints,pIntensity,samplingAlgorithm,sampling
    ,lte_only,init_lte,nThreads,nSolveIters,moldatfile,dust,gridInFile,gridOutFiles,resetRNG
    ,modelID,T,rhoc,Mstar,Rstar,Tstar,bgdens,hph,plsig1,rin,rout,sig0,mdisk,cs,h0,ab0,mdot
    ,tin,ve,mdota,mu,nu,rc,Tcloud,age,Rn)

gridInFile              = 'grid_5.ds'
moldatfile              = ['hco+@xpol.dat'] # must be a list, even when there is only 1 item.
dust                    = 'jena_thin_e6.tab'

filename                = 'image0_casa.fits'
imgres                  = 0.02
pxls                    = 101
unit                    = 0
freq                    = 0.0
rotationStyle           = 0
theta                   = 30.0
phi                     = 0.0
incl                    = 0.0
posang                  = 0.0
azimuth                 = 0.0
distance                = 140 * PC
nThreads                = 1
traceRayAlgorithm       = 1
doLine                  = True
nchan                   = 60
velres                  = 100.0
trans                   = 3
molI                    = 0
bandwidth               = 0.0
source_vel              = 0.0
doInterpolateVels       = True
polarization            = False

if doRay:
  lray.raytrace(gridInFile,moldatfile,dust,filename,imgres,pxls,unit,freq,rotationStyle
    ,theta,phi,incl,posang,azimuth,distance,nThreads,traceRayAlgorithm,doLine,nchan
    ,velres,trans,molI,bandwidth,source_vel,doInterpolateVels,polarization)

