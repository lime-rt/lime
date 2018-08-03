#!/usr/bin/python

# I've designed this as a quick-and-dirty test of the modules which are compiled by make target 'pyshared'.

import time
import math

import modellib as ml
import lime

AU = 1.49598e11    # AU to m
PC = 3.08568025e16 # PC to m
_useUserModel = True

t0 = time.time()

if _useUserModel:
  if not ml.setUserModel("model_pyshared.py"):
    raise ValueError("Could not set user model.")

else:
  # Select the Current Model:
  ml.setCurrentModel('CG97')

  # Set all parameters:
  ml.setParamDouble('Mstar', 0.5)     # Solar masses
  ml.setParamDouble('Rstar', 2.0)     # Solar radii
  ml.setParamDouble('Tstar', 4000.)   # K
  ml.setParamDouble('bgdens', 1.e-4)  # 1/cm^3
  ml.setParamDouble('hph', 4.0)
  ml.setParamDouble('plsig1', -1.0)
  ml.setParamDouble('rin', 1.0)       # AU
  ml.setParamDouble('rout', 100.0)    # AU
  ml.setParamDouble('sig0', 0.01)     # g/cm^2

  # Model CG97 does not provide the results 'abundance',
  # 'doppler', 'temperature, and 'bmag'. Hence these
  # Results have to be provided by Functions

  # Constant abundance (1.e-9):
  ml.setFunction('abundance', 'scalarConst')
  ml.setFunctionParamDouble('abundance', 'val', 1.e-9)

  # Constant doppler (100 m/s):
  ml.setFunction('doppler', 'scalarConst')
  ml.setFunctionParamDouble('doppler', 'val', 100.)

  # Zero bmag:
  ml.setFunction('bmag', 'vectorConstR')
  ml.setFunctionParamDouble('bmag', 'val', 0.)

  # Set temperature identical to t_dust:
  ml.setTempIdentTdust()

# Done setting up Model Library:
ml.finalizeConfiguration()

# Set input parameters for lime:
par = lime.createInputPars()

if _useUserModel:
  par.radius            = 2000.0*AU
  par.minScale          = 0.5*AU
  par.pIntensity        = 4000
  par.sinkPoints        = 3000
  par.dust              = "jena_thin_e6.tab"
  par.sampling          = 2 # Now only accessed if par.samplingAlgorithm==0 (the default).
  par.nSolveIters       = 14
  par.traceRayAlgorithm = 0
  par.gridOutFiles      = ['','','','',"grid_5_pylime.ds"]
  par.moldatfile        = ["hco+@xpol.dat"] # must be a list, even when there is only 1 item.

else:
  par.radius                  = 800.0 * AU
  par.minScale                = 1.0 * AU
  par.tcmb                    = 2.728
  par.pIntensity              = 1000
  par.sinkPoints              = 1000
  par.samplingAlgorithm       = 0
  par.sampling                = 2
  par.lte_only                = False
  par.init_lte                = False
  par.nThreads                = 4
  par.nSolveIters             = 12
  par.moldatfile              = ['hco+@xpol.dat'] # must be a list, even when there is only 1 item.
  par.dust                    = 'jena_thin_e6.tab'
  par.gridInFile              = ''
  par.gridOutFiles            = ['','','','',"grid_5_CG97.ds"]
  par.resetRNG                = False


# Define the images:
images = []

# Image #0: J=2-1
img = lime.createImage()

if _useUserModel:
  img.nchan             = 61             # Number of channels
  img.trans             = 3              # zero-indexed J quantum number
  img.velres            = 500.0          # Channel resolution in m/s
  img.imgres            = 0.1            # Resolution in arc seconds
  img.pxls              = 100            # Pixels per dimension
  img.unit              = 0              # 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau
  img.source_vel        = 0.0            # source velocity in m/s
  img.distance          = 140.0*PC       # source distance in m
  img.doInterpolateVels = True
  img.filename          = "image0_pylime.fits"  # Output filename

else:
  img.trans       = 2
  img.filename    = 'image0_CG97.fits'

  img.pxls        = 101
  img.imgres      = 0.02
  img.distance    = 140.0 * PC
  img.theta       = 30.0
  img.phi         = 0.0
  img.source_vel  = 0.0
  img.unit        = 0

  img.nchan       = 60
  img.velres      = 100.0

images.append(img)


# Run LIME:
lime.setSilent(False)

print "Running LIME:"
try:
    lime.runLime(par, images)
except Exception, e:
    print "Exception:", e

t1 = time.time()
print "Runtime: %ds" % (t1 - t0)
