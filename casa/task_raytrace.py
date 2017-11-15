import os

try:
  from taskinit import *
except ImportError:
  class DummyCasaLog:
    originStr=''

    def origin(self, originStr):
      self.originStr = originStr

    def post(self, errStr, priority='Message'):
      print "%s: %s" % (priority, errStr)

  casalog = DummyCasaLog()

def raytrace(gridInFile,moldatfile,dust,filename,imgres,pxls,unit,freq,rotationStyle
  ,theta,phi,incl,posang,azimuth,distance,nThreads,traceRayAlgorithm,doLine,nchan
  ,velres,trans,molI,bandwidth,source_vel,doInterpolateVels,polarization):

  # Note: rotationStyle and doLine are defined at this level purely to define some logical constraints on the subset of parameters listed; tey are not passed directly to LIME.

  casalog.post('Entering casaray version of raytrace().')

  import threading
  import time

  try:
    import casaray
  except ImportError as err:
    casalog.post("Cannot find module 'casaray'", priority='ERROR')
    return

  try:
    import par_classes as pc
  except ImportError as err:
    casalog.post("Cannot find module 'par_classes'", priority='ERROR')
    return

  class LimeExistentialLog:
    complete   = False
    errorState = False
    errorMessage = ''

  class LimeThread (threading.Thread):
    def __init__(self, par, limeLog):
      self._par = par
      self._limeLog = limeLog
      threading.Thread.__init__(self)

    def run(self):
      try:
        casaray.run(self._par)
      except Exception, e:
        self._limeLog.errorState = True
        self._limeLog.errorMessage = e

      self._limeLog.complete = True

  sleepSec = 2
  gridDoneMessageIsPrinted = False
  mainDoneMessageIsPrinted = False
  raysDoneMessageIsPrinted = False
  nItersCounted = 0

  fnDict = {'dust':dust,'gridInFile':gridInFile}
  for variableName in fnDict.keys():
    fileName = fnDict[variableName]
    if not fileName is None and fileName!='':
      if not os.path.exists(fileName):
        casalog.post("Cannot find %s file %s" % (variableName, fileName), priority='ERROR')
        return

  if not isinstance(moldatfile, list):
    moldatfile = [moldatfile]

  for i in range(len(moldatfile)):
    fileName = moldatfile[i]
    if not fileName is None and fileName!='':
      if not os.path.exists(fileName):
        casalog.post("Cannot find moldatfile[%d] file %s" % (i, fileName), priority='ERROR')
        return

  # Set input parameters for lime:
  limepar = pc.ModelParameters()

  limepar.dust              = dust
  limepar.gridInFile        = gridInFile
  limepar.polarization      = polarization
  limepar.nThreads          = nThreads
  limepar.nSolveIters       = 0
  limepar.traceRayAlgorithm = traceRayAlgorithm
  limepar.doSolveRTE        = False
  limepar.moldatfile        = moldatfile[:]

  # Define the images:
  images = []

  # Image #0:
  img = pc.ImageParameters()

  img.nchan             = nchan
  img.trans             = trans
  img.molI              = molI
  img.velres            = velres
  img.imgres            = imgres
  img.pxls              = pxls
  img.unit              = unit
  img.freq              = freq
  img.bandwidth         = bandwidth
  img.source_vel        = source_vel
  img.theta             = theta
  img.phi               = phi
  img.incl              = incl
  img.posang            = posang
  img.azimuth           = azimuth
  img.distance          = distance
  img.doInterpolateVels = doInterpolateVels
  img.filename          = filename

  images.append(img)

  limepar.img = images #**** rather images[:]??

  limeLog = LimeExistentialLog()

  # Run LIME in its own thread:
  limeThread = LimeThread(limepar, limeLog)
  limeThread.start()

  casalog.post('Starting LIME raytrace run.')

  # Loop with delay - every so often, check LIME status and print it to the casa logger.
  while not limeLog.complete:
    time.sleep(sleepSec)

  casalog.post('LIME raytrace run is complete')

  # If we get to here, LIME is finished or has crashed, so return.
  return

