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
  ,theta,phi,incl,posang,azimuth,distance,distUnit,nThreads,traceRayAlgorithm,doLine,nchan
  ,velres,trans,molI,bandwidth,source_vel,doInterpolateVels,polarization):

  # Note: rotationStyle and doLine are defined at this level purely to define some logical constraints on the subset of parameters listed; tey are not passed directly to LIME.

  _doTest = False

  casalog.post('Entering casaray version of raytrace().')
  if _doTest: print '>>> Entering casaray version of raytrace().'

  import threading
  import time

  try:
    import casaray
  except ImportError as err:
    casalog.post("Cannot find module 'casaray'", priority='ERROR')
    return
  if _doTest: print '+++ In casaray version of raytrace(). aaa'

  try:
    import par_classes as pc
  except ImportError as err:
    casalog.post("Cannot find module 'par_classes'", priority='ERROR')
    return
  if _doTest: print '+++ In casaray version of raytrace(). bbb'

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

  AU = 1.49598e11    # AU to m
  PC = 3.08568025e16 # PC to m

  sleepSec = 1
  gridDoneMessageIsPrinted = False
  mainDoneMessageIsPrinted = False
  raysDoneMessageIsPrinted = False
  nItersCounted = 0

  fnDict = {'dust':dust,'gridInFile':gridInFile,'moldatfile':moldatfile}
  for variableName in fnDict.keys():
    fileName = fnDict[variableName]
    if not fileName is None and fileName!='':
      if not os.path.exists(fileName):
        casalog.post("Cannot find %s file %s" % (variableName, fileName), priority='ERROR')
        return

  if _doTest: print '+++ In casaray version of raytrace(). ccc'

  moldatfile = [moldatfile]
  if _doTest: print '+++ In casaray version of raytrace(). ddd (there is no eee)'

#  if _doTest: print '+++ In casaray version of raytrace(). eee'

  if not(distUnit=='m' or distUnit=='pc' or distUnit=='au'):
    casalog.post("Distance unit %s not recognized." % (distUnit), priority='ERROR')
    return

  # Set input parameters for lime:
  limepar = pc.ModelParameters()
  if _doTest: print '+++ In casaray version of raytrace(). fff'

  limepar.dust              = dust
  limepar.gridInFile        = gridInFile
  limepar.polarization      = polarization
  limepar.nThreads          = nThreads
  limepar.nSolveIters       = 0
  limepar.traceRayAlgorithm = traceRayAlgorithm
  limepar.doSolveRTE        = False
  limepar.moldatfile        = moldatfile[:]
  if _doTest: print '+++ In casaray version of raytrace(). ggg'

  # Define the images:
  images = []

  # Image #0:
  img = pc.ImageParameters()
  if _doTest: print '+++ In casaray version of raytrace(). hhh'

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

  if distUnit=='pc':
    img.distance        = distance*PC
  elif distUnit=='au':
    img.distance        = distance*AU
  else:
    img.distance        = distance

  img.doInterpolateVels = doInterpolateVels
  img.filename          = filename
  if _doTest: print '+++ In casaray version of raytrace(). iii'

  images.append(img)

  limepar.img = images #**** rather images[:]??

  limeLog = LimeExistentialLog()
  if _doTest: print '+++ In casaray version of raytrace(). jjj'

  # Run LIME in its own thread:
  limeThread = LimeThread(limepar, limeLog)
  if _doTest: print '+++ In casaray version of raytrace(). kkk'
  limeThread.start()
  if _doTest: print '+++ In casaray version of raytrace(). lll'

  casalog.post('Starting LIME raytrace run.')

  # Loop with delay - every so often, check LIME status and print it to the casa logger.
  while not limeLog.complete:
    time.sleep(sleepSec)

  casalog.post('LIME raytrace run is complete')

  # If we get to here, LIME is finished or has crashed, so return.
  return

