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

  if _doTest: print '>>> Entering task_raytrace.raytrace().'

  import cPickle
  import tempfile
  import subprocess as sub
  import signal
  import time

  try:
    import limepar_classes as pc
  except ImportError as err:
    casalog.post("Could not import module 'limepar_classes'", priority='ERROR')
    casalog.post("Error message: %s" % err, priority='ERROR')
    return
  except:
    casalog.post('Import of limepar_classes failed with an unrecognized exception.', priority='ERROR')
    return

  if _doTest: print '+++ In casaray version of raytrace(). aaa'

  AU = 1.49598e11    # AU to m
  PC = 3.08568025e16 # PC to m

  if _doTest: print '+++ In casaray version of raytrace(). bbb'

  sleepSec = 1
  gridDoneMessageIsPrinted = False
  mainDoneMessageIsPrinted = False
  raysDoneMessageIsPrinted = False
  nItersCounted = 0

  if _doTest: print '+++ In casaray version of raytrace(). ccc'

  fnDict = {'dust':dust,'gridInFile':gridInFile,'moldatfile':moldatfile}
  for variableName in fnDict.keys():
    fileName = fnDict[variableName]
    if not fileName is None and fileName!='':
      if not os.path.exists(fileName):
        casalog.post("Cannot find %s file %s" % (variableName, fileName), priority='ERROR')
        return

  if _doTest: print '+++ In casaray version of raytrace(). ddd'

  moldatfile = [moldatfile]

  if not(distUnit=='m' or distUnit=='pc' or distUnit=='au'):
    casalog.post("Distance unit %s not recognized." % (distUnit), priority='ERROR')
    return

  # Set input parameters for lime:
  limepars = pc.ModelParameters()
  if _doTest: print '+++ In casaray version of raytrace(). fff'

  limepars.dust              = dust
  limepars.gridInFile        = gridInFile
  limepars.polarization      = polarization
  limepars.nThreads          = nThreads
  limepars.nSolveIters       = 0
  limepars.traceRayAlgorithm = traceRayAlgorithm
  limepars.doSolveRTE        = False
  limepars.moldatfile        = moldatfile[:]
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

  limepars.img = images #**** rather images[:]??

  tempFileObj = tempfile.NamedTemporaryFile(prefix='pars_', delete=False)
  pklFileName = tempFileObj.name

  outputFileHandle = open(pklFileName, 'wb')
  cPickle.dump((None,'',0,limepars), outputFileHandle, -1)
  outputFileHandle.close()

  if _doTest: print '+++ In casaray version of raytrace(). jjj'

  casalog.post('Starting raytrace LIME run.')

  execName = "casalime"

  sigInterp = dict((getattr(signal, n), n) for n in dir(signal) if n.startswith('SIG') and '_' not in n )

  try:
    p = sub.Popen([execName, pklFileName], bufsize=0, stdout=sub.PIPE)
  except OSError:
    casalog.post('Cannot execute %s (not in PATH, or wrong permissions?)' % (execName), priority='ERROR')
    return
  except:
    casalog.post('Execution of %s failed with an unrecognized exception.' % (execName), priority='ERROR')
    return

  while True:
    stdoutdata = p.stdout.readline()[:-1]
    if not stdoutdata is None and stdoutdata!='':
      casalog.post(stdoutdata)

    status = p.poll()
    if not status is None:
      break

  if status==0:
    casalog.post('LIME raytrace run is complete.')
  elif status<0:
    casalog.post('LIME raytrace terminated by signal: %s' % sigInterp.get(-status, "unnamed signal: %d" % -status), priority='ERROR')
  else:
    casalog.post('LIME raytrace exited with nonzero status %s' % (status), priority='ERROR')

  if os.path.exists(pklFileName):
    os.unlink(pklFileName)

  # If we get to here, raytrace is finished or has crashed, so return.
  return

