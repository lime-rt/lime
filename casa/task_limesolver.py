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

def limesolver(radius,minScale,tcmb,sinkPoints,pIntensity,samplingAlgorithm,sampling
  ,lte_only,init_lte,nThreads,nSolveIters,moldatfile,dust,gridInFile,gridOutFiles,resetRNG
  ,modelID,T,rhoc,Mstar,Rstar,Tstar,bgdens,hph,plsig1,rin,rout,sig0,mdisk,cs,h0,ab0,mdot
  ,tin,ve,mdota,mu,nu,rc,Tcloud,age,Rn):

  casalog.origin('limesolver')

  import threading
  import time

  try:
    import modellib as ml
  except ImportError as err:
    casalog.post("Cannot find module 'modellib'", priority='ERROR')
    return

  try:
    import lime
  except ImportError as err:
    casalog.post("Cannot find module 'lime'", priority='ERROR')
    return

  class LimeExistentialLog:
    complete   = False
    errorState = False
    errorMessage = ''

  class LimeThread (threading.Thread):
    def __init__(self, par, images, limeLog):
      self._par = par
      self._images = images
      self._limeLog = limeLog
      threading.Thread.__init__(self)

    def run(self):
      try:
        lime.runLime(self._par, self._images)
      except Exception, e:
        self._limeLog.errorState = True
        self._limeLog.errorMessage = e

      self._limeLog.complete = True

  sleepSec = 2#10
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

  if not ml.isRegisteredModel(modelID):
    errStr = '%s is not a registered model' % modelID
    casalog.post(errStr, priority='ERROR')
    return

  ml.setCurrentModel(modelID)

  if modelID=='BonnorEbert56':
    ml.setParamDouble('T',    T)
    ml.setParamDouble('rhoc', rhoc)

  elif modelID=='CG97':
    ml.setParamDouble('Mstar',  Mstar)
    ml.setParamDouble('Rstar',  Rstar)
    ml.setParamDouble('Tstar',  Tstar)
    ml.setParamDouble('bgdens', bgdens)
    ml.setParamDouble('hph',    hph)
    ml.setParamDouble('plsig1', plsig1)
    ml.setParamDouble('rin',    rin)
    ml.setParamDouble('rout',   rout)
    ml.setParamDouble('sig0',   sig0)

  elif modelID=='DDN01':
    ml.setParamDouble('Mstar',          Mstar)
    ml.setParamDouble('Rstar',          Rstar)
    ml.setParamDouble('Tstar',          Tstar)
    ml.setParamDouble('bgdens',         bgdens)
    ml.setParamString('dustopac_fname', dust)
    ml.setParamDouble('mdisk',          mdisk)
    ml.setParamDouble('plsig1',         plsig1)
    ml.setParamDouble('rin',            rin)
    ml.setParamDouble('rout',           rout)

  elif modelID=='LiShu96':
    ml.setParamDouble('cs', cs)
    ml.setParamDouble('h0', h0)

  elif modelID=='Mamon88':
    ml.setParamDouble('ab0',  ab0)
    ml.setParamDouble('mdot', mdot)
    ml.setParamDouble('rin',  rin)
    ml.setParamDouble('tin',  tin)
    ml.setParamDouble('ve',   ve)

  elif modelID=='Mendoza09':
    ml.setParamDouble('mdot',  mdota)
    ml.setParamDouble('mstar', Mstar)
    ml.setParamDouble('mu',    mu)
    ml.setParamDouble('nu',    nu)
    ml.setParamDouble('rc',    rc)

  elif modelID=='Shu77':
    ml.setParamDouble('Tcloud', Tcloud)
    ml.setParamDouble('time',   age)

  elif modelID=='Ulrich76':
    ml.setParamDouble('mdot',  mdota)
    ml.setParamDouble('mstar', Mstar)
    ml.setParamDouble('rc',    rc)

  elif modelID=='allen03a':
    ml.setParamDouble('Rn',  Rn)
    ml.setParamDouble('T',   Tcloud)
    ml.setParamDouble('age', age)
    ml.setParamDouble('cs',  cs)

  else:
    errStr = 'Model %s cannot yet be processed.' % modelID
    casalog.post(errStr, priority='ERROR')
    return

  # Now set some of the functions which were not included in the model:
  #
  if modelID=='Mendoza09' or modelID=='Ulrich76':
    # Constant temperature:
    ml.setFunction('temperature', 'scalarConst')
    ml.setFunctionParamDouble('temperature', 'val', 2.72)

#********* gas vs dust temps???

  if modelID=='CG97':
    # Set temperature identical to t_dust:
    ml.setTempIdentTdust()

  if modelID!='Mamon88':
    # Constant abundance:
    ml.setFunction('abundance', 'scalarConst')
    ml.setFunctionParamDouble('abundance', 'val', 1.e-9)

  # Constant doppler (100 m/s):
  ml.setFunction('doppler', 'scalarConst')
  ml.setFunctionParamDouble('doppler', 'val', 100.)

  if modelID!='LiShu96' and modelID!='allen03a':
    # Zero bmag:
    ml.setFunction('bmag', 'vectorConstR')
    ml.setFunctionParamDouble('bmag', 'val', 0.)

  # The remainder we will leave at defaults.

  # Done:
  ml.finalizeConfiguration()
  if not ml.isFinalizedConfiguration():
    casalog.post("Error in finalizing the model config.", priority='ERROR')
    return

  # Set input parameters for lime:
  par = lime.createInputPars()
  par.radius            = radius
  par.minScale          = minScale
  par.tcmb              = tcmb
  par.sinkPoints        = sinkPoints
  par.pIntensity        = pIntensity
  par.samplingAlgorithm = samplingAlgorithm
  par.sampling          = sampling
  par.lte_only          = lte_only
  par.init_lte          = init_lte
  par.nThreads          = nThreads
  par.nSolveIters       = nSolveIters
  par.moldatfile        = moldatfile[:]
  par.dust              = dust
#  par.outputfile        = "populations.pop"
#  par.binoutputfile     = "restart.pop"
#  par.gridfile          = "grid.vtk"
#  par.pregrid           = "pregrid.asc"
#  par.restart           = "restart.pop"
  par.gridInFile        = gridInFile
  par.collPartIds        = [1]#[macros["CP_H2"]] # must be a list, even when there is only 1 item.
  par.nMolWeights        = [1.0] # must be a list, even when there is only 1 item.
#  par.collPartNames     = ["phlogiston"] # must be a list, even when there is only 1 item.
#  par.collPartMolWeights = [2.0159] # must be a list, even when there is only 1 item.
#  par.gridDensMaxValues = [1.0] # must be a list, even when there is only 1 item.
#  par.gridDensMaxLoc    = [[0.0,0.0,0.0]] # must be a list, each element of which is also a list with 3 entries (1 for each spatial coordinate).
  par.gridOutFiles      = gridOutFiles[:]
  par.resetRNG          = resetRNG
  if par.moldatfile is None or par.moldatfile=='' or len(par.moldatfile)<=0:
    par.doSolveRTE = False
  else:
    par.doSolveRTE = True

  # Define an empty set of images:
  images = []

  limeLog = LimeExistentialLog()

  lime.setSilent(False)

  # Run LIME in its own thread:
  limeThread = LimeThread(par, images, limeLog)
  limeThread.start()

  casalog.post('Starting LIME run.')

  # Loop with delay - every so often, check LIME status and print it to the casa logger.
  while not limeLog.complete:
    time.sleep(sleepSec)
    limeStatus = lime.getStatus()

    if limeStatus.error!=0:
      break

    if limeStatus.statusGlobal!=0:
      casalog.post('LIME run is complete.')
      break

    if limeStatus.statusRayTracing==0:
      if nItersCounted<=0:
        if limeStatus.statusGridBuilding==0:
          pass
        elif not gridDoneMessageIsPrinted:
          gridDoneMessageIsPrinted = True
          casalog.post('Grid is complete. Starting solution.')

      if nItersCounted<par.nSolveIters:
        while nItersCounted<=limeStatus.numberIterations:
          nItersCounted += 1
          casalog.post('Iteration %d/%d' % (nItersCounted, par.nSolveIters))
          casalog.post('Min SNR %e  median %e' % (limeStatus.minsnr, limeStatus.median))

  limeStatus = lime.getStatus()
  if limeStatus.error!=0:
    casalog.post(limeStatus.message, priority='ERROR')

  elif limeLog.errorState:
    casalog.post(limeLog.errorMessage, priority='ERROR')

  # If we get to here, LIME is finished or has crashed, so return.
  return

