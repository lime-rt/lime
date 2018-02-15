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

def limesolver(radius,minScale,distUnit,tcmb,sinkPoints,pIntensity,samplingAlgorithm,sampling
  ,lte_only,init_lte,nThreads,nSolveIters,moldatfile,dust,gridInFile,gridOutFile,resetRNG
  ,modelID,T,rhoc,Mstar,Rstar,Tstar,bgdens,hph,plsig1,rin,rout,sig0,mdisk,soundspeed,h0,ab0,mdot
  ,tin,ve,mdota,mu,nu,rc,Tcloud,age,Rn,userModelPath,abundance_func,abundance_args
  ,bmag_func,bmag_args,density_func,density_args,doppler_func,doppler_args,tdust_func,tdust_args
  ,temperature_func,temperature_args,velocity_func,velocity_args):

  import threading
  import time

  try:
    import modellib as ml
  except ImportError as err:
    casalog.post("Cannot import module 'modellib'", priority='ERROR')
    return

  try:
    import lime
  except ImportError as err:
    casalog.post("Cannot import module 'lime'", priority='ERROR')
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

  moldatfile = [moldatfile]

  if not(distUnit=='m' or distUnit=='pc' or distUnit=='au'):
    casalog.post("Distance unit %s not recognized." % (distUnit), priority='ERROR')
    return

  if (not userModelPath is None) and userModelPath!='':
    casalog.post("Supply of a user model is not yet supported.", priority='Warning')

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
    ml.setParamDouble(   'cs', soundspeed)
    ml.setParamEnumIndex('h0', h0)

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
    ml.setParamDouble('T',    Tcloud)
    ml.setParamDouble('time', age)

  elif modelID=='Ulrich76':
    ml.setParamDouble('mdot',  mdota)
    ml.setParamDouble('mstar', Mstar)
    ml.setParamDouble('rc',    rc)

  elif modelID=='allen03a':
    ml.setParamDouble('Rn',  Rn)
    ml.setParamDouble('T',   Tcloud)
    ml.setParamDouble('age', age)
    ml.setParamDouble('cs',  soundspeed)

  else:
    errStr = 'Model %s cannot yet be processed.' % modelID
    casalog.post(errStr, priority='ERROR')
    return

  # Now set some of the functions which were not included in the model:
  #
  funcDict = {
     'abundance'   :{'func':abundance,    'args':abundance_args}
    ,'bmag'        :{'func':bmag,         'args':bmag_args}
    ,'density'     :{'func':density,      'args':density_args}
    ,'doppler'     :{'func':doppler,      'args':doppler_args}
    ,'tdust'       :{'func':tdust,        'args':tdust_args}
    ,'temperature' :{'func':temperature,  'args':temperature_args}
    ,'velocity'    :{'func':velocity,     'args':velocity_args}
             }

  for resultID in funcDict.keys():
    functionID = funcDict[resultID]['func']
    if functionID=='': continue

    if not ml.isRegisteredFunction(functionID):
      errStr = 'Function %s is not recognized.' % functionID
      casalog.post(errStr, priority='ERROR')
      return

    if not ml.isResultFunction(resultID, functionID):
      errStr = 'Function %s cannot be linked with result %s.' % (functionID, resultID)
      casalog.post(errStr, priority='ERROR')
      return

    if ml.isCurrentResultModel(resultID):
      errStr = "Result %s is provided by the current model. Choosing a function for this result will overwrite the model value." % (resultID)
      casalog.post(errStr, priority='Warning')

    argsRequired = ml.getFunctionParamIDs(functionID)
    numArgsRequired = len(argsRequired)
    reqArgsStr = '[%s]' % ','.join(argsRequired)

    try:
      numArgsSupplied = len(funcDict[resultID]['args'])
    except TypeError: # presumably because user has supplied a scalar for this value rather than a list.
      if not(functionID=='scalarConst' or functionID=='vectorConstR'):
        raise # we need a list.
      funcDict[resultID]['args'] = [funcDict[resultID]['args']]
      numArgsSupplied = len(funcDict[resultID]['args'])

    if numArgsSupplied!=numArgsRequired:
      errStr = 'Function %s requires %d arguments but you have supplied %d.' % (functionID, numArgsRequired, numArgsSupplied)
      casalog.post(errStr, priority='ERROR')
      return

    errStr = 'Loading function %s for result %s. Argument names:' % (functionID, resultID)
    casalog.post(errStr)
    casalog.post('  %s' % reqArgsStr)

    ml.setFunction(resultID, functionID)

    for i in range(numArgsRequired):
      ml.setFunctionParamDouble(resultID, argsRequired[i], funcDict[resultID]['args'][i])
  # End loop over results.

  if modelID=='CG97' or modelID=='DDN01':
    # Set temperature identical to t_dust:
    ml.setTempIdentTdust()

  if not ml.isConfigurationComplete():
    casalog.post("Not all results have been linked to models or functions.", priority='ERROR')
    for resultID in ml.getResultIDs():
      if not (ml.isCurrentResultModel(resultID) or ml.isCurrentResultFunction(resultID)):
        casalog.post("Result %s is provided neither by the model nor any function." % resultID, priority='ERROR')
    return

  # Done:
  ml.finalizeConfiguration()
  if not ml.isFinalizedConfiguration():
    casalog.post("Error in finalizing the model config.", priority='ERROR')
    return

  # Set input parameters for lime:
  limepars = lime.createInputPars()

  if distUnit=='pc':
    limepars.radius          = radius*PC
    limepars.minScale        = minScale*PC
  elif distUnit=='au':
    limepars.radius          = radius*AU
    limepars.minScale        = minScale*AU
  else:
    limepars.radius          = radius
    limepars.minScale        = minScale

  limepars.tcmb              = tcmb
  limepars.sinkPoints        = sinkPoints
  limepars.pIntensity        = pIntensity
  limepars.samplingAlgorithm = samplingAlgorithm
  limepars.sampling          = sampling
  limepars.lte_only          = lte_only
  limepars.init_lte          = init_lte
  limepars.nThreads          = nThreads
  limepars.nSolveIters       = nSolveIters
  limepars.moldatfile        = moldatfile[:]
  limepars.dust              = dust
#  limepars.outputfile        = "populations.pop"
#  limepars.binoutputfile     = "restart.pop"
#  limepars.gridfile          = "grid.vtk"
#  limepars.pregrid           = "pregrid.asc"
#  limepars.restart           = "restart.pop"
  limepars.gridInFile        = gridInFile
#  limepars.collPartIds        = [1]#[macros["CP_H2"]] # must be a list, even when there is only 1 item.
#  limepars.nMolWeights        = [1.0] # must be a list, even when there is only 1 item.
#  limepars.collPartNames      = ["H2"] # must be a list, even when there is only 1 item.
#  limepars.collPartMolWeights = [2.0159] # must be a list, even when there is only 1 item.
#  limepars.gridDensMaxValues = [1.0] # must be a list, even when there is only 1 item.
#  limepars.gridDensMaxLoc    = [[0.0,0.0,0.0]] # must be a list, each element of which is also a list with 3 entries (1 for each spatial coordinate).
  limepars.gridOutFiles      = ['','','','',gridOutFile]
  limepars.resetRNG          = resetRNG
  if limepars.moldatfile is None or limepars.moldatfile=='' or len(limepars.moldatfile)<=0:
    limepars.doSolveRTE = False
  else:
    limepars.doSolveRTE = True

  # Define an empty set of images:
  images = []

  limeLog = LimeExistentialLog()

  lime.setSilent(False)

  # Run LIME in its own thread:
  limeThread = LimeThread(limepars, images, limeLog)
  limeThread.start()

  casalog.post('Starting limesolver LIME run.')

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

      if nItersCounted<limepars.nSolveIters:
        while nItersCounted<=limeStatus.numberIterations:
          nItersCounted += 1
          casalog.post('Iteration %d/%d' % (nItersCounted, limepars.nSolveIters))
          casalog.post('Min SNR %e  median %e' % (limeStatus.minsnr, limeStatus.median))

  limeStatus = lime.getStatus()
  if limeStatus.error!=0:
    casalog.post(limeStatus.message, priority='ERROR')

  elif limeLog.errorState:
    casalog.post(limeLog.errorMessage, priority='ERROR')

  # If we get to here, LIME is finished or has crashed, so return.
  return

