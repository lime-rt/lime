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

  _doTest = False

  if _doTest: print '>>> Entering task_limesolver.limesolver()'

  import cPickle
  import tempfile
  import subprocess as sub
  import signal
  import time

  try:
    import modellib_classes as mc
  except ImportError as err:
    casalog.post("Could not import module 'modellib_classes'", priority='ERROR')
    casalog.post("Error message: %s" % err, priority='ERROR')
    return
  except:
    casalog.post('Import of modellib_classes failed with an unrecognized exception.', priority='ERROR')
    return

  try:
    import modellib_utils as mu
  except ImportError as err:
    casalog.post("Could not import module 'modellib_utils'", priority='ERROR')
    casalog.post("Error message: %s" % err, priority='ERROR')
    return
  except:
    casalog.post('Import of modellib_utils failed with an unrecognized exception.', priority='ERROR')
    return

  try:
    import limepar_classes as pc
  except ImportError as err:
    casalog.post("Could not import module 'limepar_classes'", priority='ERROR')
    casalog.post("Error message: %s" % err, priority='ERROR')
    return
  except:
    casalog.post('Import of limepar_classes failed with an unrecognized exception.', priority='ERROR')
    return

  if _doTest: print '+++ In task_limesolver.limesolver(). aaa'

  AU = 1.49598e11    # AU to m
  PC = 3.08568025e16 # PC to m

  if _doTest: print '+++ In task_limesolver.limesolver(). bbb'

  sleepSec = 1
  gridDoneMessageIsPrinted = False
  mainDoneMessageIsPrinted = False
  raysDoneMessageIsPrinted = False
  nItersCounted = 0

  if _doTest: print '+++ In task_limesolver.limesolver(). ccc'

  fnDict = {'dust':dust,'gridInFile':gridInFile,'moldatfile':moldatfile}
  for variableName in fnDict.keys():
    fileName = fnDict[variableName]
    if not fileName is None and fileName!='':
      if not os.path.exists(fileName):
        casalog.post("Cannot find %s file %s" % (variableName, fileName), priority='ERROR')
        return

  if _doTest: print '+++ In task_limesolver.limesolver(). ddd'

  moldatfile = [moldatfile]

  if not(distUnit=='m' or distUnit=='pc' or distUnit=='au'):
    casalog.post("Distance unit %s not recognized." % (distUnit), priority='ERROR')
    return

  if _doTest: print '+++ In task_limesolver.limesolver(). eee'

  if userModelPath!='':
    if _doTest: print '+++ In task_limesolver.limesolver(). eee-1'
    if not mu.setUserModel(userModelPath): # should test it exists, add its location to sys.path and trim off the dirs and the *.py.
      errStr = "User model %s doesn't exist or doesn't end in .py" % userModelPath
      casalog.post(errStr, priority='ERROR')
      return

  if _doTest: print '+++ In task_limesolver.limesolver(). eee-2. modelID=', modelID
  if not modelID is None:
    if not modelID in mc._modelsDict.keys():
      errStr = '%s is not a registered model' % modelID
      casalog.post(errStr, priority='ERROR')
      return

    mu.setCurrentModel(modelID)
    if _doTest: print '+++ In task_limesolver.limesolver(). eee-2 a'

    if modelID=='BonnorEbert56':
      mu.setParam('T',    T,    'double')
      mu.setParam('rhoc', rhoc, 'double')

    elif modelID=='CG97':
      mu.setParam('Mstar',  Mstar,  'double')
      mu.setParam('Rstar',  Rstar,  'double')
      mu.setParam('Tstar',  Tstar,  'double')
      mu.setParam('bgdens', bgdens, 'double')
      mu.setParam('hph',    hph,    'double')
      mu.setParam('plsig1', plsig1, 'double')
      mu.setParam('rin',    rin,    'double')
      mu.setParam('rout',   rout,   'double')
      mu.setParam('sig0',   sig0,   'double')

    elif modelID=='DDN01':
      mu.setParam('Mstar',          Mstar,  'double')
      mu.setParam('Rstar',          Rstar,  'double')
      mu.setParam('Tstar',          Tstar,  'double')
      mu.setParam('bgdens',         bgdens, 'double')
      mu.setParam('dustopac_fname', dust,   'string')
      mu.setParam('mdisk',          mdisk,  'double')
      mu.setParam('plsig1',         plsig1, 'double')
      mu.setParam('rin',            rin,    'double')
      mu.setParam('rout',           rout,   'double')

    elif modelID=='LiShu96':
      mu.setParam('cs', soundspeed, 'double')
      mu.setParam('h0', h0,         'int')

    elif modelID=='Mamon88':
      mu.setParam('ab0',  ab0,  'double')
      mu.setParam('mdot', mdot, 'double')
      mu.setParam('rin',  rin,  'double')
      mu.setParam('tin',  tin,  'double')
      mu.setParam('ve',   ve,   'double')

    elif modelID=='Mendoza09':
      mu.setParam('mdot',  mdota, 'double')
      mu.setParam('mstar', Mstar, 'double')
      mu.setParam('mu',    mu,    'double')
      mu.setParam('nu',    nu,    'double')
      mu.setParam('rc',    rc,    'double')

    elif modelID=='Shu77':
      mu.setParam('Tcloud', Tcloud, 'double')
      mu.setParam('time',   age,    'double')

    elif modelID=='Ulrich76':
      mu.setParam('mdot',  mdota, 'double')
      mu.setParam('mstar', Mstar, 'double')
      mu.setParam('rc',    rc,    'double')

    elif modelID=='allen03a':
      mu.setParam('Rn',  Rn,         'double')
      mu.setParam('T',   Tcloud,     'double')
      mu.setParam('age', age,        'double')
      mu.setParam('cs',  soundspeed, 'double')

    else:
      errStr = 'Model %s cannot yet be processed.' % modelID
      casalog.post(errStr, priority='ERROR')
      return

    if _doTest: print '+++ In task_limesolver.limesolver(). eee-2 b'

    # Now set some of the functions which were not included in the model:
    #
    funcDict = {
       'abundance'   :{'func':abundance_func,    'args':abundance_args}
      ,'bmag'        :{'func':bmag_func,         'args':bmag_args}
      ,'density'     :{'func':density_func,      'args':density_args}
      ,'doppler'     :{'func':doppler_func,      'args':doppler_args}
      ,'tdust'       :{'func':tdust_func,        'args':tdust_args}
      ,'temperature' :{'func':temperature_func,  'args':temperature_args}
      ,'velocity'    :{'func':velocity_func,     'args':velocity_args}
               }

    if _doTest: print '+++ In task_limesolver.limesolver(). eee-2 c'

    for resultID in funcDict.keys():
      if _doTest: print
      if _doTest: print '+++ In task_limesolver.limesolver(). eee-2 ca, resultID=', resultID

      functionID = funcDict[resultID]['func']
      if functionID=='': continue

      if _doTest: print '+++ In task_limesolver.limesolver(). eee-2 cb, functionID=', functionID

      if not functionID in mc._functionsDict.keys():
        errStr = 'Function %s is not recognized.' % functionID
        casalog.post(errStr, priority='ERROR')
        return

      if _doTest: print '+++ In task_limesolver.limesolver(). eee-2 cc'

      if not mu.isResultFunction(resultID, functionID):
        errStr = 'Function %s cannot be linked with result %s.' % (functionID, resultID)
        casalog.post(errStr, priority='ERROR')
        return

      if _doTest: print '+++ In task_limesolver.limesolver(). eee-2 cd'

      if resultID in mc._currentModel._listOfResultIDs:
        errStr = "Result %s is provided by the current model. Choosing a function for this result will overwrite the model value." #******** Will it? Should be tested.
        casalog.post(errStr, priority='Warning')

      if _doTest: print '+++ In task_limesolver.limesolver(). eee-2 ce'

      argsRequired = mc._functionsDict[functionID]._argsDict.keys()
      numArgsRequired = len(argsRequired)
      reqArgsStr = '[%s]' % ','.join(argsRequired)

      try:
        numArgsSupplied = len(funcDict[resultID]['args'])
      except TypeError: # presumably because user has supplied a scalar for this value rather than a list.
        if not(functionID=='scalarConst' or functionID=='vectorConstR'):
          raise # we need a list.
        funcDict[resultID]['args'] = [funcDict[resultID]['args']]
        numArgsSupplied = len(funcDict[resultID]['args'])

      if _doTest: print '+++ In task_limesolver.limesolver(). eee-2 cf'

      if numArgsSupplied!=numArgsRequired:
        errStr = 'Function %s requires %d arguments but you have supplied %d.' % (functionID, numArgsRequired, numArgsSupplied)
        casalog.post(errStr, priority='ERROR')
        return

      if _doTest: print '+++ In task_limesolver.limesolver(). eee-2 cg'

      errStr = 'Loading function %s for result %s. Argument names:' % (functionID, resultID)
      casalog.post(errStr)
      casalog.post('  %s' % reqArgsStr)

      if _doTest: print '+++ In task_limesolver.limesolver(). eee-2 ch'

      mc._currentModel.funcDict[resultID] = mc._functionsDict[functionID].copy()

      if _doTest: print '+++ In task_limesolver.limesolver(). eee-2 ci'

      for i in range(numArgsRequired):
        mc._currentModel.funcDict[resultID]._argsDict[argsRequired[i]]._value = funcDict[resultID]['args'][i]
    # End loop over results.

    if _doTest: print '+++ In task_limesolver.limesolver(). eee-2 d'

    if modelID=='CG97' or modelID=='DDN01':
      # Set temperature identical to t_dust:
      mc._copyTemp = -1

#** disable the following because in my version of modellib we can always fall back on defaults (and sometimes that is convenient).
#      if not ml.isConfigurationComplete():
#        casalog.post("Not all results have been linked to models or functions.", priority='ERROR')
#        for resultID in ml.getResultIDs():
#          if not (ml.isCurrentResultModel(resultID) or ml.isCurrentResultFunction(resultID)):
#            casalog.post("Result %s is provided neither by the model nor any function." % resultID, priority='ERROR')
#        return

    if _doTest: print '+++ In task_limesolver.limesolver(). fff'

    # Done:
    mu.finalizeConfiguration()
  # end if not modelID is None.

  if _doTest: print '+++ In task_limesolver.limesolver(). ggg'

  # Set input parameters for lime:
  limepars = pc.ModelParameters()

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
  limepars.gridInFile        = gridInFile
  limepars.gridOutFiles      = ['','','','',gridOutFile]
  limepars.resetRNG          = resetRNG
  if limepars.moldatfile is None or limepars.moldatfile=='' or len(limepars.moldatfile)<=0:
    limepars.doSolveRTE = False
  else:
    limepars.doSolveRTE = True

  # Define an empty set of images:
  limepars.img = []

  if _doTest: print '+++ In task_limesolver.limesolver(). hhh'

  if _doTest: print '+++ In task_limesolver.limesolver(). iii'

  tempFileObj = tempfile.NamedTemporaryFile(prefix='pars_', delete=False)
  pklFileName = tempFileObj.name

  outputFileHandle = open(pklFileName, 'wb')
  cPickle.dump((mc._currentModel,mc._userModuleNameNoPy,mc._copyTemp,limepars), outputFileHandle, -1)
  outputFileHandle.close()

  if _doTest: print '+++ In task_limesolver.limesolver(). jjj'

  casalog.post('Starting limesolver LIME run.')

  execName = "casalime"

  sigInterp = dict((getattr(signal, n), n) for n in dir(signal) if n.startswith('SIG') and '_' not in n )

  try:
    p = sub.Popen([execName, pklFileName], bufsize=1, stdout=sub.PIPE)
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
    casalog.post('LIME run is complete.')
  elif status<0:
    casalog.post('LIME terminated by signal: %s' % sigInterp.get(-status, "unnamed signal: %d" % -status), priority='ERROR')
  else:
    casalog.post('LIME exited with nonzero status %s' % (status), priority='ERROR')

  if os.path.exists(pklFileName):
    os.unlink(pklFileName)

  # If we get to here, LIME is finished or has crashed, so return.
  return

