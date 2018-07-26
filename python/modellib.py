#!/usr/bin/python
#
#  modellib.py
#  This file is part of LIME, the versatile line modeling engine.
#
#  The intent of this module is to allow the library of models to be accessed from within python.
#
#  Copyright (C) 2006-2014 Christian Brinch
#  Copyright (C) 2015-2017 The LIME development team
#
#TODO:
#	- Obtain as much information from the C code as possible to avoid too much duplication. (Or the other way around - i.e. pass it in via libmodellib.initialize()).
#	- Get hold of the various astronomical/physical constants somehow?
#

import libmodellib
import sys
import os
import modellib_classes as mc
import modellib_utils as mu
import limepar_classes as pc

_doTest = False

#_currentModel = None
#_copyTemp = 0 # -1 means tgas should be set to equal tdust, +1 means the other way, 0 means no action.
#userModuleNameNoPy = ''

libmodellib.initialize()#*********** could pass in all the template info to spare some duplication.

#.......................................................................
# Functions. First those that don't need to talk to the C code:

def getModelIDs():
  # return list of all recognized modelIDs
  return mc._modelsDict.keys()

def isRegisteredModel(modelID):
  if modelID in mc._modelsDict.keys():
    return True
  else:
    return False

def getModelName(modelID):
  return mc._modelsDict[modelID].name

def getModelDesc(modelID):
  return mc._modelsDict[modelID].desc

def getModelBibref(modelID):
  # returns Bibliographic reference of a Model
  return mc._modelsDict[modelID].bibref

def getParamIDs(modelID):
  return mc._modelsDict[modelID]._paramDict.keys()

def isRegisteredParam(modelID, paramID):
  if paramID in mc._modelsDict[modelID]._paramDict.keys():
    return True
  else:
    return False

def getParamName(modelID, paramID):
  return mc._modelsDict[modelID]._paramDict[paramID].name

def getParamDesc(modelID, paramID):
  returnmc. _modelsDict[modelID]._paramDict[paramID].desc

def getParamType(modelID, paramID):
  return mc._modelsDict[modelID]._paramDict[paramID].dType

def getParamUnit(modelID, paramID):
  return mc._modelsDict[modelID]._paramDict[paramID].unit

def _getParamDefVal(modelID, paramID, dTypeStr):
  if mc._modelsDict[modelID]._paramDict[paramID].dType != dTypeStr:
    raise ValueError("wrong type")
  return mc._modelsDict[modelID]._paramDict[paramID].defaultValue

def getParamDefValDouble(modelID, paramID):
  return _getParamDefVal(modelID, paramID, 'double')

def getParamDefValInt(modelID, paramID):
  return _getParamDefVal(modelID, paramID, 'int')

def getParamDefValString(modelID, paramID):
  return _getParamDefVal(modelID, paramID, 'string')

#def getParamEnumValues(modelID, paramID):
#  if _modelsDict[modelID]._paramDict[paramID].dType != 'enum':
#    raise 'wrong type'
#  return _modelsDict[modelID]._paramDict[paramID].value#*********************
#******* convert it to list and return [:]

#def getParamEnumDefIndex(modelID, paramID):
#  pass#**************************

def getResultIDs():
  return mc._resultsDict.keys()

def isRegisteredResult(resultID):
  if resultID in mc._resultsDict.keys():
    return True
  else:
    return False

def getResultName(resultID):
  return mc._resultsDict[resultID].idStr # not going to bother capitalizing the first letter to return a 'name'.

def getResultDesc(resultID):
  return mc._resultsDict[resultID].desc

def getResultUnit(resultID):
  returnmc. _resultsDict[resultID].unit

def getResultType(resultID):
  return mc._resultsDict[resultID].dType

def getModelResultIDs(modelID):
  return mc._modelsDict[modelID]._listOfResultIDs[:]

def isResultModel(resultID, modelID):
  if resultID in mc._modelsDict[modelID]._listOfResultIDs:
    return True
  else:
    return False

def getFunctionIDs():
  return mc._functionsDict.keys()

def isRegisteredFunction(functionID):
  if functionID in mc._functionsDict.keys():
    return True
  else:
    return False

def getFunctionName(functionID):
  return mc._functionsDict[functionID].name

def getFunctionDesc(functionID):
  return mc._functionsDict[functionID].desc

def getFunctionType(functionID):
  return mc._functionsDict[functionID].type

#def getFunctionResultIDs(functionID):
#  # return list of resultIDs for this function
#  pass

#def getResultFunctionIDs(resultID):
#  # The IDs of all Functions that can be linked to a Result
#  pass

def isResultFunction(resultID, functionID):
  return mu.isResultFunction(resultID, functionID)

def getFunctionParamIDs(functionID):
  return mc._functionsDict[functionID]._argsDict.keys()

def isRegisteredFunctionParam(functionID, functionParamID):
  if functionParamID in mc._functionsDict[functionID]._argsDict.keys():
    return True
  else:
    return False

#def getFunctionParamName(functionID, functionParamID):
#  # return string giving parameter name
#  pass

#def getFunctionParamDesc(functionID, functionParamID):
#  pass

#def getFunctionParamType(functionID, functionParamID):
#  pass

#### there are a variety of bespoke functions for returning values of different parameter types.

def setParamDouble(paramID, value):
  mu.setParam(paramID, value, 'double')

def setParamInt(paramID, value):
  mu.setParam(paramID, value, 'int')

def setParamString(paramID, value):
  mu.setParam(paramID, value, 'string')

def setParamEnumIndex(paramID, index):
  mu.setParam(paramID, value, 'int')####################

def _getParam(paramID, dTypeStr):
  if mc._currentModel._paramDict[paramID].dType!=dTypeStr:
    raise ValueError("parameter data type is not %s" % (dTypeStr))
  return mc._currentModel._paramDict[paramID].value ######### None if it has not been set.

def getParamDouble(paramID):
  return _getParam(paramID, 'double')

def getParamInt(paramID):
  return _getParam(paramID, 'int')

def getParamString(paramID):
  return _getParam(paramID, 'string')

def getParamEnumIndex(paramID):
  pass###################

def setFunction(resultID, functionID):
  mc._currentModel.funcDict[resultID] = mc._functionsDict[functionID].copy()

def unsetFunction(resultID):
  try:
    del mc._currentModel.funcDict[resultID]
  except KeyError:
    pass

def getFunctionID(resultID):
  return mc._currentModel.funcDict[resultID].idStr

###def _setFunctionParam(resultID, functionParamID, value, dTypeStr):
###  currentModelID = libmodellib.get_current_model_ID()
###  if _modelsDict[currentModelID].funcDict[resultID]._argsDict[functionParamID].dType!=dTypeStr:
###    raise 'parameter data type is not %s' % dTypeStr
###
###  _modelsDict[currentModelID].funcDict[resultID]._argsDict[functionParamID].value = value

def setFunctionParamDouble(resultID, functionParamID, value):
###  mu.setFunctionParam(resultID, functionParamID, value, 'double')
  mc._currentModel.funcDict[resultID]._argsDict[functionParamID]._value = value
  if _doTest: print 'In modellib.setFunctionParamDouble(). resultID=%s, functionParamID=%s, value=%e.' % (resultID, functionParamID, value)

#def setFunctionParamInt(resultID, functionParamID, value):
#  libmodellib.set_fn_param_int(functionID, functionParamID, value)

#def setFunctionParamString(resultID, functionParamID, value):
#  libmodellib.set_fn_param_str(functionID, functionParamID, value)

#def setFunctionParamEnumIndex(resultID, functionParamID, index):
#  pass

def getFunctionParamDouble(resultID, functionParamID):
  return mc._currentModel.funcDict[resultID]._argsDict[functionParamID]._value

#def getFunctionParamInt(resultID, functionParamID):
#  pass

#def getFunctionParamString(resultID, functionParamID):
#  pass

#def getFunctionParamEnumIndex(resultID, functionParamID):
#  pass

#def setTdustIdentTemp():
#  libmodellib.set_tdust_as_temp()

#def unsetTdustIdentTemp():
#  libmodellib.unset_tdust_as_temp()

#def isTdustIdentTemp():
#  pass

def setTempIdentTdust():
#  libmodellib.set_temp_as_tdust()
  mc._copyTemp = -1


#def unsetTempIdentTdust():
#  pass

#def isTempIdentTdust():
#  pass

def isCurrentResultModel(resultID):
  # Checks if a Result is provided by the Current Model
  if resultID in mc._currentModel._listOfResultIDs:
    return True
  else:
    return False

def isCurrentResultFunction(resultID):
  # Checks if a Result is linked to a Function
  pass

def isSetCurrentModel():
  if mc._currentModel is None:
    return False
  else:
    return True

def getCurrentModelID():
  return mc._currentModel.idStr

def setCurrentModel(modelID):
  mu.setCurrentModel(modelID)

def unsetCurrentModel():
  mc._currentModel = None

def setUserModel(userModelPath):
  return mu.setUserModel(userModelPath)

#.......................................................................
# Now define those functions that need to talk to the C code.

#def isConfigurationComplete():
#  # Checks if all Results are either provided by the Current Model or linked to a Function
#  if libmodellib.is_config_complete()>0:
#    return True
#  else:
#    return False

def finalizeConfiguration():
  mu.finalizeConfiguration()
  libmodellib.finalize_config((mc._currentModel,mc._userModuleNameNoPy,mc._copyTemp))

def isFinalizedConfiguration():
  # Checks if the configuration has been finalized.
  if libmodellib.is_config_finalized()>0:
    return True
  else:
    return False

def density(x, y, z):
  return libmodellib.density(x, y, z)

def doppler(x, y, z):
  return libmodellib.doppler(x, y, z)

def cleanUp():
  # This should be called after all work is finished.
  libmodellib.clean_up()

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if __name__ == '__main__':
  testI = 0

  if(testI==0):
    setCurrentModel('BonnorEbert56')

  else:
    print "Don't recognize test %d", testI

