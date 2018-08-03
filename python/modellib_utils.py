#!/usr/bin/python
#
#  modellib_utils.py
#  This file is part of LIME, the versatile line modeling engine.
#
#  This module is intended to provide sufficient direct access to 'modellib_classes' to allow the CASA task 'limesolver' to run. Note that modellib.py cannot be used for this purpose because that assumes the existence of the shared library libmodellib.so, which however is not built by the 'casa' make target.
#
#  Copyright (C) 2006-2014 Christian Brinch
#  Copyright (C) 2015-2017 The LIME development team
#
#TODO:
#

import sys
import os
import modellib_classes as mc

_doTest = False

def isResultFunction(resultID, functionID):
  # Checks if a Function can be linked to a Result
  status = True # default value

  if      mc._resultsDict[resultID].isScalar  and mc._functionsDict[functionID].typeStr!='scalar'\
  or (not mc._resultsDict[resultID].isScalar) and mc._functionsDict[functionID].typeStr=='scalar':
    status = False

  return status

def setParam(paramID, value, dTypeStr):
  if mc._currentModel._paramDict[paramID].dType!=dTypeStr:
    raise 'parameter data type is not %s' % dTypeStr
  mc._currentModel._paramDict[paramID].value = value

def setUserModel(userModelPath):
  # Test it exists, add its location to sys.path and trim off the dirs and the *.py.

  if not os.path.exists(userModelPath):
    return False

  (path, filename) = os.path.split(os.path.abspath(userModelPath))

  (mc._userModuleNameNoPy, pyExt) = os.path.splitext(filename)

  if pyExt!='.py': # file doesn't have the right suffix. 
    return False

  sys.path.append(path)

  if _doTest: print 'In modellib_utils.setUserModel(). userModelPath=%s, mc._userModuleNameNoPy set to %s.' % (userModelPath, mc._userModuleNameNoPy)

  return True

def setCurrentModel(modelID):
  if mc._currentModel is None or mc._currentModel.idStr!=modelID:
    mc._currentModel = mc._modelsDict[modelID].copy()

def finalizeConfiguration():
  if not mc._currentModel is None:
    for paramName in mc._currentModel._paramDict.keys():
      if mc._currentModel._paramDict[paramName]._value is None: 
        mc._currentModel._paramDict[paramName]._value = mc._currentModel._paramDict[paramName].defaultValue

  # The python code which imports the current module should set mc._userModuleNameNoPy as desired.
##**** tests for format/existence of mc._userModuleNameNoPy
  if _doTest: print 'In modellib_utils.finalizeConfiguration(). mc._userModuleNameNoPy is', mc._userModuleNameNoPy
  if _doTest: print 'In modellib_utils.finalizeConfiguration(). mc._currentModel is\n', mc._currentModel
  if _doTest: print 'In modellib_utils.finalizeConfiguration(). mc._copyTemp is', mc._copyTemp

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if __name__ == '__main__':
  pass

