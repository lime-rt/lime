import cPickle

def readPars(pklFileName):
  fileHandle = open(pklFileName, 'rb')
  (currentModel, userModuleNameNoPy, copyTemperatureI, limepars) = cPickle.load(fileHandle)
  # The elements in the pickle file are expected to have the following types:
  #	- currentModel should be a modellib_classes._Model instance.
  #	- userModuleNameNoPy should be a string.
  #	- copyTemperatureI should be an integer. -1 means tgas should be set to equal tdust, +1 means the other way, 0 means no action.
  #	- limepars should be a limepar_classes.ModelParameters instance.
  fileHandle.close()

  return (currentModel, userModuleNameNoPy, copyTemperatureI, limepars)

