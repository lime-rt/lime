#!/usr/bin/python
#
#  lime.py
#  This file is part of LIME, the versatile line modeling engine.
#
#  The intent of this module is to allow LIME functionality to be accessed from within python.
#
#  Copyright (C) 2006-2014 Christian Brinch
#  Copyright (C) 2015-2017 The LIME development team
#
#TODO:

import liblime as lime
import limepar_classes as pc


_doTest = False

def createInputPars():
  par = pc.ModelParameters()  
  return par

def createImage():
  img = pc.ImageParameters()
  return img

def runLime(par, images=[], debug=False):
  if _doTest: print '>>> Entering lime.runLime().'
  par.img = images #**** rather images[:]??
  if _doTest and len(images)>0: print 'xxx In lime.runLime(). par.img[0].molI=', par.img[0].molI
  lime.run_wrapper(par) # technically returns 0 when it succeeds, NULL otherwise.
  if _doTest: print '<<< Leaving lime.runLime().'

def setSilent(silent):
  lime.set_silent(silent) # technically returns 0 when it succeeds, NULL otherwise.

def getSilent():
  return lime.get_silent()

class ProgStatus(object):
  """
float progressGridBuilding
Progress of LIME's grid building (0 to 1)
float progressGridSmoothing
Progress of LIME's grid smoothing (0 to 1)
integer statusGrid
1 if LIME's grid is complete, 0 else
float progressConvergence
Progress of LIME's convergence (0 to 1)
integer numberIterations
double minsnr
double median
float progressPhotonPropagation
Progress of LIME's photon propagation (0 to 1)
float progressRayTracing
Progress of LIME's ray tracing (0 to 1)
integer statusRayTracing
1 if LIME's ray tracing is complete, 0 else
integer statusGlobal
1 if LIME's run is complete, 0 else
integer error
1 if an internal LIME error occured, 0 else. See messaeg for an error message
string message
  """
  def __init__(self):
    self.progressGridBuilding      = 0.0
    self.progressGridSmoothing     = 0.0
    self.statusGridBuilding        = 0
    self.progressConvergence       = 0.0
    self.numberIterations          = 0
    self.minsnr                    = 0.0
    self.median                    = 0.0
    self.progressPhotonPropagation = 0.0
    self.progressRayTracing        = 0.0
    self.statusRayTracing          = 0
    self.statusGlobal              = 0
    self.error                     = 0
    self.message                   = ''

  def _readStatus(self):
    statusTuple = lime.read_status()

    # Set all the attributes from the tuple.
    #
    self.progressGridBuilding      = statusTuple[0]
    self.progressGridSmoothing     = statusTuple[1]
    self.statusGridBuilding        = statusTuple[2]
    self.progressConvergence       = statusTuple[3]
    self.numberIterations          = statusTuple[4]
    self.minsnr                    = statusTuple[5]
    self.median                    = statusTuple[6]
    self.progressPhotonPropagation = statusTuple[7]
    self.progressRayTracing        = statusTuple[8]
    self.statusRayTracing          = statusTuple[9]
    self.statusGlobal              = statusTuple[10]
    self.error                     = statusTuple[11]
    self.message                   = statusTuple[12]

def getStatus():
  status = ProgStatus()
  status._readStatus()
  return status

def clearMessage():
  pass


