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

def clearMessage():
  pass


