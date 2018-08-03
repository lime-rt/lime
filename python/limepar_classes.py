#!/usr/bin/env python

#
#  limepar_classes.py
#  This file is part of LIME, the versatile line modeling engine
#
#  Copyright (C) 2006-2014 Christian Brinch
#  Copyright (C) 2015-2017 The LIME development team
#
#
# DO NOT ALTER THIS FILE UNLESS YOU KNOW WHAT YOU ARE DOING!

_NTHREADS = 1
_DEFAULT_ANGLE = -999.0
_LOCAL_CMB_TEMP = 2.72548

class ImageParameters:
  """
This is to define the complete list of 'image' parameters which the user can set for each image. Although the parameters are set up as class attributes when the module is imported, the fundamental definition is via the 'hidden' attribute _listOfAttrs. This list includes a tuple element for each parameter. Each tuple has 5 elements: the parameter name, its type (one of ['int','float','bool','str','obj']), whether it is a list, whether it is mandatory, and the default value.

***NOTE*** that the ordering of the elements in _listOfAttrs is important - the code in py_utils.c:initParImg() depends on it. Don't change the order here without also changing the relevant C code.
  """
  _listOfAttrs = []
  #                     name                type    list?  mand?  default
  _listOfAttrs.append(('nchan',            'int',   False, False, 0))
  _listOfAttrs.append(('trans',            'int',   False, False, -1))
  _listOfAttrs.append(('molI',             'int',   False, False, -1))
  _listOfAttrs.append(('velres',           'float', False, False, -1.0))
  _listOfAttrs.append(('imgres',           'float', False, True,  -1.0))
  _listOfAttrs.append(('pxls',             'int',   False, True,  -1))
  _listOfAttrs.append(('unit',             'int',   False, False,  0))
  _listOfAttrs.append(('freq',             'float', False, False, -1.0))
  _listOfAttrs.append(('bandwidth',        'float', False, False, -1.0))
  _listOfAttrs.append(('source_vel',       'float', False, False, 0.0))
  _listOfAttrs.append(('theta',            'float', False, False, 0.0))
  _listOfAttrs.append(('phi',              'float', False, False, 0.0))
  _listOfAttrs.append(('incl',             'float', False, False, _DEFAULT_ANGLE))
  _listOfAttrs.append(('posang',           'float', False, False, _DEFAULT_ANGLE))
  _listOfAttrs.append(('azimuth',          'float', False, False, _DEFAULT_ANGLE))
  _listOfAttrs.append(('distance',         'float', False, True,  -1.0))
  _listOfAttrs.append(('doInterpolateVels','bool',  False, False, False))

  _listOfAttrs.append(('filename',         'str',   False, True,  None))
  _listOfAttrs.append(('units',            'str',   False, False, None))

  def __str__(self, spaces=''):
    myStr  = spaces+'<%s instance.\n' % (self.__class__.__name__)
    for attr in self._listOfAttrs:
      paddedName = '%s%s' % (attr[0], ' '*(17-len(attr[0])))
      myStr += spaces+'  %s = %s\n' % (paddedName, str(getattr(self,attr[0])))
    myStr += spaces+'>\n'
    return myStr

for attr in ImageParameters._listOfAttrs:
  ImageParameters.__dict__[attr[0]] = attr[4]


class ModelParameters:
  """
This is to define the complete list of 'ordinary' parameters which the user can set. Although the parameters are set up as class attributes when the module is imported, the fundamental definition is via the 'hidden' attribute _listOfAttrs. This list includes a tuple element for each parameter. Each tuple has 5 elements: the parameter name, its type (one of ['int','float','bool','str','obj']), whether it is a list, whether it is mandatory, and the default value. 

***NOTE*** that the ordering of the elements in _listOfAttrs is important - the code in py_utils.c:initParImg() depends on it. Don't change the order here without also changing the relevant C code.
  """
  _listOfAttrs = []
  #                     name                 type   list?  mand?  default
  _listOfAttrs.append(('radius',            'float',False, True,  0.0))
  _listOfAttrs.append(('minScale',          'float',False, True,  0.0))
  _listOfAttrs.append(('pIntensity',        'int',  False, True,  0))
  _listOfAttrs.append(('sinkPoints',        'int',  False, True,  0))

  _listOfAttrs.append(('dust',              'str',  False, False, None))
  _listOfAttrs.append(('outputfile',        'str',  False, False, None))
  _listOfAttrs.append(('binoutputfile',     'str',  False, False, None))
  _listOfAttrs.append(('gridfile',          'str',  False, False, None))
  _listOfAttrs.append(('pregrid',           'str',  False, False, None))
  _listOfAttrs.append(('restart',           'str',  False, False, None))
  _listOfAttrs.append(('gridInFile',        'str',  False, False, None))

  _listOfAttrs.append(('collPartIds',       'int',  True,  False, []))
  _listOfAttrs.append(('nMolWeights',       'float',True,  False, []))
  _listOfAttrs.append(('dustWeights',       'float',True,  False, []))
  _listOfAttrs.append(('collPartMolWeights','float',True,  False, []))

  _listOfAttrs.append(('gridDensMaxValues', 'float',True,  False, []))
  _listOfAttrs.append(('gridDensMaxLoc',    'float',True,  False, []))

  _listOfAttrs.append(('tcmb',             'float',False, False, _LOCAL_CMB_TEMP))
  _listOfAttrs.append(('lte_only',         'bool', False, False, False))
  _listOfAttrs.append(('init_lte',         'bool', False, False, False))
  _listOfAttrs.append(('samplingAlgorithm','int',  False, False, 0))
  _listOfAttrs.append(('sampling',         'int',  False, False, 2))
  _listOfAttrs.append(('blend',            'bool', False, False, False))
  _listOfAttrs.append(('antialias',        'int',  False, False, 1))
  _listOfAttrs.append(('polarization',     'bool', False, False, False))
  _listOfAttrs.append(('nThreads',         'int',  False, False, _NTHREADS))
  _listOfAttrs.append(('nSolveIters',      'int',  False, False, 17))
  _listOfAttrs.append(('traceRayAlgorithm','int',  False, False, 0))
  _listOfAttrs.append(('resetRNG',         'bool', False, False, False))
  _listOfAttrs.append(('doSolveRTE',       'bool', False, False, False))

  _listOfAttrs.append(('gridOutFiles',     'str',  True,  False, []))
  _listOfAttrs.append(('moldatfile',       'str',  True,  False, []))
  _listOfAttrs.append(('girdatfile',       'str',  True,  False, []))
  _listOfAttrs.append(('collPartNames',    'str',  True,  False, []))

  _listOfAttrs.append(('img',              'obj',  True,  False, []))

  def __str__(self, spaces=''):
    myStr  = spaces+'<%s instance.\n' % (self.__class__.__name__)
    for attr in self._listOfAttrs:
      if attr[0]=='img':
        myStr += spaces+'  img:\n'
        for i in range(len(self.img)):
          myStr += spaces+'    image %d:\n' % i
          myStr += self.img[i].__str__(spaces+'    ')

      else:
        paddedName = '%s%s' % (attr[0], ' '*(17-len(attr[0])))
        myStr += spaces+'  %s = %s\n' % (paddedName, str(getattr(self,attr[0])))
    myStr += spaces+'>\n'
    return myStr

for attr in ModelParameters._listOfAttrs:
  ModelParameters.__dict__[attr[0]] = attr[4]

def loadTestParValues(inGridFileName, outImageName):
  limepars = ModelParameters()

  limepars.gridInFile         =  inGridFileName
  limepars.moldatfile         =  ['hco+@xpol.dat']
  limepars.dust               =  "jena_thin_e6.tab"
  limepars.traceRayAlgorithm  =  1
  limepars.polarization       =  False
  limepars.collPartIds        = [1]
  limepars.nMolWeights        = [1.0]
  limepars.nThreads           = 1
  limepars.nSolveIters        = 0
  limepars.doSolveRTE         = False

  impars = ImageParameters()

  impars.filename           =  outImageName
  impars.imgres             =  0.02
  impars.pxls               =  101
  impars.unit               =  0
  impars.freq               =  0.0
  impars.theta              =  30.0
  impars.phi                =  0.0
  impars.incl               =  0.0
  impars.posang             =  0.0
  impars.azimuth            =  0.0
  impars.distance           =  4.31995235e+18
  impars.nThreads           =  1
  impars.doLine             =  True
  impars.nchan              =  60
  impars.velres             =  100.0
  impars.trans              =  3
  impars.molI               =  0
  impars.bandwidth          =  0.0
  impars.source_vel         =  0.0
  impars.doInterpolateVels    =  True

  limepars.img = [impars]

  return limepars

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if __name__ == '__main__':

  fred = ModelParameters()

  print fred.sampling

  print fred._listOfAttrs[16][0]



