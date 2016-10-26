#!/usr/bin/env python

#
#  par_classes.py
#  This file is part of LIME, the versatile line modeling engine
#
#  Copyright (C) 2006-2014 Christian Brinch
#  Copyright (C) 2015-2016 The LIME development team
#
#
# DO NOT ALTER THIS FILE UNLESS YOU KNOW WHAT YOU ARE DOING!

_NTHREADS = 1
_DEFAULT_ANGLE = -999.0

class ImageParameters:
  """
This is to define the complete list of 'image' parameters which the user can set for each image. Although the parameters are set up as class attributes when the module is imported, the fundamental definition is via the 'hidden' attribute _listOfAttrs. This list includes a tuple element for each parameter. Each tuple has 5 elements: the parameter name, its type (one of ['int','float','bool','str','obj']), whether it is a list, whether it is mandatory, and the default value.

***NOTE*** that the ordering of the elements in _listOfAttrs is important - the code in pymodel_wrap.initParImg() depends on it. Don't change it without also changing the relevant code.
  """
  _listOfAttrs = []
  #                     name          type    list?  mand?  default
  _listOfAttrs.append(('nchan',      'int',   False, False, 0))
  _listOfAttrs.append(('trans',      'int',   False, False, -1))
  _listOfAttrs.append(('molI',       'int',   False, False, -1))
  _listOfAttrs.append(('velres',     'float', False, False, -1.0))
  _listOfAttrs.append(('imgres',     'float', False, True,  -1.0))
  _listOfAttrs.append(('pxls',       'int',   False, True,  -1))
  _listOfAttrs.append(('unit',       'int',   False, True,  -1))
  _listOfAttrs.append(('freq',       'float', False, False, -1.0))
  _listOfAttrs.append(('bandwidth',  'float', False, False, -1.0))
  _listOfAttrs.append(('source_vel', 'float', False, False, 0.0))
  _listOfAttrs.append(('theta',      'float', False, False, 0.0))
  _listOfAttrs.append(('phi',        'float', False, False, 0.0))
  _listOfAttrs.append(('incl',       'float', False, False, _DEFAULT_ANGLE))
  _listOfAttrs.append(('posang',     'float', False, False, _DEFAULT_ANGLE))
  _listOfAttrs.append(('azimuth',    'float', False, False, _DEFAULT_ANGLE))
  _listOfAttrs.append(('distance',   'float', False, True,  -1.0))

  _listOfAttrs.append(('filename',   'str',   False, True,  ''))

for attr in ImageParameters._listOfAttrs:
  ImageParameters.__dict__[attr[0]] = attr[4]

class ModelParameters:
  """
This is to define the complete list of 'ordinary' parameters which the user can set. Although the parameters are set up as class attributes when the module is imported, the fundamental definition is via the 'hidden' attribute _listOfAttrs. This list includes a tuple element for each parameter. Each tuple has 5 elements: the parameter name, its type (one of ['int','float','bool','str','obj']), whether it is a list, whether it is mandatory, and the default value. 

***NOTE*** that the ordering of the elements in _listOfAttrs is important - the code in pymodel_wrap.initParImg() depends on it. Don't change it without also changing the relevant code.
  """
  _listOfAttrs = []
  #                     name                type   list?  mand?  default
  _listOfAttrs.append(('radius',           'float',False, True,  0.0))
  _listOfAttrs.append(('minScale',         'float',False, True,  0.0))
  _listOfAttrs.append(('tcmb',             'float',False, False, 2.728))
  _listOfAttrs.append(('sinkPoints',       'int',  False, True,  0))
  _listOfAttrs.append(('pIntensity',       'int',  False, True,  0))
  _listOfAttrs.append(('blend',            'bool', False, False, False))
  _listOfAttrs.append(('traceRayAlgorithm','int',  False, False, 0))
  _listOfAttrs.append(('sampling',         'int',  False, False, 2))
  _listOfAttrs.append(('lte_only',         'bool', False, False, False))
  _listOfAttrs.append(('init_lte',         'bool', False, False, False))
  _listOfAttrs.append(('antialias',        'int',  False, False, 1))
  _listOfAttrs.append(('polarization',     'bool', False, False, False))
  _listOfAttrs.append(('nThreads',         'int',  False, False, _NTHREADS))

  _listOfAttrs.append(('outputfile',       'str',  False, False, ''))
  _listOfAttrs.append(('binoutputfile',    'str',  False, False, ''))
  _listOfAttrs.append(('gridfile',         'str',  False, False, ''))
  _listOfAttrs.append(('pregrid',          'str',  False, False, ''))
  _listOfAttrs.append(('restart',          'str',  False, False, ''))
  _listOfAttrs.append(('dust',             'str',  False, False, ''))

  _listOfAttrs.append(('nMolWeights',      'float',True,  False, []))
  _listOfAttrs.append(('dustWeights',      'float',True,  False, []))
  _listOfAttrs.append(('collPartIds',      'int',  True,  False, []))
  _listOfAttrs.append(('moldatfile',       'str',  True,  False, []))
  _listOfAttrs.append(('img',              'obj',  True,  False, []))

for attr in ModelParameters._listOfAttrs:
  ModelParameters.__dict__[attr[0]] = attr[4]


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if __name__ == '__main__':

  fred = ModelParameters()

  print fred.sampling

  print fred._listOfAttrs[16][0]



