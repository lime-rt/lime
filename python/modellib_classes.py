#!/usr/bin/python
#
#  modellib_classes.py
#  This file is part of LIME, the versatile line modeling engine
#
#  Copyright (C) 2006-2014 Christian Brinch
#  Copyright (C) 2015-2017 The LIME development team
#
#TODO:
#

_doTest = False

_currentModel = None
_copyTemp = 0 # -1 means tgas should be set to equal tdust, +1 means the other way, 0 means no action.
_userModuleNameNoPy = ''
#_userModuleNameNoPy = None

#.......................................................................
class _ModelParam:
  def __init__(self, idStr, name, desc, dType, defaultValue, unit):
    self.idStr        = idStr
    self.name         = name
    self.desc         = desc
    self.dType        = dType
    self.defaultValue = defaultValue
    self.unit         = unit

    self._value = None

  def __str__(self, spaces=''):
    myStr  = spaces+'<%s instance.\n' % (self.__class__.__name__)
    myStr += spaces+'  idStr        = %s\n' % (self.idStr)
#    myStr += spaces+'  name         = %s\n' % (self.name)
#    myStr += spaces+'  desc         = %s\n' % (self.desc)
    myStr += spaces+'  dType        = %s\n' % (self.dType)
#    myStr += spaces+'  defaultValue = %s\n' % (self.defaultValue)
#    myStr += spaces+'  unit         = %s\n' % (self.unit)

    if self._value is None:
      myStr += spaces+'  _value       = None\n'
    else:
      myStr += spaces+'  _value       = %e\n' % (self._value)

    myStr += spaces+'>\n'
    return myStr

  def copy(self):
    newParam = _ModelParam(self.idStr, self.name, self.desc, self.dType, self.defaultValue, self.unit)
    newParam._value = self._value
    return newParam

#.......................................................................
class _Result:
  def __init__(self, idStr, desc, unit, dType, numElements, isScalar=True):
    self.idStr        = idStr
    self.desc         = desc
    self.unit         = unit
    self.dType        = dType
    self.isScalar     = isScalar
    self.numElements  = numElements # None indicates >=1 but unspecified at this point.

  def copy(self):
    return _Result(self.idStr, self.desc, self.unit, self.dType, self.numElements, self.isScalar)

#.......................................................................
class _FuncParam:
  # They are all doubles.
  def __init__(self, name):
    self.name = name
    self._value = None

  def __str__(self, spaces=''):
    myStr  = spaces+'<%s instance.\n' % (self.__class__.__name__)
    myStr += spaces+'  name   = %s\n' % (self.name)

    if self._value is None:
      myStr += spaces+'  _value = None\n'
    else:
      myStr += spaces+'  _value = %e\n' % (self._value)

    myStr += spaces+'>\n'
    return myStr

  def copy(self):
    newFP = _FuncParam(self.name)
    newFP._value = self._value
    return newFP

#.......................................................................
class _Function:
  """
These are abstract procedures, each of which encodes a recipe for producing a single (either scalar or vector) output value for a given set of input arguments.
  """
  def __init__(self, idStr, name, desc, typeStr, argNamesList):
    self.idStr        = idStr
    self.name         = name
    self.desc         = desc
    self.typeStr      = typeStr # 'scalar' or 'vector'
    self.argNamesList = argNamesList

    self._argsDict = {} # argsDict[idStr] = _FuncParam
    for name in self.argNamesList:
      funcParam = _FuncParam(name)
      self._argsDict[name] = funcParam

  def __str__(self, spaces=''):
    myStr  = spaces+'<%s instance.\n' % (self.__class__.__name__)
    myStr += spaces+'  idStr   = %s\n' % (self.idStr)
#    myStr += spaces+'  name    = %s\n' % (self.name)
#    myStr += spaces+'  desc    = %s\n' % (self.desc)
    myStr += spaces+'  typeStr = %s\n' % (self.typeStr)
    myStr += spaces+'  Arguments from dict:\n'
    for name in self._argsDict.keys():
      myStr += self._argsDict[name].__str__(spaces+'  ')

    myStr += spaces+'>\n'
    return myStr

  def copy(self):
    newFunction = _Function(self.idStr, self.name, self.desc, self.typeStr, self.argNamesList[:])

    # Copy over any parameter values which have been filled:
    for key in self.argNamesList:
      newFunction._argsDict[key] = self._argsDict[key].copy()

    return newFunction

#.......................................................................
class _Model:
  """
Each model takes a list of _ModelParam values as inputs and produces a list of _Result values as output. Any _Result values missing may be supplied by _Function objects.
  """
  def __init__(self, idStr, name, desc, bibref, paramList, funcDict, resultList):
    self.idStr      = idStr
    self.name       = name
    self.desc       = desc
    self.bibref     = bibref
    self.paramList  = paramList # paramList[<i>] = _ModelParam
    self.funcDict   = funcDict # funcDict[resultName] = _Function
    self.resultList = resultList # list of _Result objects.

    self._listOfResultIDs = None
    self._paramDict  = None # _paramDict[paramName] = _ModelParam

    self._finishSetUp()

  def _finishSetUp(self):
    self._listOfResultIDs = []
    for result in self.resultList:
      self._listOfResultIDs.append(result.idStr)

    self._paramDict = {}
    for paramObj in self.paramList:
      self._paramDict[paramObj.idStr] = paramObj

  def __str__(self, spaces=''):
    myStr  = spaces+'<%s instance.\n' % (self.__class__.__name__)
    myStr += spaces+'  idStr  = %s\n' % (self.idStr)
    myStr += spaces+'  name   = %s\n' % (self.name)
#    myStr += spaces+'  desc   = %s\n' % (self.desc)
#    myStr += spaces+'  bibref = %s\n' % (self.bibref)
    myStr += spaces+'  Model params:\n'
    for i in range(len(self.paramList)):
      myStr += self.paramList[i].__str__(spaces+'  ')

#**** list of results given by the model

    myStr += spaces+'  Results linked to functions:\n'
    for resultName in self.funcDict.keys():
      myStr += spaces+'    resultName = %s\n' % (resultName)
      myStr += self.funcDict[resultName].__str__(spaces+'  ')

    myStr += spaces+'>\n'
    return myStr

  def copy(self):
    newParamList = []
    for paramObj in self.paramList:
      newParamList.append(paramObj.copy())

    newFuncDict = {}
    for key in self.funcDict.keys():
      newFuncDict[key] = self.funcDict[key].copy()

    newResultList = []
    for result in self.resultList:
      newResultList.append(result.copy())

    return _Model(self.idStr, self.name, self.desc, self.bibref, newParamList, newFuncDict, newResultList)

#.......................................................................
# Set up global dicts.

_resultsDict = {}
_resultsDict['abundance']     = _Result('abundance',   'Radiating species fractional number density', '', 'double', None)
_resultsDict['bmag']          = _Result('bmag',        'Magnetic field', 'tesla',                         'double', 3, False)
_resultsDict['density']       = _Result('density',     'Bulk gas number density', 'm^-3',                 'double', 1)
_resultsDict['doppler']       = _Result('doppler',     'Turbulent Doppler broadening', 'm/s',             'double', 1)
_resultsDict['tdust']         = _Result('tdust',       'Dust temperature', 'K',                           'double', 1)
_resultsDict['temperature']   = _Result('temperature', 'Gas temperature', 'K',                            'double', 1)
_resultsDict['velocity']      = _Result('velocity',    'Bulk gas velocity', 'm/s',                        'double', 3, False)

_functionsDict = {}
_functionsDict['scalarConst'] = _Function('scalarConst'
  , 'Const Scalar Function'
  , 'Constant scalar value (double)'
  , 'scalar'
  , ['val'])

_functionsDict['scalarPowerR'] = _Function('scalarPowerR'
  , 'Powerlaw R Scalar Function'
  , 'Scalar with power-law in spherical radius r (double)\nvalue = factor * r**exponent + offset (if r >= lowerR)\n      = factor * lowerR**exponent + offset (if r < lowerR)'
  , 'scalar'
  , ['lowerR', 'factor', 'exponent', 'offset'])

_functionsDict['scalarPowerRExpZ'] = _Function('scalarPowerRExpZ'
  , 'Powerlaw r exponential z Scalar Function'
  , 'Scalar with power-law in spherical radius r and exponential in z (double)\nvalue = valR * valZ with\nvalR = factR * r**expR + offsetR (if r >= lowerR)\n     = factR * lowerR**expR + offsetR (if r < lowerR)\nvalZ = factZ * e**(-abs(z) / scaleZ) + offsetZ'
  , 'scalar'
  , ['lowerR', 'factR', 'expR', 'offsetR', 'factZ', 'scaleZ', 'offsetZ'])

_functionsDict['scalarPowerRTheta'] = _Function('scalarPowerRTheta'
  , 'Powerlaw r theta Scalar Function'
  , 'Scalar with power-law in spherical radius r and theta (double)\nvalue = valR * valTheta with\nvalR = factR * r**expR + offsetR (if r >= lowerR)\n     = factR * lowerR**expR + offsetR (if r < lowerR)\nvalTheta = factTheta * theta**expTheta + offsetTheta (if theta >= lowerTheta)\n         =  factTheta * lowerTheta**expTheta + offsetTheta (if theta < lowerTheta)\nand with theta = acos( abs(z) / r )'
  , 'scalar'
  , ['lowerR', 'factR', 'expR', 'offsetR', 'lowerTheta', 'factTheta', 'expTheta', 'offsetTheta'])

_functionsDict['scalarPowerRZ'] = _Function('scalarPowerRZ'
  , 'Powerlaw r z Scalar Function'
  , 'Scalar with power-law in spherical radius r and z (double)\nvalue = valR * valZ with\nvalR = factR * r**expR + offsetR (if r >= lowerR)\n     = factR * lowerR**expR + offsetR (if r < lowerR)\nvalZ = factZ * abs(z)**expZ + offsetZ (if abs(z) >= lowerZ)\n     = factZ * abs(lowerZ)**expZ + offsetZ (if abs(z) < lowerZ)'
  , 'scalar'
  , ['lowerR', 'factR', 'expR', 'offsetR', 'lowerZ', 'factZ', 'expZ', 'offsetZ'])

_functionsDict['vectorConstR'] = _Function('vectorConstR'
  , 'Const radial vector Function'
  , 'Constant radial vector (double)'
  , 'vector'
  , ['val'])

_functionsDict['vectorConstXYZ'] = _Function('vectorConstXYZ'
  , 'Const cartesian vector Function'
  , 'Constant vector with x, y, and z component (double)'
  , 'vector'
  , ['valX', 'valY', 'valZ'])

_functionsDict['vectorDipole'] = _Function('vectorDipole'
  , 'Dipole vector Function'
  , 'Dipole field with dipole moment\n\nvalX = 3 * mDipole * sin(theta) * cos(theta) * cos(phi) / r**3\nvalY = 3 * mDipole * sin(theta) * cos(theta) * sin(phi) / r**3\nvalZ = 3 * mDipole * (cos(theta)**2 - 1/3) / r**3\nif r >= lowerR\n\nvalX = 3 * mDipole * sin(theta) * cos(theta) * cos(phi) / lowerR**3\nvalY = 3 * mDipole * sin(theta) * cos(theta) * sin(phi) / lowerR**3\nvalZ = 3 * mDipole * (cos(theta)**2 - 1/3) / lowerR**3\nif r < lowerR'
  , 'vector'
  , ['lowerR', 'mDipole'])

_functionsDict['vectorRadialPowerR'] = _Function('vectorRadialPowerR'
  , 'Powerlaw R radial vector Function'
  , 'Radial vector with power-law in spherical radius r (double)\nvalR = factor * r**exponent + offset (if r >= lowerR)\n     = factor * lowerR**exponent + offset (if r < lowerR)'
  , 'vector'
  , ['lowerR', 'factor', 'exponent', 'offset'])

_functionsDict['vectorRadialPowerRTheta'] = _Function('vectorRadialPowerRTheta'
  , 'Powerlaw r theta radial vector Function'
  , 'Radial vector with power-law in spherical radius r and theta (double)\nvalue = valR * valTheta with\nvalR = factR * r**expR + offsetR (if r >= lowerR)\n     = factR * lowerR**expR + offsetR (if r < lowerR)\nvalTheta = factTheta * theta**expTheta + offsetTheta (if theta >= lowerTheta)\n         =  factTheta * lowerTheta**expTheta + offsetTheta (if theta < lowerTheta)\nand with theta = acos( abs(z) / r )'
  , 'vector'
  , ['lowerR', 'factR', 'expR', 'offsetR', 'lowerTheta', 'factTheta', 'expTheta', 'offsetTheta'])

_functionsDict['vectorToroidalPowerR'] = _Function('vectorToroidalPowerR'
  , 'Powerlaw R toroidal vector Function'
  , 'Toroidal vector with power-law in cylindrical radius rho (double)\nvalPhi = factor * r**exponent + offset (if rho >= lowerRho)\n       = factor * lowerRho**exponent + offset (if rho < lowerRho)'
  , 'vector'
  , ['lowerRho', 'factor', 'exponent', 'offset'])


_modelsDict = {}

idStr = 'allen03a'
paramList = []
paramList.append(_ModelParam('Rn', 'Rn', 'Validity radius of kinematic approximation', 'double', 1.0, 'AU'))
paramList.append(_ModelParam('T', 'T', 'Temperature - constant', 'double', 30.0, 'K'))
paramList.append(_ModelParam('age', 'Age', 'Age of source', 'double', 100000.0, 'year'))
paramList.append(_ModelParam('cs', 'cs', 'Velocity of sound', 'double', 1000.0, 'm/s'))
resultList = [_resultsDict['bmag'],_resultsDict['density'],_resultsDict['temperature'],_resultsDict['velocity']]
_modelsDict[idStr] = _Model(idStr
  , 'Allen 03a'
  , 'Collapse of magnetized singular isothermal toroid (pivotal states)'
  , 'Allen et al. 2003, ApJ, 599, 351'
  , paramList, {}, resultList)

idStr = 'BonnorEbert56'
paramList = []
paramList.append(_ModelParam('T', 'T', 'Temperature of the core', 'double', 10.0, 'K'))
paramList.append(_ModelParam('rhoc', 'rhoc', 'Central volume density of the core', 'double', 1.0e6, '1/cm^3'))
resultList = [_resultsDict['density'],_resultsDict['temperature'],_resultsDict['velocity']]
_modelsDict[idStr] = _Model(idStr
  , 'Bonnor Ebert sphere'
  , 'Isothermal hydrostatic molecular cloud core by Bonnor 1956 & Ebert 1955'
  , 'Bonnor 1956, MNRAS, 116, 351 ; Ebert 1955, ZA (Zeitschrift fuer Astrophysik), 37, 217'
  , paramList, {}, resultList)

idStr = 'CG97'
paramList = []
paramList.append(_ModelParam('Mstar', 'Mstar', 'Mass of the central star', 'double', 0.5, 'Msun'))
paramList.append(_ModelParam('Rstar', 'Rstar', 'Radius of the central star', 'double', 2.0, 'Rsun'))
paramList.append(_ModelParam('Tstar', 'Tstar', 'Effective temperature of the star', 'double', 4000.0, 'K'))
paramList.append(_ModelParam('bgdens', 'bgdens', 'Background number density (floor value for the disk density)', 'double', 0.0001, '1/cm^3'))
paramList.append(_ModelParam('hph', 'hph', 'Ratio of the height of the disk atmosphere above the midplane and the pressure scale height (at all radius)', 'double', 4.0, ''))
paramList.append(_ModelParam('plsig1', 'plsig1', 'Power exponent of the radial surface density distribution', 'double', -1.0, ''))
paramList.append(_ModelParam('rin', 'rin', 'Inner Radius of the disk', 'double', 1.0, 'AU'))
paramList.append(_ModelParam('rout', 'rout', 'Outer radius of the disk', 'double', 100.0, 'AU'))
paramList.append(_ModelParam('sig0', 'sig0', 'Surface density at rout', 'double', 0.01, 'g/cm^2'))
resultList = [_resultsDict['density'],_resultsDict['tdust'],_resultsDict['velocity']]
_modelsDict[idStr] = _Model(idStr
  , 'Chiang & Goldreich 1997'
  , 'Passive irradiated disk'
  , 'Chiang & Goldreich 1997, ApJ, 490, 368'
  , paramList, {}, resultList)

idStr = 'DDN01'
paramList = []
paramList.append(_ModelParam('Mstar', 'Mstar', 'Mass of the central star', 'double', 1.0, 'Msun'))
paramList.append(_ModelParam('Rstar', 'Rstar', 'Radius of the central star', 'double', 1.0, 'Rsun'))
paramList.append(_ModelParam('Tstar', 'Mstar', 'Effective temperature of the central star', 'double', 4000.0, 'K'))
paramList.append(_ModelParam('bgdens', 'bgdens', 'Background number density (floor value for the disk density)', 'double', 0.0001, '1/cm^3'))
paramList.append(_ModelParam('dustopac_fname', 'dustopac_fname', 'File name (with full path) containing dust opacity table in two column (wavelength [um], kappa_abs [cm^2/g]) format', 'string', 'jena_thin_e6.tab', ''))
paramList.append(_ModelParam('mdisk', 'mdisk', 'Total disk mass (gas + dust ; gas-to-dust ratio of 100 is used)', 'double', 0.01, 'Msun'))
paramList.append(_ModelParam('plsig1', 'plsig1', 'Power exponent of the radial surface density distribution', 'double', -1.0, ''))
paramList.append(_ModelParam('rin', 'rin', 'Inner Radius of the disk', 'double', 1.0, 'AU'))
paramList.append(_ModelParam('rout', 'rout', 'Outer radius of the disk', 'double', 100.0, 'AU'))
resultList = [_resultsDict['density'],_resultsDict['tdust'],_resultsDict['velocity']]
_modelsDict[idStr] = _Model(idStr
  , 'DullemondDominik01'
  , 'Passive irradiated circumstellar disk with an inner hole'
  , 'Dullemond & Dominik 2001, ApJ, 560, 957'
  , paramList, {}, resultList)

idStr = 'LiShu96'
paramList = []
paramList.append(_ModelParam('cs', 'cs', 'Sound speed', 'double', 1000.0, 'm/s'))
#paramList.append(_ModelParam('h0', 'H0', 'Mass flux ratio parameter', 'enum', (0.125, 0.5, 0.75, 1.0, 1.25), 3, ''))
paramList.append(_ModelParam('h0', 'H0', 'Mass flux ratio parameter', 'int', 3, ''))
resultList = [_resultsDict['bmag'],_resultsDict['density'],_resultsDict['temperature'],_resultsDict['velocity']]
_modelsDict[idStr] = _Model(idStr
  , 'Li & Shu 96'
  , 'Axial-symetrical isothermal collapse with magnetic field (Li & Shu 1996)'
  , 'Li & Shu 1996, ApJ, 472, 211'
  , paramList, {}, resultList)

idStr = 'Mamon88'
paramList = []
paramList.append(_ModelParam('ab0', 'ab0', 'Abundance at inner edge', 'double', 0.0004, ''))
paramList.append(_ModelParam('mdot', 'Mdot', 'Mass loss rate (10^-8 .. 10^-4 Msun/yr)', 'double', 1e-05, 'Msun/year'))
paramList.append(_ModelParam('rin', 'Rin', 'Inner radius', 'double', 2.0, 'AU'))
paramList.append(_ModelParam('tin', 'Tin', 'Temperature at inner edge', 'double', 1100.0, 'K'))
paramList.append(_ModelParam('ve', 've', 'Expansion velocity (7.5 .. 30 km/s)', 'double', 15.0, 'km/s'))
resultList = [_resultsDict['abundance'],_resultsDict['density'],_resultsDict['temperature'],_resultsDict['velocity']]
_modelsDict[idStr] = _Model(idStr
  , 'Mamon et al. 1988'
  , 'Expanding circumstellar envelope'
  , 'Mamon et al. 1988, ApJ 328, 797'
  , paramList, {}, resultList)

idStr = 'Mendoza09'
paramList = []
paramList.append(_ModelParam('mdot', 'mdot', 'Accretion rate', 'double', 1e-06, 'Msun/year'))
paramList.append(_ModelParam('mstar', 'mstar', 'Mass of the central protostar', 'double', 0.5, 'Msun'))
paramList.append(_ModelParam('mu', 'mu', 'Ratio of the centrifugal radius to the cloud radius', 'double', 0.0, ''))
paramList.append(_ModelParam('nu', 'nu', 'Ratio of the radial velocity to the Kepler velocity at the cloud outer radius', 'double', 0.0, ''))
paramList.append(_ModelParam('rc', 'rc', 'Centrifugal radius', 'double', 200.0, 'AU'))
resultList = [_resultsDict['density'],_resultsDict['velocity']]
_modelsDict[idStr] = _Model(idStr
  , 'Mendoza09'
  , 'Collapse of a finite rotating cloud core  (Mendoza et al. 2009)'
  , 'Mendoza, Tejeda & Nagel, 2009, MNRAS, 393, 579'
  , paramList, {}, resultList)

idStr = 'Shu77'
paramList = []
paramList.append(_ModelParam('T', 'T', 'Temperature of the cloud', 'double', 30.0, 'K'))
paramList.append(_ModelParam('time', 'time', 'Time at which the density and velocity field should be calculated', 'double', 100000.0, 'year'))
resultList = [_resultsDict['density'],_resultsDict['temperature'],_resultsDict['tdust'],_resultsDict['velocity']]
_modelsDict[idStr] = _Model(idStr
  , 'Shu 1977'
  , 'Expansion-wave collapse solution of hydrostatic singular spheres Shu 1977'
  , 'Shu 1977, ApJ, 214, 488'
  , paramList, {}, resultList)

idStr = 'Ulrich76'
paramList = []
paramList.append(_ModelParam('mdot', 'mdot', 'Accretion rate', 'double', 1e-06, 'Msun/year'))
paramList.append(_ModelParam('mstar', 'mstar', 'Mass of the central protostar', 'double', 0.1, 'Msun'))
paramList.append(_ModelParam('rc', 'rc', 'Centrifugal radius', 'double', 60.0, 'AU'))
resultList = [_resultsDict['density'],_resultsDict['velocity']]
_modelsDict[idStr] = _Model(idStr
  , 'Ulrich76'
  , 'Cloud collapse model with rotation (Ulrich 1976)'
  , 'Ulrich 1976, ApJ, 210, 377'
  , paramList, {}, resultList)

#_modelsDict['user'] = None

