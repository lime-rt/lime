# Environment setup script for CSHELL. You should navigate to the directory this script is stored in and do
#
#	source ./pylimerc.csh

set mypwd = `pwd`

if $?PYLIMEDIR then
  if ( "$PYLIMEDIR" == "" ) then
    set basedir = $mypwd
  else
    set basedir = $PYLIMEDIR
  endif
else
  set basedir = $mypwd
endif

set py_extras = ${basedir}/python:${basedir}:${basedir}/example:${basedir}/casa
if $?PYTHONPATH then
  if ( "$PYTHONPATH" == "" ) then
    setenv PYTHONPATH $py_extras
  else
    setenv PYTHONPATH ${PYTHONPATH}:$py_extras
  endif
else
  setenv PYTHONPATH $py_extras
endif

if $?PATH then
  if ( "$PATH" == "" ) then
    setenv PATH $basedir
  else
    setenv PATH ${PATH}:$basedir
  endif
else
  setenv PATH $basedir
endif

