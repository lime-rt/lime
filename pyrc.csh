set pwd=`pwd`

if $?LIMEDIR then
  if ( "$LIMEDIR" == "" ) then
    set basedir=$pwd
  else
    set basedir=$LIMEDIR
  endif
else
  set basedir=$pwd
endif

set py_extras=${basedir}/python:${basedir}:${basedir}/example

if $?PYTHONPATH then
  if ( "$PYTHONPATH" == "" ) then
      setenv PYTHONPATH $py_extras
  else
      setenv PYTHONPATH $py_extras:$PYTHONPATH
  endif
else
    setenv PYTHONPATH $py_extras
endif

