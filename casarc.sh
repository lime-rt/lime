pwd=`pwd`

if [ -z "$CASALIMEDIR" ]; then
  basedir=$pwd
else
  basedir=$CASALIMEDIR
fi

py_extras=${basedir}/python:${basedir}:${basedir}/example
if [ -z "$PYTHONPATH" ]; then
  export PYTHONPATH=$py_extras
else
  export PYTHONPATH=${PYTHONPATH}:$py_extras
fi

ld_extras=${basedir}
if [ -z "$LD_LIBRARY_PATH" ]; then
  export LD_LIBRARY_PATH=$ld_extras
else
  export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$ld_extras
fi

