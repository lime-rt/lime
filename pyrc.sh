pwd=`pwd`

if [ -z "$LIMEDIR" ]; then
  basedir=$pwd
else
  basedir=$LIMEDIR
fi

py_extras=${basedir}/python:${basedir}:${basedir}/example
if [ -z "$PYTHONPATH" ]; then
  export PYTHONPATH=$py_extras
else
  export PYTHONPATH=${PYTHONPATH}:$py_extras
fi

