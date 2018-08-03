# Environment setup script for BASH. You should navigate to the directory this script is stored in and do
#
#	. ./pylimerc.sh

pwd=`pwd`

if [ -z "$PYLIMEDIR" ]; then
  basedir=$pwd
else
  basedir=$PYLIMEDIR
fi

py_extras=${basedir}/python:${basedir}:${basedir}/example:${basedir}/casa
if [ -z "$PYTHONPATH" ]; then
  export PYTHONPATH=$py_extras
else
  export PYTHONPATH=${PYTHONPATH}:$py_extras
fi

if [ -z "$PATH" ]; then
  export PATH=${basedir}
else
  export PATH=${PATH}:${basedir}
fi

