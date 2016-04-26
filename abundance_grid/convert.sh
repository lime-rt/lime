#!/bin/sh
#
# This small shell-script starts convert.x for all mol*.dat files
# in the directory and converts them into binary files
#
# Simon Bruderer 26.10.2008
#

for i in mol*.dat
do
  name=`basename $i .dat`
  echo convert.x $name.dat $name.bin GRIDMOL_R0.01
  ./convert.x $name.dat $name.bin GRIDMOL_R0.01
done
