#!/bin/sh

# returns the machine architecture as a string

arch=unknown

os=`uname -o`

case $os in

  *Linux*) arch=`uname -m`-linux ;;

  *Cygwin*) arch=`uname -m` ;;

  *) arch= ;;

esac


echo $arch
