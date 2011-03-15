#!/bin/sh

# returns the machine architecture as a string

arch=unknown

os=`uname`

case $os in

  Linux) arch=`uname -m`-linux ;;

  *) arch= ;;

esac


echo $arch
