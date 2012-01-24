#!/bin/bash

# returns the machine architecture as a string
# the given string has to be of type cpu-os

arch=unknown

os=$(uname -s)
cpu=$(uname -m)

if [ $# -ge 1 ]; then

  sys=( `echo $1 | sed 's/-/ /g'` )

  case ${#sys[@]} in
    1)
      cpu=
      os=${sys[0]} ;;

    *)
      cpu=${sys[0]}
      os=${sys[1]} ;;
  esac

fi


case $os in

  *Linux* | *linux*) arch=${cpu}-linux ;;

  *Darwin* | *darwin*) arch=${cpu}-darwin ;;

  *MINGW* | *mingw* ) arch=${cpu}-mingw32 ;;

  *CYGWIN* | *cygwin* ) arch=${cpu}-cygwin ;;

  *) arch= ;;

esac


echo $arch
