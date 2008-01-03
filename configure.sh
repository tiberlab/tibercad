#!/bin/sh
 
export CXX=g++-4.1.1
export CC=gcc-4.1.1

./configure --with-boost-prefix=/usr/pack/boost-1.33.1-ma \
  --with-petsc-prefix=/usr/pack/petsc-2.3.2-ma \
  --with-petsc-arch=linux-gnu-complex \
  --with-slepc-prefix=/usr/pack/slepc-2.3.2-ma \
  --disable-modules \
  --with-libmesh-prefix=/usr/pack/libmesh-0.6.0_20070607cvs-ma 
#  --with-libmesh-prefix=/home/nano0/maufder/Projects/TIBERCAD/libmesh-0.6.1 
