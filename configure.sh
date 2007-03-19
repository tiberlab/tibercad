#!/bin/sh

export CXX=g++-4.1.1
export CC=gcc-4.1.1

./configure --with-boost-prefix=/usr/pack/boost-1.33.1-ma \
  --with-libmesh-prefix=/usr/pack/libmesh-0.5.0_20070319cvs-ma \
  --with-petsc-prefix=/usr/pack/petsc-2.3.2-ma \
  --with-petsc-arch=linux-gnu-complex-opt \
  --with-slepc-prefix=/usr/pack/slepc-2.3.2-ma
