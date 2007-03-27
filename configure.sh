#!/bin/sh

export CXX=g++-4.0.0
export CC=gcc-4.0.0

./configure --with-boost-prefix=/usr/pack/boost-1.33.1-ma \
  --with-libmesh-prefix=/usr/pack/libmesh-0.5.0_20051010cvs-ma \
  --with-petsc-prefix=/usr/pack/petsc-2.3.0-ma \
  --with-petsc-arch=linux-gnu-complex-opt \
  --with-slepc-prefix=/usr/pack/slepc-2.3.0-ma
