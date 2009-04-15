#!/bin/sh
 
export CXX=g++-4.1.1
export CC=gcc-4.1.1
export FC=ifort-9.0

./configure --with-boost-prefix=/usr/pack/boost-1.33.1-ma \
  --with-petsc-prefix=/usr/pack/petsc-2.3.3-ma \
  --with-petsc-arch=linux-gnu-complex \
  --with-slepc-prefix=/usr/pack/slepc-2.3.3-ma \
  --with-mkl=/usr/pack/intel_mkl-8.0-ma \
  --with-subversion=svn-1.5.6 \
  --disable-modules \
  --enable-pardiso \
  --enable-dftb \
  --enable-uptight\
  --enable-hetero \
  --with-libmesh-prefix=/usr/pack/libmesh-0.6.1-ma 
