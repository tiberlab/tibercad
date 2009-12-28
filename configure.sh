#!/bin/sh
 
export CXX=g++-4.3.2
export CC=gcc-4.3.2
export FC=ifort-11.1
export FCFLAGS="-fPIC -fexceptions -gcc-name=${CC} -gxx-name=${CXX} -nofor-main"

ARCH=`uname -m`
if [ $ARCH == "i686" ]
then
  ARCH=""
fi


#  --enable-hetero \
#  --enable-dftb \
#  --with-mkl=/usr/pack/intel_mkl-10.2-gp \
#  --enable-pardiso \
#  --with-boost-prefix=/usr/pack/boost-1.33.1-ma \

./configure \
  --with-petsc-prefix=/usr/pack/tibercad_dev-1.3-ma/${ARCH}/petsc-3.0.0-p9 \
  --with-petsc-arch=linux-gnu-complex \
  --with-slepc-prefix=/usr/pack/tibercad_dev-1.3-ma/${ARCH}/slepc-3.0.0-p7 \
  --with-mkl=/usr/pack/intel_mkl-10.2-gp \
  --with-subversion=svn-1.6.5 \
  --disable-modules \
  --enable-uptight\
  --enable-pardiso \
  --with-boost-prefix=/usr/pack/boost-1.33.1-ma \
  --with-libmesh-prefix=/usr/pack/tibercad_dev-1.3-ma/${ARCH}/libmesh-svn
