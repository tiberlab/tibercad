#!/bin/sh
 
export CXX=g++-4.6.0
export CC=gcc-4.6.0
export FC=ifort-11.1
export FCFLAGS="-fPIC -fexceptions -gcc-name=${CC} -gxx-name=${CXX} -nofor-main"

SDKDIR=/usr/pack/tibercad_sdk-2.2-ma
ARCH=x86_64-linux



#FCFLAGS="-fPIC ${FCFLAGS}"

FORTRANDIR=/usr/pack/intel_fc-11.1-gp
export LDFLAGS="-Wl,-rpath,${FORTRANDIR}/lib/intel64"

#  --with-tao-prefix=/usr/pack/tibercad_dev-2.1-ma/tao-1.10-p1 \

./configure \
  --with-petsc-prefix=${SDKDIR}/build/petsc-3.0.0-p12 \
  --with-petsc-arch=${ARCH}-complex \
  --with-slepc-prefix=${SDKDIR}/build/slepc-3.0.0-p7 \
  --with-mkl=/usr/pack/intel_mkl-10.2-gp \
  --with-cuda=/usr/pack/cudatoolkit-5.5.11-ma \
  --with-subversion=svn-1.6.5 \
  --enable-uptight \
  --enable-pardiso \
  --with-boost-prefix=$SDKDIR --with-boost-libdir=${SDKDIR}/${ARCH}/lib \
  --with-libmesh-petsc-libdir=${SDKDIR}/${ARCH}/lib \
  --with-libmesh-prefix=${SDKDIR}/build/libmesh-svn-${ARCH}
