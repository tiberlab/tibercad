#!/bin/bash

ARCH=`uname -m`

BASEDIR=/usr/pack/tibercad_dev-3.3-ma/latest
BUILDDIR=${BASEDIR}/build
SDKDIR=${BASEDIR}/SDK
PETSC_DIR=${BUILDDIR}/petsc-3.17.4
SLEPC_DIR=${BUILDDIR}/slepc-3.17.2
LIBMESHDIR=${BUILDDIR}/libmesh-1.7.1-${ARCH}-linux

export CXX=${SDKDIR}/${ARCH}-linux/bin/mpicxx
export CC=${SDKDIR}/${ARCH}-linux/bin/mpicc
export FC=${SDKDIR}/${ARCH}-linux/bin/mpifort
export MPIEXEC=${SDKDIR}/${ARCH}-linux/bin/mpiexec
#export FCFLAGS="${FCFLAGS} -fexceptions -gcc-name=${CC} -gxx-name=${CXX} -nofor-main"
#export LDFLAGS="${LDFLAGS} `/usr/pack/tibercad_dev-2.2-ma/${ARCH}-linux/bin/mpicxx -showme:link`"
#export LDFLAGS="${LDFLAGS} `${PETSC_DIR}/externalpackages/mpich2-1.0.8/bin/mpicxx -showme:link`"


FORTRANDIR=/usr/pack/intel_oneapi-2021-ma/compiler/latest/linux/compiler
export LDFLAGS="-Wl,-rpath,${FORTRANDIR}/lib/intel64"
export FCFLAGS="-fexceptions -g -gcc-name=$(which ${CC}) -gxx-name=$(which ${CXX}) -nofor-main -qopenmp"

BOOST="${SDKDIR}"


SVN=svn-1.9.5

#  --enable-hetero \
#  --enable-dftb \
#  --with-mkl=/usr/pack/intel_mkl-10.2-gp \
#  --enable-pardiso \
#  --with-boost-prefix=/usr/pack/boost-1.33.1-ma \
#  --with-tao-prefix=/usr/pack/tibercad_dev-2.2-ma/tao-1.10-p1 \

./configure \
  --with-cuda=/usr/pack/cudatoolkit-11.4-ma \
  --with-petsc-prefix=${PETSC_DIR} \
  --with-petsc-arch=${ARCH}-linux \
  --with-slepc-prefix=${SLEPC_DIR} \
  --with-mpiexec=${MPIEXEC} \
  --with-mkl=/usr/pack/intel_oneapi-2021-ma/mkl/latest \
  --with-thread-library=intel \
  --with-subversion=${SVN} \
  --disable-license-check \
  --enable-uptight\
  --disable-pardiso \
  --with-boost-prefix=$BOOST \
  --with-boost-libdir=${BOOST}/${ARCH}-linux/lib \
  --with-libmesh-prefix=${LIBMESHDIR} \
  --with-libmesh-petsc-libdir=${SDKDIR}/${ARCH}-linux/lib
