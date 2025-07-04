#!/bin/bash

ARCH=`uname -m`

BASEDIR=/usr/pack/tibercad_dev-3.5-ma/
BUILDDIR=${BASEDIR}/build
SDKDIR=${BASEDIR}/SDK
PETSC_DIR=${BUILDDIR}/petsc-3.23.3
SLEPC_DIR=${BUILDDIR}/slepc-3.23.1
LIBMESHDIR=${BUILDDIR}/libmesh-1.8.1-${ARCH}-linux

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
