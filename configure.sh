#!/bin/bash

ARCH=`uname -m`

BASEDIR=/usr/pack/tibercad_dev-3.1-ma
BUILDDIR=${BASEDIR}/build
SDKDIR=${BASEDIR}/SDK
#PETSC_DIR=${BUILDDIR}/petsc-3.7.3
#SLEPC_DIR=${BUILDDIR}/slepc-3.7.3
#LIBMESHDIR=${BUILDDIR}/libmesh-1.0.0-${ARCH}-linux
PETSC_DIR=${BUILDDIR}/petsc-3.6.2
SLEPC_DIR=${BUILDDIR}/slepc-3.6.2
LIBMESHDIR=${BUILDDIR}/libmesh-1.0.0-${ARCH}-linux

export CXX=mpicxx-3.1.1
export CC=mpicc-3.1.1
export FC=mpifort-3.1.1
export F77=${FC}
#export CXX=/usr/pack/tibercad_dev-2.2-ma/${ARCH}-linux/bin/mpicxx
#export CC=/usr/pack/tibercad_dev-2.2-ma/${ARCH}-linux/bin/mpicc
#export FC=/usr/pack/tibercad_dev-2.2-ma/${ARCH}-linux/bin/mpif90
#export FCFLAGS="${FCFLAGS} -fexceptions -gcc-name=${CC} -gxx-name=${CXX} -nofor-main"
#export LDFLAGS="${LDFLAGS} `/usr/pack/tibercad_dev-2.2-ma/${ARCH}-linux/bin/mpicxx -showme:link`"
#export LDFLAGS="${LDFLAGS} `${PETSC_DIR}/externalpackages/mpich2-1.0.8/bin/mpicxx -showme:link`"


FORTRANDIR=/usr/pack/intel_ifort-13.0-ma
export LDFLAGS="-Wl,-rpath,${FORTRANDIR}/lib/intel64"
export FCFLAGS="-fexceptions -g -gcc-name=$(which ${CC}) -gxx-name=$(which ${CXX}) -nofor-main -openmp"

BOOST="${SDKDIR}"

#CONFIGOPTS="--with-boost-prefix=$BOOST --with-boost-libdir=${BOOST}/${sepparch}/lib"
CONFIGOPTS="--with-boost-prefix=$BOOST --with-boost-libdir=${BOOST}/${ARCH}-linux/lib"

SVN=svn-1.6.5
if test "$(cat /etc/debian_version | awk -F . '{print $1}')" = "8"
then
  SVN=svn-1.9.5
fi

#  --enable-hetero \
#  --enable-dftb \
#  --with-mkl=/usr/pack/intel_mkl-10.2-gp \
#  --enable-pardiso \
#  --with-boost-prefix=/usr/pack/boost-1.33.1-ma \
#  --with-tao-prefix=/usr/pack/tibercad_dev-2.2-ma/tao-1.10-p1 \

./configure \
  --with-petsc-prefix=${PETSC_DIR} \
  --with-petsc-arch=${ARCH}-linux \
  --with-slepc-prefix=${SLEPC_DIR} \
  --with-mpiexec=mpiexec-3.1.1 \
  --with-mkl=/usr/pack/intel_mkl-11.2-ma/mkl \
  --with-subversion=${SVN} \
  --disable-license-check \
  --enable-uptight\
  --disable-pardiso \
  ${CONFIGOPTS} \
  --with-libmesh-prefix=${LIBMESHDIR} \
  --with-libmesh-petsc-libdir=${SDKDIR}/${ARCH}-linux/lib
