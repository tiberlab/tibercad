#!/bin/sh

ARCH=`uname -m`

BUILDDIR=/usr/pack/tibercad_dev-3.0-ma/build
PETSC_DIR=${BUILDDIR}/petsc-3.6.2
SLEPC_DIR=${BUILDDIR}/slepc-3.6.2

export CXX=mpicxx-3.1.1
export CC=mpicc-3.1.1
export FC=mpif90-3.1.1 -fc=ifort-11.1
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

BOOST="/usr/pack/tibercad_dev-3.0-ma/SDK"

#CONFIGOPTS="--with-boost-prefix=$BOOST --with-boost-libdir=${BOOST}/${sepparch}/lib"
CONFIGOPTS="--with-boost-prefix=$BOOST --with-boost-libdir=${BOOST}/${ARCH}-linux/lib"


#  --enable-hetero \
#  --enable-dftb \
#  --with-mkl=/usr/pack/intel_mkl-10.2-gp \
#  --enable-pardiso \
#  --with-boost-prefix=/usr/pack/boost-1.33.1-ma \
#  --with-tao-prefix=/usr/pack/tibercad_dev-2.2-ma/tao-1.10-p1 \

./configure \
  --with-petsc-prefix=${PETSC_DIR} \
  --with-petsc-arch=${ARCH}-linux-complex \
  --with-slepc-prefix=${SLEPC_DIR} \
  --with-mkl=/usr/pack/intel_mkl-11.2-ma/mkl \
  --with-subversion=svn-1.6.5 \
  --disable-license-check \
  --enable-uptight\
  --disable-pardiso \
  ${CONFIGOPTS} \
  --with-libmesh-prefix=/usr/pack/tibercad_dev-3.0-ma/build/libmesh-0.9.5-${ARCH}-linux \
  --with-libmesh-petsc-libdir=/usr/pack/tibercad_dev-3.0-ma/SDK/${ARCH}-linux/lib
