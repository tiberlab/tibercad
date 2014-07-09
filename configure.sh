#!/bin/sh

ARCH=`uname -m`
if [ $ARCH == "i686" ]
then
  sepparch=i686-debian-linux4.0
else
  sepparch=amd64-debian-linux5.0
  FCFLAGS="-fPIC"
fi 

BUILDDIR=/usr/pack/tibercad_dev-2.2-ma/extern
PETSC_DIR=${BUILDDIR}/petsc-3.0.0-p12
SLEPC_DIR=${BUILDDIR}/slepc-3.0.0-p7

#export CXX=g++-4.6.0
#export CC=gcc-4.6.0
#export FC=ifort-11.1
export CXX=/usr/pack/tibercad_dev-2.2-ma/${ARCH}-linux/bin/mpicxx
export CC=/usr/pack/tibercad_dev-2.2-ma/${ARCH}-linux/bin/mpicc
export FC=/usr/pack/tibercad_dev-2.2-ma/${ARCH}-linux/bin/mpif90
export FCFLAGS="${FCFLAGS} -fexceptions -gcc-name=${CC} -gxx-name=${CXX} -nofor-main"
#export LDFLAGS="${LDFLAGS} `/usr/pack/tibercad_dev-2.2-ma/${ARCH}-linux/bin/mpicxx -showme:link`"
#export LDFLAGS="${LDFLAGS} `${PETSC_DIR}/externalpackages/mpich2-1.0.8/bin/mpicxx -showme:link`"


FORTRANDIR=/usr/pack/intel_ifort-13.0-ma
export LDFLAGS="-Wl,-rpath,${FORTRANDIR}/lib/intel64"

BOOST="/usr/pack/tibercad_dev-2.2-ma/SDK"

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
  --with-mkl=/usr/pack/intel_mkl-10.2-gp \
  --with-subversion=svn-1.6.5 \
  --disable-license-check \
  --enable-uptight\
  --enable-pardiso \
  ${CONFIGOPTS} \
  --with-libmesh-prefix=/usr/pack/tibercad_dev-2.2-ma/extern/libmesh-svn-${ARCH}-linux \
  --with-libmesh-petsc-libdir=/usr/pack/tibercad_dev-2.2-ma/${ARCH}-linux/lib
