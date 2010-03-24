#!/bin/sh
 
export CXX=g++-4.3.2
export CC=gcc-4.3.2
export FC=ifort-11.1
export FCFLAGS="-fexceptions -gcc-name=${CC} -gxx-name=${CXX} -nofor-main"

ARCH=`uname -m`
if [ $ARCH == "i686" ]
then
  CONFIGOPTS="--with-boost-prefix=/usr/pack/boost-1.33.1-ma"
else
  BOOST="/usr/pack/boost-1.41.0-ma"
  FCFLAGS="-fPIC ${FCFLAGS}"
  CONFIGOPTS="--with-boost-prefix=$BOOST --with-boost-libdir=${BOOST}/amd64-debian-linux5.0/lib"
fi


#  --enable-hetero \
#  --enable-dftb \
#  --with-mkl=/usr/pack/intel_mkl-10.2-gp \
#  --enable-pardiso \
#  --with-boost-prefix=/usr/pack/boost-1.33.1-ma \

./configure \
  --with-petsc-prefix=/usr/pack/tibercad_dev-2.0-ma/petsc-3.0.0-p10 \
  --with-petsc-arch=linux-gnu-${ARCH}-complex \
  --with-slepc-prefix=/usr/pack/tibercad_dev-2.0-ma/slepc-3.0.0-p7 \
  --with-mkl=/usr/pack/intel_mkl-10.2-gp \
  --with-subversion=svn-1.6.5 \
  --disable-modules \
  --enable-uptight\
  --enable-pardiso \
  ${CONFIGOPTS} \
  --with-libmesh-prefix=/usr/pack/tibercad_dev-2.0-ma/libmesh-svn-${ARCH}
