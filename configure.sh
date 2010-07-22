#!/bin/sh
 
export CXX=g++-4.3.2
export CC=gcc-4.3.2
export FC=ifort-11.1
export FCFLAGS="-fexceptions -gcc-name=${CC} -gxx-name=${CXX} -nofor-main"

BOOST="/usr/pack/boost-1.41.0-ma"

ARCH=`uname -m`
if [ $ARCH == "i686" ]
then
  sepparch=i686-debian-linux4.0
else
  sepparch=amd64-debian-linux5.0
  FCFLAGS="-fPIC ${FCFLAGS}"
fi

CONFIGOPTS="--with-boost-prefix=$BOOST --with-boost-libdir=${BOOST}/${sepparch}/lib"


#  --enable-hetero \
#  --enable-dftb \
#  --with-mkl=/usr/pack/intel_mkl-10.2-gp \
#  --enable-pardiso \
#  --with-boost-prefix=/usr/pack/boost-1.33.1-ma \

./configure \
  --with-petsc-prefix=/usr/pack/tibercad_dev-2.0-ma/petsc-3.0.0-p10 \
  --with-petsc-arch=linux-gnu-${ARCH}-complex \
  --with-slepc-prefix=/usr/pack/tibercad_dev-2.0-ma/slepc-3.0.0-p7 \
  --with-tao-prefix=/usr/pack/tibercad_dev-2.0-ma/tao-1.10-p1 \
  --with-mkl=/usr/pack/intel_mkl-10.2-gp \
  --with-subversion=svn-1.6.5 \
  --enable-uptight\
  --enable-pardiso \
  ${CONFIGOPTS} \
  --with-libmesh-prefix=/usr/pack/tibercad_dev-2.0-ma/libmesh-svn-${ARCH}
