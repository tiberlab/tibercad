#!/bin/bash

export PETSCVERSION="3.15.2"
export SLEPCVERSION="3.15.2"
export LIBMESHVERSION="1.6.2"
export GMSHVERSION="4.8.3"
export BOOSTVERSION="1_77_0"

#export BUILDDIR=/usr/pack/tibercad_dev-3.2-ma/build
#export SDKDIR=/usr/pack/tibercad_dev-3.2-ma/SDK
#export SDKDIR=/scratch/maufder/TiberCAD/3.2/SDK
ARCH=x86_64-linux

export MKL=y
export MKLDIR=/usr/pack/intel_oneapi-2021-ma/mkl/latest

export MPI=y
export MPIDIR="$(dirname `pwd`)/SDK"
#export MPIDIR=/usr/pack/tibercad_dev-3.2-ma/SDK

export FORTRANDIR=/usr/pack/intel_oneapi-2021-ma/compiler/latest/linux/compiler

export MPIEXEC=${MPIDIR}/${ARCH}/bin/mpiexec

if test -e ${MPIEXEC}
then 
  export CC=${MPIDIR}/${ARCH}/bin/mpicc
  export CXX=${MPIDIR}/${ARCH}/bin/mpicxx
  export FC=${MPIDIR}/${ARCH}/bin/mpifort
  export F77=${MPIDIR}/${ARCH}/bin/mpifort
else
  export CC=gcc-11.1
  export CXX=g++-11.1
  export FC=ifort-2021
  export F77=ifort-2021
  #export FC=gfortran-11.1
  #export F77=gfortran-11.1

  MAKE_MPICH=1
fi

export CFLAGS=
export CXXFLAGS=
export FCFLAGS="-fexceptions -gcc-name=$(which ${CC}) -gxx-name=$(which ${CXX}) -nofor-main"
#export FCFLAGS="-fexceptions"

#export FFLAGS="-fallow-argument-mismatch"

#export LINKONLY=1

./configure

if test "${MAKE_MPICH}x" != "x"
then 
  echo "MPI not yet present, I will now compile MPICH"
  echo "Re-run configure_package.sh after this has finished"
  ./build_mpich
else
#  ./build_petsc
#  ./build_slepc
#  ./build_libmesh
#  ./build_boost
  ./build_tibercad
fi
