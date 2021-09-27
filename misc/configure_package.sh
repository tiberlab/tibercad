#!/bin/bash

export BUILDDIR=/usr/pack/tibercad_dev-3.2-ma/build
export SDKDIR=/usr/pack/tibercad_dev-3.2-ma/SDK
#export SDKDIR=/work/maufder/TiberCAD/3.2/SDK
ARCH=x86_64-linux

export MKL=y
export MKLDIR=/usr/pack/intel_oneapi-2021-ma/mkl/latest

export MPI=y
export MPIDIR=${SDKDIR}

export FORTRANDIR=/usr/pack/intel_oneapi-2021-ma/compiler/latest/linux/compiler

export MPIEXEC=${SDKDIR}/${ARCH}/bin/mpiexec


if test -e ${MPIEXEC}
then 
  export CC=${SDKDIR}/${ARCH}/bin/mpicc
  export CXX=${SDKDIR}/${ARCH}/bin/mpicxx
  export FC=${SDKDIR}/${ARCH}/bin/mpifort
  export F77=${SDKDIR}/${ARCH}/bin/mpifort
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
