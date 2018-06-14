#!/bin/bash

#export BUILDDIR=/usr/pack/tibercad_sdk-2.2-ma/build
export SDKDIR=/scratch/maufder/TiberCAD/3.0.0/SDK
ARCH=x86_64-linux

export MKL=y
export MKLDIR=/usr/pack/intel_mkl-11.2-ma/mkl

export MPI=y
export MPIDIR=${SDKDIR}

#export FORTRANDIR=/usr/pack/intel_fc-11.1-gp

# for mpich
#export CC=${SDKDIR}/${ARCH}/bin/gcc
#export CXX=${SDKDIR}/${ARCH}/bin/g++
#export FC=${SDKDIR}/${ARCH}/bin/gfortran
#export F77=${SDKDIR}/${ARCH}/bin/gfortran

export CC=${SDKDIR}/bin/mpicc
export CXX=${SDKDIR}/bin/mpicxx
export FC=${SDKDIR}/bin/mpifort
export F77=${SDKDIR}/bin/mpifort

export CFLAGS=
export CXXFLAGS=
#export FCFLAGS="-fexceptions -gcc-name=$(which ${CC}) -gxx-name=$(which ${CXX}) -nofor-main"
export FCFLAGS="-fexceptions"


./configure

#./build_mpich
./build_petsc
./build_slepc
./build_libmesh
./build_boost
./build_tibercad
