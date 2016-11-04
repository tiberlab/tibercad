#!/bin/bash

#export BUILDDIR=/usr/pack/tibercad_sdk-2.2-ma/build
#export SDKDIR=/usr/pack/tibercad_sdk-2.2-ma
#export PATH=$(dirname $PWD)/SDK/x86_64-linux/bin:$PATH

export MKL=n
export LAPACKLIBS="-llapack -lblas"

export MPI=y
export MPIDIR=/usr/lib/mpich2
#/scratch/maufder/TiberCAD/sdk3/extern/petsc-3.0.0-p12/externalpackages/mpich2-1.0.8

#export FORTRANDIR=/usr/pack/intel_fc-11.1-gp

export CC=mpicc.mpich2
export CXX=mpicxx.mpich2
export FC=mpif90.mpich2
#-11.1
export F77=${FC}
#-11.1

export CFLAGS=-g
export CXXFLAGS=-g
export FCFLAGS="-fexceptions -g"
#-gcc-name=$(which ${CC}) -gxx-name=$(which ${CXX}) -nofor-main"


./configure

#./build_mpi
#./build_petsc
#./build_slepc
./build_libmesh
#./build_boost
./build_tibercad
