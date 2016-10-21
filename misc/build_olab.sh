#!/bin/bash

#export TOPDIR=/usr/pack/tibercad_dev-2.5-ma
export BUILDDIR=/usr/pack/tibercad_dev-3.0-ma/build
export SDKDIR=/usr/pack/tibercad_dev-3.0-ma/SDK
export PATH=$(dirname $PWD)/SDK/x86_64-linux/bin:$PATH

export MKL=y
export MKLDIR=/usr/pack/intel_mkl-11.2-ma/mkl

export MPI=y
export MPIDIR=/usr/pack/mpich-3.1.1-ma
#/scratch/maufder/TiberCAD/sdk3/extern/petsc-3.0.0-p12/externalpackages/mpich2-1.0.8

export FORTRANDIR=/usr/pack/intel_ifort-13.0-ma

export CC=mpicc-3.1.1
export CXX=mpicxx-3.1.1
export FC=mpif90-3.1.1
#-11.1
export F77=${FC}
#-11.1

export CFLAGS="-g -O3 -fopenmp"
export CXXFLAGS="-g O3 -fopenmp"
export FCFLAGS="-fexceptions -g -O3 -gcc-name=$(which ${CC}) -gxx-name=$(which ${CXX}) -nofor-main -openmp"

#export DEBUG=n

./configure

#./build_mpi
#./build_petsc
#./build_slepc
#./build_libmesh
#./build_boost
./build_tibercad
