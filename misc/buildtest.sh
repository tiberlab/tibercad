#!/bin/bash

export PETSCVERSION="3.16.1"
export SLEPCVERSION="3.16.1"
export LIBMESHVERSION="1.6.2"
export GMSHVERSION="4.8.3"
export BOOSTVERSION="1_77_0"

#export BUILDDIR=/usr/pack/tibercad_sdk-2.2-ma/build
#export SDKDIR=/usr/pack/tibercad_sdk-2.2-ma
#export PATH=$(dirname $PWD)/SDK/x86_64-linux/bin:$PATH

export MKL=n
export LAPACKLIBS="-llapack -lblas"

export MPI=y
export MPIDIR=/usr/lib/x86_64-linux-gnu/mpich
#/scratch/maufder/TiberCAD/sdk3/extern/petsc-3.0.0-p12/externalpackages/mpich2-1.0.8

#export FORTRANDIR=/usr/pack/intel_fc-11.1-gp

export CC=mpicc.mpich
export CXX=mpicxx.mpich
export FC=mpif90.mpich
export F77=${FC}
export MPIEXEC=/usr/bin/mpiexec.hydra

export CFLAGS=-g
export CXXFLAGS=-g
export FCFLAGS="-fexceptions -g"
#-gcc-name=$(which ${CC}) -gxx-name=$(which ${CXX}) -nofor-main"


./configure

#./build_mpi
./build_petsc
./build_slepc
./build_libmesh
./build_boost
./build_tibercad
