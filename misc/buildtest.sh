#!/bin/bash

#export BUILDDIR=/usr/pack/tibercad_sdk-2.2-ma/build
#export SDKDIR=/usr/pack/tibercad_sdk-2.2-ma
#export PATH=$(dirname $PWD)/SDK/x86_64-linux/bin:$PATH

export MKL=y

export MPI=n
export MPIDIR=
#/scratch/maufder/TiberCAD/sdk3/extern/petsc-3.0.0-p12/externalpackages/mpich2-1.0.8

export FORTRANDIR=/usr/pack/intel_fc-11.1-gp

export CC=gcc-4.6.0
export CXX=g++-4.6.0
export FC=ifort
#-11.1
export F77=ifort
#-11.1

export CFLAGS=-g
export CXXFLAGS=-g
export FCFLAGS="-fexceptions -gcc-name=$(which ${CC}) -gxx-name=$(which ${CXX}) -nofor-main"


./configure

./build_mpi
#./build_petsc
#./build_slepc
#./build_libmesh
#./build_boost
#./build_tibercad
