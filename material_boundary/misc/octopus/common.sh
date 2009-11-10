# build directory (relative)
BUILDDIR=`pwd`/../extern

# PETSc version
PETSCVERSION="2.3.3"
PATCHLEVEL="p13"

# libMesh version
LIBMESHVERSION=0.6.1

# GMSH version
GMSHVERSION=2.3.0

# whether to compile debug mode or not (y/n)
DEBUG=y

# architecture
ARCH=linux-gnu

v=4.3.2
#export CC=gcc-$v CXX=g++-$v F77=gfortran-$v F90=$F77 FC=$F77
#CC=gcc-$v CXX=g++-$v
CC=gcc CXX=g++
F77=ifort-9.1 F90=$F77 FC=$F77
export CC CXX F77 F90 FC

OTHEROPTS="--enable-uptight --disable-license-check"

#export MPICH_F90=$F90
#export MPICH_F77=$F77
#export MPIDIR=/usr/pack/mpich-1.2.7p1-ma
