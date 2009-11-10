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
DEBUG=n

# architecture
ARCH=linux-gnu

OTHEROPTS="--enable-license-check"

v=4.1.1
export CC=gcc-$v CXX=g++-$v F77=gfortran-$v F90=$F77 FC=$F77
#export CC=gcc-$v CXX=g++-$v F77=g95-0.50 F90=g95-0.50
#export F77=g95-0.91 F90=g95-0.91 FC=g95-0.91

#export MPICH_F90=$F90
#export MPICH_F77=$F77
#export MPIDIR=/usr/pack/mpich-1.2.7p1-ma
