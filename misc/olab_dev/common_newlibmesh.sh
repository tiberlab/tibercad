# build directory (relative)
BUILDDIR=`pwd`/../extern

# PETSc version
PETSCVERSION="3.0.0"
PETSCPATCHLEVEL="-p3"
SLEPCPATCHLEVEL="-p0"
#PETSCVERSION="2.3.3"
#PETSCPATCHLEVEL="-p13"
#SLEPCPATCHLEVEL=""

# libMesh version
#LIBMESHVERSION="-0.6.3-rc1"
LIBMESHVERSION="_svn"

# whether to compile debug mode or not (y/n)
DEBUG=y

# architecture
ARCH=linux-gnu


OTHEROPTS="--with-boost-prefix=/usr/pack/boost-1.33.1-ma"
#OTHEROPTS="--with-boost-prefix=/usr/pack/boost-1.33.1-ma --enable-uptight"
#OTHEROPTS="${OTHEROPTS} --enable-pardiso --with-mkl=/usr/pack/intel_mkl-8.0-ma"


v=-4.1.1
export CC=gcc$v
export CXX=g++$v #F77=gfortran-$v F90=$F77 FC=$F77
#export CC=gcc-$v CXX=g++-$v F77=g95-0.50 F90=g95-0.50
export F77=ifort-9.1
export F90=$F77
export FC=$F77

export MPICH_F90=$F90
export MPICH_F77=$F77
export MPIDIR=/usr/pack/mpich2-1.0.6-ma
