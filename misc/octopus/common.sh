# build directory
BUILDDIR=`pwd`/../extern

# PETSc version
#PETSCVERSION="2.3.3"
#PATCHLEVEL="p13"
PETSCVERSION="3.0.0"
PETSCPATCHLEVEL="-p9"
SLEPCPATCHLEVEL="-p7"

# libMesh version
#LIBMESHVERSION=0.6.1
LIBMESHVERSION="svn"

# GMSH version
GMSHVERSION=2.3.0

# whether to compile debug mode or not (y/n)
DEBUG=y

# architecture
ARCH=linux-gnu

# build shared libs
shared_libs=yes

# compilers
v=-4.3.2
CC=gcc$v
CXX=g++$v
F77=ifort-11.1
F90=$F77
FC=$F77
export CC CXX F77 F90 FC

# compiler flags
export CXXFLAGS=-fPIC
export CFLAGS="-fPIC -fexceptions"
export FFLAGS="-fPIC -fexceptions -gcc-name=${CC} -gxx-name=${CXX} -nofor-main"
export FCFLAGS="-fPIC -fexceptions -gcc-name=${CC} -gxx-name=${CXX} -nofor-main"

#OTHEROPTS="--enable-uptight --disable-license-check"
OTHEROPTS="--disable-uptight --disable-license-check --disable-hetero"

#export MPICH_F90=$F90
#export MPICH_F77=$F77
#export MPIDIR=/usr/pack/mpich-1.2.7p1-ma
