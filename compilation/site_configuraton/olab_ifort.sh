# site configuration for University of Rome Tor Vergata, OLAB

# these are the defaults
PETSCVERSION="3.23.4"
SLEPCVERSION="3.23.1"
LIBMESHVERSION="1.8.1"
GMSHVERSION="4.10.5"
BOOSTVERSION="1_88_0"

# for system-wide installation
BUILDDIR=/usr/pack/tibercad_dev-3.5-ma/build
SDKDIR=/usr/pack/tibercad_dev-3.5-ma/SDK

MPIEXEC=${SDKDIR}/x86_64-linux/bin/mpiexec
MPIDIR=${SDKDIR}

# Intel MKL, if available
MKLDIR=/usr/pack/intel_oneapi-2021-ma/mkl/latest

# Cuda, if available
CUDADIR=/usr/pack/cudatoolkit-11.4-ma


FORTRANDIR=/usr/pack/intel_oneapi-2021-ma/compiler/latest/linux/compiler
FCLAGS="-fallow-argument-mismatch"

CC=gcc-11.1
CXX=g++-11.1
FC=ifort-2021
F77=ifort-2021
