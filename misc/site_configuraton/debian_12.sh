PETSCVERSION="3.23.4"
SLEPCVERSION="3.23.1"
LIBMESHVERSION="1.8.4"
GMSHVERSION="4.13.1"
BOOSTVERSION="1_88_0"

MKL=n
LAPACKLIBS="-llapack -lblas"

MPIDIR=/usr/lib/x86_64-linux-gnu/openmpi

CC=mpicc.openmpi
CXX=mpicxx.openmpi
FC=mpif90.openmpi
F77=${FC}
MPIEXEC=/usr/bin/mpiexec.openmpi


