
# top directory
TOPDIR="$(dirname `pwd`)"

# build directory for external libraries, can be overridden
BUILDDIR="${TOPDIR}/external"

# SDK directory, can be overridden
SDKDIR="${TOPDIR}/SDK"


# load site configuration variables
[ -e site_config.sh ] && . site_config.sh


#
# Defaults
#

# PETSc version
PETSCVERSION=${PETSCVERSION:-"3.23.4"}
PETSCPATCHLEVEL=${PETSCPATCHLEVEL:-""}
SLEPCVERSION=${SLEPCVERSION:-"3.23.1"}
SLEPCPATCHLEVEL=${SLEPCPATCHLEVEL:-""}

# libMesh version
LIBMESHVERSION=${LIBMESHVERSION:-"1.8.4"}

# GMSH version
GMSHVERSION=${GMSHVERSION:-"4.13.1"}

# BOOST version
BOOSTVERSION=${BOOSTVERSION:-"1_88_0"}

# whether to compile debug mode or not (y/n)
DEBUG=${DEBUG:-'n'}



# build shared libs (yes/no)
shared_libs=${shared_libs:-'n'}

# CUDA
#CUDADIR=${CUDADIR:-'/usr/pack/cudatoolkit-11.4-ma'}

# use MKL (yes/no)
if [ -n "$MKLDIR" ]; then
    MKL="y"
else
  MKL="n"
fi
#MKLDIR=${MKLDIR:-"/usr/pack/intel_oneapi-2021-ma/mkl/latest"}

# MPI
# NOTE: now MPI is necessary, although one might try again without
#       (for Windows, e.g.)
MPI=${MPI:-'y'}
if test "${MPIDIR}x" = "x"
then
  MPIDIR=$SDKDIR
  MPIEXEC=${MPIEXEC:-${SDKDIR}/x86_64-linux/bin/mpiexec}
else
  MPIEXEC=${MPIEXEC:-${MPIDIR}/bin/mpiexec}
fi

