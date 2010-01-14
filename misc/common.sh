
# load variables
[ -e variables ] && . variables

# build directory
BUILDDIR=${BUILDDIR:-"`pwd`/../extern"}

# PETSc version
PETSCVERSION=${PETSCVERSION:-"3.0.0"}
PETSCPATCHLEVEL=${PETSCPATCHLEVEL:-"-p9"}
SLEPCPATCHLEVEL=${SLEPCPATCHLEVEL:-"-p7"}

# libMesh version
LIBMESHVERSION=${LIBMESHVERSION:-"svn"}

# GMSH version
GMSHVERSION=${GMSHVERSION:-"2.3.0"}

# whether to compile debug mode or not (y/n)
DEBUG=${DEBUG:-'n'}

# build shared libs (yes/no)
shared_libs=${shared_libs:-'y'}

# use MKL (yes/no)
MKL=${MKL:-'y'}
MKLDIR=${MKLDIR:-"/usr/pack/intel_mkl-10.2-gp"}

# MPI
MPI=${MPI:-'y'}
MPIDIR=${MPIDIR:-"/usr/pack/mpich-1.2.7p1-ma"}




# architecture
march=`uname -m`
ARCH=linux-gnu-${march}

# compilers
CC=${CC:-"gcc"}
CXX=${CXX:-"g++"}
F77=${F77:-"ifort-11.1"}
F90=${F90:-$F77}
FC=${FC:-$F77}
export CC CXX F77 F90 FC

# compiler flags
CXXFLAGS=
CFLAGS="-fexceptions"
FFLAGS="-fexceptions -gcc-name=${CC} -gxx-name=${CXX} -nofor-main"
FCFLAGS="-fexceptions -gcc-name=${CC} -gxx-name=${CXX} -nofor-main"
if [ $march == "x86_64" ]; then
  CXXFLAGS="-fPIC ${CXXFLAGS}"
  CFLAGS="-fPIC ${CFLAGS}"
  FFLAGS="-fPIC ${FFLAGS}"
  CFLAGS="-fPIC ${CFLAGS}"
  FCFLAGS="-fPIC ${FCFLAGS}"
fi
export CXXFLAGS CFLAGS FFLAGS CFLAGS FCFLAGS


if [ "x$FORTRANDIR" != "x" ]; then
  case $march in
    x86_64) FLIBDIR=${FORTRANDIR}/lib/intel64 ;;
    i?86)   FLIBDIR=${FORTRANDIR}/lib/ia32 ;;
    *)      echo "FORTRANDIR cannot be used on architecture $march"; exit 1 ;;
  esac
  LDFLAGS="${LDFLAGS} -Wl,-rpath,${FLIBDIR}"
fi


if [ $MPI == 'y' ]; then
  export MPIDIR
fi

# BLAS/LAPACK
if [ $MKL == 'y' ]; then
  if [ $march == "x86_64" ]; then
    LAPACKLIBS="--with-blas-lapack-lib=[-L${MKLDIR}/lib/em64t,-lmkl_intel_lp64,-lmkl_intel_thread,-lmkl_core,-liomp5,-lm,-lpthread]"
    LDFLAGS="${LDFLAGS} -L${MKLDIR}/lib/em64t -Wl,-rpath,${MKLDIR}/lib/em64t"
  else
    LAPACKLIBS="--with-blas-lapack-lib=[-L${MKLDIR}/lib/ia32,-lmkl_intel,-lmkl_intel_thread,-lmkl_core,-liomp5,-lm,-lpthread]"
    LDFLAGS="${LDFLAGS} -L${MKLDIR}/lib/ia32 -Wl,-rpath,${MKLDIR}/lib/ia32"
  fi
fi

export LDFLAGS

