
# load variables
[ -e variables ] && . variables

# top directory
TOPDIR=${TOPDIR:-"`pwd`/.."}

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


# system
system=`uname -s`
case $system in
  CYGWIN_NT-* | MINGW*)
    SYSTEM=mingw
    ;;

  Linux)
    SYSTEM=linux-gnu
    ;;

  *)
    echo "$system is not supported"; exit 1 ;;
esac

# target system
case $target in
  *mingw* | *cygwin* )
    target_os=win32 ;;

  * )
    target_os= ;;
esac

# architecture
march=`uname -m`

BUILD=`${TOPDIR}/share/arch.sh`

HOST=$BUILD
if [ "x$target" != "x" ]
then
  HOST=$target
fi
echo "*************************************************"
echo "* Host system :  $HOST"
echo "* Build system:  $BUILD"
echo "*************************************************"


# compilers
pre=
if [ "x$target" != "x" ]; then
  pre=${target}-
fi
if [ "x$gcc_version" != "x" ]; then
  gcc_version=-${gcc_version}
fi
CC=${CC:-"${pre}gcc${gcc_version}"}
CXX=${CXX:-"${pre}g++${gcc_version}"}
F77=${F77:-"${pre}gfortran${gcc_version}"}
F90=${F90:-$F77}
FC=${FC:-$F77}
AR=${AR:-"${pre}ar"}
LD=${LD:-"${pre}ld"}
RANLIB=${RANLIB:-"${pre}ranlib"}
CPP=${CPP:-"${pre}cpp${gcc_version}"}
export CC CXX F77 F90 FC RANLIB LD AR CPP

# compiler flags
CFLAGS="${CFLAGS} -fexceptions"
FFLAGS="${FFLAGS} -fexceptions"
if [ $march == "x86_64" ] &&  [ ${target_os} != "win32" ] ; then
  CXXFLAGS="-fPIC ${CXXFLAGS}"
  CFLAGS="-fPIC ${CFLAGS}"
  FFLAGS="-fPIC ${FFLAGS}"
fi
FCFLAGS=${FFFLAGS}
export CXXFLAGS CFLAGS FFLAGS FCFLAGS


if [ "x$FORTRANDIR" != "x" ]; then
  case $march in
    x86_64) FLIBDIR=${FORTRANDIR}/lib/intel64 ;;
    i?86)   FLIBDIR=${FORTRANDIR}/lib/ia32 ;;
    *)      echo "FORTRANDIR cannot be used on architecture $march"; exit 1 ;;
  esac
  LDFLAGS="${LDFLAGS} -Wl,-rpath-link,${FLIBDIR}"
fi


if [ $MPI == 'y' ]; then
  export MPIDIR
fi

# BLAS/LAPACK
if [ $MKL == 'y' ]; then
  if [ $march == "x86_64" ]; then
    LAPACKLIBS="--with-blas-lapack-lib=[-L${MKLDIR}/lib/em64t,-lmkl_intel_lp64,-lmkl_intel_thread,-lmkl_core,-liomp5,-lm,-lpthread]"
    LDFLAGS="${LDFLAGS} -L${MKLDIR}/lib/em64t -Wl,-rpath-link,${MKLDIR}/lib/em64t"
  else
    LAPACKLIBS="--with-blas-lapack-lib=[-L${MKLDIR}/lib/32,-lmkl_intel,-lmkl_intel_thread,-lmkl_core,-liomp5,-lm,-lpthread]"
    LDFLAGS="${LDFLAGS} -L${MKLDIR}/lib/32 -Wl,-rpath,${MKLDIR}/lib/32"
  fi
fi

export LDFLAGS

