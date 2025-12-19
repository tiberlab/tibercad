#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

dnl create verison number from git
m4_define([TCVERSION],
  [m4_esyscmd_s([ \
    (test -f VERSION && head -n 1 VERSION) || \
    (git describe --tags --dirty --always 2>/dev/null | sed 's:releases/::' | tr -d '[:space:]') || \
    echo -n "0.0.0" \
  ])]dnl
)

AC_INIT([tiberCAD],[TCVERSION()],[support@tibercad.org])
AC_PREREQ([2.71])

AC_CONFIG_MACRO_DIR([m4])

#AC_CONFIG_SRCDIR([include/tibercad/base/TiberCad.h])
AC_CONFIG_HEADERS([config.h])

AC_REVISION([$Revision: 2.0$])

# Make compilation tests using C++
AC_LANG([C++])

# The host
AC_CANONICAL_HOST


#
# Check for programs.
#

# compilers
AC_PROG_CXX
AC_PROG_CC
AC_PROG_FC
AC_PROG_RANLIB
AC_FC_LIBRARY_LDFLAGS
AC_SUBST([FCLIBS])
AC_SUBST([LDFLAGS])

AC_SUBST([cross_compiling])

# mpiexec
TC_MPIEXEC

# libtool
#AC_LIBTOOL_DLOPEN
#AC_DISABLE_STATIC
#AC_PROG_LIBTOOL

# check for libreadline
AX_LIB_READLINE
#[LIBREADLINE]=""
if test "$ax_cv_lib_readline" != "no"
then
  [LIBREADLINE]="$ax_cv_lib_readline"
fi
AC_SUBST([LIBREADLINE])

# check for perl
AC_PATH_PROG([PERL], [perl])
AC_SUBST([PERL])


# check for doxygen
AC_PATH_PROG([DOXYGEN], [doxygen])
AC_SUBST([DOXYGEN])
if test "x$DOXYGEN" != x
then
  # check for dot
  AC_PATH_PROG([DOT], [dot])
  [HAVE_DOT]="NO"
  if test "x$DOT" != x
  then
    [HAVE_DOT]="YES"
    [DOT_PATH=${DOT%/*}]
    AC_SUBST([DOT_PATH])
  fi
  AC_SUBST([HAVE_DOT])

  # check for latex
  AC_PATH_PROG([LATEX], [latex])
  [GENERATE_LATEX]="NO"
  if test "x$LATEX" != x
  then
    [GENERATE_LATEX]="YES"
    [LATEX_PATH=${LATEX%/*}]
    AC_SUBST([LATEX_PATH])
  fi
  AC_SUBST([GENERATE_LATEX])
fi

# check for gmsh
# NOTE: in principle, we do not need this anymore
AC_PATH_PROG([GMSH], [gmsh])
AC_SUBST([GMSH])



# on cygwin we check for Bytessence InstallMaker
#[BINSTALLMAKER="no"]
#case $host_os in
#  *cygwin* ) 
#	AC_PATH_PROGS([BINSTALLMAKER], [BInstallMaker.exe BInstallMakerA.exe], ["no"],
#	["/cygdrive/c/Programme/Bytessence\ InstallMaker":"/cygdrive/c/Programmi/Bytessence\ InstallMaker":"/cygdrive/c/Programs/Bytessence\ InstallMaker":"/cygdrive/c/Program\ Files/Bytessence\ InstallMaker"])
#	[WININSTALLERDIR=`cygpath --windows "${ac_pwd}/installer/windows"`]
#	AC_SUBST([WININSTALLERDIR])
#	[TIBERSHAREDIR=`cygpath --windows "${ac_pwd}/share"`]
#	AC_SUBST([TIBERSHAREDIR])
#	[BINSTALLMAKERDIR="`expr "$BINSTALLMAKER" : '\(.*\)/'`"]
#	[BINSTALLMAKERDIR=`cygpath --windows "$BINSTALLMAKERDIR"`]
#	AC_SUBST([BINSTALLMAKERDIR])
#        AC_DEFINE([CYGWIN], [1], [define if we are in cygwin]) ;;
#  * )  ;;
#esac
#AC_SUBST([BINSTALLMAKER])

#
# Checks for libraries.
#

# PARDISO
TC_PARDISO

# CUDA
TC_CUDA

# DFTB+
TC_DFTB

# UPTIGHT
TC_UPTIGHT

# HETERO
# NOTE: Obsolete
TC_HETERO


# without libmesh nothing works

TC_LIBMESH_PATH
TC_LIBMESH_CONFIG
if test "${LIBMESH_CONFIG:-no}" == "no"
then
  AC_MSG_NOTICE([])
  AC_MSG_ERROR([libMesh not found, but needed for tibercad!])
fi


TC_LIBMESH_PETSC_LIBS


# check for PETSc
TC_PETSC

# check for SLEPc
TC_SLEPC


# check for TAO
#AC_MSG_NOTICE([])
#AC_MSG_NOTICE([****************************************])
#AC_MSG_NOTICE([* Checking for tao])
#AC_MSG_NOTICE([****************************************])

#TC_TAO

#AC_MSG_NOTICE([])


# boost includes

TC_BOOST
if test "$tc_cv_boost" == no
then
  AC_MSG_NOTICE([])
  AC_MSG_ERROR([Boost is not available, but it is needed for tiberCAD!])
fi

TC_BOOST_LIBDIR

# boost regex
TC_BOOST_REGEX
if test "${tc_cv_boost_regex_lib-no}" == "no"; then
  AC_MSG_NOTICE([])
  AC_MSG_ERROR([Boost:regex is not available, but it is needed for tiberCAD!])
fi

# boost filesystem
TC_BOOST_FILESYSTEM
if test "${tc_cv_boost_filesystem_lib-no}" == "no"; then
  AC_MSG_NOTICE([])
  AC_MSG_ERROR([Boost:filesystem is not available, but it is needed for tiberCAD!])
fi

AC_MSG_NOTICE([])


# GSL GNU scientific library
# we don't use it as it is GPL
#TC_GSL
# add RPATH to LDFLAGS
#[GSL_LIBS=`echo "$GSL_LIBS" | sed -r 's/-L([^[:space:]]+)/-Wl,-rpath,\1 -L\1/'`]

  
#
# Checks for header files.
#
AC_CHECK_INCLUDES_DEFAULT
AC_PROG_EGREP

AC_CHECK_HEADERS([float.h stdlib.h dlfcn.h cassert])
ACX_TR1_UNORDERED_MAP
ACX_TR1_UNORDERED_SET
ACX_UNORDERED_MAP
ACX_UNORDERED_SET
ACX_EXT_HASH_MAP
ACX_EXT_HASH_SET
ACX_HASH_MAP
ACX_HASH_SET
AX_GCC_FUNC_ATTRIBUTE(visibility)


#
# Checks for typedefs, structures, and compiler characteristics.
#
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
TC_HAVE_CONSTEXPR


#
# Checks for library functions.
#
AC_FUNC_ERROR_AT_LINE
AC_CHECK_FUNCS([pow sqrt isinf isnan])
AC_CHECK_FUNCS([setenv putenv])

#
# check for other options
#

# enable/disable modules
TC_BUILD_MODULES


if test "x$DOXYGEN" == x
then
  AC_MSG_WARN([])
  AC_MSG_WARN([ doxygen was not found. You will not be able to create])
  AC_MSG_WARN([ the html or LaTeX documentation.])
  AC_MSG_WARN([])
fi

if test "x$GMSH" == x
then
  AC_MSG_WARN([])
  AC_MSG_WARN([ gmsh was not found. You will not be able to run all of the])
  AC_MSG_WARN([ examples. Install gmsh or provide complete path.])
  AC_MSG_WARN([])
fi


AC_CONFIG_FILES([Make.common
                 Make.defs
                 doc/Doxyfile
                 Makefile
                 doc/Makefile
                 src/apps/Makefile
                 src/solver/slepc/Makefile
                 modules/etb/libuptight/Makefile
                 modules/negf/libnegf/config.cmake
                 installer/Makefile
		 installer/windows/tibercad.bim
                 installer/linux/deb/control
                 installer/linux/tar/install
                 share/modules.conf
                 share/Doxyfile_SDK])
AC_CONFIG_FILES([installer/prepare_files.sh], [chmod +x installer/prepare_files.sh])
AC_CONFIG_FILES([bin/tibercad], [chmod +x bin/tibercad])
AC_CONFIG_FILES([bin/tiberCAD.sh], [chmod +x bin/tiberCAD.sh])

AC_OUTPUT
