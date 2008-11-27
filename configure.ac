#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT([TiberCAD], [1.0], [support@tibercad.org])
AC_PREREQ(2.59)

AC_CONFIG_SRCDIR([include/base/common/TiberCad.h])
AC_CONFIG_HEADER([include/base/common/tiber_config.h
                  include/base/common/tiber_version.h])

AC_REVISION([$Revision: 1.0$])

tc_fullversion=`cat VERSION | sed 's/\./ /g'`
tc_version_major=`echo $tc_fullversion | awk '{print $1}'`
tc_version_minor=`echo $tc_fullversion | awk '{print $2}'`
tc_version_subminor=`echo $tc_fullversion | awk '{print $3}'`
AC_DEFINE_UNQUOTED([TIBERMAJORVERSION],[$tc_version_major],
  [TiberCAD version major])
AC_DEFINE_UNQUOTED([TIBERMINORVERSION],[$tc_version_minor],
  [TiberCAD version minor])
AC_DEFINE_UNQUOTED([TIBERSUBMINORVERSION],[$tc_version_subminor],
  [TiberCAD version subminor])

[TIBERVERSION]="$tc_version_major.$tc_version_minor.$tc_version_subminor"
AC_SUBST([TIBERVERSION])

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
AC_FC_LIBRARY_LDFLAGS
AC_SUBST([FCLIBS])

# libtool
#AC_LIBTOOL_DLOPEN
#AC_DISABLE_STATIC
#AC_PROG_LIBTOOL

# check for libreadline
# we cannot use it directly, it is GPL! We are not, sigh..
# VL_LIB_READLINE
[LIBREADLINE]=""
#if test "$vl_cv_lib_readline" != "no"
#then
#  [LIBREADLINE]="$vl_cv_lib_readline"
#fi
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
AC_PATH_PROG([GMSH], [gmsh])
AC_SUBST([GMSH])


# on cygwin we check for Bytessence InstallMaker
[BINSTALLMAKER="no"]
case $host_os in
  *cygwin* ) 
	AC_PATH_PROGS([BINSTALLMAKER], [BInstallMaker.exe BInstallMakerA.exe], ["no"],
	["/cygdrive/c/Programme/Bytessence\ InstallMaker":"/cygdrive/c/Programmi/Bytessence\ InstallMaker":"/cygdrive/c/Programs/Bytessence\ InstallMaker":"/cygdrive/c/Program\ Files/Bytessence\ InstallMaker"])
	[WININSTALLERDIR=`cygpath --windows "${ac_pwd}/installer/windows"`]
	AC_SUBST([WININSTALLERDIR])
	[TIBERSHAREDIR=`cygpath --windows "${ac_pwd}/share"`]
	AC_SUBST([TIBERSHAREDIR])
	[BINSTALLMAKERDIR="`expr "$BINSTALLMAKER" : '\(.*\)/'`"]
	[BINSTALLMAKERDIR=`cygpath --windows "$BINSTALLMAKERDIR"`]
	AC_SUBST([BINSTALLMAKERDIR])
        AC_DEFINE([CYGWIN], [1], [define if we are in cygwin]) ;;
  * )  ;;
esac
AC_SUBST([BINSTALLMAKER])

#
# Checks for libraries.
#

# without libmesh nothing works
AC_MSG_NOTICE([])
AC_MSG_NOTICE([****************************************])
AC_MSG_NOTICE([* Checking for libMesh])
AC_MSG_NOTICE([****************************************])

TC_LIBMESH_PATH
TC_LIBMESH_CONFIG
if test "${LIBMESH_CONFIG:-no}" == "no"
then
  AC_MSG_NOTICE([])
  AC_MSG_ERROR([libMesh not found, but needed for tibercad!])
fi

AC_MSG_NOTICE([])


# check for SLEPc
AC_MSG_NOTICE([])
AC_MSG_NOTICE([****************************************])
AC_MSG_NOTICE([* Checking for (complex) SLEPc])
AC_MSG_NOTICE([****************************************])

TC_SLEPC

AC_MSG_NOTICE([])


# boost includes
AC_MSG_NOTICE([])
AC_MSG_NOTICE([****************************************])
AC_MSG_NOTICE([* Checking for boost...])
AC_MSG_NOTICE([****************************************])

TC_BOOST
if test "$tc_cv_boost" == no
then
  AC_MSG_NOTICE([])
  AC_MSG_ERROR([Boost is not available, but it is needed for TiberCAD!])
fi

# boost regex
TC_BOOST_REGEX
if test "${tc_cv_boost_regex_lib-no}" == "no"; then
  AC_MSG_NOTICE([])
  AC_MSG_ERROR([Boost:regex is not available, but it is needed for TiberCAD!])
fi

# boost filesystem
TC_BOOST_FILESYSTEM
if test "${tc_cv_boost_filesystem_lib-no}" == "no"; then
  AC_MSG_NOTICE([])
  AC_MSG_ERROR([Boost:filesystem is not available, but it is needed for TiberCAD!])
fi

AC_MSG_NOTICE([])


# GSL GNU scientific library
# we don't use it as it is GPL
#TC_GSL
# add RPATH to LDFLAGS
#[GSL_LIBS=`echo "$GSL_LIBS" | sed -r 's/-L([^[:space:]]+)/-Wl,-rpath,\1 -L\1/'`]


# PARDISO
TC_PARDISO

# DFTB+
TC_DFTB

# UPTIGHT
TC_UPTIGHT

# HETERO
TC_HETERO
  
#
# Checks for header files.
#
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h stdlib.h dlfcn.h cassert])


#
# Checks for typedefs, structures, and compiler characteristics.
#
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T


#
# Checks for library functions.
#
AC_FUNC_ERROR_AT_LINE
AC_CHECK_FUNCS([pow sqrt isinf isnan])

#
# check for other options
#

# enable/disable modules
TC_BUILD_MODULES

# use licensing tools
TC_LICENSE



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


AC_CONFIG_FILES([stamp-h], [echo timestamp > stamp-h])
AC_CONFIG_FILES([Make.common
                 doc/Doxyfile
                 Makefile
                 doc/Makefile
                 src/tibercad/Makefile
                 src/utils/DFTB/Makefile
                 src/utils/UPTIGHT/Makefile
                 src/utils/SLEPC/Makefile
                 installer/Makefile
		 installer/windows/tibercad.bim
                 installer/linux/deb/control
                 installer/linux/tar/install])
AC_CONFIG_FILES([installer/prepare_files.sh], [chmod +x installer/prepare_files.sh])

AC_OUTPUT
