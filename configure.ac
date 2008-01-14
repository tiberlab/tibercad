#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT([TiberCAD], [0.1.0], [developers@tibercad.org])
AC_PREREQ(2.59)

AC_CONFIG_SRCDIR([include/base/common/TiberCad.h])
AC_CONFIG_HEADER([include/base/common/tiber_config.h])

AC_REVISION([$Revision: 0.1$])

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

# libtool
#AC_LIBTOOL_DLOPEN
#AC_DISABLE_STATIC
#AC_PROG_LIBTOOL

# check for libreadline
VL_LIB_READLINE
[LIBREADLINE]=""
if test "$vl_cv_lib_readline" != "no"
then
  [LIBREADLINE]="$vl_cv_lib_readline"
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
AC_PATH_PROG([GMSH], [gmsh])
AC_SUBST([GMSH])


# on cygwin we check for Bytessence InstallMaker
[BINSTALLMAKER="no"]
case $host_os in
  *cygwin* ) 
	AC_PATH_PROG([BINSTALLMAKER], [BInstallMaker.exe], ["no"],
	["/cygdrive/c/Programme/Bytessence\ InstallMaker":"/cygdrive/c/Programmi/Bytessence\ InstallMaker":"/cygdrive/c/Programs/Bytessence\ InstallMaker":"/cygdrive/c/Program\ Files/Bytessence\ InstallMaker"])
	[WININSTALLERDIR=`cygpath --path --windows "${ac_pwd}/installer/windows"`]
	AC_SUBST([WININSTALLERDIR])
	[TIBERSHAREDIR=`cygpath --path --windows "${ac_pwd}/share"`]
	AC_SUBST([TIBERSHAREDIR])
	[BINSTALLMAKERDIR="`expr "$BINSTALLMAKER" : '\(.*\)/'`"]
	[BINSTALLMAKERDIR=`cygpath --path --windows "$BINSTALLMAKERDIR"`]
	AC_SUBST([BINSTALLMAKERDIR]) ;;
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
TC_GSL
# add RPATH to LDFLAGS
[GSL_LIBS=`echo "$GSL_LIBS" | sed -r 's/-L([^[:space:]]+)/-Wl,-rpath,\1 -L\1/'`]


# PARDISO
TC_PARDISO

# DFTB+
TC_DFTB

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

if test "${GSL_CONFIG+set}" != "set"
then
  AC_MSG_NOTICE([])
  AC_MSG_WARN([ GSL was not found. You will not be able to use])
  AC_MSG_WARN([ Fermi-Dirac statistics.])
  AC_MSG_NOTICE([])
fi


AC_CONFIG_FILES([stamp-h], [echo timestamp > stamp-h])
AC_CONFIG_FILES([Make.common
                 doc/Doxyfile
                 Makefile
                 doc/Makefile
                 src/tibercad/Makefile
                 src/utils/DFTB/Makefile
                 src/utils/SLEPC/Makefile
		 installer/windows/tibercad.bim])

AC_OUTPUT
