dnl check for boost
dnl
AC_DEFUN([TC_BOOST],
[AC_REQUIRE([AC_CXX_NAMESPACES])dnl
AC_MSG_NOTICE([])
AC_MSG_NOTICE([****************************************])
AC_MSG_NOTICE([* Checking for boost...])
AC_MSG_NOTICE([****************************************])
AC_CACHE_VAL(tc_cv_boost_prefix,
[AC_ARG_WITH([boost-prefix], AS_HELP_STRING([--with-boost-prefix],
	[specify the path to the boost installation]),
	[tc_cv_boost_prefix="$with_boost_prefix"])
])dnl
if test "${tc_cv_boost_prefix+set}" == "set"; then
 AC_SUBST([BOOST_CXXFLAGS],["-I$tc_cv_boost_prefix/include"])
fi
AC_CACHE_VAL(tc_cv_boost,
[CPPFLAGS_save=$CPPFLAGS
 CPPFLAGS=$BOOST_CXXFLAGS
 AC_CHECK_HEADER([boost/version.hpp], tc_cv_boost=yes,
 		 tc_cv_boost=no)
 CPPFLAGS=$CPPFLAGS_save
])
if test "$tc_cv_boost" == yes; then
  AC_DEFINE([HAVE_BOOST], [1], [define if boost is available.])
  AC_MSG_NOTICE([Ok. Boost is out there.])
  AC_MSG_NOTICE([])
fi
])dnl 
	

dnl check for Boost::regex in a user defined or system directory
dnl
AC_DEFUN([TC_BOOST_REGEX],
[AC_REQUIRE([TC_BOOST])dnl
AC_CACHE_CHECK([wether Boost::regex is available], tc_cv_boost_regex_lib,
[AC_ARG_WITH([boost-regex-lib], AS_HELP_STRING([--with-boost-regex-lib],
	[specify the boost regex library]),
	[tc_cv_boost_regex_lib="$with_boost_regex_lib"])
 CXXFLAGS_save=$CXXFLAGS
 CXXFLAGS=$BOOST_CXXFLAGS
 LDFLAGS_save=$LDFLAGS
 if test "x$tc_cv_boost_prefix" != "x"; then
   tc_boost_libdir="-Wl,-rpath,${tc_cv_boost_prefix}/lib -L${tc_cv_boost_prefix}/lib"
 fi
 [tc_boost_lib="boost_regex-`$CC --version | awk '{ print $1; exit}'` boost_regex"]
 if test "${tc_cv_boost_regex_lib:+set}" == "set"; then
   tc_boost_lib="boost_regex-$tc_cv_boost_regex_lib $tc_cv_boost_regex_lib $tc_boost_lib"
 fi
 AC_LANG_PUSH([C++])
 for lib in $tc_boost_lib; do
   LDFLAGS="$tc_boost_libdir -l$lib"
   AC_LINK_IFELSE(AC_LANG_PROGRAM([[#include <boost/regex.hpp>]],
			          [[boost::regex r("s/*//"); return 0;]]),
			          [tc_cv_boost_regex_lib="$LDFLAGS"; break])
 done
 AC_LANG_POP()
 CXXFLAGS=$CXXFLAGS_save
 LDFLAGS=$LDFLAGS_save])
 if test "${tc_cv_boost_regex_lib+set}" == "set"; then
   AC_SUBST([BOOST_REGEX_LIB], "$tc_cv_boost_regex_lib")
   AC_DEFINE([HAVE_BOOST_REGEX], [1], [define if Boost::regex is available])
 fi
])dnl


dnl check for GNU GSL
dnl
AC_DEFUN([TC_GSL],
[tc_gsl_config="gsl-config"
 AC_ARG_WITH([gsl-config], AS_HELP_STRING([--with-gsl-config],
	[specify the GSL config script]),
	[tc_gsl_config=$with_gsl_config])
 dnl an absolute path could have been provided
 AC_PATH_PROG([GSL_CONFIG], ["$tc_gsl_config"], [], [${PATH}:/])
 if test "${GSL_CONFIG+set}" == "set"; then
   AC_SUBST([GSL_LIBS], [`$GSL_CONFIG --libs`])
   AC_SUBST([GSL_CFLAGS], [`$GSL_CONFIG --cflags`])
   AC_DEFINE([HAVE_GSL], [1], [Define to 1 if you have the GSL library.])
 fi
])dnl


dnl check for libmesh
dnl for now we just get the installation path
AC_DEFUN([TC_LIBMESH],
[AC_MSG_NOTICE([])
 AC_MSG_NOTICE([****************************************])
 AC_MSG_NOTICE([* Checking for libMesh...])
 AC_MSG_NOTICE([****************************************])
 AC_ARG_WITH([libmesh-prefix], AS_HELP_STRING([--with-libmesh-prefix],
	[specify the libmesh installation prefix]),
	[tc_libmesh_prefix=$with_libmesh_prefix])
echo "$tc_libmesh_prefix"
  dnl check for something there
  CPPFLAGS_save="$CPPFLAGS"
  CPPFLAGS="-I$tc_libmesh_prefix/include/base"
  AC_CHECK_HEADER([libmesh_config.h], [AC_SUBST([LIBMESH_PREFIX],
  					  "$tc_libmesh_prefix")])
  CPPFLAGS="$CPPFLAGS_save"
  if test "${LIBMESH_PREFIX:+set}" == "set"; then
    AC_MSG_NOTICE([libMesh installation found in $tc_libmesh_prefix])
    AC_MSG_NOTICE([])
  else
    AC_MSG_NOTICE([No libMesh installation found!])
    AC_MSG_NOTICE([])
  fi
])dnl
