dnl check if TiberCAD models should be built as dynamic libraries
dnl
AC_DEFUN([TC_BUILD_MODULES],
[AC_ARG_ENABLE([modules], AS_HELP_STRING([--enable-modules],
	[build TiberCAD models as dynamically loadable libraries]),
	[if test $enableval != "no"; then
	  AC_DEFINE([BUILD_TIBER_MODULES], [1],
		[Define to 1 if models should be built as modules])
	  AC_SUBST([BUILD_TIBER_MODULES], ["yes"])
  cp ${topdir}/extern/${ARCH}/lib/*.dll ${files}
	 6lse 
	 AC_SUBST([BUILD_TIBER_MODULES], ["no"])
	 fi], [AC_SUBST([BUILD_TIBER_MODULES], ["yes"])dnl
	       AC_DEFINE([BUILD_TIBER_MODULES], [1],
		[Define to 1 if models should be built as modules])])
])dnl


dnl check if the compiler has constexpr
dnl
AC_DEFUN([TC_HAVE_CONSTEXPR],
[AC_MSG_CHECKING([whether compiler supports constexpr keyword])
 AC_LANG_PUSH([C++])
 CXXFLAGS_save=$CXXFLAGS
 CXXFLAGS="-std=c++0x"
 AC_TRY_COMPILE([
  class ClassWithConst {
    constexpr static double p = 3.1415;
  };], [ClassWithConst cl;],[tc_cv_have_constexpr=yes])
 CXXFLAGS=$CXXFLAGS_save
 AC_LANG_POP()
 if test "${tc_cv_have_constexpr+set}" == "set"; then
  AC_MSG_RESULT([yes])
  AC_DEFINE([HAVE_CONSTEXPR], [1],
   [Define to 1 if compiler knows about constexpr])
 else
  AC_MSG_RESULT([no])
 fi
])dnl



dnl check for boost
dnl
AC_DEFUN([TC_BOOST],
[AC_REQUIRE([AC_CXX_NAMESPACES])dnl
AC_CACHE_VAL(tc_cv_boost_prefix,
[AC_ARG_WITH([boost-prefix], AS_HELP_STRING([--with-boost-prefix=DIR],
	[specify the path to the boost installation]),
	[tc_cv_boost_prefix="$with_boost_prefix"])
])dnl
if test "${tc_cv_boost_prefix+set}" == "set"; then
 AC_SUBST([BOOST_CPPFLAGS],["-I$tc_cv_boost_prefix/include"])
fi
AC_CACHE_VAL(tc_cv_boost,
[CPPFLAGS_save=$CPPFLAGS
 CPPFLAGS=$BOOST_CPPFLAGS
 AC_CHECK_HEADER([boost/version.hpp], tc_cv_boost=yes,
 		 tc_cv_boost=no)
 CPPFLAGS=$CPPFLAGS_save
])
if test "$tc_cv_boost" == yes; then
  AC_DEFINE([HAVE_BOOST], [1], [define if boost is available.])
fi
])dnl 

dnl set boost libdir
dnl
AC_DEFUN([TC_BOOST_LIBDIR],
[AC_CACHE_VAL(tc_cv_boost_libdir,
[AC_ARG_WITH([boost-libdir], AS_HELP_STRING([--with-boost-libdir=DIR],
	[specify the path to the boost libraries]),
	[tc_cv_boost_libdir="$with_boost_libdir"])
])dnl
if test "${tc_cv_boost_libdir+set}" == "set"; then
 AC_SUBST([BOOST_CPPFLAGS],["-I$tc_cv_boost_prefix/include"])
fi
])dnl 
		

dnl check for Boost::regex in a user defined or system directory
dnl
AC_DEFUN([TC_BOOST_REGEX],
[AC_REQUIRE([TC_BOOST])dnl
AC_CACHE_CHECK([wether Boost::regex is available], tc_cv_boost_regex_lib,
[AC_ARG_WITH([boost-regex-lib], AS_HELP_STRING([--with-boost-regex-lib=LIB],
	[specify the boost regex library or library extension]),
	[tc_cv_boost_regex_lib="$with_boost_regex_lib"])
 CXXFLAGS_save=$CXXFLAGS
 CXXFLAGS=$BOOST_CPPFLAGS
 LDFLAGS_save=$LDFLAGS
 LIBS_save=$LIBS
 if test "${tc_cv_boost_libdir+set}" == "set"; then
   tc_boost_libdir="-Wl,-rpath,${tc_cv_boost_libdir} -L${tc_cv_boost_libdir}"
 elif test "x$tc_cv_boost_prefix" != "x"; then
   tc_boost_libdir="-Wl,-rpath,${tc_cv_boost_prefix}/lib -L${tc_cv_boost_prefix}/lib"
 fi
 [tc_boost_lib="boost_regex-`$CC --version | awk '{ print $1; exit}'` boost_regex-mt boost_regex"]
 if test "${tc_cv_boost_regex_lib:+set}" == "set"; then
   tc_boost_lib="boost_regex-$tc_cv_boost_regex_lib $tc_cv_boost_regex_lib $tc_boost_lib"
 fi
 AC_LANG_PUSH([C++])
 for lib in $tc_boost_lib; do
   LDFLAGS="$tc_boost_libdir"
   LIBS="-l$lib"
   AC_LINK_IFELSE(AC_LANG_PROGRAM([[#include <boost/regex.hpp>]],
			          [[boost::regex r("s/*//"); return 0;]]),
			          [tc_cv_boost_regex_lib="$LDFLAGS $LIBS"; break])
 done
 AC_LANG_POP()
 CXXFLAGS=$CXXFLAGS_save
 LDFLAGS=$LDFLAGS_save
 LIBS=$LIBS_save])
 if test "${tc_cv_boost_regex_lib+set}" == "set"; then
   AC_SUBST([BOOST_REGEX_LIB], "$tc_cv_boost_regex_lib")
   AC_DEFINE([HAVE_BOOST_REGEX], [1], [define if Boost::regex is available])
 fi
])dnl

	

dnl check for Boost::filesystem in a user defined or system directory
dnl
AC_DEFUN([TC_BOOST_FILESYSTEM],
[AC_REQUIRE([TC_BOOST])dnl
AC_CACHE_CHECK([wether Boost::filesystem is available], tc_cv_boost_filesystem_lib,
[AC_ARG_WITH([boost-filesystem-lib], AS_HELP_STRING([--with-boost-filesystem-lib=LIB],
	[specify the boost filesystem library or library extension]),
	[tc_cv_boost_filesystem_lib="$with_boost_filesystem_lib"])
 CXXFLAGS_save=$CXXFLAGS
 CXXFLAGS=$BOOST_CPPFLAGS
 LDFLAGS_save=$LDFLAGS
 LIBS_save=$LIBS
 if test "${tc_cv_boost_libdir+set}" == "set"; then
   tc_boost_libdir="-Wl,-rpath,${tc_cv_boost_libdir} -L${tc_cv_boost_libdir}"
 elif test "x$tc_cv_boost_prefix" != "x"; then
   tc_boost_libdir="-Wl,-rpath,${tc_cv_boost_prefix}/lib -L${tc_cv_boost_prefix}/lib"
 fi
 [tc_boost_lib="boost_filesystem-`$CC --version | awk '{ print $1; exit}'` boost_filesystem-mt boost_filesystem"]
 if test "${tc_cv_boost_filesystem_lib:+set}" == "set"; then
   tc_boost_lib="boost_filesystem-$tc_cv_boost_filesystem_lib $tc_cv_boost_filesystem_lib $tc_boost_lib"
 fi
 AC_LANG_PUSH([C++])
 for lib in $tc_boost_lib; do
   LDFLAGS="$tc_boost_libdir"
   LIBS="-l$lib -l`echo $lib | sed 's/boost_filesystem/boost_system/'`"
   AC_LINK_IFELSE(AC_LANG_PROGRAM([[#include <boost/filesystem/operations.hpp>]],
			          [[boost::filesystem::path p("configure");
				    boost::filesystem::exists(p); return 0;]]),
			          [tc_cv_boost_filesystem_lib="$LDFLAGS $LIBS"; break])
 done
 AC_LANG_POP()
 CXXFLAGS=$CXXFLAGS_save
 LDFLAGS=$LDFLAGS_save
 LIBS=$LIBS_save])
 if test "${tc_cv_boost_filesystem_lib+set}" == "set"; then
   AC_SUBST([BOOST_FILESYSTEM_LIB], "$tc_cv_boost_filesystem_lib")
   AC_DEFINE([HAVE_BOOST_FILESYSTEM], [1], [define if Boost::filesystem is available])
 fi
])dnl


dnl check for MKL
dnl
AC_DEFUN([TC_MKL],
[AC_CACHE_CHECK([whether MKL is available], tc_cv_mkl_dir,
 [AC_ARG_WITH([mkl], AS_HELP_STRING([--with-mkl=DIR],
 	[specify the MKL installation path]),
	[tc_cv_mkl_dir="$with_mkl"])
  HAVE_MKL="${tc_cv_mkl_dir:-no}"
  if test "$HAVE_MKL" != no
  then
    MKL_INCLUDEDIR="$tc_cv_mkl_dir/include"
    HAVE_MKL="yes"
    case $host in
      x86_64-*-*) MKL_LIBDIR="$tc_cv_mkl_dir/lib/em64t" ;;
      i?86-*-*) MKL_LIBDIR="$tc_cv_mkl_dir/lib/32" ;;
      *) tc_cv_mkl_dir="no"; HAVE_MKL="no"; MKL_LIBDIR= ; MKL_INCLUDEDIR= ;;
    esac
    AC_SUBST([MKL_LIBDIR])
    AC_SUBST([MKL_INCLUDEDIR])
    AC_SUBST([HAVE_MKL])
  fi
 ])dnl
])dnl


dnl enable or disable PARDISO
dnl
AC_DEFUN([TC_PARDISO],
[AC_REQUIRE([TC_MKL])dnl
 AC_ARG_ENABLE([pardiso], AS_HELP_STRING([--enable-pardiso],
	[enable PARDISO linear solver]),
	[if test $enableval != "no"; then
	  AC_DEFINE([ENABLE_PARDISO], [1],
	  	[Define to 1 if pardiso is enabled])
	  AC_SUBST([ENABLE_PARDISO], ["yes"])
	else
	  AC_SUBST([ENABLE_PARDISO], ["no"])
	fi], [AC_SUBST([ENABLE_PARDISO], ["no"])])
])dnl


dnl enable or disable DFTB
dnl
AC_DEFUN([TC_DFTB],
[AC_ARG_ENABLE([dftb], AS_HELP_STRING([--enable-dftb],
	[enable DFTB+ for tight-binding]),
	[if test $enableval != "no"; then
	  AC_DEFINE([ENABLE_DFTB], [1],
	  	[Define to 1 if DFTB+ is enabled])
	  AC_SUBST([ENABLE_DFTB], ["yes"])
	else
	  AC_SUBST([ENABLE_DFTB], ["no"])
	fi], [AC_SUBST([ENABLE_DFTB], ["no"])])
])dnl


dnl enable or disable UPTIGHT
dnl
AC_DEFUN([TC_UPTIGHT],
[AC_ARG_ENABLE([uptight], AS_HELP_STRING([--enable-uptight],
	[enable UPTIGHT for tight-binding]),
	[if test $enableval != "no"; then
	  AC_DEFINE([ENABLE_UPTIGHT], [1],
	  	[Define to 1 if UPTIGHT is enabled])
	  AC_SUBST([ENABLE_UPTIGHT], ["yes"])
	else
	  AC_SUBST([ENABLE_UPTIGHT], ["no"])
	fi], [AC_SUBST([ENABLE_UPTIGHT], ["no"])])
])dnl


dnl enable or disable Hetero
dnl
AC_DEFUN([TC_HETERO],
[AC_ARG_ENABLE([hetero], AS_HELP_STRING([--enable-hetero],
	[enable HETERO for tunneling current]),
	[if test $enableval != "no"; then
	  AC_DEFINE([ENABLE_HETERO], [1],
	  	[Define to 1 if HETERO is enabled])
	  AC_SUBST([ENABLE_HETERO], ["yes"])
	else
	  AC_SUBST([ENABLE_HETERO], ["no"])
	fi], [AC_SUBST([ENABLE_HETERO], ["no"])])
])dnl


dnl enable or disable license tools
dnl
AC_DEFUN([TC_LICENSE],
[AC_ARG_ENABLE([license-check], AS_HELP_STRING([--enable-license-check],
	[enable license check]),
	[if test $enableval != "no"; then
	  AC_DEFINE([LICENSE_CHECK], [1],
	  	[Define to 1 if license check is enabled])
	  AC_SUBST([LICENSE_CHECK], ["yes"])
	else
	  AC_SUBST([LICENSE_CHECK], ["no"])
	fi], [AC_SUBST([LICENSE_CHECK], ["yes"])])
])dnl




dnl check for GNU GSL
dnl
AC_DEFUN([TC_GSL],
[tc_gsl_config="gsl-config"
 AC_ARG_WITH([gsl-config], AS_HELP_STRING([--with-gsl-config=PROG],
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




dnl check for libmesh directory
dnl for now we just get the installation path
dnl
AC_DEFUN([TC_LIBMESH_PATH],
[dnl
  AC_ARG_WITH([libmesh-prefix], AS_HELP_STRING([--with-libmesh-prefix=DIR],
	[specify the libmesh installation prefix]),
	[tc_libmesh_prefix=$with_libmesh_prefix])
  dnl check for something there
  CPPFLAGS_save="$CPPFLAGS"
  CPPFLAGS="-I$tc_libmesh_prefix/include/base"
  AC_CHECK_HEADER([libmesh_config.h], [AC_SUBST([LIBMESH_PREFIX],
  					  "$tc_libmesh_prefix")])
  CPPFLAGS="$CPPFLAGS_save"
])dnl




dnl check for libmesh-config
dnl
AC_DEFUN([TC_LIBMESH_CONFIG],
[dnl
 tc_libmesh_config="libmesh-config"
 AC_ARG_WITH([libmesh-config], AS_HELP_STRING([--with-libmesh-config=PROG],
	[specify the libmesh-config tool]),
	[tc_libmesh_config=$with_libmesh_config])
  [tc_libmesh_config_name="`expr "//$tc_libmesh_config" : '.*/\([^/]*\)'`"]
  [tc_libmesh_config_path="`expr "$tc_libmesh_config" : '\(.*\)/'`"]
  if test "${LIBMESH_PREFIX:+set}" == "set"; then
    tc_libmesh_config_path="${tc_libmesh_config_path}:${LIBMESH_PREFIX}/contrib/bin"
  fi
  dnl check if it exists
  AC_PATH_PROG([LIBMESH_CONFIG], [$tc_libmesh_config_name],
  	       [""], [$tc_libmesh_config_path])
  AC_MSG_NOTICE([])
])dnl



dnl check for shared libraries in build directory
dnl
AC_DEFUN([TC_LIBMESH_PETSC_LIBS],
[AC_REQUIRE([TC_LIBMESH_CONFIG])dnl
 tc_libmesh_petsc_libdir=$(pwd)/extern/$(./share/arch.sh ${host_cpu}-${host_os})/lib
 if test "${LIBMESH_PREFIX:+set}" == "set"; then
   tc_libmesh_petsc_libdir=$(dirname ${LIBMESH_PREFIX})/$(./share/arch.sh ${host_cpu}-${host_os})/lib
 fi
 dnl
 dnl link test
 CXXFLAGS_save=$CXXFLAGS
 LDFLAGS_save=$LDFLAGS
 AC_LANG_PUSH([C++])
 LDFLAGS="-L${tc_libmesh_petsc_libdir} -Wl,-rpath,${tc_libmesh_petsc_libdir} -lmesh -lpetsc_real"
 echo "$LDFLAGS"
 AC_LINK_IFELSE(AC_LANG_PROGRAM([], []),
 			        [tc_cv_have_sharedlibs="yes"])
 echo "$tc_cv_have_sharedlibs"
 if test "${tc_cv_have_sharedlibs:+set}" != "set"; then
 echo "$LDFLAGS"
  LDFLAGS="-L${tc_libmesh_petsc_libdir} -Wl,-rpath,${tc_libmesh_petsc_libdir} -lmesh -lpetsc_real -ldl"
  AC_LINK_IFELSE(AC_LANG_PROGRAM([], []),
  			         [tc_cv_have_sharedlibs="yes"])
 fi
 AC_LANG_POP()
 if test "${tc_cv_have_sharedlibs}" == "yes"; then
   tc_libmesh_petsc_libs=${LDFLAGS}
 else
   tc_libmesh_petsc_libs=$(${LIBMESH_CONFIG} --ldflags)
 fi
 LDFLAGS=$LDFLAGS_save
 AC_SUBST([LIBMESH_PETSC_LIBS], "$tc_libmesh_petsc_libs")
])dnl



dnl check for complex petsc
dnl
AC_DEFUN([TC_COMPLEX_PETSC],
[dnl
 AC_ARG_WITH([petsc-prefix], AS_HELP_STRING([--with-petsc-prefix=DIR],
	[specify the (complex) PETSc installation prefix]),
	[tc_petsc_prefix=$with_petsc_prefix])
 AC_ARG_WITH([petsc-arch], AS_HELP_STRING([--with-petsc-arch=ARCH],
	[specify the (complex) PETSc ARCH to be used]),
	[tc_petsc_arch=$with_petsc_arch])
 AC_SUBST([COMPLEX_PETSC_DIR], "$tc_petsc_prefix")
 AC_SUBST([COMPLEX_PETSC_ARCH], "$tc_petsc_arch")
])dnl



dnl check for SLEPc
dnl
AC_DEFUN([TC_SLEPC],
[AC_REQUIRE([TC_COMPLEX_PETSC])dnl
 AC_ARG_WITH([slepc-prefix], AS_HELP_STRING([--with-slepc-prefix=DIR],
	[specify the SLEPc installation prefix]),
	[tc_slepc_prefix=$with_slepc_prefix])
 dnl
 dnl compile test
 CXXFLAGS_save=$CXXFLAGS
 LDFLAGS_save=$LDFLAGS
 AC_LANG_PUSH([C++])
 CXXFLAGS="-I${tc_slepc_prefix}/include"
 LDFLAGS="-L${tc_slepc_prefix}/lib/${tc_petsc_arch} -lslepc \
  	 -L${tc_petsc_prefix}/lib/${tc_petsc_arch} -lpetscksp"
 AC_LINK_IFELSE(AC_LANG_PROGRAM([[#include "slepceps.h"]],
 			        [[SlepcInitialize(0,0,0,0);
 				  SlepcFinalize();]]),
 			        [tc_cv_have_slepc="yes"])
  AC_LANG_POP()
  CXXFLAGS=$CXXFLAGS_save
  LDFLAGS=$LDFLAGS_save
 AC_SUBST([SLEPC_DIR], "$tc_slepc_prefix")
 if test "$tc_cv_have_slepc" == "yes"; then
   AC_DEFINE([HAVE_COMPLEX_SLEPC], [1], [Define to 1 if complex SLEPc is available])
 fi
])dnl




dnl check for TAO
dnl
AC_DEFUN([TC_TAO],
[dnl
 AC_ARG_WITH([tao-prefix], AS_HELP_STRING([--with-tao-prefix=DIR],
	[specify the TAO installation prefix]),
	[tc_tao_prefix=$with_tao_prefix])
 dnl
 dnl compile test
dnl CXXFLAGS_save=$CXXFLAGS
dnl LDFLAGS_save=$LDFLAGS
dnl AC_LANG_PUSH([C++])
dnl CXXFLAGS=-I${tc_tao_prefix}/include
dnl LDFLAGS=-L${tc_tao_prefix}/lib/${tc_petsc_arch} -ltao \
dnl  	 -L${tc_tao_prefix}/lib/${tc_petsc_arch} -ltaopetsc
dnl AC_LINK_IFELSE(AC_LANG_PROGRAM([[#include "tao.h"]],
dnl 			        [[TaoInitialize(0,0,0,0);
dnl 				  TaoFinalize();]]);
dnl 			        [tc_cv_have_tao="yes"])
 dnl AC_LANG_POP()
 dnl CXXFLAGS=$CXXFLAGS_save
 dnl LDFLAGS=$LDFLAGS_save
 AC_SUBST([TAO_DIR], "$tc_tao_prefix")
 AC_SUBST([TAO_INCLUDES], "-I${tc_tao_prefix} -I${tc_tao_prefix}/include")
 if test "$tc_cv_have_tao" == "yes"; then
   AC_DEFINE([HAVE_TAO], [1], [Define to 1 if tao is available])
 fi
])dnl


dnl which SVN version to use
dnl
AC_DEFUN([TC_SUBVERSION],
[tc_subversion="svn"
 AC_ARG_WITH([subversion], AS_HELP_STRING([--with-subversion=PROG],
 	[specify the subversion executable]),
	[tc_subversion=$with_subversion])
 AC_PATH_PROG([SVN], [$tc_subversion])])dnl



dnl ----------------------------------------------------------------------------
dnl Check to see if the compiler can compile a test program using
dnl std::tr1::unordered_map
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_TR1_UNORDERED_MAP],
[AC_CACHE_CHECK(whether the compiler supports std::tr1::unordered_map,
ac_cv_cxx_tr1_unordered_map,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <tr1/unordered_map>],
[
  std::tr1::unordered_map<int, int> m;
  m.insert(std::make_pair(1, 2));
],
 ac_cv_cxx_tr1_unordered_map=yes, ac_cv_cxx_tr1_unordered_map=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_tr1_unordered_map" = yes; then
  AC_DEFINE(HAVE_TR1_UNORDERED_MAP,1,
            [define if the compiler supports std::tr1::unordered_map])
fi
])



dnl ----------------------------------------------------------------------------
dnl Check to see if the compiler can compile a test program using
dnl std::tr1::unordered_set
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_TR1_UNORDERED_SET],
[AC_CACHE_CHECK(whether the compiler supports std::tr1::unordered_set,
ac_cv_cxx_tr1_unordered_set,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <tr1/unordered_set>],
[
  std::tr1::unordered_set<int> m;
  m.insert(1);  m.insert(2);
],
 ac_cv_cxx_tr1_unordered_set=yes, ac_cv_cxx_tr1_unordered_set=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_tr1_unordered_set" = yes; then
  AC_DEFINE(HAVE_TR1_UNORDERED_SET,1,
            [define if the compiler supports std::tr1::unordered_set])
fi
])



dnl ----------------------------------------------------------------------------
dnl Check to see if the compiler can compile a test program using
dnl std::unordered_map
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_UNORDERED_MAP],
[AC_CACHE_CHECK(whether the compiler supports std::unordered_map,
ac_cv_cxx_unordered_map,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <unordered_map>],
[
  std::unordered_map<int, int> m;
  m.insert(std::make_pair(1, 2));
],
 ac_cv_cxx_unordered_map=yes, ac_cv_cxx_unordered_map=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_unordered_map" = yes; then
  AC_DEFINE(HAVE_UNORDERED_MAP,1,
            [define if the compiler supports std::unordered_map])
fi
])



dnl ----------------------------------------------------------------------------
dnl Check to see if the compiler can compile a test program using
dnl std::unordered_set
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_UNORDERED_SET],
[AC_CACHE_CHECK(whether the compiler supports std::unordered_set,
ac_cv_cxx_unordered_set,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <unordered_set>],
[
  std::unordered_set<int, int> m;
  m.insert(std::make_pair(1, 2));
],
 ac_cv_cxx_unordered_set=yes, ac_cv_cxx_unordered_set=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_unordered_set" = yes; then
  AC_DEFINE(HAVE_UNORDERED_SET,1,
            [define if the compiler supports std::unordered_set])
fi
])




dnl ----------------------------------------------------------------------------
dnl Check to see if the compiler can compile a test program using
dnl __gnu_cxx::hash_map
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_EXT_HASH_MAP],
[AC_CACHE_CHECK(whether the compiler supports __gnu_cxx::hash_map,
ac_cv_cxx_ext_hash_map,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <ext/hash_map>],
[
  __gnu_cxx::hash_map<int, int> m;
  m.insert(std::make_pair(1, 2));
],
 ac_cv_cxx_ext_hash_map=yes, ac_cv_cxx_ext_hash_map=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_ext_hash_map" = yes; then
  AC_DEFINE(HAVE_EXT_HASH_MAP,1,
            [define if the compiler supports __gnu_cxx::hash_map])
fi
])



dnl ----------------------------------------------------------------------------
dnl Check to see if the compiler can compile a test program using
dnl __gnu_cxx::hash_set
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_EXT_HASH_SET],
[AC_CACHE_CHECK(whether the compiler supports __gnu_cxx::hash_set,
ac_cv_cxx_ext_hash_set,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <ext/hash_set>],
[
  __gnu_cxx::hash_set<int> m;
  m.insert(1);  m.insert(2);
],
 ac_cv_cxx_ext_hash_set=yes, ac_cv_cxx_ext_hash_set=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_ext_hash_set" = yes; then
  AC_DEFINE(HAVE_EXT_HASH_SET,1,
            [define if the compiler supports __gnu_cxx::hash_set])
fi
])



dnl ----------------------------------------------------------------------------
dnl Check to see if the compiler can compile a test program using
dnl std::hash_map (This is unlikely...)
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_HASH_MAP],
[AC_CACHE_CHECK(whether the compiler supports std::hash_map,
ac_cv_cxx_hash_map,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <hash_map>],
[
  std::hash_map<int, int> m;
  m.insert(std::make_pair(1, 2));
],
 ac_cv_cxx_hash_map=yes, ac_cv_cxx_hash_map=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_hash_map" = yes; then
  AC_DEFINE(HAVE_HASH_MAP,1,
            [define if the compiler supports std::hash_map])
fi
])



dnl ----------------------------------------------------------------------------
dnl Check to see if the compiler can compile a test program using
dnl std::hash_set (This is unlikely...)
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_HASH_SET],
[AC_CACHE_CHECK(whether the compiler supports std::hash_set,
ac_cv_cxx_hash_set,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <hash_set>],
[
  std::hash_set<int> m;
  m.insert(1);  m.insert(2);
],
 ac_cv_cxx_hash_set=yes, ac_cv_cxx_hash_set=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_hash_set" = yes; then
  AC_DEFINE(HAVE_HASH_SET,1,
            [define if the compiler supports std::hash_set])
fi
])
