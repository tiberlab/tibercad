#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT([TiberCAD], [1.0], [support@tibercad.org])
AC_PREREQ(2.59)

AC_CONFIG_SRCDIR([.])

AC_REVISION([$Revision: 1.0$])

# Make compilation tests using C++
AC_LANG([C++])


# The host
AC_CANONICAL_HOST


#
# Check for programs.
#

# compilers
AC_PROG_CPP
AC_PROG_CC
AC_PROG_CXX
AC_PROG_FC
AC_PROG_RANLIB
AC_FC_LIBRARY_LDFLAGS
AC_SUBST([FCLIBS])
AC_SUBST([LDFLAGS])

# libtool
#AC_LIBTOOL_DLOPEN
#AC_DISABLE_STATIC
#AC_PROG_LIBTOOL

AC_CONFIG_FILES([common.sh])
#AC_CONFIG_FILES([installer/prepare_files.sh], [chmod +x installer/prepare_files.sh])

AC_OUTPUT
