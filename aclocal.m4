# generated automatically by aclocal 1.10 -*- Autoconf -*-

# Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004,
# 2005, 2006  Free Software Foundation, Inc.
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.

# ===========================================================================
#     http://www.gnu.org/software/autoconf-archive/ax_cxx_namespaces.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_CXX_NAMESPACES
#
# DESCRIPTION
#
#   If the compiler can prevent names clashes using namespaces, define
#   HAVE_NAMESPACES.
#
# LAST MODIFICATION
#
#   2004-02-04
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 6

AU_ALIAS([AC_CXX_NAMESPACES], [AX_CXX_NAMESPACES])
AC_DEFUN([AX_CXX_NAMESPACES],
[AC_CACHE_CHECK(whether the compiler implements namespaces,
ax_cv_cxx_namespaces,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([namespace Outer { namespace Inner { int i = 0; }}],
                [using namespace Outer::Inner; return i;],
 ax_cv_cxx_namespaces=yes, ax_cv_cxx_namespaces=no)
 AC_LANG_RESTORE
])
if test "$ax_cv_cxx_namespaces" = yes; then
  AC_DEFINE(HAVE_NAMESPACES,,[define if the compiler implements namespaces])
fi
])

m4_include([acinclude.m4])
