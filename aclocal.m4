# generated automatically by aclocal 1.16.1 -*- Autoconf -*-

# Copyright (C) 1996-2018 Free Software Foundation, Inc.

# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.

m4_ifndef([AC_CONFIG_MACRO_DIRS], [m4_defun([_AM_CONFIG_MACRO_DIRS], [])m4_defun([AC_CONFIG_MACRO_DIRS], [_AM_CONFIG_MACRO_DIRS($@)])])
# ===========================================================================
#    https://www.gnu.org/software/autoconf-archive/ax_cxx_namespaces.html
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
# LICENSE
#
#   Copyright (c) 2008 Todd Veldhuizen
#   Copyright (c) 2008 Luc Maisonobe <luc@spaceroots.org>
#   Copyright (c) 2013 Bastien Roucaries <roucaries.bastien+autoconf@gmail.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 9

AU_ALIAS([AC_CXX_NAMESPACES], [AX_CXX_NAMESPACES])
AC_DEFUN([AX_CXX_NAMESPACES],
[AC_CACHE_CHECK(whether the compiler implements namespaces,
ax_cv_cxx_namespaces,
[AC_LANG_PUSH([C++])
 AC_COMPILE_IFELSE([AC_LANG_SOURCE([namespace Outer { namespace Inner { int i = 0; }}
                                   using namespace Outer::Inner; int foo(void) { return i;} ])],
                   ax_cv_cxx_namespaces=yes, ax_cv_cxx_namespaces=no)
 AC_LANG_POP
])
if test "$ax_cv_cxx_namespaces" = yes; then
  AC_DEFINE(HAVE_NAMESPACES,,[define if the compiler implements namespaces])
fi
])

# ===========================================================================
#     https://www.gnu.org/software/autoconf-archive/ax_lib_readline.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_LIB_READLINE
#
# DESCRIPTION
#
#   Searches for a readline compatible library. If found, defines
#   `HAVE_LIBREADLINE'. If the found library has the `add_history' function,
#   sets also `HAVE_READLINE_HISTORY'. Also checks for the locations of the
#   necessary include files and sets `HAVE_READLINE_H' or
#   `HAVE_READLINE_READLINE_H' and `HAVE_READLINE_HISTORY_H' or
#   'HAVE_HISTORY_H' if the corresponding include files exists.
#
#   The libraries that may be readline compatible are `libedit',
#   `libeditline' and `libreadline'. Sometimes we need to link a termcap
#   library for readline to work, this macro tests these cases too by trying
#   to link with `libtermcap', `libcurses' or `libncurses' before giving up.
#
#   Here is an example of how to use the information provided by this macro
#   to perform the necessary includes or declarations in a C file:
#
#     #ifdef HAVE_LIBREADLINE
#     #  if defined(HAVE_READLINE_READLINE_H)
#     #    include <readline/readline.h>
#     #  elif defined(HAVE_READLINE_H)
#     #    include <readline.h>
#     #  else /* !defined(HAVE_READLINE_H) */
#     extern char *readline ();
#     #  endif /* !defined(HAVE_READLINE_H) */
#     char *cmdline = NULL;
#     #else /* !defined(HAVE_READLINE_READLINE_H) */
#       /* no readline */
#     #endif /* HAVE_LIBREADLINE */
#
#     #ifdef HAVE_READLINE_HISTORY
#     #  if defined(HAVE_READLINE_HISTORY_H)
#     #    include <readline/history.h>
#     #  elif defined(HAVE_HISTORY_H)
#     #    include <history.h>
#     #  else /* !defined(HAVE_HISTORY_H) */
#     extern void add_history ();
#     extern int write_history ();
#     extern int read_history ();
#     #  endif /* defined(HAVE_READLINE_HISTORY_H) */
#       /* no history */
#     #endif /* HAVE_READLINE_HISTORY */
#
# LICENSE
#
#   Copyright (c) 2008 Ville Laurikari <vl@iki.fi>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 7

AU_ALIAS([VL_LIB_READLINE], [AX_LIB_READLINE])
AC_DEFUN([AX_LIB_READLINE], [
  AC_CACHE_CHECK([for a readline compatible library],
                 ax_cv_lib_readline, [
    ORIG_LIBS="$LIBS"
    for readline_lib in readline edit editline; do
      for termcap_lib in "" termcap curses ncurses; do
        if test -z "$termcap_lib"; then
          TRY_LIB="-l$readline_lib"
        else
          TRY_LIB="-l$readline_lib -l$termcap_lib"
        fi
        LIBS="$ORIG_LIBS $TRY_LIB"
        AC_TRY_LINK_FUNC(readline, ax_cv_lib_readline="$TRY_LIB")
        if test -n "$ax_cv_lib_readline"; then
          break
        fi
      done
      if test -n "$ax_cv_lib_readline"; then
        break
      fi
    done
    if test -z "$ax_cv_lib_readline"; then
      ax_cv_lib_readline="no"
    fi
    LIBS="$ORIG_LIBS"
  ])

  if test "$ax_cv_lib_readline" != "no"; then
    LIBS="$LIBS $ax_cv_lib_readline"
    AC_DEFINE(HAVE_LIBREADLINE, 1,
              [Define if you have a readline compatible library])
    AC_CHECK_HEADERS(readline.h readline/readline.h)
    AC_CACHE_CHECK([whether readline supports history],
                   ax_cv_lib_readline_history, [
      ax_cv_lib_readline_history="no"
      AC_TRY_LINK_FUNC(add_history, ax_cv_lib_readline_history="yes")
    ])
    if test "$ax_cv_lib_readline_history" = "yes"; then
      AC_DEFINE(HAVE_READLINE_HISTORY, 1,
                [Define if your readline library has \`add_history'])
      AC_CHECK_HEADERS(history.h readline/history.h)
    fi
  fi
])dnl

m4_include([acinclude.m4])
