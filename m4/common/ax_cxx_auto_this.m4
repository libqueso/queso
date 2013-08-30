# ============================================================================
#  http://www.gnu.org/software/autoconf-archive/ax_cxx_auto_this.html
# ============================================================================
#
# SYNOPSIS
#
#   AX_CXX_AUTO_THIS([mandatory|optional])
#
# DESCRIPTION
#
#   This macro is based on AX_CXX_COMPILE_STDCXX_11, which can
#   usefully precede it to configure CXXFLAGS for C++11 support.
#
#   Check for coverage in the compiler for the ability to reference
#   "this" in C++11 trailing return types in class member functions.
#   CXXFLAGS should already include any flags necessary to enable
#   C++11 or better language support.
#
#   The argument, if specified 'mandatory' or if left unspecified,
#   indicates that support for "this" in trailing return types is
#   required and that the macro should error out if no mode with that
#   support is found.  If specified 'optional', then configuration
#   proceeds regardless, after defining HAVE_CXX_AUTO_THIS if and only
#   if a supporting mode is found.
#
# LICENSE
#
#   Copyright (c) 2008 Benjamin Kosnik <bkoz@redhat.com>
#   Copyright (c) 2012 Zack Weinberg <zackw@panix.com>
#   Copyright (c) 2013 Roy Stogner <roystgnr@ices.utexas.edu>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 3

m4_define([_AX_CXX_AUTO_THIS_testbody], [
    struct check
    {
      int f() { return 0; }
      auto g() -> decltype(this->f()) { return f(); }
    };

    check c;
    int gval = c->g();
])

AC_DEFUN([AX_CXX_AUTO_THIS], [dnl
  m4_if([$1], [], [ax_cxx_auto_this_required=true],
        [$1], [mandatory], [ax_cxx_auto_this_required=true],
        [$1], [optional], [ax_cxx_auto_this_required=false],
        [m4_fatal([invalid second argument `$1' to AX_CXX_AUTO_THIS])])dnl
  AC_LANG_PUSH([C++])dnl
  ac_success=no
  AC_CACHE_CHECK([whether $CXX supports C++11 -> decltype(this->f())],
  ax_cv_cxx_auto_this,
  [AC_COMPILE_IFELSE([AC_LANG_SOURCE([_AX_CXX_AUTO_THIS_testbody])],
    [ax_cv_cxx_auto_this=yes],
    [ax_cv_cxx_auto_this=no])])
  if test x$ax_cv_cxx_auto_this = xyes; then
    ac_success=yes
  fi

  AC_LANG_POP([C++])
  if test x$ax_cxx_auto_this_required = xtrue; then
    if test x$ac_success = xno; then
      AC_MSG_ERROR([*** A compiler supporting trailing decltype(this->f()) is required.])
    fi
  else
    if test x$ac_success = xno; then
      HAVE_AUTO_THIS=0
      AC_MSG_NOTICE([No compiler supporting trailing decltype(this->f()) was found])
    else
      HAVE_AUTO_THIS=1
      AC_DEFINE(HAVE_AUTO_THIS,1,
                [define if the compiler supports trailing decltype(this->f())])
    fi

    AC_SUBST(HAVE_AUTO_THIS)
  fi
])
