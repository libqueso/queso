# SYNOPSIS
#
#   Test for Antioch
#
#   AX_PATH_ANTIOCH( <Minimum Required Version>, <package-required=yes/no> )
#
# DESCRIPTION
#
#   Provides a --with-antioch=DIR option. Searches --with-antioch,
#   $ANTIOCH_DIR, and the usual places for ANTIOCH headers and libraries.
#
#   On success, sets ANTIOCH_CPPFLAGS, ANTIOCH_LIBS, and #defines HAVE_ANTIOCH.
#   Also defines automake conditional ANTIOCH_ENABLED.  Assumes package
#   is optional unless overridden with $2=yes.
#
# COPYLEFT
#
#   Copyright (c) 2013 Roy H. Stogner <roystgnr@ices.utexas.edu>
#   Copyright (c) 2013 Paul T. Bauman <pbauman@ices.utexas.edu>
#   Copyright (c) 2010 Karl W. Schulz <karl@ices.utexas.edu>
#   Copyright (c) 2009 Rhys Ulerich <rhys.ulerich@gmail.com>
#   Copyright (c) 2008 Thomas Porschberg <thomas@randspringer.de>
#   Copyright (c) 2008 Caolan McNamara <caolan@skynet.ie>
#   Copyright (c) 2008 Alexandre Duret-Lutz <adl@gnu.org>
#   Copyright (c) 2008 Matthew Mueller <donut@azstarnet.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_PATH_ANTIOCH],
[

AC_ARG_VAR(ANTIOCH_DIR,[root directory of Antioch installation])

AC_ARG_WITH(antioch,
  [AS_HELP_STRING([--with-antioch[=DIR]],[root directory of Antioch installation (default = ANTIOCH_DIR)])],
  [with_antioch=$withval
if test "${with_antioch}" != yes; then
    ANTIOCH_PREFIX=$withval
fi
],[
with_antioch=$withval
if test "x${ANTIOCH_DIR}" != "x"; then
   ANTIOCH_PREFIX=${ANTIOCH_DIR}
fi
])

# package requirement; if not specified, the default is to assume that
# the package is optional

is_package_required=ifelse([$2], ,no, $2 )

HAVE_ANTIOCH=0

# logic change: if the user called the macro, check for package,
# decide what to do based on whether the package is required or not.

# if test "${with_antioch}" != no ; then

    if test -d "${ANTIOCH_PREFIX}/lib" ; then
       ANTIOCH_LDFLAGS="-L${ANTIOCH_PREFIX}/lib -Wl,-rpath,${ANTIOCH_PREFIX}/lib"
       ANTIOCH_LIBS="-lantioch"
    fi

    if test -d "${ANTIOCH_PREFIX}/include" ; then
       ANTIOCH_CPPFLAGS="-I${ANTIOCH_PREFIX}/include -I${ANTIOCH_PREFIX}/src"
    fi

    ac_ANTIOCH_save_CPPFLAGS="$CPPFLAGS"
    ac_ANTIOCH_save_LDFLAGS="$LDFLAGS"
    ac_ANTIOCH_save_LIBS="$LIBS"

    CPPFLAGS="${ANTIOCH_CPPFLAGS} ${CPPFLAGS}"
    LDFLAGS="${ANTIOCH_LDFLAGS} ${LDFLAGS}"
    LIBS="${ANTIOCH_LIBS} ${LIBS}"

    AC_LANG_PUSH([C++])
    AC_CHECK_HEADER([antioch/antioch_version.h],[found_header=yes],[found_header=no])
    AC_LANG_POP([C++])

    #-----------------------
    # Minimum version check
    #----------------------

    min_antioch_version=ifelse([$1], ,0.0.0, $1)

    # looking for major.minor.micro style versioning

    MAJOR_VER=`echo $min_antioch_version | sed 's/^\([[0-9]]*\).*/\1/'`
    if test "x${MAJOR_VER}" = "x" ; then
       MAJOR_VER=0
    fi

    MINOR_VER=`echo $min_antioch_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\2/'`
    if test "x${MINOR_VER}" = "x" ; then
       MINOR_VER=0
    fi

    MICRO_VER=`echo $min_antioch_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\3/'`
    if test "x${MICRO_VER}" = "x" ; then
       MICRO_VER=0
    fi

    # begin additional test(s) if header if available

    if test "x${found_header}" = "xyes" ; then

        AC_MSG_CHECKING(for antioch - version >= $min_antioch_version)
        version_succeeded=no

	AC_LANG_PUSH([C++])
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
        @%:@include "antioch/antioch_version.h"
            ]], [[
            #if ANTIOCH_MAJOR_VERSION > $MAJOR_VER
            /* Sweet nibblets */
            #elif (ANTIOCH_MAJOR_VERSION >= $MAJOR_VER) && (ANTIOCH_MINOR_VERSION > $MINOR_VER)
            /* Winner winner, chicken dinner */
            #elif (ANTIOCH_MAJOR_VERSION >= $MAJOR_VER) && (ANTIOCH_MINOR_VERSION >= $MINOR_VER) && (ANTIOCH_MICRO_VERSION >= $MICRO_VER)
            /* I feel like chicken tonight, like chicken tonight? */
            #else
            #  error version is too old
            #endif
        ]])],[
            AC_MSG_RESULT(yes)
            version_succeeded=yes
        ],[
            AC_MSG_RESULT(no)
        ])
	AC_LANG_POP([C++])

    if test "$version_succeeded" != "yes";then
       if test "$is_package_required" = yes; then
          AC_MSG_ERROR([

   Your Antioch version does not meet the minimum versioning
   requirements ($min_antioch_version).  Please use --with-antioch to specify the location
   of an updated installation or consider upgrading the system version.

          ])
       fi
    fi

    # Library availability

    AC_MSG_CHECKING([for -lantioch linkage])

    AC_LANG_PUSH([C++])

    AC_LINK_IFELSE(
                  [AC_LANG_PROGRAM([#include "antioch/antioch_version.h"],
                                   [Antioch::get_antioch_version()])],
                  [AC_MSG_RESULT(yes)
                   found_library=yes],
                  [AC_MSG_RESULT(no) 
                   found_library=no])

    fi   dnl end test if header if available

    AC_LANG_POP([C++])

    CPPFLAGS="$ac_ANTIOCH_save_CPPFLAGS"
    LDFLAGS="$ac_ANTIOCH_save_LDFLAGS"
    LIBS="$ac_ANTIOCH_save_LIBS"

    succeeded=no
    if test "$found_header" = yes; then
        if test "$version_succeeded" = yes; then
           if test "$found_library" = yes; then
              succeeded=yes
           fi
        fi
    fi

    if test "$succeeded" = no; then
       if test "$is_package_required" = yes; then
          AC_MSG_ERROR([Antioch not found.  Try either --with-antioch or setting ANTIOCH_DIR.])
       else
          AC_MSG_NOTICE([optional Antioch library not found])
          ANTIOCH_CPPFLAGS=""   # ANTIOCH_CFLAGS empty on failure
          ANTIOCH_LDFLAGS=""    # ANTIOCH_LDFLAGS empty on failure
          ANTIOCH_LIBS=""       # ANTIOCH_LIBS empty on failure
       fi
    else
        HAVE_ANTIOCH=1
        AC_DEFINE(HAVE_ANTIOCH,1,[Define if Antioch is available])
        AC_SUBST(ANTIOCH_CPPFLAGS)
        AC_SUBST(ANTIOCH_LDFLAGS)
        AC_SUBST(ANTIOCH_LIBS)
        AC_SUBST(ANTIOCH_PREFIX)
    fi

    AC_SUBST(HAVE_ANTIOCH)

# fi

AM_CONDITIONAL(ANTIOCH_ENABLED,test x$HAVE_ANTIOCH = x1)

])
