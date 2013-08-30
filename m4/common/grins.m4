# SYNOPSIS
#
#   Test for GRINS
#
#   AX_PATH_GRINS( <Minimum Required Version>, <package-required=yes/no> )
#
# DESCRIPTION
#
#   Provides a --with-grins=DIR option. Searches --with-grins,
#   $GRINS_DIR, and the usual places for GRINS headers and libraries.
#
#   On success, sets GRINS_CPPFLAGS, GRINS_LIBS, and #defines HAVE_GRINS.
#   Also defines automake conditional GRINS_ENABLED.  Assumes package
#   is optional unless overridden with $2=yes.
#
# LAST MODIFICATION
#
#   $Id: grins.m4 -1   $
#
# COPYLEFT
#
#   Copyright (c) 2012 Paul T. Bauman <pbauman@ices.utexas.edu>
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

AC_DEFUN([AX_PATH_GRINS],
[

AC_ARG_VAR(GRINS_DIR,[root directory of GRINS installation])

AC_ARG_WITH(grins,
  [AS_HELP_STRING([--with-grins[=DIR]],[root directory of GRINS installation (default = GRINS_DIR)])],
  [with_grins=$withval
if test "${with_grins}" != yes; then
    GRINS_PREFIX=$withval
fi
],[
with_grins=$withval
if test "x${GRINS_DIR}" != "x"; then
   GRINS_PREFIX=${GRINS_DIR}
fi
])

# package requirement; if not specified, the default is to assume that
# the package is optional

is_package_required=ifelse([$2], ,no, $2 )

HAVE_GRINS=0

# logic change: if the user called the macro, check for package,
# decide what to do based on whether the package is required or not.

# if test "${with_grins}" != no ; then

    AC_REQUIRE([AX_PATH_LIBMESH_NEW(0.9.1, yes)])

    if test -d "${GRINS_PREFIX}/lib" ; then
       GRINS_LDFLAGS="-L${GRINS_PREFIX}/lib -Wl,-rpath,${GRINS_PREFIX}/lib ${LIBMESH_LDFLAGS}"
       GRINS_LIBS="-lgrins ${LIBMESH_LIBS}"
    fi

    if test -d "${GRINS_PREFIX}/include" ; then
       GRINS_CPPFLAGS="-I${GRINS_PREFIX}/include -I${GRINS_PREFIX}/src ${LIBMESH_CPPFLAGS}"
    fi

    ac_GRINS_save_CPPFLAGS="$CPPFLAGS"
    ac_GRINS_save_LDFLAGS="$LDFLAGS"
    ac_GRINS_save_LIBS="$LIBS"

    CPPFLAGS="${GRINS_CPPFLAGS} ${CPPFLAGS}"
    LDFLAGS="${GRINS_LDFLAGS} ${LDFLAGS}"
    LIBS="${GRINS_LIBS} ${LIBS}"

    AC_LANG_PUSH([C++])
    AC_CHECK_HEADER([grins/grins_version.h],[found_header=yes],[found_header=no])
    AC_LANG_POP([C++])

    #-----------------------
    # Minimum version check
    #----------------------

    min_grins_version=ifelse([$1], ,0.2.0, $1)

    # looking for major.minor.micro style versioning

    MAJOR_VER=`echo $min_grins_version | sed 's/^\([[0-9]]*\).*/\1/'`
    if test "x${MAJOR_VER}" = "x" ; then
       MAJOR_VER=0
    fi

    MINOR_VER=`echo $min_grins_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\2/'`
    if test "x${MINOR_VER}" = "x" ; then
       MINOR_VER=0
    fi

    MICRO_VER=`echo $min_grins_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\3/'`
    if test "x${MICRO_VER}" = "x" ; then
       MICRO_VER=0
    fi

    # begin additional test(s) if header if available

    if test "x${found_header}" = "xyes" ; then

        AC_MSG_CHECKING(for grins - version >= $min_grins_version)
        version_succeeded=no

	AC_LANG_PUSH([C++])
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
        @%:@include "grins/grins_version.h"
            ]], [[
            #if GRINS_MAJOR_VERSION > $MAJOR_VER
            /* Sweet nibblets */
            #elif (GRINS_MAJOR_VERSION >= $MAJOR_VER) && (GRINS_MINOR_VERSION > $MINOR_VER)
            /* Winner winner, chicken dinner */
            #elif (GRINS_MAJOR_VERSION >= $MAJOR_VER) && (GRINS_MINOR_VERSION >= $MINOR_VER) && (GRINS_MICRO_VERSION >= $MICRO_VER)
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

   Your GRINS version does not meet the minimum versioning
   requirements ($min_grins_version).  Please use --with-grins to specify the location
   of an updated installation or consider upgrading the system version.

          ])
       fi
    fi

    # Library availability

    AC_MSG_CHECKING([for -lgrins linkage])

    AC_LANG_PUSH([C++])

    AC_LINK_IFELSE(
                  [AC_LANG_PROGRAM([#include "grins/grins_version.h"],
                                   [GRINS::get_grins_version()])],
                  [AC_MSG_RESULT(yes)
                   found_library=yes],
                  [AC_MSG_RESULT(no) 
                   found_library=no])

    fi   dnl end test if header if available

    AC_LANG_POP([C++])

    CPPFLAGS="$ac_GRINS_save_CPPFLAGS"
    LDFLAGS="$ac_GRINS_save_LDFLAGS"
    LIBS="$ac_GRINS_save_LIBS"

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
          AC_MSG_ERROR([GRINS not found.  Try either --with-grins or setting GRINS_DIR.])
       else
          AC_MSG_NOTICE([optional GRINS library not found])
          GRINS_CPPFLAGS="" # GRINS_CFLAGS empty on failure
          GRINS_LDFLAGS=""  # GRINS_LDFLAGS empty on failure
          GRINS_LIBS=""     # GRINS_LIBS empty on failure
       fi
    else
        HAVE_GRINS=1
        AC_DEFINE(HAVE_GRINS,1,[Define if GRINS is available])
        AC_SUBST(GRINS_CPPFLAGS)
        AC_SUBST(GRINS_LDFLAGS)
        AC_SUBST(GRINS_LIBS)
        AC_SUBST(GRINS_PREFIX)
    fi

    AC_SUBST(HAVE_GRINS)

# fi

AM_CONDITIONAL(GRINS_ENABLED,test x$HAVE_GRINS = x1)

])
