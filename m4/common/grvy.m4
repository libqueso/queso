# SYNOPSIS
#
#   Test for libGRVY
#
#   AM_PATH_grvy([, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
#
# DESCRIPTION
#
#   Provides a --with-grvy=DIR option. Searches --with-grvy,
#   $GRVY_DIR, and the usual places for GRVY headers and libraries.
#
#   On success, sets GRVY_CFLAGS, GRVY_LIBS, and
#   #defines HAVE_GRVY.  When ACTION-IF-NOT-FOUND is not specified,
#   the default behavior is for configure to fail.
#
# LAST MODIFICATION
#
#   $Id: grvy.m4 13380 2010-09-14 20:23:07Z karl $
#
# COPYLEFT
#
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

AC_DEFUN([AX_PATH_GRVY],
[

AC_ARG_VAR(GRVY_DIR,[root directory of GRVY installation])

AC_ARG_WITH(grvy, 
  [AS_HELP_STRING([--with-grvy[=DIR]],[root directory of GRVY installation (default = GRVY_DIR)])],
  [with_grvy=$withval
if test "${with_grvy}" != yes; then
    GRVY_PREFIX=$withval
fi
],[
with_grvy=$withval
if test "x${GRVY_DIR}" != "x"; then
   GRVY_PREFIX=${GRVY_DIR}
fi
])

if test "${with_grvy}" != no ; then

    if test -d "${GRVY_PREFIX}/lib" ; then
       GRVY_LIBS="-L${GRVY_PREFIX}/lib -lgrvy -Wl,-rpath,${GRVY_PREFIX}/lib"
       GRVY_FCFLAGS="-I${GRVY_PREFIX}/lib"
    fi

    if test -d "${GRVY_PREFIX}/include" ; then
        GRVY_CFLAGS="-I${GRVY_PREFIX}/include"
    fi

    ac_GRVY_save_CFLAGS="$CFLAGS"
    ac_GRVY_save_CPPFLAGS="$CPPFLAGS"
    ac_GRVY_save_LDFLAGS="$LDFLAGS"
    ac_GRVY_save_LIBS="$LIBS"

    CFLAGS="${GRVY_CFLAGS} ${CFLAGS}"
    CPPFLAGS="${GRVY_CFLAGS} ${CPPFLAGS}"
    LDFLAGS="${GRVY_LIBS} ${LDFLAGS}"
    AC_LANG_PUSH([C])
    AC_CHECK_HEADER([grvy.h],[found_header=yes],[found_header=no])
    AC_CHECK_LIB([grvy],grvy_input_fopen,[found_library=yes],[found_library=no])
    AC_LANG_POP([C])

    CFLAGS="$ac_GRVY_save_CFLAGS"
    CPPFLAGS="$ac_GRVY_save_CPPFLAGS"
    LDFLAGS="$ac_GRVY_save_LDFLAGS"
    LIBS="$ac_GRVY_save_LIBS"

    succeeded=no
    if test "$found_header" = yes; then
        if test "$found_library" = yes; then
            succeeded=yes
        fi
    fi

    if test "$succeeded" = no; then
        ifelse([$2],,AC_MSG_ERROR([libGRVY not found.  Try either --with-grvy or setting GRVY_DIR.]),
            [$2])
    else
        AC_DEFINE(HAVE_GRVY,1,[Define if GRVY is available])
        AC_SUBST(GRVY_CFLAGS)
        AC_SUBST(GRVY_FCFLAGS)
        AC_SUBST(GRVY_LIBS)
	AC_SUBST(GRVY_PREFIX)
        ifelse([$1],,,[$1])
    fi

fi

])
