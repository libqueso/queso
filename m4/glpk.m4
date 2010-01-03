# SYNOPSIS
#
#   Test for HPCT 
#
#   AM_PATH_HPCT([, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
#
# DESCRIPTION
#
#   Provides a --with-hpct=DIR option. Searches --with-hpct,
#   $HPCT_DIR, and the usual places for HPCT headers and libraries.
#
#   On success, sets HPCT_CFLAGS, HPCT_LIBS, and
#   #defines HAVE_HPCT.  When ACTION-IF-NOT-FOUND is not specified,
#   the default behavior is for configure to fail.
#
# LAST MODIFICATION
#
#   2009-07-14
#
# COPYLEFT
#
#   Copyright (c) 2009 Rhys Ulerich <rhys.ulerich@gmail.com>
#   Copyright (c) 2008 Thomas Porschberg <thomas@randspringer.de>
#   Copyright (c) 2008 Caolan McNamara <caolan@skynet.ie>
#   Copyright (c) 2008 Alexandre Duret-Lutz <adl@gnu.org>
#   Copyright (c) 2008 Matthew Mueller <donut@azstarnet.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_PATH_GLPK],
[

AC_ARG_VAR(GLPK_HOME,[root directory of GLPK installation])

AC_ARG_WITH(glpk, 
  [AS_HELP_STRING([--with-glpk[=DIR]],[root directory of GLPK installation (default = GLPK_DIR)])],
  [with_glpk=$withval
if test "${with_glpk}" != yes; then
    GLPK_PREFIX=$withval
fi
],[
with_glpk=$withval
if test "x${GLPK_DIR}" != "x"; then
   GLPK_PREFIX=${GLPK_DIR}
fi
])

if test "${with_glpk}" != no ; then

    if test -d "${GLPK_PREFIX}/lib" ; then
       GLPK_LIBS="-L${GLPK_PREFIX}/lib -lglpk_for_queso"
    fi

    if test -d "${GLPK_PREFIX}/include" ; then
        GLPK_CFLAGS="-I${GLPK_PREFIX}/include"
    fi

    ac_GLPK_save_CFLAGS="$CFLAGS"
    ac_GLPK_save_CPPFLAGS="$CPPFLAGS"
    ac_GLPK_save_LDFLAGS="$LDFLAGS"
    ac_GLPK_save_LIBS="$LIBS"

    CFLAGS="${GLPK_CFLAGS} ${CFLAGS}"
    CPPFLAGS="${GLPK_CFLAGS} ${CPPFLAGS}"
    LDFLAGS="${GLPK_LIBS} ${LDFLAGS}"
    AC_LANG_PUSH([C])
    AC_CHECK_HEADER([glpk.h],[found_header=yes],[found_header=no])
#    AC_CHECK_LIB([glpk],[found_library=yes],[found_library=no])
    AC_LANG_POP([C])

    CFLAGS="$ac_GLPK_save_CFLAGS"
    CPPFLAGS="$ac_GLPK_save_CPPFLAGS"
    LDFLAGS="$ac_GLPK_save_LDFLAGS"
    LIBS="$ac_GLPK_save_LIBS"

    succeeded=no
    if test "$found_header" = yes; then
#        if test "$found_library" = yes; then
            succeeded=yes
#        fi
    fi

    if test "$succeeded" = no; then
        ifelse([$2],,AC_MSG_ERROR([GLPK not found.  Try either --with-glpk or setting GLPK_DIR.]),
            [$2])
    else
        AC_DEFINE(HAVE_GLPK,1,[Define if GLPK is available])
        AC_SUBST(GLPK_CFLAGS)
        AC_SUBST(GLPK_LIBS)
	AC_SUBST(GLPK_PREFIX)
        ifelse([$1],,,[$1])
    fi

fi

])
