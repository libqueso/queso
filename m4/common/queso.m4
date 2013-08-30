# SYNOPSIS
#
#   Test for QUESO 
#
#   AM_PATH_QUESO([, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
#
# DESCRIPTION
#
#   Provides a --with-queso=DIR option. Searches --with-queso,
#   $QUESO_DIR, and the usual places for QUESO headers and libraries.
#
#   On success, sets QUESO_CFLAGS, QUESO_LIBS, and
#   #defines HAVE_QUESO.  When ACTION-IF-NOT-FOUND is not specified,
#   the default behavior is for configure to fail.
#
# LAST MODIFICATION
#
#   2010-03-22 by Paul T. Bauman
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

AC_DEFUN([AX_PATH_QUESO],
[

AC_ARG_VAR(QUESO_DIR,[root directory of QUESO installation])

AC_ARG_WITH(queso, 
  [AS_HELP_STRING([--with-queso[=DIR]],[root directory of QUESO installation (default = QUESO_DIR)])],
  [with_queso=$withval
if test "${with_queso}" != yes; then
    QUESO_PREFIX=$withval
fi
],[
with_queso=$withval
if test "x${QUESO_DIR}" != "x"; then
   QUESO_PREFIX=${QUESO_DIR}
fi
])

if test "${with_queso}" != no ; then

    if test -d "${QUESO_PREFIX}/lib" ; then
       QUESO_LIBS="-L${QUESO_PREFIX}/lib -lqueso"
    fi

    if test -d "${QUESO_PREFIX}/include" ; then
        dnl FIXME: CFLAGS should be reserved for C-compiler flags
        dnl FIXME: conflicts with use of CFLAGS below.
        QUESO_CFLAGS="-I${QUESO_PREFIX}/include"

	dnl PB: Added QUESO_CPPFLAGS
        QUESO_CPPFLAGS="-I${QUESO_PREFIX}/include"
    fi

    ac_QUESO_save_CFLAGS="$CFLAGS"
    ac_QUESO_save_CPPFLAGS="$CPPFLAGS"
    ac_QUESO_save_LDFLAGS="$LDFLAGS"
    ac_QUESO_save_LIBS="$LIBS"

    CFLAGS="${QUESO_CFLAGS} ${CFLAGS}"
    CPPFLAGS="${QUESO_CPPFLAGS} ${CPPFLAGS}"
    LDFLAGS="${QUESO_LIBS} ${LDFLAGS}"

    AC_LANG_PUSH([C++])

    
    dnl AC_CHECK_HEADER([uqDefines.h], [found_header=yes], [found_header=no])

    dnl PB: Replaced the previous CHECK_HEADER test with COMPILE_IFELSE test becaue
    dnl PB: uqDefines.h depends on MPI and the user will always get a CPP warning
    dnl PB: since the c-preprocessor test will not use MPI.
    dnl PB: However, I'm still not satisfied because no message is emitted and
    dnl PB: if somehow the library was present, but not the includes, the user
    dnl PB: will be a little lost.
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([#include "uqDefines.h"])],[found_header=yes],[found_header=no])

    dnl PB: Check below is "newer" over deprecated AC_HAVE_LIBRARY. 
    AC_CHECK_LIB([queso], [main], [found_library=yes], [found_library=no],[ ${GSL_LIBS} 
                                                                            ${EPETRA_LDFLAGS} -lepetra 
                                                                            ${GRVY_LIBS} 
                                                                            ${BOOST_PROGRAM_OPTIONS_LIBS}
                                                                            ${BOOST_PROGRAM_OPTIONS_LDFLAGS}
									    ${GLPK_LIBS}
									    ${HDF5_LIBS} ])

    dnl PB: LDFLAGS and LIBS do not expand as I would have expected so using above test instead to resolve
    dnl PB: library dependencies.
    dnl FIXME: Question: Why do we even have to resolve those at all???
    dnl AC_CHECK_LIB([queso], [main], [found_library=yes], [found_library=no],[ ${LDFLAGS} ${LIBS} ] )									   

    dnl PB: http://nerdland.net/2009/07/detecting-c-libraries-with-autotools
    dnl PB: Gives a good guide to rolling up a test for C++. Will take a peek in the future.
    dnl AC_CHECK_LIB([queso], [uqBaseEnvironmentClass::worldRank], [found_library=yes], [found_library=no])

    AC_LANG_POP([C++])

    CFLAGS="$ac_QUESO_save_CFLAGS"
    CPPFLAGS="$ac_QUESO_save_CPPFLAGS"
    LDFLAGS="$ac_QUESO_save_LDFLAGS"
    LIBS="$ac_QUESO_save_LIBS"

    succeeded=no
    if test "$found_header" = yes; then
        if test "$found_library" = yes; then
            succeeded=yes
        fi
    fi

    if test "$succeeded" = no; then
        ifelse([$2],,AC_MSG_ERROR([QUESO not found.  Try either --with-queso or setting QUESO_DIR.]),
            [$2])
    else
        AC_DEFINE(HAVE_QUESO,1,[Define if QUESO is available])
        AC_SUBST(QUESO_CFLAGS)
        AC_SUBST(QUESO_CPPFLAGS)
        AC_SUBST(QUESO_LIBS)
	AC_SUBST(QUESO_PREFIX)
        ifelse([$1],,,[$1])
    fi

fi

])
