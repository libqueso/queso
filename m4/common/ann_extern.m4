# SYNOPSIS
#
#   Test for  
#
#   AM_PATH_ANN([, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
#
# DESCRIPTION
#
#   Provides a --with-ann=DIR option. Searches --with-ann,
#   $ANN_DIR, and the usual places for ANN headers and libraries.
#
#   On success, sets ANN_CFLAGS, ANN_LIBS, and
#   #defines HAVE_ANN.  When ACTION-IF-NOT-FOUND is not specified,
#   the default behavior is for configure to fail.
#
# LAST MODIFICATION
#
#   2011-02-08 by Gabriel Terejanu
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

AC_DEFUN([AX_PATH_ANN],
[

AC_ARG_VAR(ANN_DIR,[root directory of ANN installation])

AC_ARG_WITH(ann, 
  [AS_HELP_STRING([--with-ann[=DIR]],[root directory of ANN installation (default = ANN_DIR)])],
  [with_ann=$withval
if test "${with_ann}" != yes; then
    ANN_PREFIX=$withval
fi
],[
with_ann=$withval
if test "x${ANN_DIR}" != "x"; then
   ANN_PREFIX=${ANN_DIR}
fi
])

if test "${with_ann}" != no ; then

    if test -d "${ANN_PREFIX}/lib" ; then
       ANN_LIBS="-L${ANN_PREFIX}/lib -lANN"
    fi

    if test -d "${ANN_PREFIX}/include" ; then
        dnl FIXME: CFLAGS should be reserved for C-compiler flags
        dnl FIXME: conflicts with use of CFLAGS below.
        ANN_CFLAGS="-I${ANN_PREFIX}/include"

	dnl PB: Added ANN_CPPFLAGS
        ANN_CPPFLAGS="-I${ANN_PREFIX}/include"
    fi

    ac_ANN_save_CFLAGS="$CFLAGS"
    ac_ANN_save_CPPFLAGS="$CPPFLAGS"
    ac_ANN_save_LDFLAGS="$LDFLAGS"
    ac_ANN_save_LIBS="$LIBS"

    CFLAGS="${ANN_CFLAGS} ${CFLAGS}"
    CPPFLAGS="${ANN_CPPFLAGS} ${CPPFLAGS}"
    LDFLAGS="${ANN_LIBS} ${LDFLAGS}"


    AC_LANG_PUSH([C++])

    # Header check
    
    AC_CHECK_HEADER([ANN/ANN.h], [found_header=yes], [found_header=no])

    # Library availability

    AC_MSG_CHECKING([for -lANN linkage])

    AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([#include <ANN/ANN.h>],[annClose])],
    [TEST_LIBS="$TEST_LIBS -lANN"] [
    AC_MSG_RESULT(yes)
    found_library=yes ],[AC_MSG_RESULT(no)])

    AC_LANG_POP([C++])


    CFLAGS="$ac_ANN_save_CFLAGS"
    CPPFLAGS="$ac_ANN_save_CPPFLAGS"
    LDFLAGS="$ac_ANN_save_LDFLAGS"
    LIBS="$ac_ANN_save_LIBS"

    succeeded=no
    if test "$found_header" = yes; then
       if test "$found_library" = yes; then
            succeeded=yes
       fi
    fi

    if test "$succeeded" = yes; then
        AC_DEFINE(HAVE_ANN,1,[Define if ANN is available])
        AC_SUBST(ANN_CFLAGS)
        AC_SUBST(ANN_CPPFLAGS)
        AC_SUBST(ANN_LIBS)
	AC_SUBST(ANN_PREFIX)
        ifelse([$1],,,[$1])
    else
	ANN_CFLAGS=""
	ANN_CPPFLAGS=""
	ANN_LIBS=""
	ANN_PREFIX=""
        AC_SUBST(ANN_CFLAGS)
        AC_SUBST(ANN_CPPFLAGS)
        AC_SUBST(ANN_LIBS)
	AC_SUBST(ANN_PREFIX)
    fi

fi

])
