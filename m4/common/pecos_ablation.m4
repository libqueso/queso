# SYNOPSIS
#
#   Test for PECOS 1D Ablation library
#
#   AM_PATH_ABLATION([, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
#
# DESCRIPTION
#
#   Provides a --with-ablation=DIR option. Searches --with-ablation,
#   $ABLATION_DIR, and the usual places for ablation headers and libraries.
#
#   On success, sets ABLATION_CFLAGS, ABLATION_FCFLAGS, ABLATION_LIBS, and
#   #defines HAVE_ABLATION.  When ACTION-IF-NOT-FOUND is not specified,
#   the default behavior is for configure to fail.
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

AC_DEFUN([AM_PATH_ABLATION],
[

AC_ARG_VAR(ABLATION_DIR,[root directory of PECOS ablation installation])

AC_ARG_WITH(ablation, 
  [AS_HELP_STRING([--with-ablation[=DIR]],[root directory of PECOS ablation installation (default = ABLATION_DIR)])],
  [with_ablation=$withval
if test "${with_ablation}" != yes; then
    ABLATION_PREFIX=$withval
fi
],[
# assume a sensible default of --with-ablation=yes
with_ablation=yes
if test "x${ABLATION_DIR}" != "x"; then
   ABLATION_PREFIX=${ABLATION_DIR}
fi
])

if test "${with_ablation}" != no -a "x${ABLATION_PREFIX}" != "x" ; then
    if test -d "${ABLATION_PREFIX}/conf" ; then
       ABLATION_LIBS=`{ echo "include ${ABLATION_PREFIX}/conf/Makefile.export.ablation"; echo "echo:"; echo "	@echo \\$(PECOS_ABLATION_LIBS)"; } | make -s -f -`
    elif test -d "${ABLATION_PREFIX}/lib" ; then
       # This only works for ablation libraries without PETSc enabled, and requires AM_PATH_GRVY
       ABLATION_LIBS="-L${ABLATION_PREFIX}/lib -lablation -Wl,-rpath,${ABLATION_PREFIX}/lib ${GRVY_LIBS}"
    fi

    if test -d "${ABLATION_PREFIX}/include" ; then
        ABLATION_CFLAGS="-I${ABLATION_PREFIX}/include"
    fi

    ac_ABLATION_save_CFLAGS="$CFLAGS"
    ac_ABLATION_save_CPPFLAGS="$CPPFLAGS"
    ac_ABLATION_save_LDFLAGS="$LDFLAGS"
    ac_ABLATION_save_LIBS="$LIBS"

    CFLAGS="${ABLATION_CFLAGS} ${CFLAGS}"
    CPPFLAGS="${ABLATION_CFLAGS} ${CPPFLAGS}"
    LDFLAGS="${LDFLAGS}"
    AC_LANG_PUSH([C])
    AC_CHECK_HEADER([ablation1d.h],[found_header=yes],[found_header=no])
    AC_CHECK_LIB([ablation],[ablationC1d_version],[found_library=yes],[found_library=no], [ ${ABLATION_LIBS} ])
    AC_LANG_POP([C])

    CFLAGS="$ac_ABLATION_save_CFLAGS"
    CPPFLAGS="$ac_ABLATION_save_CPPFLAGS"
    LDFLAGS="$ac_ABLATION_save_LDFLAGS"
    LIBS="$ac_ABLATION_save_LIBS"

    succeeded=no
    if test "$found_header" = yes; then
        if test "$found_library" = yes; then
            succeeded=yes
        fi
    fi

    if test "$succeeded" = no; then
        ENABLE_ABLATION=0
        ifelse([$2],,AC_MSG_ERROR([libablation not found.  Try either --with-ablation or setting ABLATION_DIR.]),
            [$2])
    else
        ENABLE_ABLATION=1
        AC_DEFINE(HAVE_ABLATION,1,[Define if PECOS libablation is available])
        AC_SUBST(ABLATION_CFLAGS)
        AC_SUBST(ABLATION_LIBS)
	AC_SUBST(ABLATION_PREFIX)
        ifelse([$1],,,[$1])
    fi
else
    ENABLE_ABLATION=0
    AC_MSG_NOTICE([Disabling optional ablation library])
fi

])
