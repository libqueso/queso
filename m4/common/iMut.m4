# SYNOPSIS
#
#   Test for Mutation C++ interface
#
#   AM_PATH_IMUT( <Minimum Required Version>, <package-required=yes/no> )
#
# DESCRIPTION
#
#   Provides a --with-imut=DIR option. Searches --with-imut,
#   $IMUT_DIR, and the usual places for IMUT headers and libraries.
#
#   On success, sets IMUT_CFLAGS, IMUT_LIBS, and #defines HAVE_IMUT.
#   Also defines automake conditional IMUT_ENABLED.  Assumes package
#   is optional unless overridden with $2=yes.
#
# LAST MODIFICATION
#
#   $Id: iMut.m4 17817 2011-02-18 17:13:16Z karl $
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

AC_DEFUN([AX_PATH_IMUT],
[

AC_ARG_VAR(IMUT_DIR,[root directory of IMUT installation])

AC_ARG_WITH(imut,
  [AS_HELP_STRING([--with-imut[=DIR]],[root directory of IMUT installation (default = IMUT_DIR)])],
  [with_imut=$withval
if test "${with_imut}" != yes; then
    IMUT_PREFIX=$withval
fi
],[
with_imut=$withval
if test "x${IMUT_DIR}" != "x"; then
   IMUT_PREFIX=${IMUT_DIR}
fi
])

# package requirement; if not specified, the default is to assume that
# the package is optional

is_package_required=ifelse([$2], ,no, $2 )

HAVE_IMUT=0

    if test -d "${IMUT_PREFIX}/lib" ; then
       IMUT_LIBS="-L${IMUT_PREFIX}/lib -limut -Wl,-rpath,${IMUT_PREFIX}/lib"
       IMUT_FCFLAGS="-I${IMUT_PREFIX}/lib"
    fi

    if test -d "${IMUT_PREFIX}/include" ; then
        IMUT_CFLAGS="-I${IMUT_PREFIX}/include"
    fi

    ac_IMUT_save_CFLAGS="$CFLAGS"
    ac_IMUT_save_CPPFLAGS="$CPPFLAGS"
    ac_IMUT_save_LDFLAGS="$LDFLAGS"
    ac_IMUT_save_LIBS="$LIBS"

    CFLAGS="${IMUT_CFLAGS} ${CFLAGS}"
    CPPFLAGS="${IMUT_CFLAGS} ${CPPFLAGS}"
    LDFLAGS="${IMUT_LIBS} ${LDFLAGS}"
    AC_LANG_PUSH([C])
    AC_CHECK_HEADER([imut.h],[found_header=yes],[found_header=no])

    #-----------------------
    # Minimum version check
    #----------------------

    min_imut_version=ifelse([$1], ,0.29, $1)

    # looking for major.minor.micro style versioning

    MAJOR_VER=`echo $min_imut_version | sed 's/^\([[0-9]]*\).*/\1/'`
    if test "x${MAJOR_VER}" = "x" ; then
       MAJOR_VER=0
    fi

    MINOR_VER=`echo $min_imut_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\2/'`
    if test "x${MINOR_VER}" = "x" ; then
       MINOR_VER=0
    fi

    MICRO_VER=`echo $min_imut_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\3/'`
    if test "x${MICRO_VER}" = "x" ; then
       MICRO_VER=0
    fi

    # begin additional test(s) if header if available

    if test "x${found_header}" = "xyes" ; then

        AC_MSG_CHECKING(for imut - version >= $min_imut_version)
        version_succeeded=no

        AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
        @%:@include <imut.h>
            ]], [[
            #if IMUT_MAJOR_VERSION > $MAJOR_VER
            /* Sweet nibblets */
            #elif (IMUT_MAJOR_VERSION >= $MAJOR_VER) && (IMUT_MINOR_VERSION > $MINOR_VER)
            /* Winner winner, chicken dinner */
            #elif (IMUT_MAJOR_VERSION >= $MAJOR_VER) && (IMUT_MINOR_VERSION >= $MINOR_VER) && (IMUT_MICRO_VERSION >= $MICRO_VER)
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

    if test "$version_succeeded" != "yes";then
       if test "$is_package_required" = yes; then
          AC_MSG_ERROR([

   Your IMUT library version does not meet the minimum versioning
   requirements ($min_imut_version).  Please use --with-imut to specify the location
   of an updated installation or consider upgrading the system version.

          ])
       fi
    fi

    # Library availability

dnl    AC_MSG_CHECKING([for -limut linkage])

dnl    AC_CHECK_LIB([imut],imut_input_fopen,[found_library=yes],[found_library=no])

    fi   dnl end test if header if available

    AC_LANG_POP([C])

    CFLAGS="$ac_IMUT_save_CFLAGS"
    CPPFLAGS="$ac_IMUT_save_CPPFLAGS"
    LDFLAGS="$ac_IMUT_save_LDFLAGS"
    LIBS="$ac_IMUT_save_LIBS"

    succeeded=no
    if test "$found_header" = yes; then
        if test "$version_succeeded" = yes; then
dnl           if test "$found_library" = yes; then
              succeeded=yes
dnl           fi
        fi
    fi

    if test "$succeeded" = no; then
       if test "$is_package_required" = yes; then
          AC_MSG_ERROR([libIMUT not found.  Try either --with-imut or setting IMUT_DIR.])
       else
          AC_MSG_NOTICE([optional IMUT library not found])
          IMUT_CFLAGS=""   # IMUT_CFLAGS empty on failure
          IMUT_LIBS=""     # IMUT_LIBS empty on failure
          IMUT_FLIBS=""    # IMUT_FLIBS empty on failure
       fi
    else
        HAVE_IMUT=1
        AC_DEFINE(HAVE_IMUT,1,[Define if IMUT is available])
        AC_SUBST(IMUT_CFLAGS)
        AC_SUBST(IMUT_LIBS)
        AC_SUBST(IMUT_PREFIX)
        AC_SUBST(IMUT_FLIBS)
    fi

    AC_SUBST(HAVE_IMUT)

# fi

AM_CONDITIONAL(IMUT_ENABLED,test x$HAVE_IMUT = x1)

])
