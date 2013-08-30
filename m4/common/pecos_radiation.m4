# SYNOPSIS
#
#   Test for pecos libradiation
#
#   AM_PATH_RADIATION( <Minimum Required Version>, <package-required=yes/no> )
#
# DESCRIPTION
#
#   Provides a --with-radiation=DIR option. Searches --with-radiation,
#   $RADIATION_DIR, and the usual places for libradiation headers and
#   libraries.
#
#   On success, sets RADIATION_CFLAGS, RADIATION_LIBS, and #defines
#   HAVE_RADIATION.
#   Also defines automake conditional RADIATION_ENABLED.  Assumes
#   package is optional unless overridden with $2=yes.
#
# LAST MODIFICATION
#
#   $Id: pecos_radiation.m4 33353 2012-09-25 22:19:37Z roystgnr $
#
# COPYLEFT
#
#   Copyright (c) 2012 Roy H. Stogner <roystgnr@ices.utexas.edu>
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

AC_DEFUN([AX_PATH_RADIATION],
[

AC_ARG_VAR(RADIATION_DIR,[root directory of libradiation installation])

AC_ARG_WITH(radiation,
  [AS_HELP_STRING([--with-radiation[=DIR]],[root directory of libradiation installation (default = RADIATION_DIR)])],
  [with_radiation=$withval
if test "${with_radiation}" != yes; then
    RADIATION_PREFIX=$withval
fi
],[
# assume a sensible default of --with-radiation=yes
with_radiation=yes
if test "x${RADIATION_DIR}" != "x"; then
   RADIATION_PREFIX=${RADIATION_DIR}
fi
])

# package requirement; if not specified, the default is to assume that
# the package is optional

is_package_required=ifelse([$2], ,no, $2 )

HAVE_RADIATION=0

# logic change: if the user called the macro, check for package,
# decide what to do based on whether the package is required or not.

# if test "${with_radiation}" != no ; then

    if test -d "${RADIATION_PREFIX}/lib" ; then
       RADIATION_LIBS="-L${RADIATION_PREFIX}/lib -lradiation -Wl,-rpath,${RADIATION_PREFIX}/lib ${GRVY_FLIBS}"
       RADIATION_FLIBS="$RADIATION_LIBS"
       RADIATION_FCFLAGS="-I${RADIATION_PREFIX}/lib"
    fi

    if test -d "${RADIATION_PREFIX}/include" ; then
        RADIATION_CFLAGS="-I${RADIATION_PREFIX}/include"
    fi

    ac_RADIATION_save_CFLAGS="$CFLAGS"
    ac_RADIATION_save_CPPFLAGS="$CPPFLAGS"
    ac_RADIATION_save_LDFLAGS="$LDFLAGS"
    ac_RADIATION_save_LIBS="$LIBS"

    CFLAGS="${RADIATION_CFLAGS} ${CFLAGS}"
    CPPFLAGS="${RADIATION_CFLAGS} ${CPPFLAGS}"
    LDFLAGS="${RADIATION_LIBS} ${LDFLAGS}"
    AC_LANG_PUSH([C])
    AC_CHECK_HEADER([radiation.h],[found_header=yes],[found_header=no])

    #-----------------------
    # Minimum version check
    #----------------------

    min_radiation_version=ifelse([$1], ,0.20, $1)

    # looking for major.minor.micro style versioning

    MAJOR_VER=`echo $min_radiation_version | sed 's/^\([[0-9]]*\).*/\1/'`
    if test "x${MAJOR_VER}" = "x" ; then
       MAJOR_VER=0
    fi

    MINOR_VER=`echo $min_radiation_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\2/'`
    if test "x${MINOR_VER}" = "x" ; then
       MINOR_VER=0
    fi

    MICRO_VER=`echo $min_radiation_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\3/'`
    if test "x${MICRO_VER}" = "x" ; then
       MICRO_VER=0
    fi

    # begin additional test(s) if header if available

    if test "x${found_header}" = "xyes" ; then

        AC_MSG_CHECKING(for radiation - version >= $min_radiation_version)
        version_succeeded=no

        AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
        @%:@include <radiation.h>
            ]], [[
            #if RADIATION_MAJOR_VERSION > $MAJOR_VER
            /* Sweet nibblets */
            #elif (RADIATION_MAJOR_VERSION >= $MAJOR_VER) && (RADIATION_MINOR_VERSION > $MINOR_VER)
            /* Winner winner, chicken dinner */
            #elif (RADIATION_MAJOR_VERSION >= $MAJOR_VER) && (RADIATION_MINOR_VERSION >= $MINOR_VER) && (RADIATION_MICRO_VERSION >= $MICRO_VER)
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

   Your libradiation version does not meet the minimum versioning
   requirements ($min_radiation_version).  Please use --with-radiation
   to specify the location of an updated installation or consider
   upgrading the system version.

          ])
       fi
    fi

    # Library availability

    AC_MSG_CHECKING([for -lradiation linkage])

    AC_CHECK_LIB([radiation],absorption_coefficient_,[found_library=yes],[found_library=no], [ ${RADIATION_LIBS} ])

    fi   dnl end test if header if available

    AC_LANG_POP([C])

    CFLAGS="$ac_RADIATION_save_CFLAGS"
    CPPFLAGS="$ac_RADIATION_save_CPPFLAGS"
    LDFLAGS="$ac_RADIATION_save_LDFLAGS"
    LIBS="$ac_RADIATION_save_LIBS"

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
          AC_MSG_ERROR([libradiation not found.  Try either --with-radiation or setting RADIATION_DIR.])
       else
          AC_MSG_NOTICE([optional radiation library not found])
          RADIATION_CFLAGS=""   # RADIATION_CFLAGS empty on failure
          RADIATION_FCFLAGS=""  # RADIATION_FCFLAGS empty on failure
          RADIATION_LIBS=""     # RADIATION_LIBS empty on failure
          RADIATION_FLIBS=""    # RADIATION_FLIBS empty on failure
       fi
    else
        HAVE_RADIATION=1
        AC_DEFINE(HAVE_RADIATION,1,[Define if RADIATION is available])
        AC_SUBST(RADIATION_CFLAGS)
        AC_SUBST(RADIATION_FCFLAGS)
        AC_SUBST(RADIATION_LIBS)
        AC_SUBST(RADIATION_PREFIX)
        AC_SUBST(RADIATION_FLIBS)
    fi

    AC_SUBST(HAVE_RADIATION)

# fi

])
