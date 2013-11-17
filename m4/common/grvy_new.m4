# SYNOPSIS
#
#   Test for libGRVY
#
#   AM_PATH_GRVY_NEW( <Minimum Required Version>, <package-required=yes/no> )
#
# DESCRIPTION
#
#   Provides a --with-grvy=DIR option. Searches --with-grvy,
#   $GRVY_DIR, and the usual places for GRVY headers and libraries.
#
#   On success, sets GRVY_CFLAGS, GRVY_LIBS, and #defines HAVE_GRVY.
#   Also defines automake conditional GRVY_ENABLED.  Assumes package
#   is optional unless overridden with $2=yes.
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

AC_DEFUN([AX_PATH_GRVY_NEW],
[

AC_ARG_VAR(GRVY_DIR,[root directory of GRVY installation])

AC_ARG_WITH(grvy,
  [AS_HELP_STRING([--with-grvy[=DIR]],[root directory of GRVY installation (default = GRVY_DIR)])],
  [with_grvy=$withval
if test "${with_grvy}" != yes; then
    GRVY_PREFIX=$withval
fi
],[
# assume a sensible default of --with-grvy=yes
with_grvy=yes
if test "x${GRVY_DIR}" != "x"; then
   GRVY_PREFIX=${GRVY_DIR}
fi
])

# package requirement; if not specified, the default is to assume that
# the package is optional

is_package_required=ifelse([$2], ,no, $2 )

HAVE_GRVY=0

# logic change: if the user called the macro, check for package,
# decide what to do based on whether the package is required or not.

# if test "${with_grvy}" != no ; then

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

    #-----------------------
    # Minimum version check
    #----------------------

    min_grvy_version=ifelse([$1], ,0.29, $1)

    # looking for major.minor.micro style versioning

    MAJOR_VER=`echo $min_grvy_version | sed 's/^\([[0-9]]*\).*/\1/'`
    if test "x${MAJOR_VER}" = "x" ; then
       MAJOR_VER=0
    fi

    MINOR_VER=`echo $min_grvy_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\2/'`
    if test "x${MINOR_VER}" = "x" ; then
       MINOR_VER=0
    fi

    MICRO_VER=`echo $min_grvy_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\3/'`
    if test "x${MICRO_VER}" = "x" ; then
       MICRO_VER=0
    fi

    # begin additional test(s) if header if available

    if test "x${found_header}" = "xyes" ; then

        AC_MSG_CHECKING(for grvy - version >= $min_grvy_version)
        version_succeeded=no

        AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
        @%:@include <grvy.h>
            ]], [[
            #if GRVY_MAJOR_VERSION > $MAJOR_VER
            /* Sweet nibblets */
            #elif (GRVY_MAJOR_VERSION >= $MAJOR_VER) && (GRVY_MINOR_VERSION > $MINOR_VER)
            /* Winner winner, chicken dinner */
            #elif (GRVY_MAJOR_VERSION >= $MAJOR_VER) && (GRVY_MINOR_VERSION >= $MINOR_VER) && (GRVY_MICRO_VERSION >= $MICRO_VER)
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

	# do we need separate fortran linkage?

        AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
        @%:@include <grvy.h>
            ]], [[
            #if GRVY_MINOR_VERSION > 29
            #else
            #  error version is too old
            #endif
        ]])],[
            fortran_separate=yes
        ],[
            fortran_separate=no
        ])

    if test "$version_succeeded" != "yes";then
       if test "$is_package_required" = yes; then
          AC_MSG_ERROR([

   Your GRVY library version does not meet the minimum versioning
   requirements ($min_grvy_version).  Please use --with-grvy to specify the location
   of an updated installation or consider upgrading the system version.

          ])
       fi
    fi

    # Fortran support - version 0.30.0 introduced an additional link
    # requirement for Fortran support.  

    if test "$fortran_separate" = "yes" ; then
       GRVY_FLIBS="-L${GRVY_PREFIX}/lib -lgrvy -lgrvyf -Wl,-rpath,${GRVY_PREFIX}/lib"
    else
       GRVY_FLIBS="${GRVY_LIBS}"
    fi

    # Library availability

    AC_MSG_CHECKING([for -lgrvy linkage])

    AC_CHECK_LIB([grvy],grvy_input_fopen,[found_library=yes],[found_library=no])

    fi   dnl end test if header if available

    AC_LANG_POP([C])

    CFLAGS="$ac_GRVY_save_CFLAGS"
    CPPFLAGS="$ac_GRVY_save_CPPFLAGS"
    LDFLAGS="$ac_GRVY_save_LDFLAGS"
    LIBS="$ac_GRVY_save_LIBS"

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
          AC_MSG_ERROR([libGRVY not found.  Try either --with-grvy or setting GRVY_DIR.])
       else
          AC_MSG_NOTICE([optional GRVY library not found])
          GRVY_CFLAGS=""   # GRVY_CFLAGS empty on failure
          GRVY_FCFLAGS=""  # GRVY_FCFLAGS empty on failure
          GRVY_LIBS=""     # GRVY_LIBS empty on failure
          GRVY_FLIBS=""    # GRVY_FLIBS empty on failure
       fi
    else
        HAVE_GRVY=1
        AC_DEFINE(HAVE_GRVY,1,[Define if GRVY is available])
        AC_SUBST(GRVY_CFLAGS)
        AC_SUBST(GRVY_FCFLAGS)
        AC_SUBST(GRVY_LIBS)
        AC_SUBST(GRVY_PREFIX)
        AC_SUBST(GRVY_FLIBS)
    fi

    AC_SUBST(HAVE_GRVY)

# fi

AM_CONDITIONAL(GRVY_ENABLED,test x$HAVE_GRVY = x1)

])
