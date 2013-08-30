# SYNOPSIS
#
#   Test for libPPPH
#
#   AX_PATH_PPPH( <Minimum Required Version>, <package-required=yes/no> )
#
# DESCRIPTION
#
#   Provides a --with-3PH=DIR option. Searches --with-3PH,
#   $PPPH_DIR, and the usual places for 3PH headers and libraries.
#
#   On success, sets PPPH_CFLAGS, PPPH_LIBS, and #defines HAVE_PPPH.
#   Also defines automake conditional PPPH_ENABLED.  Assumes package
#   is optional unless overridden with $2=yes.
#
# COPYLEFT
#
#   Copyright (c) 2012 Sylvain <splessis@ices.utexas.edu>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_PATH_3PH],
[

AC_ARG_VAR(PPPH_DIR,[root directory of 3PH installation])

AC_ARG_WITH(3PH,
  [AS_HELP_STRING([--with-3PH[=DIR]],[root directory of 3PH installation (default = PPPH_DIR)])],
  [with_3PH=$withval
if test "${with_3PH}" != yes; then
    PPPH_PREFIX=$withval
fi
],[
# assume a sensible default of --with-3PH=yes
with_3PH=yes
if test "x${PPPH_DIR}" != "x"; then
   PPPH_PREFIX=${PPPH_DIR}
fi
])

# package requirement; if not specified, the default is to assume that
# the package is optional

is_package_required=ifelse([$2], ,no, $2 )

HAVE_PPPH=0

# logic change: if the user called the macro, check for package,
# decide what to do based on whether the package is required or not.

# if test "${with_3PH}" != no ; then

    if test -d "${PPPH_PREFIX}/lib" ; then
       PPPH_LDFLAGS="-L${PPPH_PREFIX}/lib -Wl,-rpath,${PPPH_PREFIX}/lib"
       PPPH_LIBS="-l3PH"
       PPPH_FCFLAGS="-I${PPPH_PREFIX}/lib"
    fi

    if test -d "${PPPH_PREFIX}/include" ; then
        PPPH_CPPFLAGS="-I${PPPH_PREFIX}/include"
    fi

    ac_PPPH_save_CFLAGS="$CFLAGS"
    ac_PPPH_save_CPPFLAGS="$CPPFLAGS"
    ac_PPPH_save_LDFLAGS="$LDFLAGS"
    ac_PPPH_save_LIBS="$LIBS"

    CFLAGS="${PPPH_CFLAGS} ${CFLAGS}"
    CPPFLAGS="${PPPH_CPPFLAGS} ${CPPFLAGS}"
    LDFLAGS="${PPPH_LDFLAGS} ${LDFLAGS}"
    LIBS="${PPPH_LIBS} ${LIBS}"
    AC_LANG_PUSH([C])
    AC_CHECK_HEADER([3PH.h],[found_header=yes],[found_header=no])

    #-----------------------
    # Minimum version check
    #----------------------

    min_PPPH_version=ifelse([$1], ,0, $1)

    # looking for major.minor.micro style versioning

    MAJOR_VER=`echo $min_PPPH_version | sed 's/^\([[0-9]]*\).*/\1/'`
    if test "x${MAJOR_VER}" = "x" ; then
       MAJOR_VER=0
    fi

    MINOR_VER=`echo $min_PPPH_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\2/'`
    if test "x${MINOR_VER}" = "x" ; then
       MINOR_VER=0
    fi

    MICRO_VER=`echo $min_PPPH_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\3/'`
    if test "x${MICRO_VER}" = "x" ; then
       MICRO_VER=0
    fi

    # begin additional test(s) if header if available

    if test "x${found_header}" = "xyes" ; then

        AC_MSG_CHECKING(for 3PH - version >= $min_PPPH_version)
        version_succeeded=no

        AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
        @%:@include <3PH.h>
            ]], [[
            #if PPPH_MAJOR_VERSION > $MAJOR_VER
            /* Sweet nibblets */
            #elif (PPPH_MAJOR_VERSION >= $MAJOR_VER) && (PPPH_MINOR_VERSION > $MINOR_VER)
            /* Winner winner, chicken dinner */
            #elif (PPPH_MAJOR_VERSION >= $MAJOR_VER) && (PPPH_MINOR_VERSION >= $MINOR_VER) && (PPPH_MICRO_VERSION >= $MICRO_VER)
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
        @%:@include <3PH.h>
            ]], [[
            #if PPPH_MINOR_VERSION > 29
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

   Your 3PH library version does not meet the minimum versioning
   requirements ($min_PPPH_version).  Please use --with-3PH to specify the location
   of an updated installation or consider upgrading the system version.

          ])
       fi
    fi

    # Fortran support - version 0.30.0 introduced an additional link
    # requirement for Fortran support.  

    if test "$fortran_separate" = "yes" ; then
       PPPH_FLIBS="-L${PPPH_PREFIX}/lib -l3PH -l3PHf -Wl,-rpath,${PPPH_PREFIX}/lib"
    else
       PPPH_FLIBS="${PPPH_LIBS}"
    fi

    # Library availability

    AC_MSG_CHECKING([for -l3PH linkage])

    AC_LANG_PUSH([C++])                                                                                                                                                                                                                      
    AC_LINK_IFELSE(
                  [AC_LANG_PROGRAM([#include "3PH_version.hpp"],
                                   [PPPH::print_3PH_version()])],
                  [AC_MSG_RESULT(yes)
                   found_library=yes],
                  [AC_MSG_RESULT(no) 
                   found_library=no])

    fi   dnl end test if header if available

    AC_LANG_POP([C++])

    CFLAGS="$ac_PPPH_save_CFLAGS"
    CPPFLAGS="$ac_PPPH_save_CPPFLAGS"
    LDFLAGS="$ac_PPPH_save_LDFLAGS"
    LIBS="$ac_PPPH_save_LIBS"

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
          AC_MSG_ERROR([lib3PH not found.  Try either --with-3PH or setting PPPH_DIR.])
       else
          AC_MSG_NOTICE([optional 3PH library not found])
          PPPH_CFLAGS=""   # PPPH_CFLAGS empty on failure
          PPPH_CPPFLAGS="" # PPPH_CPPFLAGS empty on failure
          PPPH_FCFLAGS=""  # PPPH_FCFLAGS empty on failure
          PPPH_LIBS=""     # PPPH_LIBS empty on failure
          PPPH_FLIBS=""    # PPPH_FLIBS empty on failure
       fi
    else
        HAVE_PPPH=1
        AC_DEFINE(HAVE_PPPH,1,[Define if 3PH is available])
        AC_SUBST(PPPH_CFLAGS)
        AC_SUBST(PPPH_CPPFLAGS)
        AC_SUBST(PPPH_FCFLAGS)
        AC_SUBST(PPPH_LIBS)
        AC_SUBST(PPPH_PREFIX)
        AC_SUBST(PPPH_FLIBS)
    fi

    AC_SUBST(HAVE_PPPH)

# fi

AM_CONDITIONAL(PPPH_ENABLED,test x$HAVE_PPPH = x1)

])
