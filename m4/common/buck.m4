# SYNOPSIS
#
#   Test for libBUCK
#
#   AX_PATH_BUCK( <Minimum Required Version>, <package-required=yes/no> )
#
# DESCRIPTION
#
#   Provides a --with-BUCK=DIR option. Searches --with-BUCK,
#   $BUCK_DIR, and the usual places for BUCK headers and libraries.
#
#   On success, sets BUCK_CFLAGS, BUCK_LIBS, and #defines HAVE_BUCK.
#   Also defines automake conditional BUCK_ENABLED.  Assumes package
#   is optional unless overridden with $2=yes.
#
# COPYLEFT
#
#   Copyright (c) 2012 Sylvain <splessis@ices.utexas.edu>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_PATH_BUCK],
[

AC_ARG_VAR(BUCK_DIR,[root directory of BUCK installation])

AC_ARG_WITH(BUCK,
  [AS_HELP_STRING([--with-BUCK[=DIR]],[root directory of BUCK installation (default = BUCK_DIR)])],
  [with_BUCK=$withval
if test "${with_BUCK}" != yes; then
    BUCK_PREFIX=$withval
fi
],[
# assume a sensible default of --with-BUCK=yes
with_BUCK=yes
if test "x${BUCK_DIR}" != "x"; then
   BUCK_PREFIX=${BUCK_DIR}
fi
])

# package requirement; if not specified, the default is to assume that
# the package is optional

is_package_required=ifelse([$2], ,no, $2 )

HAVE_BUCK=0

# logic change: if the user called the macro, check for package,
# decide what to do based on whether the package is required or not.

# if test "${with_BUCK}" != no ; then

    if test -d "${BUCK_PREFIX}/lib" ; then
       BUCK_LDFLAGS="-L${BUCK_PREFIX}/lib -Wl,-rpath,${BUCK_PREFIX}/lib"
       BUCK_LIBS="-lBUCK"
       BUCK_FCFLAGS="-I${BUCK_PREFIX}/lib"
    fi

    if test -d "${BUCK_PREFIX}/include" ; then
        BUCK_CPPFLAGS="-I${BUCK_PREFIX}/include"
    fi

    ac_BUCK_save_CFLAGS="$CFLAGS"
    ac_BUCK_save_CPPFLAGS="$CPPFLAGS"
    ac_BUCK_save_LDFLAGS="$LDFLAGS"
    ac_BUCK_save_LIBS="$LIBS"

    CFLAGS="${BUCK_CFLAGS} ${CFLAGS}"
    CPPFLAGS="${BUCK_CPPFLAGS} ${CPPFLAGS}"
    LDFLAGS="${BUCK_LDFLAGS} ${LDFLAGS}"
    LIBS="${BUCK_LIBS} ${LIBS}"
    AC_LANG_PUSH([C])
    AC_CHECK_HEADER([BUCK.h],[found_header=yes],[found_header=no])

    #-----------------------
    # Minimum version check
    #----------------------

    min_BUCK_version=ifelse([$1], ,0, $1)

    # looking for major.minor.micro style versioning

    MAJOR_VER=`echo $min_BUCK_version | sed 's/^\([[0-9]]*\).*/\1/'`
    if test "x${MAJOR_VER}" = "x" ; then
       MAJOR_VER=0
    fi

    MINOR_VER=`echo $min_BUCK_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\2/'`
    if test "x${MINOR_VER}" = "x" ; then
       MINOR_VER=0
    fi

    MICRO_VER=`echo $min_BUCK_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\3/'`
    if test "x${MICRO_VER}" = "x" ; then
       MICRO_VER=0
    fi

    # begin additional test(s) if header if available

    if test "x${found_header}" = "xyes" ; then

        AC_MSG_CHECKING(for BUCK - version >= $min_BUCK_version)
        version_succeeded=no

        AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
        @%:@include <BUCK.h>
            ]], [[
            #if BUCK_MAJOR_VERSION > $MAJOR_VER
            /* Sweet nibblets */
            #elif (BUCK_MAJOR_VERSION >= $MAJOR_VER) && (BUCK_MINOR_VERSION > $MINOR_VER)
            /* Winner winner, chicken dinner */
            #elif (BUCK_MAJOR_VERSION >= $MAJOR_VER) && (BUCK_MINOR_VERSION >= $MINOR_VER) && (BUCK_MICRO_VERSION >= $MICRO_VER)
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
        @%:@include <BUCK.h>
            ]], [[
            #if BUCK_MINOR_VERSION > 29
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

   Your BUCK library version does not meet the minimum versioning
   requirements ($min_BUCK_version).  Please use --with-BUCK to specify the location
   of an updated installation or consider upgrading the system version.

          ])
       fi
    fi

    # Fortran support - version 0.30.0 introduced an additional link
    # requirement for Fortran support.  

    if test "$fortran_separate" = "yes" ; then
       BUCK_FLIBS="-L${BUCK_PREFIX}/lib -lBUCK -lBUCKf -Wl,-rpath,${BUCK_PREFIX}/lib"
    else
       BUCK_FLIBS="${BUCK_LIBS}"
    fi

    # Library availability

    AC_MSG_CHECKING([for -lBUCK linkage])

    AC_LANG_PUSH([C++])                                                                                                                                                                                                                      
    AC_LINK_IFELSE(
                  [AC_LANG_PROGRAM([#include "BUCK_version.hpp"],
                                   [BUCK::print_BUCK_version()])],
                  [AC_MSG_RESULT(yes)
                   found_library=yes],
                  [AC_MSG_RESULT(no) 
                   found_library=no])

    fi   dnl end test if header if available

    AC_LANG_POP([C++])

    CFLAGS="$ac_BUCK_save_CFLAGS"
    CPPFLAGS="$ac_BUCK_save_CPPFLAGS"
    LDFLAGS="$ac_BUCK_save_LDFLAGS"
    LIBS="$ac_BUCK_save_LIBS"

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
          AC_MSG_ERROR([libBUCK not found.  Try either --with-BUCK or setting BUCK_DIR.])
       else
          AC_MSG_NOTICE([optional BUCK library not found])
          BUCK_CFLAGS=""   # BUCK_CFLAGS empty on failure
          BUCK_CPPFLAGS="" # BUCK_CPPFLAGS empty on failure          
          BUCK_FCFLAGS=""  # BUCK_FCFLAGS empty on failure
          BUCK_LDFLAGS=""  # BUCK_LDFLAGS empty on failure
          BUCK_LIBS=""     # BUCK_LIBS empty on failure
          BUCK_FLIBS=""    # BUCK_FLIBS empty on failure
       fi
    else
        HAVE_BUCK=1
        AC_DEFINE(HAVE_BUCK,1,[Define if BUCK is available])
        AC_SUBST(BUCK_CFLAGS)
        AC_SUBST(BUCK_CPPFLAGS)
        AC_SUBST(BUCK_FCFLAGS)
        AC_SUBST(BUCK_LDFLAGS)
        AC_SUBST(BUCK_LIBS)
        AC_SUBST(BUCK_PREFIX)
        AC_SUBST(BUCK_FLIBS)
    fi

    AC_SUBST(HAVE_BUCK)

# fi

AM_CONDITIONAL(BUCK_ENABLED,test x$HAVE_BUCK = x1)

])
