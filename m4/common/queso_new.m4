# SYNOPSIS - (ks note: work in progress.....)
#
#   Test for QUESO 
#
#   AM_PATH_QUESO_NEW( <Minimum Required Version>, <package-required=yes/no> )
#
# DESCRIPTION
#
#   Provides a --with-queso=DIR option. Searches --with-queso,
#   $QUESO_DIR, and the usual places for QUESO headers and libraries.
#
#   On success, sets QUESO_CFLAGS, QUESO_LIBS, and #defines HAVE_QUESO.  
#   Also defines automake conditional GRVY_ENABLED.  Assumes package
#   is optional unless overridden with $2=yes.
#
# COPYLEFT
#
#   Copyright (c) 2011 Karl W. Schulz <karl@ices.utexas.edu>
#   Copyright (c) 2010 Paul T. Bauman <pbauman@ices.utexas.edu>
#   Copyright (c) 2009 Rhys Ulerich <rhys.ulerich@gmail.com>
#   Copyright (c) 2008 Thomas Porschberg <thomas@randspringer.de>
#   Copyright (c) 2008 Caolan McNamara <caolan@skynet.ie>
#   Copyright (c) 2008 Alexandre Duret-Lutz <adl@gnu.org>
#   Copyright (c) 2008 Matthew Mueller <donut@azstarnet.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_PATH_QUESO_NEW],
[

AC_ARG_VAR(QUESO_DIR,[root directory of QUESO installation])

AC_ARG_WITH(queso, 
  [AS_HELP_STRING([--with-queso[=DIR]],[root directory of QUESO installation (default = QUESO_DIR)])],
  [with_queso=$withval
if test "${with_queso}" != yes; then
    QUESO_PREFIX=$withval
fi
],[
# assume a sensible default of --with-queso=yes
with_queso=yes
if test "x${QUESO_DIR}" != "x"; then
   QUESO_PREFIX=${QUESO_DIR}
fi
])

# package requirement; if not specified, the default is to assume that
# the package is optional

is_package_required=ifelse([$2], ,no, $2 )

HAVE_QUESO=0

# logic change: if the user called the macro, check for package,
# decide what to do based on whether the package is required or not.

# if test "${with_queso}" != no ; then

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

    AC_CHECK_HEADER([queso/queso.h], [found_header=yes], [found_header=no])

    #-----------------------
    # Minimum version check
    #----------------------

    min_queso_version=ifelse([$1], ,0.29, $1)

    # looking for major.minor.micro style versioning

    MAJOR_VER=`echo $min_queso_version | sed 's/^\([[0-9]]*\).*/\1/'`
    if test "x${MAJOR_VER}" = "x" ; then
       MAJOR_VER=0
    fi

    MINOR_VER=`echo $min_queso_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\2/'`
    if test "x${MINOR_VER}" = "x" ; then
       MINOR_VER=0
    fi

    MICRO_VER=`echo $min_queso_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\3/'`
    if test "x${MICRO_VER}" = "x" ; then
       MICRO_VER=0
    fi

    # begin additional test(s) if header if available

    if test "x${found_header}" = "xyes" ; then

        AC_MSG_CHECKING(for QUESO - version >= $min_queso_version)
        version_succeeded=no

        AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
        @%:@include <queso/queso.h>
            ]], [[
            #if QUESO_MAJOR_VERSION > $MAJOR_VER
            /* Sweet nibblets */
            #elif (QUESO_MAJOR_VERSION >= $MAJOR_VER) && (QUESO_MINOR_VERSION > $MINOR_VER)
            /* Winner winner, chicken dinner */
            #elif (QUESO_MAJOR_VERSION >= $MAJOR_VER) && (QUESO_MINOR_VERSION >= $MINOR_VER) && (QUESO_MICRO_VERSION >= $MICRO_VER)
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

    fi   dnl end test if header if available


#####     AC_COMPILE_IFELSE([#include "queso/queso.h"],[found_header=yes],[found_header=no])
##### 
#####     ac_QUESO_BOOST_PROGRAM_OPTIONS_LDFLAGS_COMPILER=''
#####     
#####     for ldflag in $BOOST_PROGRAM_OPTIONS_LDFLAGS; do
#####         ac_QUESO_BOOST_PROGRAM_OPTIONS_LDFLAGS_COMPILER="$ac_QUESO_BOOST_PROGRAM_OPTIONS_LDFLAGS_COMPILER -Wl,$ldflag"
#####     done
##### 
#####     dnl PB: Check below is "newer" over deprecated AC_HAVE_LIBRARY. 
#####     AC_CHECK_LIB([queso], [main], [found_library=yes], [found_library=no],[ ${GSL_LIBS} 
#####                                                                             ${EPETRA_LDFLAGS} -lepetra 
#####                                                                             ${GRVY_LIBS} 
#####                                                                             ${ac_QUESO_BOOST_PROGRAM_OPTIONS_LDFLAGS_COMPILER} ${BOOST_PROGRAM_OPTIONS_LIBS}
##### 									    ${GLPK_LIBS}
##### 									    ${HDF5_LIBS} ])
##### 
    dnl PB: LDFLAGS and LIBS do not expand as I would have expected so using above test instead to resolve
    dnl PB: library dependencies.
    dnl FIXME: Question: Why do we even have to resolve those at all???
    dnl AC_CHECK_LIB([queso], [main], [found_library=yes], [found_library=no],[ ${LDFLAGS} ${LIBS} ] )									   

    dnl PB: http://nerdland.net/2009/07/detecting-c-libraries-with-autotools
    dnl PB: Gives a good guide to rolling up a test for C++. Will take a peek in the future.
    dnl AC_CHECK_LIB([queso], [uqBaseEnvironmentClass::worldRank], [found_library=yes], [found_library=no])

    if test "$version_succeeded" != "yes";then
       if test "$is_package_required" = yes; then
          AC_MSG_ERROR([

   Your QUESO library version does not meet the minimum versioning
   requirements ($min_queso_version).  Please use --with-queso to specify the location
   of an updated installation or consider upgrading the system version.

          ])
       fi
    fi

    AC_LANG_POP([C++])

    CFLAGS="$ac_QUESO_save_CFLAGS"
    CPPFLAGS="$ac_QUESO_save_CPPFLAGS"
    LDFLAGS="$ac_QUESO_save_LDFLAGS"
    LIBS="$ac_QUESO_save_LIBS"

    succeeded=no
    if test "$found_header" = yes; then
            succeeded=yes
    fi

    if test "$succeeded" = no; then
       if test "$is_package_required" = yes; then
          AC_MSG_ERROR([QUESO not found.  Try either --with-queso or setting QUESO_DIR.])
       else
          AC_MSG_NOTICE([optional QUESO library not found])
          QUESO_CFLAGS=""   # empty on failure
          QUESO_CPPFLAGS="" # empty on failure
          QUESO_LIBS=""     # empty on failure
       fi
    else
        HAVE_QUESO=1
        AC_DEFINE(HAVE_QUESO,1,[Define if QUESO is available])
        AC_SUBST(QUESO_CFLAGS)
        AC_SUBST(QUESO_CPPFLAGS)
        AC_SUBST(QUESO_LIBS)
	AC_SUBST(QUESO_PREFIX)
    fi

  AC_SUBST(HAVE_QUESO)

# fi

AM_CONDITIONAL(QUESO_ENABLED,test x$HAVE_QUESO = x1)

])
