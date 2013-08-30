# SYNOPSIS
#
#   Test for MASA Library
#
#   AM_PATH_MASA( <Minimum Required Version>, <package-required=yes/no> )
#
# DESCRIPTION
#
#   Provides a --with-esio=DIR option. Searches --with-esio,
#   $ESIO_DIR, and the usual places for ESIO headers and libraries.
#
#   On success, sets ESIO_CXXFLAGS, ESIO_LIBS, ESIO_FC_LIBS (for
#   Fortran) and #defines HAVE_ESIO.  Assumes package is optional
#   unless overridden with $2=yes
#
# LAST MODIFICATION
#
#   $Id: 
#
# COPYLEFT
#
#   Copyright (c) 2010 Nicholas Malaya   <nick@ices.utexas.edu>
#   Copyright (c) 2010 Karl W. Schulz    <karl@ices.utexas.edu>
#   Copyright (c) 2009 Rhys Ulerich      <rhys.ulerich@gmail.com>
#   Copyright (c) 2008 Thomas Porschberg <thomas@randspringer.de>
#   Copyright (c) 2008 Caolan McNamara   <caolan@skynet.ie>
#   Copyright (c) 2008 Alexandre Duret-Lutz <adl@gnu.org>
#   Copyright (c) 2008 Matthew Mueller   <donut@azstarnet.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.
#

AC_DEFUN([AX_PATH_ESIO],
[

AC_ARG_VAR(ESIO_DIR,[root directory of ESIO installation])

PKG_CHECK_MODULES([ESIO],
    [esio],
    AC_DEFINE([HAVE_ESIO],1,[Define ESIO available]),AC_MSG_WARN([Could not find ESIO pkg-config file; Continuing...]))

AC_ARG_WITH(esio, 
  [AS_HELP_STRING([--with-esio[=DIR]],[root directory of ESIO installation (default = ESIO_DIR)])],
  [with_esio=$withval
if test "${with_esio}" != yes; then
    ESIO_PREFIX=$withval
fi
],[
with_esio=$withval
if test "x${ESIO_DIR}" != "x"; then
   ESIO_PREFIX=${ESIO_DIR}
fi
])

# package requirement; if not specified, the default is to assume that
# the package is optional

is_package_required=ifelse([$2], ,no, $2 )

HAVE_ESIO=0

if test "${with_esio}" != no ; then

    if test -d "${ESIO_PREFIX}/lib" ; then
       ESIO_LIBS="-L${ESIO_PREFIX}/lib -lesio -Wl,-rpath,${ESIO_PREFIX}/lib"
       ESIO_FC_LIBS="-L${ESIO_PREFIX}/lib -lfesio -lesio -Wl,-rpath,${ESIO_PREFIX}/lib"
       ESIO_FCFLAGS="-I${ESIO_PREFIX}/lib -DHAVE_ESIO"
    fi

    if test -d "${ESIO_PREFIX}/include" ; then
        ESIO_CXXFLAGS="-I${ESIO_PREFIX}/include"
    fi

    ac_ESIO_save_CXXFLAGS="$CXXFLAGS"
    ac_ESIO_save_CPPFLAGS="$CPPFLAGS"
    ac_ESIO_save_LDFLAGS="$LDFLAGS"
    ac_ESIO_save_LIBS="$LIBS"

    CXXFLAGS="${ESIO_CXXFLAGS} ${CXXFLAGS}"
    CPPFLAGS="${ESIO_CXXFLAGS} ${CPPFLAGS}"
    LDFLAGS="${ESIO_LIBS} ${LDFLAGS}"

    AC_LANG_PUSH([C++])
    AC_CHECK_HEADER([esio/esio.h],[found_header=yes],[found_header=no])

    #-----------------------
    # Minimum version check
    #----------------------

    min_esio_version=ifelse([$1], ,0.10, $1)	

    # looking for major.minor.micro style versioning

    MAJOR_VER=`echo $min_esio_version | sed 's/^\([[0-9]]*\).*/\1/'`
    if test "x${MAJOR_VER}" = "x" ; then
       MAJOR_VER=0
    fi

    MINOR_VER=`echo $min_esio_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\2/'`
    if test "x${MINOR_VER}" = "x" ; then
       MINOR_VER=0
    fi

    MICRO_VER=`echo $min_esio_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\3/'`
    if test "x${MICRO_VER}" = "x" ; then
       MICRO_VER=0
    fi

    dnl begin additional test(s) if header if available

    if test "x${found_header}" = "xyes" ; then

        AC_MSG_CHECKING(for esio - version >= $min_esio_version)
        version_succeeded=no

    	AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
       	@%:@include <esio/esio.h>
            ]], [[
            #if ESIO_MAJOR_VERSION > $MAJOR_VER
            /* Sweet nibblets */
            #elif (ESIO_MAJOR_VERSION >= $MAJOR_VER) && (ESIO_MINOR_VERSION >= $MINOR_VER) && (ESIO_MICRO_VERSION >= $MICRO_VER)
            /* Winner winner, chicken dinner */
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

   Your ESIO library version does not meet the minimum versioning
   requirements ($min_esio_version).  Please use --with-esio to specify the location
   of an updated installation or consider upgrading the system version.

          ]) 
       fi 
    fi     

    # Library availability

    AC_MSG_CHECKING([for -lesio linkage])

    AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([#include <esio/esio.h>],[ESIO::esio_version_stdout])],
    [TEST_LIBS="$TEST_LIBS -lesio"] [
    AC_MSG_RESULT(yes)
    found_library=yes ],[AC_MSG_RESULT(no)])

    fi   dnl end test if header if available

    AC_LANG_POP([C++])

    CXXFLAGS="$ac_ESIO_save_CXXFLAGS"
    CPPFLAGS="$ac_ESIO_save_CPPFLAGS"
    LDFLAGS="$ac_ESIO_save_LDFLAGS"
    LIBS="$ac_ESIO_save_LIBS"

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
       	  AC_MSG_ERROR([ESIO not found.  Try either --with-esio or setting ESIO_DIR.])
       else
          AC_MSG_NOTICE([optional ESIO library not found])
       fi
    else
        HAVE_ESIO=1
        AC_DEFINE(HAVE_ESIO,1,[Define if ESIO is available])
        AC_SUBST(ESIO_CXXFLAGS)
	AC_SUBST(ESIO_FCFLAGS)
        AC_SUBST(ESIO_LIBS)
        AC_SUBST(ESIO_FC_LIBS)
	AC_SUBST(ESIO_PREFIX)
    fi

    AC_SUBST(HAVE_ESIO)
fi

AM_CONDITIONAL(ESIO_ENABLED,test x$HAVE_ESIO = x1)

])



