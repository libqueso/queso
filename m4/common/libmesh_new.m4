# SYNOPSIS
#
#   Test for LIBMESH
#
#   AM_PATH_LIBMESH_NEW( <Minimum Required Version>, <package-required=yes/no> )
#
# DESCRIPTION
#
#   Provides a --with-libmesh=DIR option. Searches --with-libmesh,
#   $LIBMESH_DIR, and the usual places for LIBMESH headers and libraries.
#
#   On success, sets LIBMESH_CFLAGS, LIBMESH_LIBS, and
#   #defines HAVE_LIBMESH.  When ACTION-IF-NOT-FOUND is not specified,
#   the default behavior is for configure to fail.
#
# LAST MODIFICATION
#
#   2009-08-28
#
# COPYLEFT#

#   Copyright (c) 2012 Paul T. Bauman <pbauman@ices.utexas.edu>
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

AC_DEFUN([AX_PATH_LIBMESH_NEW],
[

AC_ARG_VAR(LIBMESH_DIR,[root directory of LIBMESH installation])

AC_ARG_WITH(libmesh, 
  [AS_HELP_STRING([--with-libmesh[=DIR]],[root directory of LIBMESH installation (default = LIBMESH_DIR)])],
  [with_libmesh=$withval
if test "${with_libmesh}" != yes; then
    LIBMESH_PREFIX=$withval
fi
],[
with_libmesh=$withval
if test "x${LIBMESH_DIR}" != "x"; then
   LIBMESH_PREFIX=${LIBMESH_DIR}
fi
if test "x${LIBMESH_ROOT}" != "x"; then
   LIBMESH_PREFIX=${LIBMESH_ROOT}
fi
])

# --------------------------------------------------------------
# package requirement: if not specified, the default is to 
# assume that the package is optional.
# --------------------------------------------------------------
is_package_required=ifelse([$2], ,no, $2 )

# --------------------------------------------------------------
# libMesh installs libs and binaries based on the host and 
# $METHOD so we need to check there for relevant stuff. 
# If $METHOD is not set, we default to "opt" as libMesh does.
# --------------------------------------------------------------
AC_REQUIRE([AC_CANONICAL_TARGET])

libmesh_host=$host

if test "x${METHOD}" = "x"; then
   METHOD=opt
fi

LIBMESH_LIBDIR=${LIBMESH_PREFIX}/lib
libmesh_bindir=${LIBMESH_PREFIX}/bin

#--------------------------------------------------------------
# First check for libmesh-config since we need that for getting 
# CPPFLAGS among other things.
#--------------------------------------------------------------
AC_CHECK_FILE(${libmesh_bindir}/libmesh-config,
              [LIBMESH_CONFIG=${libmesh_bindir}/libmesh-config
               libmesh_config_found=yes],
              [libmesh_config_found=no] )

# If we didn't find a recent version of libMesh with autotools-based build system
# check for 0.8.0-type version
if test "$libmesh_config_found" != "yes";then
  LIBMESH_LIBDIR=${LIBMESH_PREFIX}/lib/${libmesh_host}_${METHOD}
  libmesh_bindir=${LIBMESH_PREFIX}/bin/${libmesh_host}_${METHOD}
  AC_CHECK_FILE(${libmesh_bindir}/libmesh-config,
              [LIBMESH_CONFIG=${libmesh_bindir}/libmesh-config
               libmesh_config_found=yes],
              [libmesh_config_found=no] )
fi
  
if test "$libmesh_config_found" != "yes";then
  if test "$is_package_required" = "yes";then
     AC_MSG_ERROR([Cannot find libmesh-config! Please use --with-libmesh to specify
                 the location of a valid libmesh installation.])
  fi
  found_header=no
else

  LIBMESH_CPPFLAGS=`$LIBMESH_CONFIG --cppflags --include | tr '\n' ' '`
  LIBMESH_INCLUDE=$LIBMESH_CPPFLAGS
  LIBMESH_CXXFLAGS=`$LIBMESH_CONFIG --cxxflags`

  # Not all libmesh-config versions support all arguments
  LIBMESH_LDFLAGS=`$LIBMESH_CONFIG --ldflags | grep -v "Unknown argument" | grep -v libmesh-config || $LIBMESH_CONFIG --libs`
  LIBMESH_LIBS=`$LIBMESH_CONFIG --libs | grep -v "Unknown argument" | grep -v libmesh-config || $LIBMESH_CONFIG --ldflags`

  LIBMESH_CXX=`$LIBMESH_CONFIG --cxx`
  LIBMESH_CC=`$LIBMESH_CONFIG --cc`
  LIBMESH_FC=`$LIBMESH_CONFIG --fc`

  ac_LIBMESH_save_CPPFLAGS="$CPPFLAGS"
  ac_LIBMESH_save_LDFLAGS="$LDFLAGS"
  ac_LIBMESH_save_LIBS="$LIBS"

  CPPFLAGS="${LIBMESH_CPPFLAGS} ${CPPFLAGS}"
  LDFLAGS="${LIBMESH_LDFLAGS} ${LDFLAGS}"
  LIBS="${LIBMESH_LIBS} ${LIBS}"

  #--------------------------------------------------------------
  # Now check for the libmesh_version.h header
  #--------------------------------------------------------------
  AC_LANG_PUSH([C++])
  AC_CHECK_HEADER([libmesh/libmesh_version.h],[found_header=yes],[found_header=no])
  AC_LANG_POP([C++])


  #--------------------------------------------------------------
  # Minimum version check: looking for major.minor.micro style 
  #                        versioning
  #--------------------------------------------------------------
  min_libmesh_version=ifelse([$1], ,0.8.0, $1)


  MAJOR_VER=`echo $min_libmesh_version | sed 's/^\([[0-9]]*\).*/\1/'`
  if test "x${MAJOR_VER}" = "x" ; then
     MAJOR_VER=0
  fi

  MINOR_VER=`echo $min_libmesh_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\2/'`
  if test "x${MINOR_VER}" = "x" ; then
     MINOR_VER=0
  fi

  MICRO_VER=`echo $min_libmesh_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\3/'`
  if test "x${MICRO_VER}" = "x" ; then
     MICRO_VER=0
  fi

  if test "x${found_header}" = "xyes" ; then

    AC_MSG_CHECKING(for libMesh version >= $min_libmesh_version)
    version_succeeded=no

    AC_LANG_PUSH([C++])
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
                @%:@include "libmesh/libmesh_version.h"
                ]], [[
                #if LIBMESH_MAJOR_VERSION > $MAJOR_VER
                /* Sweet nibblets */
                #elif (LIBMESH_MAJOR_VERSION >= $MAJOR_VER) && (LIBMESH_MINOR_VERSION > $MINOR_VER)
                /* Winner winner, chicken dinner */
                #elif (LIBMESH_MAJOR_VERSION >= $MAJOR_VER) && (LIBMESH_MINOR_VERSION >= $MINOR_VER) && (LIBMESH_MICRO_VERSION >= $MICRO_VER)
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
    AC_LANG_POP([C++])

    if test "$version_succeeded" != "yes";then
      if test "$is_package_required" = yes; then
        AC_MSG_ERROR([Your LIBMESH version does not meet the minimum versioning
                      requirements ($min_libmesh_version). Please use 
                      --with-libmesh to specify the location of an updated 
                      installation or consider upgrading the system version. ])
      fi
    fi


    #--------------------------------------------------------------
    # Check for libMesh library linking
    #--------------------------------------------------------------
    AC_MSG_CHECKING([for libMesh linkage])

    AC_LANG_PUSH([C++])
    AC_LINK_IFELSE( [AC_LANG_PROGRAM([#include "libmesh/libmesh_version.h"],
                                     [libMesh::get_libmesh_version()])],
                                     [AC_MSG_RESULT(yes)
                                      found_library=yes],
                                     [AC_MSG_RESULT(no) 
                                      found_library=no] )
    AC_LANG_POP([C++])
  fi   dnl end test if header available
fi   dnl end test if libmesh-config available

CPPFLAGS="$ac_LIBMESH_save_CPPFLAGS"
LDFLAGS="$ac_LIBMESH_save_LDFLAGS"
LIBS="$ac_LIBMESH_save_LIBS"

succeeded=no
if test "x$found_header" = xyes; then
   if test "x$version_succeeded" = xyes; then
      if test "x$found_library" = xyes; then
         succeeded=yes
      fi
   fi
fi

HAVE_LIBMESH=0
LIBMESH_HAVE_LIBTOOL=0

if test "$succeeded" = no; then
   if test "$is_package_required" = yes; then
      AC_MSG_ERROR([libMesh not found.  Try either --with-libmesh or setting LIBMESH_DIR.])
   else
      AC_MSG_NOTICE([optional libMesh library not found])
      LIBMESH_CPPFLAGS=""   # LIBMESH_ variables empty on failure
      LIBMESH_LDFLAGS=""
      LIBMESH_LIBS=""
      LIBMESH_CXXFLAGS=""
      LIBMESH_PREFIX=""
      LIBMESH_INCLUDE=""
      LIBMESH_LIBDIR=""
      LIBMESH_CXX=""
      LIBMESH_CC=""
      LIBMESH_FC=""
   fi
else
   HAVE_LIBMESH=1

   AC_CHECK_FILE(${LIBMESH_LIBDIR}/libmesh_${METHOD}.la,
              [LIBMESH_HAVE_LIBTOOL=1],
              [LIBMESH_HAVE_LIBTOOL=0] )

   AC_DEFINE(HAVE_LIBMESH,1,[Define if LIBMESH is available])
   AC_SUBST(LIBMESH_CPPFLAGS)
   AC_SUBST(LIBMESH_INCLUDE)
   AC_SUBST(LIBMESH_LDFLAGS)
   AC_SUBST(LIBMESH_LIBS)
   AC_SUBST(LIBMESH_PREFIX)
   AC_SUBST(LIBMESH_CXXFLAGS)
   AC_SUBST(LIBMESH_LIBDIR)
   AC_SUBST(LIBMESH_CXX)
   AC_SUBST(LIBMESH_CC)
   AC_SUBST(LIBMESH_FC)
fi

LIBMESH_METHOD=${METHOD}
AC_SUBST(LIBMESH_METHOD)
AC_SUBST(HAVE_LIBMESH)
AC_SUBST(LIBMESH_HAVE_LIBTOOL)

AM_CONDITIONAL(LIBMESH_ENABLED,test x$HAVE_LIBMESH = x1)
AM_CONDITIONAL(LIBMESH_LIBTOOL,test x$LIBMESH_HAVE_LIBTOOL = x1)

])
