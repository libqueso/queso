# SYNOPSIS
#
#   Test for GSL
#
#   AX_PATH_GSL_NEW( <Minimum Required Version>, <package-required=yes/no> )
#
# DESCRIPTION
#
#   Provides a --with-gsl=DIR option. Searches --with-gsl,
#   $GSL_DIR, and the usual places for GSL headers and libraries.
#
#   On success, sets GSL_CPPFLAGS, GSL_LIBS, and #defines HAVE_GSL.
#   Also defines automake conditional GSL_ENABLED.  Assumes package
#   is optional unless overridden with $2=yes.
#
# LAST MODIFICATION
#
#   $Id: gsl_new.m4 37680 2013-03-08 01:29:27Z roystgnr $
#
# COPYLEFT
#
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

AC_DEFUN([AX_PATH_GSL_NEW],
[

AC_ARG_VAR(GSL_DIR,[root directory of GSL installation])

AC_ARG_WITH(gsl,
  [AS_HELP_STRING([--with-gsl[=DIR]],[root directory of GSL installation (default = GSL_DIR)])],
  [with_gsl=$withval
if test "${with_gsl}" != yes; then
    GSL_PREFIX=$withval
fi
],[
with_gsl=$withval
if test "x${GSL_DIR}" != "x"; then
   GSL_PREFIX=${GSL_DIR}
fi
])

# package requirement; if not specified, the default is to assume that
# the package is optional

is_package_required=ifelse([$2], ,no, $2 )

HAVE_GSL=0

# logic change: if the user called the macro, check for package,
# decide what to do based on whether the package is required or not.

# if test "${with_gsl}" != no ; then

    if test -d "${GSL_PREFIX}/lib" ; then
       GSL_LDFLAGS="-L${GSL_PREFIX}/lib -Wl,-rpath,${GSL_PREFIX}/lib"
       GSL_LIBS="-lgsl -lgslcblas"
    fi

    if test -d "${GSL_PREFIX}/include" ; then
       GSL_CPPFLAGS="-I${GSL_PREFIX}/include"
    fi

    ac_GSL_save_CPPFLAGS="$CPPFLAGS"
    ac_GSL_save_LDFLAGS="$LDFLAGS"
    ac_GSL_save_LIBS="$LIBS"

    CPPFLAGS="${GSL_CPPFLAGS} ${CPPFLAGS}"
    LDFLAGS="${GSL_LDFLAGS} ${LDFLAGS}"
    LIBS="${GSL_LIBS} ${LIBS}"

    #-----------------------
    # Minimum version check
    #----------------------

    AC_PATH_PROG([GSL_CONFIG], 
                 [gsl-config], 
                 [])

    if test "x$GSL_CONFIG" = "x"; then
       AC_MSG_ERROR([ Could not find gsl-config problem. 
                   Please use --with-gsl to specify the location
                   of the GSL library. ])
    fi

    min_gsl_version=ifelse([$1], ,1.0, $1)

    # looking for major.minor.micro style versioning

    MIN_MAJOR_VER=`echo $min_gsl_version | sed 's/^\([[0-9]]*\).*/\1/'`
    if test "x${MIN_MAJOR_VER}" = "x" ; then
       MIN_MAJOR_VER=0
    fi

    MIN_MINOR_VER=`echo $min_gsl_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\2/'`
    if test "x${MIN_MINOR_VER}" = "x" ; then
       MIN_MINOR_VER=0
    fi

    MIN_MICRO_VER=`echo $min_gsl_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\3/'`
    if test "x${MIN_MICRO_VER}" = "x" ; then
       MIN_MICRO_VER=0
    fi

    GSL_MAJOR_VER=`$GSL_CONFIG --version | \
           sed 's/^\([[0-9]]*\).*/\1/'`
    if test "x${GSL_MAJOR_VER}" = "x" ; then
       GSL_MAJOR_VER=0
    fi

    GSL_MINOR_VER=`$GSL_CONFIG --version | \
           sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\2/'`
    if test "x${GSL_MINOR_VER}" = "x" ; then
       GSL_MINOR_VER=0
    fi

    GSL_MICRO_VER=`$GSL_CONFIG --version | \
           sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\3/'`
    if test "x${GSL_MICRO_VER}" = "x" ; then
       GSL_MICRO_VER=0
    fi

    AC_MSG_CHECKING(for gsl - version >= $min_gsl_version)
    version_succeeded=no

    AC_LANG_PUSH([C++])
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[]], [[
            #if $GSL_MAJOR_VER > $MIN_MAJOR_VER
            /* Sweet nibblets */
            #elif ($GSL_MAJOR_VER >= $MIN_MAJOR_VER) && ($GSL_MINOR_VER > $MIN_MINOR_VER)
            /* Winner winner, chicken dinner */
            #elif ($GSL_MAJOR_VER >= $MIN_MAJOR_VER) && ($GSL_MINOR_VER >= $MIN_MINOR_VER) && ($GSL_MICRO_VER >= $MIN_MICRO_VER)
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
          AC_MSG_ERROR([

   Your GSL version does not meet the minimum versioning
   requirements ($min_gsl_version).  Please use --with-gsl to specify the location
   of an updated installation or consider upgrading the system version.

          ])
       fi
    fi

    #-------------------------
    # Check for a basic header
    #-------------------------

    AC_LANG_PUSH([C])
    AC_CHECK_HEADER([gsl/gsl_math.h],[found_header=yes],[found_header=no])
    AC_LANG_POP([C])


    #-------------------------------
    # Check for a common function
    # gamma function has been around 
    # since before v1.0
    #-------------------------------

    AC_MSG_CHECKING([for -lgsl linkage])

    AC_LANG_PUSH([C])

    AC_CHECK_LIB([gsl],
                 [gsl_sf_gamma],
                 [found_library=yes],
                 [found_library=no])

    AC_LANG_POP([C])

    CPPFLAGS="$ac_GSL_save_CPPFLAGS"
    LDFLAGS="$ac_GSL_save_LDFLAGS"
    LIBS="$ac_GSL_save_LIBS"

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
          AC_MSG_ERROR([GSL not found.  Try either --with-gsl or setting GSL_DIR.])
       else
          AC_MSG_NOTICE([optional GSL library not found])
          GSL_CPPFLAGS=""   # GSL_CFLAGS empty on failure
          GSL_LDFLAGS=""    # GSL_LDFLAGS empty on failure
          GSL_LIBS=""       # GSL_LIBS empty on failure
       fi
    else
        HAVE_GSL=1
        AC_DEFINE(HAVE_GSL,1,[Define if GSL is available])
        AC_SUBST(GSL_CPPFLAGS)
        AC_SUBST(GSL_LDFLAGS)
        AC_SUBST(GSL_LIBS)
        AC_SUBST(GSL_PREFIX)
    fi

    AC_SUBST(HAVE_GSL)

# fi

AM_CONDITIONAL(GSL_ENABLED,test x$HAVE_GSL = x1)

])
