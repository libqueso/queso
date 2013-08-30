# SYNOPSIS
#
#   Test for PECOS 1D Ablation library
#
#   AM_PATH_ABLATION1D( <Minimum Required Version>, <package-required=yes/no> )
#
# DESCRIPTION
#
#   Provides a --with-ablation=DIR option. Searches --with-ablation,
#   $ABLATION_DIR, and the usual places for ablation headers and libraries.
#
#   The macro performs a minimum version check and on success, sets
#   ABLATION_CFLAGS, ABLATION_FCFLAGS, ABLATION_LIBS, and #defines
#   HAVE_ABLATION.
#
#   Assumes package is optional unless overridden with $2=yes
#
# LAST MODIFICATION
#
#   $Id: pecos_ablation_new.m4 32126 2012-07-26 18:38:36Z benkirk $
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

AC_DEFUN([AM_PATH_ABLATION1D],
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

# package requirement; if not specified, the default is to assume that
# the package is optional

is_package_required=ifelse([$2], ,no, $2 )

HAVE_ABLATION=0

# if test "${with_ablation}" != no ; then

    if test -d "${ABLATION_PREFIX}/lib" ; then
       ABLATION_LIBS="-L${ABLATION_PREFIX}/lib -lablation -Wl,-rpath,${ABLATION_PREFIX}/lib"
       ABLATION_FCFLAGS="-I${ABLATION_PREFIX}/lib"
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

    #-----------------------
    # Minimum version check
    #----------------------

    min_ablation_version=ifelse([$1], ,0.25.1, $1)	

    # looking for major.minor.micro style versioning

    MAJOR_VER=`echo $min_ablation_version | sed 's/^\([[0-9]]*\).*/\1/'`
    if test "x${MAJOR_VER}" = "x" ; then
       MAJOR_VER=0
    fi

    MINOR_VER=`echo $min_ablation_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\2/'`
    if test "x${MINOR_VER}" = "x" ; then
       MINOR_VER=0
    fi

    MICRO_VER=`echo $min_ablation_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\3/'`
    if test "x${MICRO_VER}" = "x" ; then
       MICRO_VER=0
    fi

    dnl begin additional test(s) if header if available

    if test "x${found_header}" = "xyes" ; then

        AC_MSG_CHECKING(for ablation - version >= $min_ablation_version)
        version_succeeded=no

    	AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
       	@%:@include <ablation1d.h>
            ]], [[
            #if ABLATION_MAJOR_VERSION > $MAJOR_VER
            /* Sweet nibblets */
            #elif (ABLATION_MAJOR_VERSION >= $MAJOR_VER) && (ABLATION_MINOR_VERSION >  $MINOR_VER)
            /* Winner winner, chicken dinner */
            #elif (ABLATION_MAJOR_VERSION >= $MAJOR_VER) && (ABLATION_MINOR_VERSION >= $MINOR_VER) && (ABLATION_MICRO_VERSION >= $MICRO_VER)
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

   Your Ablation library version does not meet the minimum versioning
   requirements ($min_ablation_version).  Please use --with-ablation to specify the location
   of an updated installation or consider upgrading the system version.

          ]) 
       fi 
    fi     

    # Library availability

    AC_CHECK_LIB([ablation],[ablationC1d_version],[found_library=yes],[found_library=no], [ ${ABLATION_LIBS} ])

    fi   dnl end test if header if available

    AC_LANG_POP([C])

    CFLAGS="$ac_ABLATION_save_CFLAGS"
    CPPFLAGS="$ac_ABLATION_save_CPPFLAGS"
    LDFLAGS="$ac_ABLATION_save_LDFLAGS"
    LIBS="$ac_ABLATION_save_LIBS"

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
       	  AC_MSG_ERROR([Ablation library not found.  Try either --with-ablation or setting ABLATION_DIR.])
       else
          AC_MSG_NOTICE([optional Ablation library not found])
          ABLATION_CFLAGS=""   # GRVY_CFLAGS empty on failure
          ABLATION_FCFLAGS=""  # GRVY_FCFLAGS empty on failure
          ABLATION_LIBS=""     # GRVY_LIBS empty on failure
       fi
    else
        HAVE_ABLATION=1
        AC_DEFINE(HAVE_ABLATION,1,[Define if PECOS libablation is available])
        AC_SUBST(ABLATION_CFLAGS)
        AC_SUBST(ABLATION_LIBS)
	AC_SUBST(ABLATION_PREFIX)
    fi

    AC_SUBST(HAVE_ABLATION)

# fi

AM_CONDITIONAL(ABLATION_ENABLED,test x$HAVE_ABLATION = x1)

])
