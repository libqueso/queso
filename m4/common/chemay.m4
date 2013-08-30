# SYNOPSIS
#
#   Test for libchemay
#
#   AM_PATH_CHEMAY( <Minimum Required Version>, <package-required=yes/no> )
#
# DESCRIPTION
#
#   Provides a --with-chemay=DIR option. Searches --with-chemay,
#   $CHEMAY_DIR, and the usual places for Chemay headers and libraries.
#
#   On success, sets CHEMAY_CFLAGS, CHEMAY_LIBS, and #defines HAVE_CHEMAY.
#   Also defines automake conditional CHEMAY_ENABLED.  Assumes package
#   is optional unless overridden with $2=yes.
#
# LAST MODIFICATION
#
#   $Id: chemay.m4 18587 2011-03-14 23:29:42Z roystgnr $
#
# COPYLEFT
#
#   Copyright (c) 2010 Roy H. Stogner <roystgnr@ices.utexas.edu>
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

AC_DEFUN([AX_PATH_CHEMAY],
[

AC_ARG_VAR(CHEMAY_DIR,[root directory of Chemay installation])

AC_ARG_WITH(chemay, 
  [AS_HELP_STRING([--with-chemay[=DIR]],[root directory of CHEMAY installation (default = CHEMAY_DIR)])],
  [with_chemay=$withval
if test "${with_chemay}" != yes; then
    CHEMAY_PREFIX=$withval
fi
],[
with_chemay=$withval
if test "x${CHEMAY_DIR}" != "x"; then
   CHEMAY_PREFIX=${CHEMAY_DIR}
fi
])

# package requirement; if not specified, the default is to assume that
# the package is optional

is_package_required=ifelse([$2], ,no, $2 )

HAVE_CHEMAY=0

if test "${with_chemay}" != no ; then

    BOOST_SMART_PTR

    AX_PATH_GSL

    if test ! -d "${CHEMAY_PREFIX}/" ; then
       AC_MSG_NOTICE([Chemay directory ${CHEMAY_PREFIX} not found])	
    fi

    if test -d "${CHEMAY_PREFIX}/lib" ; then
       CHEMAY_LIBS="-L${CHEMAY_PREFIX}/lib -lchemay -lcantera_compare -ltinyxml -Wl,-rpath,${CHEMAY_PREFIX}/lib ${GSL_LIBS}"
       CHEMAY_FCFLAGS="-I${CHEMAY_PREFIX}/lib"
    else
       AC_MSG_NOTICE([Chemay library directory ${CHEMAY_PREFIX}/lib not found])	
    fi

    if test -d "${CHEMAY_PREFIX}/include" ; then
        CHEMAY_CPPFLAGS="-I${CHEMAY_PREFIX}/include ${BOOST_CPPFLAGS}"
    else
       AC_MSG_NOTICE([Chemay include directory ${CHEMAY_PREFIX}/include not found])	
    fi

    ac_CHEMAY_save_CFLAGS="$CFLAGS"
    ac_CHEMAY_save_CPPFLAGS="$CPPFLAGS"
    ac_CHEMAY_save_LDFLAGS="$LDFLAGS"
    ac_CHEMAY_save_LIBS="$LIBS"

    CPPFLAGS="${CHEMAY_CPPFLAGS} ${CPPFLAGS}"
    LDFLAGS="${CHEMAY_LIBS} ${LDFLAGS}"
    AC_LANG_PUSH([C++])
    AC_CHECK_HEADER([chemay/Arrhenius.hh],[found_header=yes],[found_header=no])

dnl    #-----------------------
dnl    # Minimum version check
dnl    #----------------------
dnl
dnl    min_chemay_version=ifelse([$1], ,0.29, $1)	
dnl
dnl    # looking for major.minor.micro style versioning
dnl
dnl    MAJOR_VER=`echo $min_chemay_version | sed 's/^\([[0-9]]*\).*/\1/'`
dnl    if test "x${MAJOR_VER}" = "x" ; then
dnl       MAJOR_VER=0
dnl    fi
dnl
dnl    MINOR_VER=`echo $min_chemay_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\2/'`
dnl    if test "x${MINOR_VER}" = "x" ; then
dnl       MINOR_VER=0
dnl    fi
dnl
dnl    MICRO_VER=`echo $min_chemay_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\3/'`
dnl    if test "x${MICRO_VER}" = "x" ; then
dnl       MICRO_VER=0
dnl    fi
dnl
dnl    # begin additional test(s) if header if available
dnl
dnl    if test "x${found_header}" = "xyes" ; then
dnl
dnl        AC_MSG_CHECKING(for Chemay - version >= $min_chemay_version)
dnl        version_succeeded=no
dnl
dnl    	AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
dnl       	@%:@include <chemay/Arrhenius.hh>
dnl            ]], [[
dnl            #if CHEMAY_MAJOR_VERSION > $MAJOR_VER
dnl            /* Sweet nibblets */
dnl            #elif (CHEMAY_MAJOR_VERSION >= $MAJOR_VER) && (CHEMAY_MINOR_VERSION > $MINOR_VER)
dnl            /* Winner winner, chicken dinner */
dnl            #elif (CHEMAY_MAJOR_VERSION >= $MAJOR_VER) && (CHEMAY_MINOR_VERSION >= $MINOR_VER) && (CHEMAY_MICRO_VERSION >= $MICRO_VER)
dnl            /* I feel like chicken tonight, like chicken tonight? */
dnl            #else
dnl            #  error version is too old
dnl            #endif
dnl        ]])],[
dnl            AC_MSG_RESULT(yes)
dnl            version_succeeded=yes
dnl        ],[
dnl            AC_MSG_RESULT(no)
dnl        ])
dnl
dnl    if test "$version_succeeded" != "yes";then
dnl       if test "$is_package_required" = yes; then	
dnl       	  AC_MSG_ERROR([
dnl
dnl   Your Chemay library version does not meet the minimum versioning
dnl   requirements ($min_chemay_version).  Please use --with-chemay to specify the location
dnl   of an updated installation or consider upgrading the system version.
dnl
dnl          ]) 
dnl       fi 
dnl    fi     

    # Library availability

    AC_MSG_CHECKING([for -lchemay linkage])

    AC_LINK_IFELSE(
      [AC_LANG_PROGRAM([#include <chemay/Arrhenius.hh>],
        [Chemay::Arrhenius dummy])],
      [AC_MSG_RESULT(yes)]
      [found_library=yes],
      [AC_MSG_RESULT(no)]
      [AC_MSG_WARN([Failed to link to libchemay.])]
      [found_library=no])

dnl    AC_CHECK_LIB([chemay],chemay_input_fopen,[found_library=yes],[found_library=no])

dnl fi   dnl end test if header if available

    AC_LANG_POP([C++])

    CFLAGS="$ac_CHEMAY_save_CFLAGS"
    CPPFLAGS="$ac_CHEMAY_save_CPPFLAGS"
    LDFLAGS="$ac_CHEMAY_save_LDFLAGS"
    LIBS="$ac_CHEMAY_save_LIBS"

    succeeded=no
    if test "$found_header" = yes; then
dnl        if test "$version_succeeded" = yes; then
	   if test "$found_library" = yes; then
              succeeded=yes
           fi	
dnl        fi
    fi

    if test "$succeeded" = no; then
        CHEMAY_CPPFLAGS=""
        CHEMAY_FCFLAGS=""
        CHEMAY_LIBS=""
        CHEMAY_PREFIX=""
        if test "$is_package_required" = yes; then
       	    AC_MSG_ERROR([libchemay not found.  Try either --with-chemay or setting CHEMAY_DIR.])
        else			
       	    AC_MSG_NOTICE([optional Chemay library not found])	
        fi
    else
        HAVE_CHEMAY=1
        AC_DEFINE(HAVE_CHEMAY,1,[Define if Chemay is available])
    fi

    AC_SUBST(CHEMAY_CPPFLAGS)
    AC_SUBST(CHEMAY_FCFLAGS)
    AC_SUBST(CHEMAY_LIBS)
    AC_SUBST(CHEMAY_PREFIX)
    AC_SUBST(HAVE_CHEMAY)
fi

AM_CONDITIONAL(CHEMAY_ENABLED,test x$HAVE_CHEMAY = x1)

])
