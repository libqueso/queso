# SYNOPSIS
#
#   Test for the GNU Linear Programming Kit (GLPK)
#
#   AM_PATH_GLPK([, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
#
# DESCRIPTION
#
#   Provides a --with-glpk=DIR option. Searches --with-glpk,
#   $GLPK_DIR, and the usual places for GLPK headers and libraries.
#
#   On success, sets GLPK_CFLAGS, GLPK_LIBS, and
#   #defines HAVE_GLPK.  When ACTION-IF-NOT-FOUND is not specified,
#   the default behavior is for configure to fail.
#
# LAST MODIFICATION
#
#   2010-02-24
#
# COPYLEFT
#
#   Copyright (c) 2010 Karl W. Schulz <karl@ices.utexas.edu>
#   Copyright (c) 2010 Ernesto Prudenci <prudenci@ices.utexas.edu>
#   Copyright (c) 2009 Rhys Ulerich <rhys.ulerich@gmail.com>
#   Copyright (c) 2008 Thomas Porschberg <thomas@randspringer.de>
#   Copyright (c) 2008 Caolan McNamara <caolan@skynet.ie>
#   Copyright (c) 2008 Alexandre Duret-Lutz <adl@gnu.org>
#   Copyright (c) 2008 Matthew Mueller <donut@azstarnet.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_PATH_GLPK_OLD],
[

AC_ARG_VAR(GLPK_DIR,[root directory of GLPK installation])

AC_ARG_WITH(glpk, 
  [AS_HELP_STRING([--with-glpk[=DIR]],[root directory of GLPK installation (default = GLPK_DIR)])],
  [with_glpk=$withval
if test "${with_glpk}" != yes; then
    GLPK_PREFIX=$withval
fi
],[
with_glpk=$withval
if test "x${GLPK_DIR}" != "x"; then
   GLPK_PREFIX=${GLPK_DIR}
fi
])

if test "${with_glpk}" != no ; then

    if test -d "${GLPK_PREFIX}/lib" ; then
       GLPK_LIBS="-L${GLPK_PREFIX}/lib -lglpk"
    fi

    if test -d "${GLPK_PREFIX}/include" ; then
        GLPK_CFLAGS="-I${GLPK_PREFIX}/include"
    fi

    ac_GLPK_save_CFLAGS="$CFLAGS"
    ac_GLPK_save_CPPFLAGS="$CPPFLAGS"
    ac_GLPK_save_LDFLAGS="$LDFLAGS"
    ac_GLPK_save_LIBS="$LIBS"

    CFLAGS="${GLPK_CFLAGS} ${CFLAGS}"
    CPPFLAGS="${GLPK_CFLAGS} ${CPPFLAGS}"
    LDFLAGS="${GLPK_LIBS} ${LDFLAGS}"
    AC_LANG_PUSH([C])
    AC_CHECK_HEADER([glpk.h],[found_header=yes],[found_header=no])


    #-----------------------
    # Minimum version check
    #----------------------

    min_glpk_version=ifelse([$1], ,4.35,$1)	

    # koomie note: although GLPK appears to follow a major.minor
    # versioning scheme, we will also look for the possibility of a
    # micro version from the autoconf query to keep the sed queries
    # reasonably protected.

    ver_form=`echo $min_glpk_version | sed 's/[[0-9]]*\.[[0-9]]*/major-minor/'`	

    if test $ver_form = "major-minor" ; then
        MIN_GLPK_MAJOR_VERSION=`echo $min_glpk_version | sed 's/\([[0-9]]*\)\.\([[0-9]]*\)/\1/'`
        MIN_GLPK_MINOR_VERSION=`echo $min_glpk_version | sed 's/\([[0-9]]*\)\.\([[0-9]]*\)/\2/'`
    else
        ver_form=`echo $min_glpk_version | sed 's/[0-9]*\.[0-9]*\.[0-9]*/major-minor-micro/'`

	if test $ver_form = "major-minor-micro" ; then
           MIN_GLPK_MAJOR_VERSION=`echo $min_glpk_version | sed 's/\([[0-9]]*\)\.\([[0-9]]*\)\.[[0-9]]*/\1/'`	
           MIN_GLPK_MINOR_VERSION=`echo $min_glpk_version | sed 's/\([[0-9]]*\)\.\([[0-9]]*\)\.[[0-9]]*/\2/'`	
        else				
           AC_MSG_ERROR([Unable to parse desired GLPK version string correctly.])
        fi			
    fi	

    AC_MSG_CHECKING(for glpk - version >= $MIN_GLPK_MAJOR_VERSION.$MIN_GLPK_MINOR_VERSION)

    succeeded=no
    AC_LANG_PUSH([C])

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
       @%:@include <glpk.h>
            ]], [[
	    #if GLP_MAJOR_VERSION > $MIN_GLPK_MAJOR_VERSION
	    /* Sweet nibblets */
       	    #elif (GLP_MAJOR_VERSION >= $MIN_GLPK_MAJOR_VERSION) && (GLP_MINOR_VERSION >= $MIN_GLPK_MINOR_VERSION)
            /* Winner winner, chicken dinner */
            #else
            #  error GLPK version is too old
            #endif
        ]])],[
            AC_MSG_RESULT(yes)
            succeeded=yes
        ],[
            AC_MSG_RESULT(no)
        ])

    AC_LANG_POP([C])

    if test "$succeeded" != "yes";then
       AC_MSG_ERROR([

       Your GNU Linear Programming Kit (GLPK) version does not meet the minimum 
       versioning requirements.  Please use --with-glpk to specify the location of
       an updated installation or consider upgrading the system version.

       ]) 
    fi     

    # Test library linkage 

    AC_CHECK_LIB([glpk],glp_create_prob,[found_library=yes],[found_library=no])


    CFLAGS="$ac_GLPK_save_CFLAGS"
    CPPFLAGS="$ac_GLPK_save_CPPFLAGS"
    LDFLAGS="$ac_GLPK_save_LDFLAGS"
    LIBS="$ac_GLPK_save_LIBS"

    succeeded=no
    if test "$found_header" = yes; then
        if test "$found_library" = yes; then
            succeeded=yes
        fi
    fi

    if test "$succeeded" = no; then
        ifelse([$2],,AC_MSG_ERROR([GLPK not found.  Try either --with-glpk or setting GLPK_DIR.]),
            [$2])
    else
        AC_DEFINE(HAVE_GLPK,1,[Define if GLPK is available])
        AC_SUBST(GLPK_CFLAGS)
        AC_SUBST(GLPK_LIBS)
	AC_SUBST(GLPK_PREFIX)
#####        ifelse([$1],,,[$1])
    fi

fi

])
