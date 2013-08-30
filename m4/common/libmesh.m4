# SYNOPSIS
#
#   Test for LIBMESH
#
#   AM_PATH_LIBMESH([, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
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
# COPYLEFT
#
#   Copyright (c) 2009 Karl W. Schulz <karl@ices.utexas.edu>
#   Copyright (c) 2009 Rhys Ulerich <rhys.ulerich@gmail.com>
#   Copyright (c) 2008 Thomas Porschberg <thomas@randspringer.de>
#   Copyright (c) 2008 Caolan McNamara <caolan@skynet.ie>
#   Copyright (c) 2008 Alexandre Duret-Lutz <adl@gnu.org>
#   Copyright (c) 2008 Matthew Mueller <donut@azstarnet.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_PATH_LIBMESH],
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


if test "${with_libmesh}" != no ; then

   DEFAULT_SYS_LIB=/usr/lib


   AC_CHECK_FILE($LIBMESH_PREFIX/contrib/bin/libmesh-config,
                 LIBMESH_CONFIG=$LIBMESH_PREFIX/contrib/bin/libmesh-config,
                 [AC_CHECK_FILE($LIBMESH_PREFIX/bin/libmesh-config,
                                LIBMESH_CONFIG=$LIBMESH_PREFIX/bin/libmesh-config,
                                AC_MSG_ERROR([Cannot find libmesh-config!]))
                 ])


    if test -x "${LIBMESH_CONFIG}" ; then
	LIBMESH_CXXFLAGS=`${LIBMESH_CONFIG} --include`
        LIBMESH_INCLUDE=`${LIBMESH_CONFIG} --include`

	# A tiny bit of tom'foolery here. Pop the default system lib directory to the end to avoid
	# problems when non-default libs are used (eg. C++ libs).  We will add it at the end after the fact.
	LIBMESH_LIBS=`${LIBMESH_CONFIG} --ldflags | sed s#"\-Wl,\-rpath,${DEFAULT_SYS_LIB}\ "##g `
#        echo libmesh libs : $LIBMESH_LIBS
	LIBMESH_LIBS="${LIBMESH_LIBS} -Wl,-rpath,${DEFAULT_SYS_LIB}"
    fi

    ac_LIBMESH_save_CXXFLAGS="$CXXFLAGS"
    ac_LIBMESH_save_CPPFLAGS="$CPPFLAGS"
    ac_LIBMESH_save_LDFLAGS="$LDFLAGS"
#    ac_LIBMESH_save_LIBS="$LIBS"

    CFLAGS="${LIBMESH_CXXFLAGS} ${CXXFLAGS}"
    CPPFLAGS="${LIBMESH_CXXFLAGS} ${CPPFLAGS}"
    LDFLAGS="${LIBMESH_LIBS} ${LDFLAGS}"

     AC_LANG_PUSH([C++])
#     AC_CHECK_HEADER([libmesh.h],[found_header=yes],[found_header=no])
#     AC_CHECK_LIB([libmesh],libmesh_input_fopen,[found_library=yes],[found_library=no])
     AC_LANG_POP([C++])

     CXXFLAGS="$ac_LIBMESH_save_CXXFLAGS"
     CPPFLAGS="$ac_LIBMESH_save_CPPFLAGS"
     LDFLAGS="$ac_LIBMESH_save_LDFLAGS"
#     LIBS="$ac_LIBMESH_save_LIBS"

#     succeeded=no
#     if test "$found_header" = yes; then
#         if test "$found_library" = yes; then
#             succeeded=yes
#         fi
#     fi

#     if test "$succeeded" = no; then
#         ifelse([$2],,AC_MSG_ERROR([LIBMESH not found.  Try either --with-libmesh or setting LIBMESH_DIR.]),
#             [$2])
#     else
        AC_DEFINE(HAVE_LIBMESH,1,[Define if LIBMESH is available])
	AC_SUBST(LIBMESH_CONFIG)
        AC_SUBST(LIBMESH_CXXFLAGS)
        AC_SUBST(LIBMESH_INCLUDE)
        AC_SUBST(LIBMESH_LIBS)
	AC_SUBST(LIBMESH_PREFIX)
        ifelse([$1],,,[$1])
#    fi

fi

])
