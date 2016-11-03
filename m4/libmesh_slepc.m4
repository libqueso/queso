# SYNOPSIS
#
#   Test for LIBMESH
#
#   AM_LIBMESH_SLEPC( <Minimum Required Version>, <package-required=yes/no> )
#
# DESCRIPTION
#
#   Uses AM_PATH_LIBMESH_NEW to provides a --with-libmesh=DIR option.
#   Searches --with-libmesh, $LIBMESH_DIR, and the usual places for
#   LIBMESH headers and libraries.
#
#   Fails not only if no libMesh installation is found, but also if
#   the first libMesh installation discovered has not been configured
#   with SLEPc enabled.
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

#   Copyright (c) 2016 Roy H. Stogner <roystgnr@ices.utexas.edu>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_LIBMESH_SLEPC],
[

AX_PATH_LIBMESH_NEW($1, $2)

AC_ARG_VAR(LIBMESH_DIR,[root directory of LIBMESH installation])

HAVE_LIBMESH_SLEPC=0

if test "x${HAVE_LIBMESH}" = "x1"; then

  CPPFLAGS="${LIBMESH_CPPFLAGS} ${CPPFLAGS}"
  LDFLAGS="${LIBMESH_LDFLAGS} ${LDFLAGS}"
  LIBS="${LIBMESH_LIBS} ${LIBS}"

  AC_MSG_CHECKING(for libMesh SLEPc support)

  AC_LANG_PUSH([C++])
  AC_LINK_IFELSE( [AC_LANG_PROGRAM([#include "libmesh/slepc_eigen_solver.h"],
                                   [libMesh::LibMeshInit init (0, 0);
                                    libMesh::SlepcEigenSolver<libMesh::Number>(init.comm());])],
                                   [AC_MSG_RESULT(yes)
                                    HAVE_LIBMESH_SLEPC=1],
                                   [AC_MSG_RESULT(no) 
                                    HAVE_LIBMESH_SLEPC=0] )
  AC_LANG_POP([C++])
fi
if test "x${HAVE_LIBMESH_SLEPC}" = "x1"; then
  AC_DEFINE(HAVE_LIBMESH_SLEPC,1,[Define if libMesh with SLEPc support is available])
fi
AC_SUBST(HAVE_LIBMESH_SLEPC)
AM_CONDITIONAL(LIBMESH_SLEPC_ENABLED,test x$HAVE_LIBMESH_SLEPC = x1)
])
