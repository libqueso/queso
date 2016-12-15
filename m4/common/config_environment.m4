# SYNOPSIS
#
#   Queries configuration environment.
#
#   AX_SUMMARIZE_ENV([, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
#
# DESCRIPTION
#
#   Queries compile environment and git revision for use in configure summary 
#   and pre-processing macros.
#
# COPYLEFT
#
#   Copyright (c) 2010 Karl W. Schulz <karl@ices.utexas.edu>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_SUMMARIZE_ENV],
[

AC_CANONICAL_HOST

BUILD_USER=${USER}
BUILD_ARCH=${host}
BUILD_HOST=${ac_hostname}
BUILD_DATE=`date +'%F %H:%M'`

# Determine method for querying Source code revisioning (assumes git)

AC_PATH_PROG(gitbin,git)

if test "x${gitbin}" != "x" && test "`${gitbin} rev-parse --is-inside-work-tree`" = "true"; then
  GIT_REVISION=`${gitbin} rev-parse --short HEAD`
  BUILD_DEVSTATUS="Development Build"
else
  GIT_REVISION="N/A"
  BUILD_DEVSTATUS="External Release"
fi


AC_SUBST(GIT_REVISION)
AC_SUBST(BUILD_DEVSTATUS)

# Query current version.

BUILD_VERSION=${GIT_REVISION}

# Versioning info - check local developer version (if cloned)

AC_DEFINE_UNQUOTED([BUILD_USER],     "${BUILD_USER}",     [The fine user who built the package])
AC_DEFINE_UNQUOTED([BUILD_ARCH],     "${BUILD_ARCH}",     [Architecture of the build host])
AC_DEFINE_UNQUOTED([BUILD_HOST],     "${BUILD_HOST}",     [Build host name])
AC_DEFINE_UNQUOTED([BUILD_VERSION],  "${BUILD_VERSION}",  [git revision])
AC_DEFINE_UNQUOTED([BUILD_DEVSTATUS],"${BUILD_DEVSTATUS}",[Dev/Release build])
AC_DEFINE(         [BUILD_DATE],     __DATE__ " " __TIME__, [Build date])

# Defines needed by Environment.C (formerly in queso.h)
AC_DEFINE_UNQUOTED([CXX],          "${CXX}",           [C++ compiler])
AC_DEFINE_UNQUOTED([CXXFLAGS],     "${CXXFLAGS}",      [C++ compiler flags])
# The following must be unconditionally defined per Environment.C
# TODO: Should these go in the respective <tpl>.m4?
AC_DEFINE_UNQUOTED([TRILINOS_DIR], "${TRILINOS_HOME}", [Trilinos home directory])
AC_DEFINE_UNQUOTED([GSL_DIR],      "${GSL_LIBS}",      [GSL libraries])
AC_DEFINE_UNQUOTED([GRVY_DIR],     "${GRVY_PREFIX}",   [GRVY directory])
AC_DEFINE_UNQUOTED([GLPK_DIR],     "${GLPK_PREFIX}",   [GLPK directory])
AC_DEFINE_UNQUOTED([HDF5_DIR],     "${HDF5_PREFIX}",   [HDF directory])

AC_SUBST(BUILD_USER)
AC_SUBST(BUILD_ARCH)
AC_SUBST(BUILD_HOST)
AC_SUBST(BUILD_DATE)
AC_SUBST(BUILD_VERSION)

# Defines needed by Environment.C (formerly in queso.h)
AC_SUBST(CXX)
AC_SUBST(CXX_FLAGS)

])
