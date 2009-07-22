# SYNOPSIS
#
#   Queries configuration environment.
#
#   AX_SUMMARIZE_ENV([, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
#
# DESCRIPTION
#
#   Queries compile environment for use in configure summary 
#   and versioning information.
#
# LAST MODIFICATION
#
#   2009-07-16
#

AC_DEFUN([AX_SUMMARIZE_ENV],
[

AC_CANONICAL_HOST

BUILD_USER=${USER}
BUILD_ARCH=${host}
BUILD_HOST=${ac_hostname}
BUILD_DATE=`date +'%F %H:%M'`
BUILD_REV=`svn info | grep "^Revision" | cut -d ' ' -f 2`
BUILD_UP2DATE=`svn status -q | wc -l`

AC_DEFINE_UNQUOTED([BUILD_USER],"${BUILD_USER}",[The fine user who built the library])
AC_DEFINE_UNQUOTED([BUILD_ARCH],"${BUILD_ARCH}",[Architecture of the build host])
AC_DEFINE_UNQUOTED([BUILD_HOST],"${BUILD_HOST}",[Build host])
AC_DEFINE_UNQUOTED([BUILD_DATE],"${BUILD_DATE}",[Build date])
AC_DEFINE_UNQUOTED([BUILD_REV], "${BUILD_REV}", [Build revision])

if test \( $BUILD_UP2DATE -eq 0 \);then
   BUILD_STATUS="current"
else
   BUILD_STATUS="local modifications made"
fi

AC_DEFINE_UNQUOTED([BUILD_STATUS],"${BUILD_STATUS}",[Build update status])

])
