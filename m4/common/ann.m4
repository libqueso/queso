# SYNOPSIS
#
#   Test for  
#
#   AM_ENABLE_ANN([, ACTION-IF-GIVEN [, ACTION-IF-NOT-GIVEN]]])
#
# DESCRIPTION
#
#   Provides a --enable-ann=[yes|no] option. 
#
#   On success, sets ANN_CFLAGS, ANN_LIBS, and
#   #defines HAVE_ANN.
#
# LAST MODIFICATION
#
#   2011-02-08 by Gabriel Terejanu
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_ENABLE_ANN],
[

AC_MSG_CHECKING([whether to enable ANN])

HAVE_ANN=0

AC_ARG_ENABLE(ann, 
  [AS_HELP_STRING([--enable-ann[=yes|no]],[enable ANN installation (default = no)])],
  [if test "$enableval" = "yes" ; then
      enable_ann=yes
  else
      enable_ann=no
  fi],
  [enable_ann=no]
)

# force auto-make to build ANN
AM_CONDITIONAL(HAVE_ANN, test x$enable_ann = xyes)

# define it in config.h
if test "${enable_ann}" == yes ; then

      AC_DEFINE(HAVE_ANN,1,[Define if ANN is available])

      HAVE_ANN=1

      ANN_PREFIX="\$(top_srcdir)/src/contrib/ANN"
      ANN_LIBS="-L${ANN_PREFIX}/lib -lANN"
      ANN_CFLAGS="-I${ANN_PREFIX}/include"
      ANN_CPPFLAGS="-I${ANN_PREFIX}/include"

else

      ANN_CFLAGS=""
      ANN_CPPFLAGS=""
      ANN_LIBS=""
      ANN_PREFIX=""

fi

# create flags
AC_SUBST(ANN_CFLAGS)
AC_SUBST(ANN_CPPFLAGS)
AC_SUBST(ANN_LIBS)
AC_SUBST(ANN_PREFIX)

# end
AC_MSG_RESULT([$enable_ann])

])
