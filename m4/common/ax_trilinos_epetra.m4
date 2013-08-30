# ===========================================================================
#         http://autoconf-archive.cryp.to/ax_trilinos_epetra.html
# ===========================================================================
#
# SYNOPSIS
#
#   Test for the Trilinos Epetra
#   (http://trilinos.sandia.gov/packages/epetra) library.
#
#   AX_TRILINOS_EPETRA([, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   On success, adds "include Makefile.export.epetra" statements to
#   every Automake file containing @INC_AMINCLUDE@.  Requires that
#   Trilinos was configured with the --enable-export-makefiles option.
#   When ACTION-IF-NOT-FOUND is not specified, the default behavior is
#   for configure to fail.
#
# LAST MODIFICATION
#
#   2008-08-13
#
# COPYLEFT
#
#   Copyright (c) 2008 Rhys Ulerich <rhys.ulerich@gmail.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_TRILINOS_EPETRA],[
    AC_REQUIRE([AX_TRILINOS_BASE])
    ax_trilinos_epetra=yes

    AC_LANG_SAVE
    AC_LANG_CPLUSPLUS

    AC_HAVE_LIBRARY([epetra],[:],[ax_trilinos_epetra=no])

    AC_CHECK_HEADER([Epetra_Vector.h],[:],[ax_trilinos_epetra=no])

    AC_LANG_RESTORE

    ax_epetra=no
    ax_trilinos=no # kemelli: not convinced that we need this variable and/or its test bellow

    if test "$ax_trilinos_epetra" = yes; then
dnl      ax_trilinos_epetra_9=yes
dnl      ax_trilinos_epetra_10=yes
dnl      AX_ADD_AM_TRILINOS_MAKEFILE_EXPORT([epetra.macros],[ax_trilinos_epetra_9=no])
dnl      AX_ADD_AM_TRILINOS_MAKEFILE_EXPORT([Epetra],[ax_trilinos_epetra_9=no])
dnl      AX_ADD_AM_TRILINOS_MAKEFILE_EXPORT([Trilinos],[ax_trilinos_epetra_10=no])

      AX_ADD_AM_TRILINOS_MAKEFILE_EXPORT([Trilinos],[ax_trilinos=no])
      AX_ADD_AM_TRILINOS_MAKEFILE_EXPORT([Epetra],[ax_epetra=no])
    fi

dnl    if test "$ax_trilinos_epetra_9" = yes -o "$ax_trilinos_epetra_10" = yes; then
dnl        : # NOP
dnl		ifelse([$1],,,
dnl            [$1])
dnl    else
dnl        : # NOP
dnl		ifelse([$2],,AC_MSG_ERROR([Trilinos Epetra not usable.]),
dnl            [$2])
dnl    fi
 # Epetra is required case Trilinos is enabled. Throws an error message and exit otherwise.
    if test "$ax_epetra" = "yes" -a "$ax_trilinos" = "yes"; then 
      ifelse([$1],,,[$1])
    else
      ifelse([$2],,AC_MSG_ERROR([Trilinos Epetra not usable. QUESO will not be build with Trilinos]),[$2])      
    fi
])
