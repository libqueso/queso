# ===========================================================================
#                http://autoconf-archive.cryp.to/acx_mkl.html
# ===========================================================================
#
# SYNOPSIS
#
#   Test for Intel(R) Math Kernel Library
#   (http://software.intel.com/en-us/intel-mkl/)
#
#   ACX_MKL([USE_THREADED,[ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]]])
#
# DESCRIPTION
#
#   Provides a --with-mkl=DIR option. Searches --with-mkl, $MKLROOT, and the
#   usual places for MKL headers and libraries.  Attempts to autodetect 32-bit,
#   64-bit, and em64t architecture linking.
#
#   Supports separately specifying --with-mkl-include or --with-mkl-libdir to
#   override default locations underneath either --with-mkl or $MKLROOT.
#   Allows specifying whether or not threaded MKL routines should be used,
#   with a default to use non-threaded MKL routines.
#
#   On success, sets MKL_CFLAGS, MKL_LIBS, and #defines HAVE_MKL.  When
#   ACTION-IF-NOT-FOUND is not specified, the default behavior is for configure
#   to fail.
#
# LAST MODIFICATION
#
#   2009-07-14
#
# COPYLEFT
#
#   Copyright (c) 2009 Rhys Ulerich <rhys.ulerich@gmail.com>
#   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
#   Copyright (c) 2008 Thomas Porschberg <thomas@randspringer.de>
#   Copyright (c) 2008 Caolan McNamara <caolan@skynet.ie>
#   Copyright (c) 2008 Alexandre Duret-Lutz <adl@gnu.org>
#   Copyright (c) 2008 Matthew Mueller <donut@azstarnet.com>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Macro Archive. When you make and
#   distribute a modified version of the Autoconf Macro, you may extend this
#   special exception to the GPL to apply to your modified version as well.

AC_DEFUN([ACX_MKL], [
AC_PREREQ(2.60)
AC_REQUIRE([AC_CANONICAL_TARGET])
AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
AC_REQUIRE([AX_COMPILER_VENDOR])
AC_ARG_VAR(MKLROOT,[root directory of MKL installation])

dnl Note the assumption that lp64 and not ilp64 should be used
dnl Please add entries to the case statement as required
case $target_cpu in
    x86_64)
             acx_mkl_libdirsuffix="em64t"
             acx_mkl_libsuffix="_lp64"
             ;;
    i686)
             acx_mkl_libdirsuffix="32"
             acx_mkl_libsuffix=""
             ;;
    unknown)
             AC_MSG_WARN([Unknown target_cpu; defaulting to 32-bit MKL unless --with-mkl-lib=<DIR> supplied])
             acx_mkl_libdirsuffix="32"
             acx_mkl_libsuffix=""
             ;;
    *)
             AC_MSG_ERROR([m4 macro [$0][ unable to handle target_cpu: $target_cpu]])
             ;;
esac

dnl TODO Add a [default=something] message to AS_HELP_STRING below
AC_MSG_CHECKING(whether to link threaded MKL routines)
acx_mkl_enable_threads_default="m4_tolower(m4_normalize(ifelse([$1],,[no],[$1])))"
AC_ARG_ENABLE(mkl-threads,
    [AS_HELP_STRING([--enable-mkl-threads],
        [enable threaded MKL routines])],
    [acx_mkl_enable_threads=yes],acx_mkl_enable_threads=$acx_mkl_enable_threads_default)
AC_MSG_RESULT($acx_mkl_enable_threads)

dnl Please add entries to the case statement as required
case $ax_cv_c_compiler_vendor in
    intel)
           acx_mkl_interfacelayer="-lmkl_intel${acx_mkl_libsuffix}"
           if test "${acx_mkl_enable_threads}" = "yes"; then
               acx_mkl_threadinglayer="-lmkl_intel_thread"
               acx_mkl_rtllayer="-liomp5 -lpthread"
           else
               acx_mkl_threadinglayer="-lmkl_sequential"
               acx_mkl_rtllayer=""
           fi
           ;;
    gnu)
           acx_mkl_interfacelayer="-lmkl_gf${acx_mkl_libsuffix}"
           if test "${acx_mkl_enable_threads}" = "yes"; then
               acx_mkl_threadinglayer="-lmkl_gnu_thread"
               acx_mkl_rtllayer="-liomp5 -lpthread"
           else
               acx_mkl_threadinglayer="-lmkl_sequential"
               acx_mkl_rtllayer=""
           fi
           ;;
    *)
           AC_MSG_ERROR([m4 macro [$0][ unable to handle ax_cv_c_compiler_vendor: $ax_cv_c_compiler_vendor]])
           ;;
esac

AC_ARG_WITH(mkl, [AS_HELP_STRING([--with-mkl[=DIR]],[root directory of MKL installation])],[
with_mkl=$withval
if test "${with_mkl}" != yes; then
    MKLROOT=$withval
    acx_mkl_include="$withval/include"
    acx_mkl_libdir="$withval/lib/$acx_mkl_libdirsuffix"
fi
],[
with_mkl=$withval
if test "x${MKLROOT}" != "x"; then
    acx_mkl_include="${MKLROOT}/include"
    acx_mkl_libdir="${MKLROOT}/lib/$acx_mkl_libdirsuffix"
fi
])

AC_ARG_WITH(mkl-include,
[AS_HELP_STRING([--with-mkl-include=DIR],[specify exact directory for MKL headers])],[
if test -d "$withval"; then
    acx_mkl_include="$withval"
else
    AC_MSG_ERROR([--with-mkl-include expected directory name])
fi
])

AC_ARG_WITH(mkl-libdir, [AS_HELP_STRING([--with-mkl-libdir=DIR],[specify exact directory for MKL libraries])],[
if test -d "$withval"; then
    acx_mkl_libdir="$withval"
else
    AC_MSG_ERROR([--with-mkl-libdir expected directory name])
fi
])

if test "${with_mkl}" != no ; then
    MKL_LIBS=""
    MKL_LIBS="${MKL_LIBS} ${acx_mkl_interfacelayer}"
    MKL_LIBS="${MKL_LIBS} ${acx_mkl_threadinglayer}"
    MKL_LIBS="${MKL_LIBS} -lmkl_lapack -lmkl_core"
    MKL_LIBS="${MKL_LIBS} ${acx_mkl_rtllayer}"
    MKL_LIBS="${MKL_LIBS} -lm"

    if test -d "${acx_mkl_libdir}" ; then
        MKL_LIBS="-L${acx_mkl_libdir} ${MKL_LIBS}"
    fi
    MKL_CFLAGS=""
    if test -d "${acx_mkl_include}" ; then
        MKL_CFLAGS="-I${acx_mkl_include} ${MKL_CFLAGS}"
    fi

    acx_mkl_save_CFLAGS="$CFLAGS"
    acx_mkl_save_LDFLAGS="$LDFLAGS"
    acx_mkl_save_LIBS="$LIBS"
    CFLAGS="${MKL_CFLAGS} ${CFLAGS}"
    LDFLAGS="${MKL_LIBS} ${LDFLAGS}"
    AC_LANG_PUSH([C])
    AC_MSG_NOTICE([Ensuring we can use Intel MKL routines using a known link line])
    AC_CHECK_HEADER([mkl.h],[acx_mkl_found_header=yes],[acx_mkl_found_header=no])
    AC_SEARCH_LIBS(MKLGetVersion,[],[acx_mkl_found_library=yes],[acx_mkl_found_library=no])
    AC_SEARCH_LIBS(dgemm,[],[acx_mkl_found_blas=yes],[acx_mkl_found_blas=no])
    AC_SEARCH_LIBS(dgbtrf,[],[acx_mkl_found_lapack=yes],[acx_mkl_found_lapack=no])
    AC_LANG_POP([C])
    LDFLAGS="$acx_mkl_save_LDFLAGS"
    CFLAGS="$acx_mkl_save_CFLAGS"
    LIBS="$acx_mkl_save_LIBS"

    acx_mkl_succeeded=no
    if test "$acx_mkl_found_header" = yes; then
        if test "$acx_mkl_found_library" = yes; then
            if test "$acx_mkl_found_blas" = yes; then
                if test "$acx_mkl_found_lapack" = yes; then
                    acx_mkl_succeeded=yes
                fi
            fi
        fi
    fi

    if test "$acx_mkl_succeeded" = no; then
        ifelse([$2],,AC_MSG_ERROR([Intel MKL not found.  Try either --with-mkl or setting MKLROOT.]),
            [$2])
    else
        AC_DEFINE(HAVE_MKL,1,[Define if MKL is available])
        AC_SUBST(MKL_CFLAGS)
        AC_SUBST(MKL_LIBS)
        ifelse([$1],,,[$1])
    fi

fi
])dnl ACX_MKL
