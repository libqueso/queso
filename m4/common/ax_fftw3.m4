#
#   AX_FFTW3([FFTW3_MODULES])
#   where FFTW3_MODULES is a PKG_CHECK_MODULES string like [fftw3 > 3.2].
#
AC_DEFUN([AX_FFTW3], [
AC_PREREQ(2.60)
AC_REQUIRE([AC_PROG_SED])


AC_ARG_WITH(fftw3, 
  [AS_HELP_STRING([--with-fftw3[=DIR]],[root directory of fftw3 installation (default = FFTW3_DIR)])],
  [with_fftw3=$withval
if test "${with_fftw3}" != yes; then
    FFTW3_PREFIX=$withval
fi
],[
with_fftw3=$withval
if test "x${FFTW3_DIR}" != "x"; then
   FFTW3_PREFIX=${FFTW3_DIR}
fi
])

dnl echo "incoming config path = $PKG_CONFIG_PATH"
ac_PKG_CONFIG_PATH_save="$PKG_CONFIG_PATH"
export PKG_CONFIG_PATH="${FFTW3_PREFIX}/lib/pkgconfig:$PKG_CONFIG_PATH"
dnl echo "updated config path = $PKG_CONFIG_PATH"

AC_LANG_PUSH([C])

# Search for FFTW itself using pkg-config
PKG_CHECK_MODULES(FFTW3,m4_default([$1],[fftw3 > 3.2]),,AC_MSG_ERROR([Valid FFTW3 library not found. Try either --with-fftw3 or setting FFTW3_DIR.]))
AC_SUBST(FFTW3_CFLAGS)
AC_SUBST(FFTW3_LIBS)

# Prepare for FFTW threading search
ax_fftw3_libs_pre=`echo  $FFTW3_LIBS | $SED -e 's/-lfftw3.*$//'`
ax_fftw3_libs_post=`echo $FFTW3_LIBS | $SED -e 's/^.*-lfftw3/-lfftw3/'`
# Default values overridden on success
FFTW3_THREADS_CFLAGS=$FFTW3_CFLAGS
FFTW3_THREADS_LIBS=$FFTW3_LIBS
ax_fftw3_found_threads=no

# Check if we can link in threading without any additional help
AC_CHECK_LIB([fftw3_threads],[fftw_init_threads],[
    ax_fftw3_found_threads=yes 
    FFTW3_THREADS_CFLAGS="$OPENMP_CFLAGS $FFTW3_CFLAGS"
    FFTW3_THREADS_LIBS="$ax_fftw3_libs_pre -lfftw3_threads $ax_fftw3_libs_post"
    AC_MSG_NOTICE([Detected FFTW threading not requiring any external dependencies])
])

# Check if OpenMP is the way to go...
if test "${ax_fftw3_found_threads}" != yes; then
    ax_fftw3_save_LDFLAGS=$LDFLAGS
    LDFLAGS="$LDFLAGS $ax_fftw3_libs_pre"
    AX_OPENMP([
        AC_DEFINE(HAVE_OPENMP,[1],[Define if OpenMP is enabled])
        # Note the new function name to thwart caching from earlier AC_CHECK_LIB
        AC_CHECK_LIB([fftw3_threads],[fftw_cleanup_threads], [
                ax_fftw3_found_threads=yes
                FFTW3_THREADS_CFLAGS="$OPENMP_CFLAGS $FFTW3_CFLAGS"
                FFTW3_THREADS_LIBS="$ax_fftw3_libs_pre -lfftw3_threads $ax_fftw3_libs_post $OPENMP_CFLAGS"
                AC_MSG_NOTICE([Detected OpenMP-enabled FFTW threading])
            ], [
                AC_MSG_NOTICE([Detected OpenMP but not OpenMP-enabled FFTW threading])
            ], [$ax_fftw3_libs_post $OPENMP_CFLAGS]
        )
    ])
    LDFLAGS=$ax_fftw3_save_LDFLAGS
fi

# Check if POSIX Threads is the way to go...
if test "${ax_fftw3_found_threads}" != yes; then
    ax_fftw3_save_LDFLAGS=$LDFLAGS
    LDFLAGS="$LDFLAGS $ax_fftw3_libs_pre"
    ACX_PTHREAD([
        AC_DEFINE(HAVE_PTHREAD,[1],
                [Define if you have POSIX threads libraries and header files])
        # Note the new function name to thwart caching from earlier AC_CHECK_LIB
        AC_CHECK_LIB([fftw3_threads],[fftw_plan_with_nthreads], [
            ax_fftw3_found_threads=yes
            FFTW3_THREADS_CFLAGS="$PTHREAD_CFLAGS $FFTW3_CFLAGS"
            FFTW3_THREADS_LIBS="$ax_fftw3_libs_pre -lfftw3_threads $ax_fftw3_libs_post $PTHREAD_LIBS"
            AC_MSG_NOTICE([Detected pthread-enabled FFTW threading])
        ], [
            AC_MSG_NOTICE([Detected pthread but not pthread-enabled FFTW threading])
        ], [$ax_fftw3_libs_post $PTHREAD_LIBS])
    ])
    LDFLAGS=$ax_fftw3_save_LDFLAGS
fi

# Finish up FFTW thread-related processing
if test "${ax_fftw3_found_threads}" = yes; then
    AC_DEFINE([HAVE_FFTW3_THREADS],[1],[Defined if FFTW3 threads available])
fi
AC_SUBST(FFTW3_THREADS_CFLAGS)
AC_SUBST(FFTW3_THREADS_LIBS)

#
# 8/17/10 - karl disabling MPI checks for use with compDNS. Would need ability to define
# as optional to use with both suzerain/compDNS.
#

# Search for FFTW MPI independent of FFTW threading result

# ax_fftw3_found_mpi=no
# AC_REQUIRE([ACX_MPI])
# ax_fftw3_save_CC=$CC
# CC=$MPICC
# ax_fftw3_save_LDFLAGS=$LDFLAGS
# LDFLAGS="$LDFLAGS $ax_fftw3_libs_pre"
# AC_CHECK_LIB([fftw3_mpi],[fftw_mpi_init], [
#     ax_fftw3_found_mpi=yes
#     AC_DEFINE([HAVE_FFTW3_MPI],[1],[Defined if FFTW3 MPI available])

#     # Prepare libraries for FFTW MPI only
#     FFTW3_MPI_LIBS="$ax_fftw3_libs_pre -lfftw3_mpi $ax_fftw3_libs_post"
#     AC_SUBST(FFTW3_MPI_LIBS)

#     # Prepare libraries for FFTW MPI and threading
#     if test "${ax_fftw3_found_threads}" = yes; then
#         ax_fftw3_threads_libs_pre=`echo  $FFTW3_THREADS_LIBS | $SED -e 's/-lfftw3_threads.*$//'`
#         ax_fftw3_threads_libs_post=`echo $FFTW3_THREADS_LIBS | $SED -e 's/^.*-lfftw3_threads/-lfftw3_threads/'`
#         FFTW3_MPI_THREADS_LIBS="$ax_fftw3_threads_libs_pre -lfftw3_mpi $ax_fftw3_threads_libs_post"
#     else
#         FFTW3_MPI_THREADS_LIBS="$ax_fftw3_libs_pre -lfftw3_mpi $ax_fftw3_libs_post"
#     fi
#     AC_SUBST(FFTW3_MPI_THREADS_LIBS)
# ], [
#     AC_MSG_ERROR([Could not find FFTW MPI functionality])
# ], [$ax_fftw3_libs_post])
# LDFLAGS=$ax_fftw3_save_LDFLAGS
# CC=$ax_fftw3_save_CC

export PKG_CONFIG_PATH="$ac_PKG_CONFIG_PATH_save"

AC_LANG_POP([C])
])dnl AX_FFTW3
