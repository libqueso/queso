dnl Init autoconf w/ Package name, version number, and developer email
AC_INIT([queso], [0.4.0], [pecos_dev@ices.utexas.edu])

dnl configuration header
AM_CONFIG_HEADER(config.h)

AC_CONFIG_AUX_DIR([build-aux])

dnl Init automake
AM_INIT_AUTOMAKE

dnl Generating a libtool library
dnl AC_PROG_LIBTOOL

dnl Required for some macros
AX_AM_MACROS

dnl Get compilers
AC_PROG_CC
echo "configure found CC = $CC"

AC_PROG_CXX
echo "configure found CXX = $CXX"

AC_PROG_RANLIB

dnl Check for mpi c and c++ compilers
AC_LANG([C])
ACX_MPI([CC="$MPICC"], [AC_MSG_ERROR([Could not find MPI.])])

AC_LANG([C++])
ACX_MPI([CXX="$MPICXX"], [AC_MSG_ERROR([Could not find MPI.])])

dnl Check for trilinos\
dnl TO DO: Add AM_CONDITIONAL here to determine whether to build trilinos versions of libs
AX_TRILINOS_BASE(8.0.7,,AC_MSG_ERROR([Could not find required Trilinos version.]))
# AX_TRILINOS_TEUCHOS
AX_TRILINOS_EPETRA(AC_DEFINE(HAVE_LIBEPETRA,1,[Define if Epetra is available]),AC_MSG_ERROR([Could not find Epetra.]))
AC_CACHE_SAVE

dnl Check for gsl
AC_LANG([C])
AX_PATH_GSL(1.10,AM_CONDITIONAL([UQBT_GSL], [test 'TRUE']),AC_MSG_ERROR([Could not find required GSL version.]))
AC_CACHE_SAVE

dnl Check for boost
AC_LANG([C++])
AX_BOOST_BASE(1.35)
AX_BOOST_PROGRAM_OPTIONS
AC_CACHE_SAVE

dnl HPCT Config

HPCT_PREFIX=${HPCT_DIR}
AC_ARG_WITH(hpct,
        [  --with-hpct=<path>      prefix of HPCT installation (default = HPCT_DIR) ], [HPCT_PREFIX="$withval"])

HPCT_LIBS="-L${HPCT_PREFIX}/lib -lhpct"
HPCT_LIBSTATIC="${HPCT_PREFIX}/lib/libhpct.a"
HPCT_CFLAGS="-I${HPCT_PREFIX}/include"
AC_SUBST(HPCT_LIBS)
AC_SUBST(HPCT_LIBSTATIC)
AC_SUBST(HPCT_CFLAGS)

dnl Generate makefiles
AC_CONFIG_FILES([
	Makefile
	libs/Makefile
	examples/Makefile
	examples/statisticalInverseProblem1/Makefile
	examples/statisticalInverseProblem1/src/Makefile	
	examples/validationPyramid/Makefile
	examples/validationPyramid/src/Makefile
	examples/validationPyramid/src/physics_1/Makefile
	examples/validationPyramid/src/physics_2/Makefile
	examples/validationPyramid/src/physics_1_2/Makefile
	examples/validationCycle/Makefile
	examples/validationCycle/src/Makefile
	examples/statisticalForwardProblem/Makefile
	examples/statisticalForwardProblem/src/Makefile
	examples/inverse_basic_api/Makefile
	test/Makefile
        test/valid_cycle/Makefile
])

AC_OUTPUT()
