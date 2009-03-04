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

dnl Check for mpi (FIX ME: this can't be the best way)
ACX_MPI
MPI_CXXFLAGS=`mpicxx --showme:compile` 
AC_SUBST(MPI_CXXFLAGS)
MPI_CXXLIBS=`mpicxx --showme:link`
AC_SUBST(MPI_CXXLIBS)
AC_CACHE_SAVE

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

dnl Generate makefiles
AC_CONFIG_FILES([
	Makefile
	libs/Makefile
	examples/Makefile
	examples/markovChain1/Makefile
	examples/markovChain1/src/Makefile	
	examples/pyramid/Makefile
	examples/pyramid/src/Makefile
	examples/pyramid/src/physics_1/Makefile
	examples/pyramid/src/physics_2/Makefile
	examples/pyramid/src/physics_1_2/Makefile
	examples/tga/Makefile
	examples/tga/src/Makefile
	examples/monteCarlo/Makefile
	examples/monteCarlo/src/Makefile
])

AC_OUTPUT()