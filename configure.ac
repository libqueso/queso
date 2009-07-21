dnl Init autoconf w/ Package name, version number, and developer email
AC_INIT([queso], [0.4.0], [pecos_dev@ices.utexas.edu])

dnl configuration header
AM_CONFIG_HEADER(config.h)

AC_CONFIG_AUX_DIR([build-aux])

dnl Init automake
AM_INIT_AUTOMAKE

dnl Generating a libtool library

dnl Release versioning

GENERIC_MAJOR_VERSION=0
GENERIC_MINOR_VERSION=4
GENERIC_MICRO_VERSION=0

GENERIC_RELEASE=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION
AC_SUBST(GENERIC_RELEASE)
AC_PROG_LIBTOOL

dnl Required for some macros
AX_AM_MACROS

#-----------------------
# Check for debug mode?
#-----------------------

AC_MSG_CHECKING(whether to enable debugging)
debug_default="no"
AC_ARG_ENABLE(debug, [  --enable-debug=[no/yes] turn on debugging
                     [default=$debug_default]],, enable_debug=$debug_default)

if test "x$enable_debug" = "xyes"; then
   CFLAGS="-g -O0 -DDEBUG"
   CXXFLAGS="-g -O0 -DDEBUG"
   AC_MSG_RESULT(yes)
else
   CFLAGS="-g -O3"
   CXXFLAGS="-g -O3"
   AC_MSG_RESULT(no)
fi

dnl Get compilers
AC_PROG_CC
###echo "configure found CC = $CC"

AC_PROG_CXX
###echo "configure found CXX = $CXX"

#AC_PROG_RANLIB

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

# Check for boost
AC_LANG([C++])
AX_BOOST_BASE(1.35)
AX_BOOST_PROGRAM_OPTIONS
AC_CACHE_SAVE

# HPCT Config

AX_PATH_HPCT




#PKG_CHECK_MODULES(hpct, hpct >= 0.23)
#AC_SUBST(HPCT_CFLAGS)
#AC_SUBST(HPCT_LIBS)

# Query configuration environment

AX_SUMMARIZE_ENV

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
	examples/basic_interface/Makefile
	test/Makefile
        test/t01_valid_cycle/valid_cycle/Makefile
        test/t02_sip_sfp/sip_sfp/Makefile
])


AC_OUTPUT()

# Final summary

AX_SUMMARIZE_CONFIG