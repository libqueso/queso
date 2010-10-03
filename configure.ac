# Init autoconf w/ Package name, version number, and developer email
AC_INIT([queso], [0.42.0], [pecos_dev@ices.utexas.edu])

# Set default CFLAGS & CXXFLAGS IF the user hasn't specified anything
AC_MSG_CHECKING(if CFLAGS is set)
if test "x$CFLAGS" = "x"; then
   CFLAGS="-g -O2 -Wall"
   AC_MSG_RESULT(no... setting to -g -O2 -Wall)
else
   AC_MSG_RESULT(yes)
fi

AC_MSG_CHECKING(if CXXFLAGS is set)
if test "x$CXXFLAGS" = "x"; then
   CXXFLAGS="-g -O2 -Wall"
   AC_MSG_RESULT(no... setting to -g -O2 -Wall)
else
   AC_MSG_RESULT(yes)
fi

# configuration header
AM_CONFIG_HEADER(config.h)

AC_CONFIG_AUX_DIR([build-aux])

# Init automake
AM_INIT_AUTOMAKE

# Generating a libtool library

# Release versioning

GENERIC_MAJOR_VERSION=0
GENERIC_MINOR_VERSION=40
GENERIC_MICRO_VERSION=2

GENERIC_RELEASE=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION
AC_SUBST(GENERIC_RELEASE)
AC_PROG_LIBTOOL

# Required for some macros
# AX_AM_MACROS

#-----------------------
# Check for debug mode?
#-----------------------

AC_MSG_CHECKING(whether to enable debugging)
debug_default="no"
AC_ARG_ENABLE(debug, [  --enable-debug=[no/yes] turn on debugging
                     [default=$debug_default]],, enable_debug=$debug_default)

if test "x$enable_debug" = "xyes"; then
   CFLAGS="$CFLAGS -DDEBUG"
   CXXFLAGS="$CXXFLAGS -DDEBUG"
   AC_MSG_RESULT(yes)
else
   AC_MSG_RESULT(no)
fi

# Get compilers
AC_PROG_CC
###echo "configure found CC = $CC"

AC_PROG_CXX
###echo "configure found CXX = $CXX"

#AC_PROG_RANLIB

# Check for mpi c and c++ compilers
AC_LANG([C])
ACX_MPI([CC="$MPICC"], [AC_MSG_ERROR([Could not find MPI.])])

AC_LANG([C++])
ACX_MPI([CXX="$MPICXX"], [AC_MSG_ERROR([Could not find MPI.])])


# Check for math libraries. Makes sure we are compatible for builds on Mac OS X.
# Here we force the check to only include cmath OR math.h. Note that Trilinos
# checks for both, and if both are present, will include both.
# 
# Commenting out this check below until we decide how to move forward with
# this stuff.
#
# AC_CHECK_HEADERS([cmath], [break],
#			  [AC_CHECK_HEADERS( [math.h], [break], 
#			  		     [AC_MSG_ERROR([Cannot find cmath or math.h. One or the other is required.])],[])
#			  ], [] )

# Check for trilinos\
# TO DO: Add AM_CONDITIONAL here to determine whether to build trilinos versions of libs
AX_TRILINOS_BASE(8.0.7,,AC_MSG_ERROR([Could not find required Trilinos version.]))
# AX_TRILINOS_TEUCHOS
AX_TRILINOS_EPETRA(AC_DEFINE(HAVE_LIBEPETRA,1,[Define if Epetra is available]),AC_MSG_ERROR([Could not find Epetra.]))
AC_CACHE_SAVE

# Check for gsl
AC_LANG([C])
AX_PATH_GSL(1.10,AM_CONDITIONAL([UQBT_GSL], [test 'TRUE']),AC_MSG_ERROR([Could not find required GSL version.]))
AC_CACHE_SAVE

# Check for boost
AC_LANG([C++])
BOOST_REQUIRE(1.35)
BOOST_MATH
BOOST_PROGRAM_OPTIONS
BOOST_TEST
BOOST_LDFLAGS="$BOOST_PROGRAM_OPTIONS_LDFLAGS $BOOST_PROGRAM_OPTIONS_LIBS"
AC_SUBST(BOOST_LDFLAGS)
AC_CACHE_SAVE

# Check for glpk
# AC_LANG([C])
# AX_PATH_GLPK(4.41,AC_MSG_ERROR([Could not find required GLPK version.]))
# AX_GLPK_BASE(4.41)
# AX_GLPK_PROGRAM_OPTIONS
# AC_CACHE_SAVE

AX_PATH_GLPK

# HPCT Config

AX_PATH_HPCT

#PKG_CHECK_MODULES(hpct, hpct >= 0.24)
#AC_SUBST(HPCT_CFLAGS)
#AC_SUBST(HPCT_LIBS)

#----------------
# Doxygen support
#----------------

DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(ON)
DX_PS_FEATURE(OFF)

DX_INIT_DOXYGEN([queso],[doxygen/queso.dox],docs)

# Query configuration environment

AX_SUMMARIZE_ENV

# Generate makefiles

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
	examples/validationCycle2/Makefile
	examples/validationCycle2/src/Makefile
	examples/statisticalForwardProblem/Makefile
	examples/statisticalForwardProblem/src/Makefile
	examples/basic_interface/Makefile
	test/Makefile
	test/t01_valid_cycle/Makefile
	test/t02_sip_sfp/Makefile
	test/t03_sequence/Makefile
        test/t01_valid_cycle/valid_cycle/Makefile
        test/t02_sip_sfp/sip_sfp/Makefile
        test/t03_sequence/sequence/Makefile
        test/gsl_tests/Makefile
	doxygen/Makefile
])


AC_OUTPUT()

# Final summary

AX_SUMMARIZE_CONFIG
