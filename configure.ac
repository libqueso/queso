# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script. 

AC_PREREQ(2.61)
AC_INIT([queso], [0.43.0], [pecos_dev@ices.utexas.edu])
AM_CONFIG_HEADER(config.h)
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE

# Release versioning

GENERIC_MAJOR_VERSION=0
GENERIC_MINOR_VERSION=43
GENERIC_MICRO_VERSION=0

# Generating a libtool library

GENERIC_RELEASE=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION

AC_SUBST(GENERIC_RELEASE)
AC_SUBST(GENERIC_MAJOR_VERSION)
AC_SUBST(GENERIC_MINOR_VERSION)
AC_SUBST(GENERIC_MICRO_VERSION)

AC_PROG_LIBTOOL

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

#-------------------
# Compilers and MPI
#-------------------

AC_PROG_CC
AC_PROG_CXX
AC_PROG_FC
AC_LANG([C])
ACX_MPI([CC="$MPICC"], [AC_MSG_ERROR([Could not find MPI.])])
AC_LANG([C++])
ACX_MPI([CXX="$MPICXX"], [AC_MSG_ERROR([Could not find MPI.])])

#-------------------------
# External Library Checks
#-------------------------

# Check for Trilinos
# TODO: Add AM_CONDITIONAL here to determine whether to build trilinos versions of libs

AX_TRILINOS_BASE(8.0.7,,AC_MSG_ERROR([Could not find required Trilinos version.]))
# AX_TRILINOS_TEUCHOS
AX_TRILINOS_EPETRA(AC_DEFINE(HAVE_LIBEPETRA,1,[Define if Epetra is available]),AC_MSG_ERROR([Could not find Epetra.]))
AC_CACHE_SAVE

# Check for GSL

AC_LANG([C])
AX_PATH_GSL(1.10,AM_CONDITIONAL([UQBT_GSL], [test 'TRUE']),AC_MSG_ERROR([Could not find required GSL version.]))
AC_CACHE_SAVE

# Check for boost

AC_LANG([C++])
BOOST_REQUIRE([1.35])
BOOST_PROGRAM_OPTIONS
AC_CACHE_SAVE

# Check for GNU Linear Programming kit

AX_PATH_GLPK(4.35,,AC_MSG_ERROR([Could not find required GLPK error]))

# Check for HDF5 library

AX_PATH_HDF5(1.8.0,,AC_MSG_ERROR([Could not find required HDF5 error]))

# libGRVY Config

AX_PATH_GRVY

AC_CACHE_SAVE

#PKG_CHECK_MODULES(grvy, grvy >= 0.24)

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

#--------------------------
# Checks for code coverage
#-------------------------

AX_CODE_COVERAGE

# Query configuration environment

AX_SUMMARIZE_ENV

# Generate makefiles

AC_CONFIG_FILES([
	Makefile
	src/Makefile
	src/core/inc/queso.h
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
	test/t04_bimodal/Makefile
        test/t02_sip_sfp/sip_sfp/Makefile
        test/t03_sequence/sequence/Makefile
        test/t04_bimodal/bimodal/Makefile
        test/gsl_tests/Makefile
	doxygen/Makefile
])


AC_OUTPUT()

#---------------
# Final summary
#---------------

AX_SUMMARIZE_CONFIG
