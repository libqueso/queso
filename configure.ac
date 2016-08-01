# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.65)
AC_INIT([queso], [0.55.0], [queso-users@googlegroups.com])
PACKAGE_DESCRIPTION="The parallel C++ statistical library QUESO: Quantification of uncertainty for estimation, simulation and optimization"
AC_SUBST([PACKAGE_DESCRIPTION])
PACKAGE_URL="https://github.com/libqueso/queso"
AC_SUBST([PACKAGE_URL])

AC_CONFIG_MACRO_DIR([m4])

AC_CONFIG_HEADERS(config_queso.h.tmp)
dnl Prefixes everything in the temporary file with QUESO_
AX_PREFIX_CONFIG_H(config_queso.h,QUESO,config_queso.h.tmp)

AC_CONFIG_AUX_DIR([build-aux])
AC_CANONICAL_TARGET()
AM_INIT_AUTOMAKE([color-tests subdir-objects])
AM_SILENT_RULES(yes)  # use silent rules where available - automake 1.11

# Release versioning

AX_SPLIT_VERSION
GENERIC_MAJOR_VERSION=$AX_MAJOR_VERSION
GENERIC_MINOR_VERSION=$AX_MINOR_VERSION
GENERIC_MICRO_VERSION=$AX_POINT_VERSION

# Generating a libtool library

GENERIC_RELEASE=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION

AC_SUBST(GENERIC_RELEASE)
AC_SUBST(GENERIC_MAJOR_VERSION)
AC_SUBST(GENERIC_MINOR_VERSION)
AC_SUBST(GENERIC_MICRO_VERSION)

LT_INIT

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
   QUESO_CPPFLAGS="$CPPFLAGS -DDEBUG"
   AC_SUBST(QUESO_CPPFLAGS)
   AC_MSG_RESULT(yes)
else
   AC_MSG_RESULT(no)
fi

#############################################
# Check if the user explicitly disabled mpi #
#############################################
AC_ARG_ENABLE(mpi,
              AS_HELP_STRING([--disable-mpi],
                             [build without message passing support]),
              [enable_mpi=no],
              [enable_mpi=yes])

##############################################
# Check if the user disabled program options #
##############################################
AC_ARG_ENABLE([boost-program-options],
              AS_HELP_STRING([--disable-boost-program-options], [Disables use
                              of boost program options for parsing the QUESO
                              input file.  GetPot will be used instead.]),
              [enable_program_options="no"],
              [enable_program_options="yes"])

#-------------------
# Compilers and MPI
#-------------------
AC_PROG_CC
AC_PROG_CXX
AC_PROG_FC
AC_LANG([C])

HAVE_MPI=0
if test "x$enable_mpi" = "xyes"; then
    ACX_MPI(CC="$MPICC", [AC_MSG_ERROR([Could not find MPI.])])
    HAVE_MPI=1
fi

AC_LANG([C++])

if test "x$enable_mpi" = "xyes"; then
    ACX_MPI(CXX="$MPICXX", [AC_MSG_ERROR([Could not find MPI.])])
    HAVE_MPI=1
fi
AC_SUBST(HAVE_MPI)
AM_CONDITIONAL(MPI_ENABLED, test x$HAVE_MPI = x1)


# --------------------------------------------------------------
# Allow the user to explicitly turn off C++11 support by saying
# --disable-cxx11.  Otherwise, C++11 support is enabled by default
# although libmesh does not (currently) use any C++11 features.
# --------------------------------------------------------------
AC_ARG_ENABLE(cxx11,
              AC_HELP_STRING([--disable-cxx11],
                             [build without C++11 support]),
                             [case "${enableval}" in
                   yes) enablecxx11=yes ;;
                   no)  enablecxx11=no ;;
                   *)   AC_MSG_ERROR(bad value ${enableval} for --disable-cxx11) ;;
                   esac],
                   [enablecxx11=yes])

if (test "x$enablecxx11" = "xyes"); then
  # --------------------------------------------------------------
  # Standard autoconf macro for determining the proper -std=c++11
  # flag for the current compiler.  Adds the result to CXXFLAGS if
  # one is found.  See ax_cxx_compile_stdcxx_11.m4 for details.
  # --------------------------------------------------------------
  AX_CXX_COMPILE_STDCXX_11([noext],[optional])
fi

QUESO_TEST_CXX11_ISNAN
QUESO_TEST_CXX11_ISFINITE
QUESO_TEST_CXX11_UNIQUE_PTR
QUESO_TEST_CXX11_SHARED_PTR

#-------------------------
# External Library Checks
#-------------------------

# Check for Trilinos (optional)
#Kemelli: The lowest version used to be 8.0.7, now is 11.0.2, to accomodate the changes in -lteuchos*.

HAVE_TRILINOS=0

AX_TRILINOS_BASE(11.0.2,[HAVE_TRILINOS=1],[HAVE_TRILINOS=0])

if test "$HAVE_TRILINOS" = "1" ;then
   #AX_TRILINOS_EPETRA([HAVE_TRILINOS=1],[HAVE_TRILINOS=0])
   AX_TRILINOS_EPETRA([HAVE_TRILINOS=1],AC_MSG_ERROR([Could not find Epetra - which is required if Trilinos is enabled.]))
fi

AC_SUBST(HAVE_TRILINOS)
AM_CONDITIONAL(TRILINOS_ENABLED,test x$HAVE_TRILINOS = x1)

AC_CACHE_SAVE

# Check for libmesh (optional)
AX_PATH_LIBMESH_NEW(0.9.1,no)
AC_CACHE_SAVE

# Check for GSL (required)

AC_LANG([C])
AX_PATH_GSL(1.10,AM_CONDITIONAL([UQBT_GSL], [test 'TRUE']),AC_MSG_ERROR([Could not find required GSL version.]))
AC_CACHE_SAVE

# Check for boost (required)

AC_LANG([C++])
BOOST_REQUIRE([1.35])
AS_IF([test "x$enable_program_options" != "xno"],
      [
       BOOST_PROGRAM_OPTIONS
      ])

AS_IF([test "x$enable_program_options" = "xno"],
      [
       QUESO_CPPFLAGS="$QUESO_CPPFLAGS -DDISABLE_BOOST_PROGRAM_OPTIONS"
       AC_SUBST(QUESO_CPPFLAGS)
      ])

BOOST_FIND_HEADER([boost/scoped_ptr.hpp])
BOOST_FIND_HEADER([boost/shared_ptr.hpp])
BOOST_FIND_HEADER([boost/math/special_functions.hpp])
AC_CACHE_SAVE

# Check for GLPK (optional)

AX_PATH_GLPK([4.35],[no])

# Check for HDF5 (optional)

AX_PATH_HDF5_NEW([1.8.0],[no])

# Check for ANN (external library)
# AX_PATH_ANN
#### TODO: Make sure that the ANN uses L-infinity (Max) norm

# Check for ANN feature
AX_ENABLE_ANN

# Check for libGRVY (optional as of QUESO version 0.46.0)

AX_PATH_GRVY_NEW([0.29],[no])

# Check for slepc
#AX_PATH_SLEPC_NEW([3.3],[no])

AC_CACHE_SAVE

#----------------
# Doxygen support
#----------------

DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(OFF)
DX_PS_FEATURE(OFF)

DX_INIT_DOXYGEN([queso],[doxygen/queso.dox],[docs])

#--------------------------
# Checks for code coverage
#-------------------------

AX_CODE_COVERAGE

# Query configuration environment

AX_SUMMARIZE_ENV

# Generate makefiles

AC_CONFIG_FILES([
  queso.pc
  queso-deps.pc
  Makefile
  inc/queso/Makefile
  src/Makefile
  src/contrib/ANN/Makefile
  src/contrib/ANN/test/Makefile
  src/core/inc/queso.h
  examples/Makefile
  test/Makefile
  test/t01_valid_cycle/Makefile
  test/t02_sip_sfp/Makefile
  test/t03_sequence/Makefile
  test/t04_bimodal/Makefile
  test/gsl_tests/Makefile
  test/gsl_tests/input
  doxygen/Makefile
  doxygen/queso.dox
  doxygen/txt_common/about_vpath.page
])

AC_CONFIG_FILES(src/apps/queso-config, [chmod +x src/apps/queso-config])

dnl ----------------------------------------------
dnl Collect files for licence header stamping here
dnl ----------------------------------------------
_AC_SRCDIRS(.)

# We have to do this by subdirectory because otherwise distcheck breaks as we
# start picking up files in the directories that it uses
QUESO_STAMPED_FILES=$(find $ac_abs_top_srcdir/{src,examples,test} -name "*.h" -or -name "*.C" | tr "\n" " ")
AC_SUBST(QUESO_STAMPED_FILES)

AC_OUTPUT()

#---------------
# Final summary
#---------------

AX_SUMMARIZE_CONFIG
