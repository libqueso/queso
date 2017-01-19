//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2017 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-

#ifndef UQ_DEFINES_H
#define UQ_DEFINES_H

#include <queso/config_queso.h>
#include <queso/asserts.h>

#ifdef QUESO_HAVE_GLPK
#define QUESO_HAS_GLPK
#endif

#ifdef QUESO_HAVE_HDF5
#define QUESO_HAS_HDF5
#endif

#ifdef QUESO_HAVE_TRILINOS
#define QUESO_HAS_TRILINOS
#endif

#ifdef QUESO_HAVE_ANN
#define QUESO_HAS_ANN
#endif

#ifdef QUESO_HAVE_MPI
#define QUESO_HAS_MPI
#endif

//! This define is deprecated.  Remove any #ifdef statements in user code.
#define QUESO_EXPECTS_LN_LIKELIHOOD_INSTEAD_OF_MINUS_2_LN

// Use GSL inline functions
#define HAVE_INLINE

// So we don't clash with other getpots
#define GETPOT_NAMESPACE QUESO

// And only do GSL range-checking if we're really debugging
#ifndef DEBUG
#define GSL_RANGE_CHECK_OFF
#endif

#include <iostream>
#include <stdlib.h> // For exit()
#include <set>
#include <vector>

#ifdef QUESO_HAS_MPI
#include <mpi.h>
#endif

#include <queso/asserts.h> // for queso_error handler

namespace QUESO {

//! Returns the rank of the calling process in the communicator.
int MyWorldfullRank();

#define ML_CODE_HAS_NEW_RESTART_CAPABILITY
#undef  QUESO_MEMORY_DEBUGGING
#undef  UQ_DEBUG_PARALLEL_RUNS_IN_DETAIL
#undef  UQ_ALSO_COMPUTE_MDFS_WITHOUT_KDE
#define QUESO_CLASSES_INSTANTIATE_NEW_MAPS
#undef  QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
#undef  QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
#define QUESO_USES_NEW_RNG_CLASS

const int UQ_UNAVAILABLE_RANK = -1;

const int UQ_OK_RC                         =  0;
const int UQ_INCOMPLETE_IMPLEMENTATION_RC  = -1;
const int UQ_INVALID_PARAMETER_SPEC_RC     = -2;
const int UQ_INVALID_OBSERVABLE_SPEC_RC    = -3;
const int UQ_INVALID_QOI_SPEC_RC           = -4;
const int UQ_INVALID_INTERNAL_RESULT_RC    = -5;
const int UQ_INVALID_INTERNAL_STATE_RC     = -6;
const int UQ_FAILED_TO_OPEN_FILE_RC        = -7;
const int UQ_MATRIX_IS_NOT_POS_DEFINITE_RC = -8;
const int UQ_FAILED_READING_FILE_RC        = -9;
const int UQ_INVALID_SPACE_COMPONENT_ID_RC = -10;
const int UQ_MATRIX_SVD_FAILED_RC          = -11;

#define UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT "m"
#define UQ_FILE_EXTENSION_FOR_TXT_FORMAT    "txt"
#define UQ_FILE_EXTENSION_FOR_HDF_FORMAT    "h5"


/*! \file Defines.h
    \brief Definitions and a class to provide default options to  pass to a QUESO environment.
*/

// Macros

// The following code is a copy-pasta from libmesh.  The same licence applies,
// so we're good here.
//
// The queso_do_once macro helps us avoid redundant repeated
// repetitions of the same warning messages
#undef queso_do_once
#define queso_do_once(do_this)             \
  do {                                     \
    static bool did_this_already = false;  \
    if (!did_this_already) {               \
      did_this_already = true;             \
      do_this;                             \
    } } while (0)

// The queso_warning macro outputs a file/line/time stamped warning
// message, if warnings are enabled.
#define queso_warning(message)         \
  queso_do_once(std::cerr << message  \
                << __FILE__ << ", line " << __LINE__ << ", compiled " << __DATE__ << " at " << __TIME__ << " ***" << std::endl;)

// The queso_deprecated macro warns that you are using obsoleted code
#define queso_deprecated()  \
  queso_warning("*** Warning, this code is deprecated and likely to be removed in future library versions:  ");

#define UQ_RC_MACRO(macroIRc,givenRank,where,what,retValue) \
  queso_deprecated();                                       \
  if (macroIRc) {                                           \
    int macroRank = givenRank;                              \
    if (macroRank < 0) {                                    \
      macroRank = QUESO::MyWorldfullRank();                 \
    }                                                       \
    std::cerr << "UQ RC ERROR"                              \
              << ", rank "  << macroRank                    \
              << ", in "    << where                        \
              << ": "       << what                         \
              << ", iRc = " << macroIRc                     \
              << std::endl;                                 \
    return retValue;                                        \
  }

#define UQ_TEST_MACRO(test,givenRank,where,what,retValue) \
  queso_deprecated();                                     \
  if (test) {                                             \
    int macroRank = givenRank;                            \
    if (macroRank < 0) {                                  \
      macroRank = QUESO::MyWorldfullRank();               \
    }                                                     \
    std::cerr << "UQ TEST ERROR"                          \
              << ", rank " << macroRank                   \
              << ", in "   << where                       \
              << ": "      << what                        \
              << std::endl;                               \
    return retValue;                                      \
  }

#define UQ_FATAL_RC_MACRO(macroIRc,givenRank,where,what) \
  queso_deprecated();                                    \
  if (macroIRc) {                                        \
    int macroRank = givenRank;                           \
    if (macroRank < 0) {                                 \
      macroRank = QUESO::MyWorldfullRank();              \
    }                                                    \
    std::cerr << "UQ RC FATAL ERROR"                     \
              << ", rank "  << macroRank                 \
              << ", in "    << where                     \
              << ": "       << what                      \
              << ", iRC = " << macroIRc                  \
              << ". Exiting..."                          \
              << std::endl;                              \
    queso_error(); \
  }

#define UQ_FATAL_TEST_MACRO(test,givenRank,where,what)  \
  queso_deprecated();                                   \
  if (test) {                                           \
    int macroRank = givenRank;                          \
    if (macroRank < 0) {                                \
      macroRank = QUESO::MyWorldfullRank();             \
    }                                                   \
    std::cerr << "UQ TEST FATAL ERROR"                  \
              << ", rank "  << macroRank                \
              << ", in "    << where                    \
              << ": "       << what                     \
              << ". Exiting..."                         \
              << std::endl;                             \
    queso_error();	                   		\
  }

}  // End namespace QUESO

#endif // UQ_DEFINES_H
