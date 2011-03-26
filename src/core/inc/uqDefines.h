//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011 The PECOS Development Team
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
// 
// $Id$
//
//--------------------------------------------------------------------------

#ifndef __UQ_DEFINES_H__
#define __UQ_DEFINES_H__

#include <config_queso.h>

// Define available optional libraries (GLPK, HDF5, Trilinos)

#ifdef HAVE_GLPK
#define QUESO_HAS_GLPK
#endif

#ifdef HAVE_HDF5
#define QUESO_HAS_HDF5
#endif

#ifdef HAVE_TRILINOS
#define QUESO_HAS_TRILINOS
#endif

#ifdef HAVE_ANN
#define QUESO_HAS_ANN
#endif

#define QUESO_HAS_MPI

#ifdef QUESO_HAS_TRILINOS
#ifdef QUESO_HAS_MPI
// Ok
#else
Incompatible combination of defines in QUESO 'uqDefines.h': QUESO_HAS_TRILINOS is defined but QUESO_HAS_MPI is undefined
#endif
#endif

#include <iostream>
#include <stdlib.h> // For exit()
#include <set>
#include <vector>

int uqMyWorldfullRank();

#define ML_CODE_HAS_NEW_RESTART_CAPABILITY
#undef QUESO_MEMORY_DEBUGGING
#undef UQ_DEBUG_PARALLEL_RUNS_IN_DETAIL
#undef UQ_ALSO_COMPUTE_MDFS_WITHOUT_KDE
#define QUESO_EXPECTS_LN_LIKELIHOOD_INSTEAD_OF_MINUS_2_LN

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

#define UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT "m"
#define UQ_FILE_EXTENSION_FOR_HDF_FORMAT    "h5"

class uqEnvOptionsValuesClass
{
public:
  uqEnvOptionsValuesClass            ();
  uqEnvOptionsValuesClass            (const uqEnvOptionsValuesClass& src);
  uqEnvOptionsValuesClass& operator= (const uqEnvOptionsValuesClass& rhs);
 ~uqEnvOptionsValuesClass            ();

  unsigned int           m_numSubEnvironments;
  std::string            m_subDisplayFileName;
  bool                   m_subDisplayAllowAll;
  std::set<unsigned int> m_subDisplayAllowedSet;
  unsigned int           m_displayVerbosity;
  unsigned int           m_syncVerbosity;
  int                    m_seed;
  std::string            m_identifyingString;
  unsigned int           m_numDebugParams;
  std::vector<double>    m_debugParams;

private:
  void copy(const uqEnvOptionsValuesClass& src);
};

#define UQ_RC_MACRO(macroIRc,givenRank,where,what,retValue) \
  if (macroIRc) {                                           \
    int macroRank = givenRank;                              \
    if (macroRank < 0) {                                    \
      macroRank = uqMyWorldfullRank();                      \
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
  if (test) {                                             \
    int macroRank = givenRank;                            \
    if (macroRank < 0) {                                  \
      macroRank = uqMyWorldfullRank();                    \
    }                                                     \
    std::cerr << "UQ TEST ERROR"                          \
              << ", rank " << macroRank                   \
              << ", in "   << where                       \
              << ": "      << what                        \
              << std::endl;                               \
    return retValue;                                      \
  }

#define UQ_FATAL_RC_MACRO(macroIRc,givenRank,where,what) \
  if (macroIRc) {                                        \
    int macroRank = givenRank;                           \
    if (macroRank < 0) {                                 \
      macroRank = uqMyWorldfullRank();                   \
    }                                                    \
    std::cerr << "UQ RC FATAL ERROR"                     \
              << ", rank "  << macroRank                 \
              << ", in "    << where                     \
              << ": "       << what                      \
              << ", iRC = " << macroIRc                  \
              << ". Exiting..."                          \
              << std::endl;                              \
    exit(1);                                             \
  }

#ifdef QUESO_HAS_MPI
#define UQ_FATAL_TEST_MACRO(test,givenRank,where,what) \
  if (test) {                                          \
    int macroRank = givenRank;                         \
    if (macroRank < 0) {                               \
      macroRank = uqMyWorldfullRank();                 \
    }                                                  \
    std::cerr << "UQ TEST FATAL ERROR"                 \
              << ", rank "  << macroRank               \
              << ", in "    << where                   \
              << ": "       << what                    \
              << ". Exiting..."                        \
              << std::endl;                            \
    /*int macroMpiRC = 0;*/                            \
    /*macroMpiRC = */MPI_Abort(MPI_COMM_WORLD,-999);   \
    exit(1);                                           \
  }
#else
#define UQ_FATAL_TEST_MACRO(test,givenRank,where,what) \
  if (test) {                                          \
    int macroRank = givenRank;                         \
    if (macroRank < 0) {                               \
      macroRank = uqMyWorldfullRank();                 \
    }                                                  \
    std::cerr << "UQ TEST FATAL ERROR"                 \
              << ", rank "  << macroRank               \
              << ", in "    << where                   \
              << ": "       << what                    \
              << ". Exiting..."                        \
              << std::endl;                            \
    exit(1);                                           \
  }
#endif

#endif // __UQ_DEFINES_H__
