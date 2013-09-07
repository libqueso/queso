//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
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

//! Defines available optional libraries (GLPK, HDF5, Trilinos)

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

namespace QUESO {

//! Returns the rank of the calling process in the communicator.
int MyWorldfullRank();

#define ML_CODE_HAS_NEW_RESTART_CAPABILITY
#undef  QUESO_MEMORY_DEBUGGING
#undef  UQ_DEBUG_PARALLEL_RUNS_IN_DETAIL
#undef  UQ_ALSO_COMPUTE_MDFS_WITHOUT_KDE
#define QUESO_EXPECTS_LN_LIKELIHOOD_INSTEAD_OF_MINUS_2_LN
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
#define UQ_FILE_EXTENSION_FOR_HDF_FORMAT    "h5"


/*! \file uqDefines.h
    \brief Definitions and a class to provide default options to  pass to a QUESO environment.
*/

/*! \class EnvOptionsValuesClass
 *  \brief This class provides a suite options one can pass to a QUESO environment.
 * 
 *  QUESO expects the user to provide an input file with environment options for the library variables.
 *  If no input file, a collection of default values is assigned to some of the variables. The class
 *  EnvOptionsValuesClass is responsible for this task.
 */

class EnvOptionsValuesClass
{
public:
    //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor
  EnvOptionsValuesClass            ();
  
  //! Copy constructor
  EnvOptionsValuesClass            (const EnvOptionsValuesClass& src);
 
  //! Destructor
  ~EnvOptionsValuesClass            ();
  //@}
  
   //! @name Set methods
  //@{    
  //! Operator for copying the options of an environment.
  EnvOptionsValuesClass& operator= (const EnvOptionsValuesClass& rhs);
 //@}

  //! @name Attributes
 //! Number of sub-environments.
  unsigned int           m_numSubEnvironments;
  
  //! Output filename for sub-screen writing.
  std::string            m_subDisplayFileName;
  
  //! Allows (or not) all sub-environments to write to output file.
  bool                   m_subDisplayAllowAll;
  
  //! Allows (or not) all inter0 nodes to write to output file
  bool                   m_subDisplayAllowInter0;
  
  //! Sub-environments that will write to output.
  std::set<unsigned int> m_subDisplayAllowedSet;
  
  //! Verbosity.
  unsigned int           m_displayVerbosity;
  
  //! Synchronized verbosity.
  unsigned int           m_syncVerbosity;
  
  //! Checking level
  unsigned int           m_checkingLevel;

  //! Type of the random number generator.
  std::string            m_rngType;
  
  //! Seed of the random number generator.
  /*! If env_seed = -z, with z>=1, then each processor sets the seed to value MPI_RANK + z.
   It is crucial that \verb+env_seed+ takes a \underline{negative} value, otherwise all chain samples are going to be the same.*/
  int                    m_seed;
  
  //! Platform name.
  std::string            m_platformName;
  
  //! Identifying string.
  std::string            m_identifyingString;
  
  //! Number of debug parameters.
  unsigned int           m_numDebugParams;
  
  //! Debug parameters
  std::vector<double>    m_debugParams;
  //@}
  
private:
  //! Makes an exact copy of an existing EnvOptionsValues instance.
  void copy(const EnvOptionsValuesClass& src);
};

//! Macros
#define UQ_RC_MACRO(macroIRc,givenRank,where,what,retValue) \
  if (macroIRc) {                                           \
    int macroRank = givenRank;                              \
    if (macroRank < 0) {                                    \
      macroRank = QUESO::MyWorldfullRank();                      \
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
      macroRank = QUESO::MyWorldfullRank();                    \
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
      macroRank = QUESO::MyWorldfullRank();                   \
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
      macroRank = QUESO::MyWorldfullRank();                 \
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
      macroRank = QUESO::MyWorldfullRank();                 \
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

}  // End namespace QUESO

#endif // __UQ_DEFINES_H__
