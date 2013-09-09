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

#ifndef __UQ_MOC_SG_OPTIONS_H__
#define __UQ_MOC_SG_OPTIONS_H__

#include <uqEnvironment.h>
#include <uqSequenceStatisticalOptions.h>

#define UQ_MOC_SG_FILENAME_FOR_NO_FILE "."

// _ODV = option default value
#define UQ_MOC_SG_DATA_OUTPUT_FILE_NAME_ODV        UQ_MOC_SG_FILENAME_FOR_NO_FILE
#define UQ_MOC_SG_DATA_OUTPUT_ALLOWED_SET_ODV      ""

#define UQ_MOC_SG_PSEQ_DATA_OUTPUT_PERIOD_ODV      0
#define UQ_MOC_SG_PSEQ_DATA_OUTPUT_FILE_NAME_ODV   UQ_MOC_SG_FILENAME_FOR_NO_FILE
#define UQ_MOC_SG_PSEQ_DATA_OUTPUT_FILE_TYPE_ODV   UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT
#define UQ_MOC_SG_PSEQ_DATA_OUTPUT_ALLOWED_SET_ODV ""
#define UQ_MOC_SG_PSEQ_COMPUTE_STATS_ODV           0

#define UQ_MOC_SG_QSEQ_DATA_INPUT_FILE_NAME_ODV    UQ_MOC_SG_FILENAME_FOR_NO_FILE
#define UQ_MOC_SG_QSEQ_DATA_INPUT_FILE_TYPE_ODV    UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT
#define UQ_MOC_SG_QSEQ_SIZE_ODV                    100
#define UQ_MOC_SG_QSEQ_DISPLAY_PERIOD_ODV          500
#define UQ_MOC_SG_QSEQ_MEASURE_RUN_TIMES_ODV       0
#define UQ_MOC_SG_QSEQ_DATA_OUTPUT_PERIOD_ODV      0
#define UQ_MOC_SG_QSEQ_DATA_OUTPUT_FILE_NAME_ODV   UQ_MOC_SG_FILENAME_FOR_NO_FILE
#define UQ_MOC_SG_QSEQ_DATA_OUTPUT_FILE_TYPE_ODV   UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT
#define UQ_MOC_SG_QSEQ_DATA_OUTPUT_ALLOWED_SET_ODV ""
#define UQ_MOC_SG_QSEQ_COMPUTE_STATS_ODV           0

namespace QUESO {

/*! \file uqMonteCarloSGOptions.h
    \brief Classes to allow options to be passed to a Monte Carlo sequence generator.
*/

/*! \class McOptionsValuesClass
 *  \brief This class provides options for the  Monte Carlo sequence generator if no input file is available.
 * 
 *  Monte Carlo sequence generator expects options for its methods. This class provides default
 * values for such options if no input file is available. */

class McOptionsValuesClass
{
public:
  //! @name Constructor/Destructor methods
  //@{ 
  
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  McOptionsValuesClass            (const SsOptionsValuesClass* alternativePSsOptionsValues,
                                     const SsOptionsValuesClass* alternativeQSsOptionsValues);
#else
  //! Default constructor.
  /*! Assigns the default suite of options to the Monte Carlo sequence generator.*/
  McOptionsValuesClass            ();
#endif
  //! Copy constructor.
  /*! It assigns the same options values from  \c src to \c this.*/
  McOptionsValuesClass            (const McOptionsValuesClass& src);
  
  //! Destructor
  ~McOptionsValuesClass            ();
  //@}
  
  //! @name Set methods
  //@{ 
  //! Assignment operator; it copies \c rhs to \c this. 
  McOptionsValuesClass& operator= (const McOptionsValuesClass& rhs);
  //@}

  std::string                        m_dataOutputFileName;
  std::set<unsigned int>             m_dataOutputAllowedSet;

  unsigned int                       m_pseqDataOutputPeriod;
  std::string                        m_pseqDataOutputFileName;
  std::string                        m_pseqDataOutputFileType;
  std::set<unsigned int>             m_pseqDataOutputAllowedSet;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  bool                               m_pseqComputeStats;
#endif

  std::string                        m_qseqDataInputFileName;
  std::string                        m_qseqDataInputFileType;
  unsigned int                       m_qseqSize;
  unsigned int                       m_qseqDisplayPeriod;
  bool                               m_qseqMeasureRunTimes;
  unsigned int                       m_qseqDataOutputPeriod;
  std::string                        m_qseqDataOutputFileName;
  std::string                        m_qseqDataOutputFileType;
  std::set<unsigned int>             m_qseqDataOutputAllowedSet;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  bool                               m_qseqComputeStats;
#endif

private:
  //! Copies the option values from \c src to \c this.
  void copy(const McOptionsValuesClass& src);

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  friend class MonteCarloSGOptionsClass;
  SsOptionsValuesClass             m_alternativePSsOptionsValues;
  SsOptionsValuesClass             m_alternativeQSsOptionsValues;
#endif
};

// --------------------------------------------------
// --------------------------------------------------
// --------------------------------------------------

/*! \class MonteCarloSGOptionsClass
 *  \brief This class reads the options for the  Monte Carlo sequence generator from  an input file.
 * 
 * Monte Carlo sequence generator expects options for its methods. This class reads the 
 * options for the Monte Carlo sequence generator from an input file provided by the user. 
 * The class expects the prefix '\<prefix\>_mc_'. For instance, if 'prefix' is 'foo_775_fp_', 
 * then the constructor will read all options that begin with 'foo_775_fp_mc_'. */

class MonteCarloSGOptionsClass
{
public:
  
  //! @name Constructor/Destructor methods
  //@{ 
  //! Constructor: reads options from the input file.
  MonteCarloSGOptionsClass(const BaseEnvironmentClass& env, const char* prefix);
  
  //! Constructor: with alternative option values.
  /*! In this constructor, the input options are given by \c alternativeOptionsValues, thus, they
   * are not read from an input file.*/
  MonteCarloSGOptionsClass(const BaseEnvironmentClass& env, const char* prefix, const McOptionsValuesClass& alternativeOptionsValues);
 
  //! Destructor
  ~MonteCarloSGOptionsClass();
  //@}
  
  //! @name I/O methods
  //@{
  //! It scans the option values from the options input file.
  void scanOptionsValues();
  
  //!  It prints the option values.
  void print            (std::ostream& os) const;
  //@}
  
  McOptionsValuesClass             m_ov;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  SequenceStatisticalOptionsClass* m_pseqStatisticalOptionsObj;
  SequenceStatisticalOptionsClass* m_qseqStatisticalOptionsObj;
#endif
  std::string                        m_prefix;

private:
  //! Defines the options for the Monte Carlo sequence generator as the default options.
  void   defineMyOptions  (po::options_description& optionsDesc) const;
  
  //! Gets the sequence options.
  void   getMyOptionValues(po::options_description& optionsDesc);

  const BaseEnvironmentClass& m_env;
  po::options_description*      m_optionsDesc;

  std::string                   m_option_help;
  std::string                   m_option_dataOutputFileName;
  std::string                   m_option_dataOutputAllowedSet;

  std::string                   m_option_pseq_dataOutputPeriod;
  std::string                   m_option_pseq_dataOutputFileName;
  std::string                   m_option_pseq_dataOutputFileType;
  std::string                   m_option_pseq_dataOutputAllowedSet;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  std::string                   m_option_pseq_computeStats;
#endif

  std::string                   m_option_qseq_dataInputFileName;
  std::string                   m_option_qseq_dataInputFileType;
  std::string                   m_option_qseq_size;
  std::string                   m_option_qseq_displayPeriod;
  std::string                   m_option_qseq_measureRunTimes;
  std::string                   m_option_qseq_dataOutputPeriod;
  std::string                   m_option_qseq_dataOutputFileName;
  std::string                   m_option_qseq_dataOutputFileType;
  std::string                   m_option_qseq_dataOutputAllowedSet;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  std::string                   m_option_qseq_computeStats;
#endif
};

//! Prints the object \c obj, overloading an operator.
std::ostream& operator<<(std::ostream& os, const MonteCarloSGOptionsClass& obj);

}  // End namespace QUESO

#endif // __UQ_MOC_SG_OPTIONS_H__
