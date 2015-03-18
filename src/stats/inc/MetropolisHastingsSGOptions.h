//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2015 The PECOS Development Team
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

#ifndef UQ_MH_SG_OPTIONS_H
#define UQ_MH_SG_OPTIONS_H

#include <queso/Environment.h>
#include <queso/MLSamplingLevelOptions.h>
#include <queso/SequenceStatisticalOptions.h>

#undef  UQ_MH_SG_REQUIRES_INVERTED_COV_MATRICES
#define UQ_NOTHING_JUST_FOR_TEST_OF_SVN_ID 1

#define UQ_MH_SG_FILENAME_FOR_NO_FILE   "."

// _ODV = option default value
#define UQ_MH_SG_DATA_OUTPUT_FILE_NAME_ODV                            UQ_MH_SG_FILENAME_FOR_NO_FILE
#define UQ_MH_SG_DATA_OUTPUT_ALLOW_ALL_ODV                            0
#define UQ_MH_SG_DATA_OUTPUT_ALLOWED_SET_ODV                          ""

#define UQ_MH_SG_TOTALLY_MUTE_ODV                                     0
#define UQ_MH_SG_INITIAL_POSITION_DATA_INPUT_FILE_NAME_ODV            UQ_MH_SG_FILENAME_FOR_NO_FILE
#define UQ_MH_SG_INITIAL_POSITION_DATA_INPUT_FILE_TYPE_ODV            UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT
#define UQ_MH_SG_INITIAL_PROPOSAL_COV_MATRIX_DATA_INPUT_FILE_NAME_ODV UQ_MH_SG_FILENAME_FOR_NO_FILE
#define UQ_MH_SG_INITIAL_PROPOSAL_COV_MATRIX_DATA_INPUT_FILE_TYPE_ODV UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT
#define UQ_MH_SG_LIST_OF_DISABLED_PARAMETERS_ODV                      ""
#define UQ_MH_SG_RAW_CHAIN_DATA_INPUT_FILE_NAME_ODV                   UQ_MH_SG_FILENAME_FOR_NO_FILE
#define UQ_MH_SG_RAW_CHAIN_DATA_INPUT_FILE_TYPE_ODV                   UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT
#define UQ_MH_SG_RAW_CHAIN_SIZE_ODV                                   100
#define UQ_MH_SG_RAW_CHAIN_GENERATE_EXTRA_ODV                         0
#define UQ_MH_SG_RAW_CHAIN_DISPLAY_PERIOD_ODV                         500
#define UQ_MH_SG_RAW_CHAIN_MEASURE_RUN_TIMES_ODV                      1
#define UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_PERIOD_ODV                     0
#define UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_FILE_NAME_ODV                  UQ_MH_SG_FILENAME_FOR_NO_FILE
#define UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_FILE_TYPE_ODV                  UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT
#define UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_ALLOW_ALL_ODV                  0
#define UQ_MH_SG_RAW_CHAIN_DATA_OUTPUT_ALLOWED_SET_ODV                ""
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
#define UQ_MH_SG_RAW_CHAIN_COMPUTE_STATS_ODV                          0
#endif
#define UQ_MH_SG_FILTERED_CHAIN_GENERATE_ODV                          0
#define UQ_MH_SG_FILTERED_CHAIN_DISCARDED_PORTION_ODV                 0.
#define UQ_MH_SG_FILTERED_CHAIN_LAG_ODV                               1
#define UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_FILE_NAME_ODV             UQ_MH_SG_FILENAME_FOR_NO_FILE
#define UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_FILE_TYPE_ODV             UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT
#define UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_ALLOW_ALL_ODV             0
#define UQ_MH_SG_FILTERED_CHAIN_DATA_OUTPUT_ALLOWED_SET_ODV           ""
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
#define UQ_MH_SG_FILTERED_CHAIN_COMPUTE_STATS_ODV                     0
#endif
#define UQ_MH_SG_DISPLAY_CANDIDATES_ODV                               0
#define UQ_MH_SG_PUT_OUT_OF_BOUNDS_IN_CHAIN_ODV                       1
#define UQ_MH_SG_TK_USE_LOCAL_HESSIAN_ODV                             0
#define UQ_MH_SG_TK_USE_NEWTON_COMPONENT_ODV                          1
#define UQ_MH_SG_DR_MAX_NUM_EXTRA_STAGES_ODV                          0
#define UQ_MH_SG_DR_LIST_OF_SCALES_FOR_EXTRA_STAGES_ODV               ""
#define UQ_MH_SG_DR_DURING_AM_NON_ADAPTIVE_INT_ODV                    1
#define UQ_MH_SG_AM_KEEP_INITIAL_MATRIX_ODV                           0
#define UQ_MH_SG_AM_INIT_NON_ADAPT_INT_ODV                            0
#define UQ_MH_SG_AM_ADAPT_INTERVAL_ODV                                0
#define UQ_MH_SG_AM_ADAPTED_MATRICES_DATA_OUTPUT_PERIOD_ODV           0
#define UQ_MH_SG_AM_ADAPTED_MATRICES_DATA_OUTPUT_FILE_NAME_ODV        UQ_MH_SG_FILENAME_FOR_NO_FILE
#define UQ_MH_SG_AM_ADAPTED_MATRICES_DATA_OUTPUT_FILE_TYPE_ODV        UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT
#define UQ_MH_SG_AM_ADAPTED_MATRICES_DATA_OUTPUT_ALLOW_ALL_ODV        0
#define UQ_MH_SG_AM_ADAPTED_MATRICES_DATA_OUTPUT_ALLOWED_SET_ODV      ""
#define UQ_MH_SG_AM_ETA_ODV                                           1.
#define UQ_MH_SG_AM_EPSILON_ODV                                       1.e-5
#define UQ_MH_SG_ENABLE_BROOKS_GELMAN_CONV_MONITOR                    0
#define UQ_MH_SG_BROOKS_GELMAN_LAG                                    100
#define UQ_MH_SG_OUTPUT_LOG_LIKELIHOOD                                1
#define UQ_MH_SG_OUTPUT_LOG_TARGET                                    1
#define UQ_MH_SG_DO_LOGIT_TRANSFORM                                   1

namespace QUESO {

/*! \file MetropolisHastingsSGOptions.h
    \brief Classes to allow options to be passed to a Metropolis-Hastings algorithm.
*/

/*!
 * \class MhOptionsValues
 * \brief This class provides options for the Metropolis-Hastings generator of samples if no input file is available.
 *
 * Metropolis-Hastings generator of samples expects some options for its
 * methods to be fully defined.  This class provides default values for such
 * options if no input file is available.
 */

class MhOptionsValues
{
public:

  //! @name Constructor/Destructor methods
  //@{

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  MhOptionsValues            (const SsOptionsValues* alternativeRawSsOptionsValues,
                                     const SsOptionsValues* alternativeFilteredSsOptionsValues);
#else

  //! Default constructor.
  /*! Assigns the default suite of options to the Metropolis-Hastings generator of samples.*/
  MhOptionsValues            ();
#endif

  //! Copy constructor.
  /*! It assigns the same options values from  \c src to \c this.*/
  MhOptionsValues            (const MhOptionsValues& src);

  //! Destructor
  ~MhOptionsValues            ();
  //@}

  //! @name Set methods
  //@{
  //! Assignment operator; it copies \c rhs to \c this.
  MhOptionsValues& operator= (const MhOptionsValues& rhs);
  //@}

  std::string                        m_dataOutputFileName;
  bool                               m_dataOutputAllowAll;
  std::set<unsigned int>             m_dataOutputAllowedSet;

  bool                               m_totallyMute;
  std::string                        m_initialPositionDataInputFileName;
  std::string                        m_initialPositionDataInputFileType;
  std::string                        m_initialProposalCovMatrixDataInputFileName;
  std::string                        m_initialProposalCovMatrixDataInputFileType;
  std::set<unsigned int>             m_parameterDisabledSet;  // gpmsa2
  std::string                        m_rawChainDataInputFileName;
  std::string                        m_rawChainDataInputFileType;
  unsigned int                       m_rawChainSize;
  bool                               m_rawChainGenerateExtra;
  unsigned int                       m_rawChainDisplayPeriod;
  bool                               m_rawChainMeasureRunTimes;
  unsigned int                       m_rawChainDataOutputPeriod;
  std::string                        m_rawChainDataOutputFileName;
  std::string                        m_rawChainDataOutputFileType;
  bool                               m_rawChainDataOutputAllowAll;
  std::set<unsigned int>             m_rawChainDataOutputAllowedSet;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  bool                               m_rawChainComputeStats;
#endif

  //! Toggle the option to save a filtered chain.  Default is 0 (off).
  /*!
   * A filtered chain is one where only every k-th sample is saved.  Here k
   * is called the 'lag' and can be set through m_filteredChainLag.
   */
  bool                               m_filteredChainGenerate;

  //! What initial fraction of the filtered chain is discarded.  Default is 0.
  /*!
   * For example, if set to 0.2 then the first 20% of the filtered chain is
   * discarded.  Useful for discarding 'burn-in'.
   */
  double                             m_filteredChainDiscardedPortion; // input or set during run time

  //! Set the lag for the filtered chain.  Default is 1.
  /*!
   * The lag is the number of samples to compute before saving.  For example,
   * if the lag is set to 20, then only every 20th sample is saved.
   */
  unsigned int                       m_filteredChainLag;              // input or set during run time

  //! File name to save the filtered chain to.  Default is ".".
  std::string                        m_filteredChainDataOutputFileName;
  std::string                        m_filteredChainDataOutputFileType;
  bool                               m_filteredChainDataOutputAllowAll;
  std::set<unsigned int>             m_filteredChainDataOutputAllowedSet;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  bool                               m_filteredChainComputeStats;
#endif

  bool                               m_displayCandidates;
  bool                               m_putOutOfBoundsInChain;
  bool                               m_tkUseLocalHessian;
  bool                               m_tkUseNewtonComponent;
  unsigned int                       m_drMaxNumExtraStages;
  std::vector<double>                m_drScalesForExtraStages;
  bool                               m_drDuringAmNonAdaptiveInt;
  bool                               m_amKeepInitialMatrix;
  unsigned int                       m_amInitialNonAdaptInterval;
  unsigned int                       m_amAdaptInterval;
  unsigned int                       m_amAdaptedMatricesDataOutputPeriod;
  std::string                        m_amAdaptedMatricesDataOutputFileName;
  std::string                        m_amAdaptedMatricesDataOutputFileType;
  bool                               m_amAdaptedMatricesDataOutputAllowAll;
  std::set<unsigned int>             m_amAdaptedMatricesDataOutputAllowedSet;

  /*! \brief Proposal covariance scaling factor, usually 2.4 * 2.4 / d
   *
   * This is a parameter in the DRAM algorithm and can be found in Haario et al
   * (2006).
   *
   * The parameter defines how much the proposal covariance matrix is to be
   * scaled by, and should usually be set to 2.4 * 2.4 / d, where d is the
   * dimension of the state space being sampled.
   */
  double                             m_amEta;

  /*! \brief Regularisation parameter for the DRAM covariance matrix
   *
   * This is a parameter in the DRAM algorithm that regularises the proposal
   * covariance matrix.  Details can be found in Haario et al (2006).
   *
   * The parameter defines how much the diagonal of the proposal covariance
   * matrix is perturbed.  Usually this is small, of order 1e-5.
   */
  double                             m_amEpsilon;

  unsigned int                       m_enableBrooksGelmanConvMonitor;
  unsigned int                       m_BrooksGelmanLag;

  //! Flag for deciding whether or not to dump log likelihood values in output
  bool m_outputLogLikelihood;

  //! Flag for deciding whether or not to dump log target values in output
  bool m_outputLogTarget;

  //! Flag for deciding whether or not to do logit transform of bounded domains
  bool m_doLogitTransform;

private:
  //! Copies the option values from \c src to \c this.
  void copy(const MhOptionsValues& src);

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  friend class MetropolisHastingsSGOptions;
  SsOptionsValues             m_alternativeRawSsOptionsValues;
  SsOptionsValues             m_alternativeFilteredSsOptionsValues;
#endif
};

/*! \class MetropolisHastingsSGOptions
 *  \brief This class reads the options for the Metropolis-Hastings generator of samples from an input file.
 *
 * This class implements a Metropolis-Hastings generator of samples.  'SG'
 * stands for 'Sequence Generator'.  Metropolis-Hastings generator of samples
 * expects some options to be fully defined.  This class reads the options for
 * the Metropolis-Hastings generator of samples from an input file provided by
 * the user.  The class expects the prefix '\<prefix\>_mh_'.  For instance, if
 * 'prefix' is 'foo_775_fp_', then the constructor will read all options that
 * begin with 'foo_775_fp_mh_'.  Options reading is hpandled by class
 * 'MetropolisHastingsOptions'.  To set options by hand, use the \c
 * MhOptionsValues class.
 */

class MetropolisHastingsSGOptions
{
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Constructor: reads options from the input file.
  MetropolisHastingsSGOptions(const BaseEnvironment& env, const char* prefix);

  //! Constructor: with alternative option values.
  /*! In this constructor, the input options are given by \c alternativeOptionsValues.*/
  MetropolisHastingsSGOptions(const BaseEnvironment& env, const char* prefix, const MhOptionsValues& alternativeOptionsValues);

  //! Copy constructor
  MetropolisHastingsSGOptions(const MLSamplingLevelOptions& mlOptions);

  //! Destructor
  ~MetropolisHastingsSGOptions();
  //@}

  //! @name I/O methods
  //@{
  //! It scans the option values from the options input file.
  void scanOptionsValues();

  //!  It prints the option values.
  void print            (std::ostream& os) const;
  //@}

  //! This class is where the actual options are stored
  MhOptionsValues             m_ov;

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  SequenceStatisticalOptions* m_rawChainStatisticalOptionsObj;
  bool                               m_rawChainStatOptsInstantiated;
  SequenceStatisticalOptions* m_filteredChainStatisticalOptionsObj;
  bool                               m_filteredChainStatOptsInstantiated;
#endif
  std::string                        m_prefix;

private:
  //! Defines the options for the Metropolis-Hastings generator of samples as the default options.
  void   defineMyOptions  (po::options_description& optionsDesc) const;

  //! Gets the sequence options defined to the  Metropolis-Hastings algorithm.
  void   getMyOptionValues(po::options_description& optionsDesc);

  const BaseEnvironment& m_env;
  po::options_description*      m_optionsDesc;

  std::string                   m_option_help;

  std::string                   m_option_dataOutputFileName;
  std::string                   m_option_dataOutputAllowAll;
  std::string                   m_option_dataOutputAllowedSet;

  std::string                   m_option_totallyMute;
  std::string                   m_option_initialPosition_dataInputFileName;
  std::string                   m_option_initialPosition_dataInputFileType;
  std::string                   m_option_initialProposalCovMatrix_dataInputFileName;
  std::string                   m_option_initialProposalCovMatrix_dataInputFileType;
  std::string                   m_option_listOfDisabledParameters;  // gpmsa2
  std::string                   m_option_rawChain_dataInputFileName;
  std::string                   m_option_rawChain_dataInputFileType;
  std::string                   m_option_rawChain_size;
  std::string                   m_option_rawChain_generateExtra;
  std::string                   m_option_rawChain_displayPeriod;
  std::string                   m_option_rawChain_measureRunTimes;
  std::string                   m_option_rawChain_dataOutputPeriod;
  std::string                   m_option_rawChain_dataOutputFileName;
  std::string                   m_option_rawChain_dataOutputFileType;
  std::string                   m_option_rawChain_dataOutputAllowAll;
  std::string                   m_option_rawChain_dataOutputAllowedSet;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  std::string                   m_option_rawChain_computeStats;
#endif
  std::string                   m_option_filteredChain_generate;
  std::string                   m_option_filteredChain_discardedPortion;
  std::string                   m_option_filteredChain_lag;
  std::string                   m_option_filteredChain_dataOutputFileName;
  std::string                   m_option_filteredChain_dataOutputFileType;
  std::string                   m_option_filteredChain_dataOutputAllowAll;
  std::string                   m_option_filteredChain_dataOutputAllowedSet;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  std::string                   m_option_filteredChain_computeStats;
#endif
  std::string                   m_option_displayCandidates;
  std::string                   m_option_putOutOfBoundsInChain;
  std::string                   m_option_tk_useLocalHessian;
  std::string                   m_option_tk_useNewtonComponent;
  std::string                   m_option_dr_maxNumExtraStages;
  std::string                   m_option_dr_listOfScalesForExtraStages;
  std::string                   m_option_dr_duringAmNonAdaptiveInt;
  std::string                   m_option_am_keepInitialMatrix;
  std::string                   m_option_am_initialNonAdaptInterval;
  std::string                   m_option_am_adaptInterval;
  std::string                   m_option_am_adaptedMatrices_dataOutputPeriod;
  std::string                   m_option_am_adaptedMatrices_dataOutputFileName;
  std::string                   m_option_am_adaptedMatrices_dataOutputFileType;
  std::string                   m_option_am_adaptedMatrices_dataOutputAllowAll;
  std::string                   m_option_am_adaptedMatrices_dataOutputAllowedSet;

  //! See MhOptionsValues::m_amEta
  std::string                   m_option_am_eta;
  //! See MhOptionsValues::m_amEpsilon
  std::string                   m_option_am_epsilon;

  std::string                   m_option_enableBrooksGelmanConvMonitor;
  std::string                   m_option_BrooksGelmanLag;

  std::string                   m_option_outputLogLikelihood;
  std::string                   m_option_outputLogTarget;
  std::string                   m_option_doLogitTransform;
};

std::ostream& operator<<(std::ostream& os, const MetropolisHastingsSGOptions& obj);

}  // End namespace QUESO

#endif // UQ_MH_SG_OPTIONS_H
