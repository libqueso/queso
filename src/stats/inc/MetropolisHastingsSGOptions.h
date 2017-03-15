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

#ifndef UQ_MH_SG_OPTIONS_H
#define UQ_MH_SG_OPTIONS_H

#include <queso/Environment.h>
#include <queso/MLSamplingLevelOptions.h>
#include <queso/SequenceStatisticalOptions.h>

#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
#include <queso/BoostInputOptionsParser.h>
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

#undef  UQ_MH_SG_REQUIRES_INVERTED_COV_MATRICES
#define UQ_NOTHING_JUST_FOR_TEST_OF_SVN_ID 1

#define UQ_MH_SG_FILENAME_FOR_NO_FILE   "."

// _ODV = option default value
#define UQ_MH_SG_HELP                                                 ""
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
#define UQ_MH_SG_ALGORITHM                                            "logit_random_walk"
#define UQ_MH_SG_TK                                                   "logit_random_walk"
#define UQ_MH_SG_UPDATE_INTERVAL                                      1

#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
namespace boost {
  namespace program_options {
    class options_description;
  }
}
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

namespace QUESO {

class BaseEnvironment;

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
  MhOptionsValues            (const SsOptionsValues* alternativeRawSsOptionsValues,
                                     const SsOptionsValues* alternativeFilteredSsOptionsValues,
                                     const BaseEnvironment * env, const char * prefix);
#else

  //! Default constructor.
  /*! Assigns the default suite of options to the Metropolis-Hastings generator of samples.*/
  MhOptionsValues            ();

  //! Prefix constructor for input file parsing purposes
  MhOptionsValues(const BaseEnvironment * env, const char * prefix);
#endif

  //! Copy constructor.
  /*! It assigns the same options values from  \c src to \c this.*/
  MhOptionsValues            (const MhOptionsValues& src);

  //! Destructor
  virtual ~MhOptionsValues            ();
  //@}

  //! @name Set methods
  //@{
  //! Assignment operator; it copies \c rhs to \c this.
  MhOptionsValues& operator= (const MhOptionsValues& rhs);
  //@}

  //! Prefix for input file option names.  Prepends all options for this class.
  std::string                        m_prefix;

  //! If non-empty string, print options and values to the output file
  /*!
   * Default is empty string
   */
  std::string m_help;

  //! The base name of output files where the chain (and related information) will be written.
  /*!
   * For multiple environments, the respective chains will append "_subN"
   * where N is the environment number.
   *
   * For log likelihood, log target, and log prior, a further "_loglikelihood",
   * "_logtarget" and "_logprior" is appended respectively
   *
   * Default value is "."
   */
  std::string                        m_dataOutputFileName;

  //! If true, all processes write output and m_dataOutputAllowedSet is ignored
  /*!
   * If false, m_dataOutputAllowedSet determines the set of MPI ranks that
   * can write output
   *
   * Default is false
   */
  bool                               m_dataOutputAllowAll;

  //! The set of MPI ranks that can write output.  See m_dataOutputAllowAll
  /*!
   * Default value is the empty set
   */
  std::set<unsigned int>             m_dataOutputAllowedSet;

  //! If true, zero output is written to files.  Default is false.
  bool                               m_totallyMute;

  //! If not ".", reads the contents of the file and uses that to start the MCMC.  Default is "."
  /*!
   * If ".", the input file is not read.
   */
  std::string                        m_initialPositionDataInputFileName;

  //! The filetype of m_initialPositionDataInputFileName.  Only "m" (matlab) is currently supported.  Default is "m"
  std::string                        m_initialPositionDataInputFileType;

  //! If not ".", reads the contents of the file as the initial proposal covariance matrix.
  /*!
   * To use this, m_tkUseLocalHessian must be false.  If it is true, and this
   * option is ".", no input file is read and the user-provided initial
   * proposal covariance matrix is used.
   *
   * If this option is ".", no input file is read.
   *
   * Default is "."
   */
  std::string                        m_initialProposalCovMatrixDataInputFileName;

  //! The filetype of m_initialProposalCovMatrixDataInputFileName.  Only "m" (matlab) is currently supported.  Default is "m"
  std::string                        m_initialProposalCovMatrixDataInputFileType;

  //! Set of parameters that don't get sampled
  /*!
   * Default is empty set
   */
  std::set<unsigned int>             m_parameterDisabledSet;  // gpmsa2

  //! Filename for reading an already-produced Markov chain
  /*!
   * If ".", generates the Markov chain without reading anything.
   * If not ".", this is the name of the file to read a Markov chain from
   * instead of generating a chain.  Filtering and MLE/MAP computations are
   * still done.
   *
   * Default is "."
   */
  std::string                        m_rawChainDataInputFileName;

  //! The filetype of m_rawChainDataInputFileName.  Only "m" (matlab) is currently supported.  Default is "m"
  std::string                        m_rawChainDataInputFileType;

  //! The size of the chain (number of posterior samples) to generate.  Default is 100
  unsigned int                       m_rawChainSize;

  //! If true, extra chain information is computed/stored
  /*!
   * The extra information written is: log target, acceptance ratio, number
   * of rejected samples.
   *
   * Default is false
   */
  bool                               m_rawChainGenerateExtra;

  //! The frequency with which to output diagnostic information
  /*!
   * Diagnostics are current sample iteration and current rejection percentage.
   *
   * Must be an integer.
   *
   * Default is 500
   */
  unsigned int                       m_rawChainDisplayPeriod;

  //! If true, measures timings spent in various chain computions and writes them to the output file
  /*!
   * The measurements are:
   *
   * time spent computing proposal
   * time spent computing target
   * time spent computing metropolis
   * time spent computing hastings acceptance ratio
   * time spent computing delayed rejection
   * average time spent computing target
   * total delayed rejection run time
   * total adaptive metropolis run time
   *
   * Default is true
   */
  bool                               m_rawChainMeasureRunTimes;

  //! The frequency with which to write chain output.  Defaults to 0.
  unsigned int                       m_rawChainDataOutputPeriod;

  //! If not ".", filename to write the Markov chain to
  /*!
   * If ".", the chain is not written to a file.
   *
   * If m_totallyMute is true, the chain is not written regardless of the value
   * of this option.
   *
   * Default is "."
   */
  std::string                        m_rawChainDataOutputFileName;

  //! The filetype of m_rawChainDataOutputFileName
  /*!
   * Only "m" (matlab) and "txt" are currently supported.
   *
   * If "m", the raw chain data that is written will be matlab-friendly.  I.e.,
   * arrays will be set up with the correct dimensions.  The first dimension
   * is the sample dimension and the second dimension (for high-dimensional
   * statistical inverse problems) is the parameter dimension.
   *
   * If "txt", none of the matlab-specific information is written to the file.
   * Instead, space-delimited column data is written.  The first line
   * specifies the dimensions of the saved array data.  The first dimension
   * is the sample dimension and the second dimension is the parameter
   * dimension.  For example, if the first line written is "12 34" then the
   * following 12 lines are samples that live in a 34-dimensional space.
   *
   * Default is "m"
   */
  std::string                        m_rawChainDataOutputFileType;

  //! Toggle for whether or not to allow all processes to write Markov chain output to a file
  /*!
   * If true, all processes write Markov chain output and
   * m_rawChainDataOutputAllowedSet is ignored.
   *
   * If false, m_rawChainDataOutputAllowedSet determines the set of
   * processor ranks that will write Markov chain output to a file
   *
   * Default is false
   */
  bool                               m_rawChainDataOutputAllowAll;

  //! The set of MPI ranks that will write Markov chain output to a file.  See also m_rawChainDataOutputAllowAll.  Default is empty set.
  std::set<unsigned int>             m_rawChainDataOutputAllowedSet;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  //! Flag to tell QUESO whether or not to compute chain statistics.  Default is false.
  bool                               m_rawChainComputeStats;
#endif

  //! Toggle the option to save a filtered chain.
  /*!
   * A filtered chain is one where only every k-th sample is saved.  Here k
   * is called the 'lag' and can be set through m_filteredChainLag.
   *
   * Default is false
   */
  bool                               m_filteredChainGenerate;

  //! What initial fraction of the filtered chain is discarded.
  /*!
   * For example, if set to 0.2 then the first 20% of the filtered chain is
   * discarded.  Useful for discarding 'burn-in'.
   *
   * Default is 0.0
   */
  double                             m_filteredChainDiscardedPortion; // input or set during run time

  //! Set the lag for the filtered chain.  Default is 1.
  /*!
   * The lag is the number of samples to compute before saving.  For example,
   * if the lag is set to 20, then only every 20th sample is saved.
   */
  unsigned int                       m_filteredChainLag;              // input or set during run time

  //! If not ".", file name to save the filtered chain to.  Default is ".".
  /*!
   * If ".", the filtered chain is not written to a file.
   *
   * If m_totallyMute is true, then the chain is not written to a file
   * regardless of the value of this option.
   */
  std::string                        m_filteredChainDataOutputFileName;

  //! The filetype of m_filteredChainDataOutputFileName.  Only "m" (matlab) is currently supported.  Default is "m"
  std::string                        m_filteredChainDataOutputFileType;

  //! Toggle for whether or not to allow all processes to write *filtered* Markov chain output to a file
  /*!
   * If true, all processes write filtered Markov chain output and
   * m_filteredChainDataOutputAllowedSet is ignored.
   *
   * If false, m_filteredChainDataOutputAllowedSet determines the set of
   * processor ranks that will write filtered Markov chain output to a file
   *
   * Default is false
   */
  bool                               m_filteredChainDataOutputAllowAll;

  //! The set of MPI ranks that will write filtered Markov chain output to a file.  See also m_filteredChainDataOutputAllowAll.  Default is empty set.
  std::set<unsigned int>             m_filteredChainDataOutputAllowedSet;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  //! Toggle to tell QUESO whether or not to compute statistics on the filtered chain.  Default is false
  bool                               m_filteredChainComputeStats;
#endif

  //! Toggle to tell QUESO whether or not to write proposal (candidate) state to output file.  Default is false
  bool                               m_displayCandidates;

  //! Flag to tell QUESO how chains should be upon generating a proposal that is out of the problem domain.
  /*!
   * If true, the chain will reject any proposal states outside of the
   * problem domain and the chain will *advance* to the next iteration, staying
   * where it is.
   *
   * If false, the chain will not do an accept-reject on a proposed state that
   * is out of bounds.  Instead, QUESO will regenerate proposal until it is in
   * the domain.  Only when a proposal is in the problem domain does QUESO do
   * an accept-reject step.
   *
   * The default is true
   */
  bool                               m_putOutOfBoundsInChain;

  //! Flag to tell QUESO whether or not to use Hessian information for the proposal covariance matrix
  /*!
   * The interaction with delayed rejection:
   *
   * If m_drDuringAmNonAdaptiveInt is false and m_tkUseLocalHessian is false
   * and m_amInitialNonAdaptInterval > 0 and m_amAdaptInterval is > 0 and
   * the current sampler iteration is <= m_amInitialNonAdaptInterval, then
   * no delayed rejection is done.
   *
   * The above is hard to parse, it essentially means that, unless you ask for
   * it, no delayed rejection is dont in the initial phase of sampling when no
   * adaptive Metropolis is happening.
   *
   * The interaction with adaptive Metropolis:
   *
   * If m_tkUseLocalHessian is true, no adaptive Metropolis is done.  If it's
   * false, then we also need that m_amInitialNonAdaptInterval > 0 and
   * m_amAdaptInterval > 0 otherwise no adaptive Metropolis is done.
   *
   * The default is false.
   *
   * DM: I think m_tkUseLocalHessian == true is not supported, since we have no
   * algorithms that leverage Hessian information right now.
   */
  bool                               m_tkUseLocalHessian;

  //! This option is a no-op.  Default is true.
  bool                               m_tkUseNewtonComponent;

  //! The number of delayed rejection stages to do.  Default is 0
  /*!
   * Delayed rejection happens when a proposal is about to be rejected.
   *
   * If m_drMaxNumExtraStages is zero, no delayed rejection is done.
   *
   * Default is 0
   */
  unsigned int                       m_drMaxNumExtraStages;

  //! The vector of scale factors for the proposal covariance matrix to use for delayed rejection.
  /*!
   * Example:
   *
   * Assume the user provides "2.0 3.0 5.0" as their list of scales.  This
   * means QUESO will multiply the proposal covariance matrix for delayed
   * rejection stage 1 by 1.0 / (2.0 * 2.0).  For delayed rejection stage 2
   * the proposal covariance matrix will be multiplied by 1.0 / (3.0 * 3.0).
   * And for stage 3 it will be multiplied by 1.0 / (5.0 * 5.0).
   *
   * Generally, the user-provided list of scales should be *increasing*, but
   * this is not checked for.
   *
   * QUESO will prepend a scale of 1.0 to the user-provided list.
   *
   * The default is the empty list.
   */
  std::vector<double>                m_drScalesForExtraStages;

  //! Do delayed rejection during the initial non-adaptive part of sampling?
  /*!
   * This option interacts with the following other options:
   *
   * m_tkUseLocalHessian
   * m_amInitialNonAdaptInterval
   * m_amAdaptInterval
   *
   *
   * The interaction is as follows.  If this option is false, then as long as
   * m_tkUseLocalHessian is false, m_amInitialNonAdaptInterval is > 0,
   * m_amAdaptInterval is > 0, and the current sampler iteration is <=
   * m_amInitialNonAdaptInterval *NO* delayed rejection is done.
   *
   * If this option is true (regardless of the other above conditions) then
   * delayed rejection is done (including the initial non-adaptive part of
   * sampling).
   *
   * The default is true.
   */
  bool                               m_drDuringAmNonAdaptiveInt;

  //! This option is a no-op.  The default is false.
  bool                               m_amKeepInitialMatrix;

  //! The number of initial samples to do without adapting the proposal covariance matrix
  /*!
   * If positive and the current sampler iteration is <=
   * m_amInitialNonAdaptInterval, then no delayed rejection is done if
   * m_drDuringAmNonAdaptiveInt is false and m_tkUseLocalHessian is false
   * and m_amAdaptInterval is > 0.
   *
   * The default is 0.
   */
  unsigned int                       m_amInitialNonAdaptInterval;

  //! The frequency at which to adapt the proposal covariance matrix.
  /*!
   * If zero, no adaptive metropolis is done.  If positive, the proposal
   * covariance matrix will be adapted every m_amAdaptInterval samples
   * after m_amInitialNonAdaptInterval have been done.
   *
   * The default is 0
   */
  unsigned int                       m_amAdaptInterval;

  //! The frequency (after m_amInitialNonAdaptInterval samples are done) of printing the last adapted proposal covariance matrix.
  /*!
   * Provided m_amAdaptedMatricesDataOutputFileName is not ".", the last
   * adapted proposal covariance matrix will be written regardless of the value
   * of this option if the current sampler iteration is exactly
   * m_amInitialNonAdaptInterval.
   *
   * Beyond m_amInitialNonAdaptInterval, the last adapted proposal covariance
   * matrix will be printed to the output file every
   * m_amAdaptedMatricesDataOutputPeriod samples.
   *
   * If zero, no output is written.
   *
   * The default is 0
   */
  unsigned int                       m_amAdaptedMatricesDataOutputPeriod;

  //! If not ".", this is the file to write adapted proposal covariance matrices to.  Default is "."
  /*!
   * If ".", the adapted proposal covariance matrices are not written to a file
   */
  std::string                        m_amAdaptedMatricesDataOutputFileName;

  //! The filetype of m_amAdaptedMatricesDataOutputFileName.  Only "m" (matlab) is currently supported.  Default is "m"
  std::string                        m_amAdaptedMatricesDataOutputFileType;

  //! This option is a no-op.  The default is false.
  bool                               m_amAdaptedMatricesDataOutputAllowAll;

  //! This option is a no-op.  The default is the empty set.
  std::set<unsigned int>             m_amAdaptedMatricesDataOutputAllowedSet;

  /*! \brief Proposal covariance scaling factor, usually 2.4 * 2.4 / d
   *
   * This is a parameter in the DRAM algorithm and can be found in Haario et al
   * (2006).
   *
   * The parameter defines how much the proposal covariance matrix is to be
   * scaled by, and should usually be set to 2.4 * 2.4 / d, where d is the
   * dimension of the state space being sampled.
   *
   * The default is 1.0.
   */
  double                             m_amEta;

  /*! \brief Regularisation parameter for the DRAM covariance matrix
   *
   * This is a parameter in the DRAM algorithm that regularises the proposal
   * covariance matrix.  Details can be found in Haario et al (2006).
   *
   * The parameter defines how much the diagonal of the proposal covariance
   * matrix is perturbed.  Usually this is small, of order 1e-5.
   *
   * The default is 1.e-5
   */
  double                             m_amEpsilon;

  //! The frequency with which to compute the Brooks-Gelman convergence statistic.
  /*!
   * If zero, it is not computed.
   *
   * If positive, the Brooks-Gelman convergence statistic is computed every
   * m_enableBrooksGelmanConvMonitor samples as long as the current sampler
   * iteration is larger than m_BrooksGelmanLag + 1.  The "+ 1" is to ensure
   * there are at least two samples with which to compute the convergence
   * statistic.
   *
   * Needs at least 2 chains (sub environments) to be computed.
   *
   * The default is 0
   */
  unsigned int                       m_enableBrooksGelmanConvMonitor;

  //! The lag with which to compute the Brooks-Gelman convergence statistic
  /*!
   * The convergence statistic will be computed from sampler iteration
   * m_BrooksGelmanLag to sampler iteration
   * (current_sampler_iteration - m_BrooksGelmanLag).
   *
   * The default is 100.
   */
  unsigned int                       m_BrooksGelmanLag;

  //! Flag for deciding whether or not to dump log likelihood values in output.  Default is true.
  bool m_outputLogLikelihood;

  //! Flag for deciding whether or not to dump log target values in output Default is true.
  bool m_outputLogTarget;

  //! Flag for deciding whether or not to do logit transform of bounded domains Default is true.
  bool m_doLogitTransform;

  //! Which algorithm to use for the MCMC.  Default is "random_walk"
  std::string m_algorithm;

  //! Which transition kernel to use for MCMC.  Default is "random_walk"
  std::string m_tk;

  //! How often to call the TK's updateTK method.  Default is 1.
  unsigned int m_updateInterval;

private:
  // Cache a pointer to the environment.
  const BaseEnvironment * m_env;

#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  ScopedPtr<BoostInputOptionsParser>::Type m_parser;
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

  //! Option name for MhOptionsValues::m_help.  Option name is m_prefix + "mh_help"
  std::string                   m_option_help;

  //! Option name for MhOptionsValues::m_dataOutputFileName.  Option name is m_prefix + "mh_dataOutputFileName"
  std::string                   m_option_dataOutputFileName;
  //! Option name for MhOptionsValues::m_dataOutputAllowAll.  Option name is m_prefix + "mh_dataOutputAllowAll"
  std::string                   m_option_dataOutputAllowAll;
  //! Option name for MhOptionsValues::m_dataOutputAllowedSet.  Option name is m_prefix + "mh_dataOutputAllowedSet"
  std::string                   m_option_dataOutputAllowedSet;

  //! Option name for MhOptionsValues::m_totallyMute.  Option name is m_prefix + "mh_totallyMute"
  std::string                   m_option_totallyMute;
  //! Option name for MhOptionsValues::m_initialPositionDataInputFileName.  Option name is m_prefix + "mh_initialPosition_dataInputFileName"
  std::string                   m_option_initialPosition_dataInputFileName;
  //! Option name for MhOptionsValues::m_initialPositionDataInputFileType.  Option name is m_prefix + "mh_initialPosition_dataInputFileType"
  std::string                   m_option_initialPosition_dataInputFileType;
  //! Option name for MhOptionsValues::m_initialProposalCovMatrixDataInputFileName.  Option name is m_prefix + "mh_initialProposalCovMatrix_dataInputFileName"
  std::string                   m_option_initialProposalCovMatrix_dataInputFileName;
  //! Option name for MhOptionsValues::m_initialProposalCovMatrixDataInputFileType.  Option name is m_prefix + "mh_initialProposalCovMatrix_dataInputFileType"
  std::string                   m_option_initialProposalCovMatrix_dataInputFileType;
  //! Option name for MhOptionsValues::m_parameterDisabledSet.  Option name is m_prefix + "mh_listOfDisabledParameters"
  std::string                   m_option_listOfDisabledParameters;  // gpmsa2
  //! Option name for MhOptionsValues::m_rawChainDataInputFileName.  Option name is m_prefix + "mh_rawChain_dataInputFileName"
  std::string                   m_option_rawChain_dataInputFileName;
  //! Option name for MhOptionsValues::m_rawChainDataInputFileType.  Option name is m_prefix + "mh_rawChain_dataInputFileType"
  std::string                   m_option_rawChain_dataInputFileType;
  //! Option name for MhOptionsValues::m_rawChainSize.  Option name is m_prefix + "mh_rawChain_size"
  std::string                   m_option_rawChain_size;
  //! Option name for MhOptionsValues::m_rawChainGenerateExtra.  Option name is m_prefix + "mh_rawChain_generateExtra"
  std::string                   m_option_rawChain_generateExtra;
  //! Option name for MhOptionsValues::m_rawChainDisplayPeriod.  Option name is m_prefix + "mh_rawChain_displayPeriod"
  std::string                   m_option_rawChain_displayPeriod;
  //! Option name for MhOptionsValues::m_rawChainMeasureRunTimes.  Option name is m_prefix + "mh_rawChain_measureRunTimes"
  std::string                   m_option_rawChain_measureRunTimes;
  //! Option name for MhOptionsValues::m_rawChainDataOutputPeriod.  Option name is m_prefix + "mh_rawChain_dataOutputPeriod"
  std::string                   m_option_rawChain_dataOutputPeriod;
  //! Option name for MhOptionsValues::m_rawChainDataOutputFileName.  Option name is m_prefix + "mh_rawChain_dataOutputFileName"
  std::string                   m_option_rawChain_dataOutputFileName;
  //! Option name for MhOptionsValues::m_rawChainDataOutputFileType.  Option name is m_prefix + "mh_rawChain_dataOutputFileType"
  std::string                   m_option_rawChain_dataOutputFileType;
  //! Option name for MhOptionsValues::m_rawChainDataOutputAllowAll.  Option name is m_prefix + "mh_rawChain_dataOutputAllowAll"
  std::string                   m_option_rawChain_dataOutputAllowAll;
  //! Option name for MhOptionsValues::m_rawChainDataOutputAllowedSet.  Option name is m_prefix + "mh_rawChain_dataOutputAllowedSet"
  std::string                   m_option_rawChain_dataOutputAllowedSet;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  //! Option name for MhOptionsValues::m_rawChainComputeStats.  Option name is m_prefix + "mh_rawChain_computeStats"
  std::string                   m_option_rawChain_computeStats;
#endif
  //! Option name for MhOptionsValues::m_filteredChainGenerate.  Option name is m_prefix + "mh_filteredChain_generate"
  std::string                   m_option_filteredChain_generate;
  //! Option name for MhOptionsValues::m_filteredChainDiscardedPortion.  Option name is m_prefix + "mh_filteredChain_discardedPortion"
  std::string                   m_option_filteredChain_discardedPortion;
  //! Option name for MhOptionsValues::m_filteredChainLag.  Option name is m_prefix + "mh_filteredChain_lag"
  std::string                   m_option_filteredChain_lag;
  //! Option name for MhOptionsValues::m_filteredChainDataOutputFileName.  Option name is m_prefix + "mh_filteredChain_dataOutputFileName"
  std::string                   m_option_filteredChain_dataOutputFileName;
  //! Option name for MhOptionsValues::m_filteredChainDataOutputFileType.  Option name is m_prefix + "mh_filteredChain_dataOutputFileType"
  std::string                   m_option_filteredChain_dataOutputFileType;
  //! Option name for MhOptionsValues::m_filteredChainDataOutputAllowAll.  Option name is m_prefix + "mh_filteredChain_dataOutputAllowAll"
  std::string                   m_option_filteredChain_dataOutputAllowAll;
  //! Option name for MhOptionsValues::m_filteredChainDataOutputAllowedSet.  Option name is m_prefix + "mh_filteredChain_dataOutputAllowedSet"
  std::string                   m_option_filteredChain_dataOutputAllowedSet;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  //! Option name for MhOptionsValues::m_filteredChainComputeStats.  Option name is m_prefix + "mh_filteredChain_computeStats"
  std::string                   m_option_filteredChain_computeStats;
#endif
  //! Option name for MhOptionsValues::m_displayCandidates.  Option name is m_prefix + "mh_displayCandidates"
  std::string                   m_option_displayCandidates;
  //! Option name for MhOptionsValues::m_putOutOfBoundsInChain.  Option name is m_prefix + "mh_putOutOfBoundsInChain"
  std::string                   m_option_putOutOfBoundsInChain;
  //! Option name for MhOptionsValues::m_tkUseLocalHessian.  Option name is m_prefix + "mh_tk_useLocalHessian"
  std::string                   m_option_tk_useLocalHessian;
  //! Option name for MhOptionsValues::m_tkUseNewtonComponent.  Option name is m_prefix + "mh_tk_useNewtonComponent"
  std::string                   m_option_tk_useNewtonComponent;
  //! Option name for MhOptionsValues::m_drMaxNumExtraStages.  Option name is m_prefix + "mh_dr_maxNumExtraStages"
  std::string                   m_option_dr_maxNumExtraStages;
  //! Option name for MhOptionsValues::m_drScalesForExtraStages.  Option name is m_prefix + "mh_dr_listOfScalesForExtraStages"
  std::string                   m_option_dr_listOfScalesForExtraStages;
  //! Option name for MhOptionsValues::m_drDuringAmNonAdaptiveInt.  Option name is m_prefix + "mh_dr_duringAmNonAdaptiveInt"
  std::string                   m_option_dr_duringAmNonAdaptiveInt;
  //! Option name for MhOptionsValues::m_amKeepInitialMatrix.  Option name is m_prefix + "mh_am_keepInitialMatrix"
  std::string                   m_option_am_keepInitialMatrix;
  //! Option name for MhOptionsValues::m_amInitialNonAdaptInterval.  Option name is m_prefix + "mh_am_initialNonAdaptInterval"
  std::string                   m_option_am_initialNonAdaptInterval;
  //! Option name for MhOptionsValues::m_amAdaptInterval.  Option name is m_prefix + "mh_am_adaptInterval"
  std::string                   m_option_am_adaptInterval;
  //! Option name for MhOptionsValues::m_amAdaptedMatricesDataOutputPeriod.  Option name is m_prefix + "mh_am_adaptedMatrices_dataOutputPeriod"
  std::string                   m_option_am_adaptedMatrices_dataOutputPeriod;
  //! Option name for MhOptionsValues::m_amAdaptedMatricesDataOutputFileName.  Option name is m_prefix + "mh_am_adaptedMatrices_dataOutputFileName"
  std::string                   m_option_am_adaptedMatrices_dataOutputFileName;
  //! Option name for MhOptionsValues::m_amAdaptedMatricesDataOutputFileType.  Option name is m_prefix + "mh_am_adaptedMatrices_dataOutputFileType"
  std::string                   m_option_am_adaptedMatrices_dataOutputFileType;
  //! Option name for MhOptionsValues::m_amAdaptedMatricesDataOutputAllowAll.  Option name is m_prefix + "mh_am_adaptedMatrices_dataOutputAllowAll"
  std::string                   m_option_am_adaptedMatrices_dataOutputAllowAll;
  //! Option name for MhOptionsValues::m_amAdaptedMatricesDataOutputAllowedSet.  Option name is m_prefix + "mh_am_adaptedMatrices_dataOutputAllowedSet"
  std::string                   m_option_am_adaptedMatrices_dataOutputAllowedSet;

  //! Option name for MhOptionsValues::m_amEta.  Option name is m_prefix + "mh_am_eta"
  std::string                   m_option_am_eta;
  //! Option name for MhOptionsValues::m_amEpsilon.  Option name is m_prefix + "mh_am_epsilon"
  std::string                   m_option_am_epsilon;

  //! Option name for MhOptionsValues::m_enableBrooksGelmanConvMonitor.  Option name is m_prefix + "mh_enableBrooksGelmanConvMonitor"
  std::string                   m_option_enableBrooksGelmanConvMonitor;
  //! Option name for MhOptionsValues::m_BrooksGelmanLag.  Option name is m_prefix + "mh_BrooksGelmanLag"
  std::string                   m_option_BrooksGelmanLag;

  //! Option name for MhOptionsValues::m_outputLogLikelihood.  Option name is m_prefix + "mh_outputLogLikelihood"
  std::string                   m_option_outputLogLikelihood;
  //! Option name for MhOptionsValues::m_outputLogTarget.  Option name is m_prefix + "mh_outputLogTarget"
  std::string                   m_option_outputLogTarget;
  //! Option name for MhOptionsValues::m_doLogitTransform.  Option name is m_prefix + "mh_doLogitTransform"
  std::string                   m_option_doLogitTransform;
  //! Option name for MhOptionsValues::m_algorithm.  Option name is m_prefix + "mh_algorithm"
  std::string                   m_option_algorithm;
  //! Option name for MhOptionsValues::m_tk.  Option name is m_prefix + "mh_tk"
  std::string                   m_option_tk;
  //! Option name for MhOptionsValues::m_updateInterval.  Option name is m_prefix + "mh_updateInterval"
  std::string                   m_option_updateInterval;

  //! Copies the option values from \c src to \c this.
  void copy(const MhOptionsValues& src);

  // We pass the the passed environment to get access to the MPI ranks etc for
  // sanity checks
  void checkOptions(const BaseEnvironment * env);

  friend std::ostream & operator<<(std::ostream & os,
      const MhOptionsValues & obj);

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
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  //! Defines the options for the Metropolis-Hastings generator of samples as the default options.
  void   defineMyOptions  (boost::program_options::options_description& optionsDesc) const;

  //! Gets the sequence options defined to the  Metropolis-Hastings algorithm.
  void   getMyOptionValues(boost::program_options::options_description& optionsDesc);
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

  const BaseEnvironment& m_env;
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  ScopedPtr<boost::program_options::options_description>::Type m_optionsDesc;
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

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
  std::string                   m_option_algorithm;
  std::string                   m_option_tk;
  std::string                   m_option_updateInterval;
};

std::ostream& operator<<(std::ostream& os, const MetropolisHastingsSGOptions& obj);

}  // End namespace QUESO

#endif // UQ_MH_SG_OPTIONS_H
