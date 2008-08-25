/* uq/libs/mcmc/inc/uqDRAM_MarkovChainGenerator.h
 *
 * Copyright (C) 2008 The PECOS Team, http://queso.ices.utexas.edu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __UQ_DRAM_MCG_H__
#define __UQ_DRAM_MCG_H__

#undef UQ_DRAM_MCG_REQUIRES_INVERTED_COV_MATRICES

// _ODV = option default value
#define UQ_MCMC_MH_CHAIN_SIZE_ODV              "100"
#define UQ_MCMC_DR_MAX_NUM_EXTRA_STAGES_ODV    0
#define UQ_MCMC_DR_SCALES_FOR_EXTRA_STAGES_ODV "1."
#define UQ_MCMC_AM_INIT_NON_ADAPT_INT_ODV      0
#define UQ_MCMC_AM_ADAPT_INTERVAL_ODV          0
#define UQ_MCMC_AM_ETA_ODV                     1.
#define UQ_MCMC_AM_EPSILON_ODV                 1.e-5
#define UQ_MCMC_MH_USE_CHAIN2                  0
#define UQ_MCMC_MH_GENERATE_UNIQUE_CHAIN_ODV   0
#define UQ_MCMC_MH_GENERATE_EXTRA_CHAINS_ODV   0
#define UQ_MCMC_MH_OUTPUT_FILE_NAME_ODV        "."
#define UQ_MCMC_MH_CHAIN_DISPLAY_PERIOD_ODV    500
#define UQ_MCMC_MH_MEASURE_RUN_TIMES_ODV       0
#define UQ_MCMC_MH_GENERATE_WHITE_NOISE_ODV    0
#define UQ_MCMC_FINAL_PERCENTS_FOR_STATS_ODV   "100."
#define UQ_MCMC_RUN_BMM_ODV                    0
#define UQ_MCMC_LENGTHS_FOR_BMM_ODV            "0"
#define UQ_MCMC_COMPUTE_PSDS_ODV               0
#define UQ_MCMC_COMPUTE_GEWEKE_COEFS_ODV       0
#define UQ_MCMC_COMPUTE_CORRELATIONS_ODV       0
#define UQ_MCMC_SECOND_LAG_FOR_CORRS_ODV       0
#define UQ_MCMC_LAG_SPACING_FOR_CORRS_ODV      0
#define UQ_MCMC_PRINT_CORRS_ODV                0
#define UQ_MCMC_WRITE_CORRS_ODV                0
#define UQ_MCMC_NUMBER_OF_LAGS_FOR_CORRS_ODV   0
#define UQ_MCMC_COMPUTE_HISTOGRAMS_ODV         0
#define UQ_MCMC_COMPUTE_KDES_ODV               0

#include <uqProbDensity.h>
#include <uqLikelihoodFunction.h>
#include <uqParamSpace.h>
#include <uqObservableSpace.h>
#include <uqChainPosition.h>
#include <uqMiscellaneous.h>
#include <uqArrayOfSequences.h>
#include <uqSequenceStatistics.h>
#include <uq2dArrayOfStuff.h>
#include <sys/time.h>
#include <fstream>

/*! A templated class that generates a Markov chain using the DRAM algorithm.
 */
template <class V, class M>
class uqDRAM_MarkovChainGeneratorClass
{
public:
  typedef typename std::vector<const V*>::iterator chainPositionIteratorTypedef;
  uqDRAM_MarkovChainGeneratorClass(const uqEnvironmentClass&                   env,                        /*! The QUESO toolkit environment.   */
                                   const char*                                 prefix,                     /*! Prefix for the validation phase. */
                                   const uqParamSpaceClass<V,M>&               paramSpace,                 /*! The parameter space.             */
                                   const uqObservableSpaceClass<V,M>&          observableSpace,            /*! The observable space.            */
                                   const uq_ProbDensity_BaseClass<V,M>&        m2lPriorProbDensity_Obj,    /*! -2*ln(prior())                   */
                                   const uq_LikelihoodFunction_BaseClass<V,M>& m2lLikelihoodFunction_Obj); /*! -2*ln(likelihood())              */
 ~uqDRAM_MarkovChainGeneratorClass();

  void generateChains             (const M* proposalCovMatrix,
                                   //const M* proposalPrecMatrix,
                                   const M* mahalanobisMatrix = NULL,
                                   bool     applyMahalanobisInvert = true);

  void print                      (std::ostream& os) const;

  //const std::vector<const V*>& chain              () const;
  //const std::vector<const V*>& misfitChain        () const;
  //const std::vector<const V*>& misfitVarianceChain() const;
  //const std::string&           outputFileName     () const;

private:
  void   resetChainAndRelatedInfo();
  void   defineMyOptions         (po::options_description& optionsDesc) const;
  void   getMyOptionValues       (po::options_description& optionsDesc);

  int    prepareForNextChain     (const M* proposalCovMatrix);
                                  //const M* proposalPrecMatrix,

  void   generateChains1         (const M* proposalCovMatrix,
                                  //const M* proposalPrecMatrix,
                                  const M* mahalanobisMatrix = NULL,
                                  bool     applyMahalanobisInvert = true);
  int    generateChain1          (unsigned int chainId,
                                  const V&     valuesOf1stPosition,
                                  const M*     proposalCovMatrix,
                                  const M*     mahalanobisMatrix = NULL,
                                  bool         applyMahalanobisInvert = true);
  int    generateWhiteNoise1     (unsigned int chainId);
  void   computeStatistics1      (const std::vector<const V*>& chain1, std::ofstream* passedOfs);
  void   updateCovMatrix1        (const std::vector<V*>& subChain1,
                                  unsigned int           idOfFirstPositionInSubChain,
                                  double&                lastChainSize,
                                  V&                     lastMean,
                                  M&                     lastAdaptedCovMatrix);
  int    writeChain1             (std::ofstream& ofs,
                                  const M*       mahalanobisMatrix = NULL,
                                  bool           applyMahalanobisInvert = true) const;

  void   generateChains2         (const M* proposalCovMatrix,
                                  //const M* proposalPrecMatrix,
                                  const M* mahalanobisMatrix = NULL,
                                  bool     applyMahalanobisInvert = true);
  int    generateChain2          (unsigned int chainId,
                                  const V&     valuesOf1stPosition,
                                  const M*     proposalCovMatrix,
                                  const M*     mahalanobisMatrix = NULL,
                                  bool         applyMahalanobisInvert = true);
  int    generateWhiteNoise2     (unsigned int chainId);
  void   computeStatistics2      (const uqArrayOfSequencesClass<V>& chain2, std::ofstream* passedOfs);
  void   updateCovMatrix2        (const uqArrayOfSequencesClass<V>& subChain2,
                                  unsigned int                      idOfFirstPositionInSubChain,
                                  double&                           lastChainSize,
                                  V&                                lastMean,
                                  M&                                lastAdaptedCovMatrix);
  int    writeChain2             (std::ofstream& ofs,
                                  const M*       mahalanobisMatrix = NULL,
                                  bool           applyMahalanobisInvert = true) const;

  double logProposal             (const uqChainPositionClass<V>& x,
                                  const uqChainPositionClass<V>& y,
                                  unsigned int                   idOfProposalCovMatrix);
  double logProposal             (const std::vector<uqChainPositionClass<V>*>& inputPositions);
  double alpha                   (const uqChainPositionClass<V>& x,
                                  const uqChainPositionClass<V>& y,
                                  double* alphaQuotientPtr = NULL);
  double alpha                   (const std::vector<uqChainPositionClass<V>*>& inputPositions);
  bool   acceptAlpha             (double alpha);
  void   updateCovMatrices       ();
  //void   gammar                  (double a,
  //                                double b,
  //                                M&     mat);

  const uqEnvironmentClass&                   m_env;
        std::string                           m_prefix;
  const uqParamSpaceClass<V,M>&               m_paramSpace;
  const uqObservableSpaceClass<V,M>&          m_observableSpace;
  const uq_ProbDensity_BaseClass<V,M>&        m_m2lPriorProbDensity_Obj;
  const uq_LikelihoodFunction_BaseClass<V,M>& m_m2lLikelihoodFunction_Obj;

  std::string m_option_help;
  std::string m_option_mh_sizesOfChains;
  std::string m_option_dr_maxNumberOfExtraStages;
  std::string m_option_dr_scalesForExtraStages;
  std::string m_option_am_initialNonAdaptInterval;
  std::string m_option_am_adaptInterval;
  std::string m_option_am_eta;
  std::string m_option_am_epsilon;
  std::string m_option_mh_useChain2;
  std::string m_option_mh_generateUniqueChain;
  std::string m_option_mh_generateExtraChains;
  std::string m_option_mh_namesOfOutputFiles;
  std::string m_option_mh_chainDisplayPeriod;
  std::string m_option_mh_measureRunTimes;
  std::string m_option_mh_generateWhiteNoise;
  std::string m_option_finalPercentsForStats;
  std::string m_option_runBMM;
  std::string m_option_lengthsForBMM;
  std::string m_option_computePSDs;
  std::string m_option_numBlocksForPSD;
  std::string m_option_hopSizeRatioForPSD;
  std::string m_option_computeGewekeCoefs;
  std::string m_option_ratioNaForGeweke;
  std::string m_option_ratioNbForGeweke;
  std::string m_option_computeCorrelations;
  std::string m_option_secondLagForCorrs;
  std::string m_option_lagSpacingForCorrs;
  std::string m_option_numberOfLagsForCorrs;
  std::string m_option_printCorrs;
  std::string m_option_writeCorrs;
  std::string m_option_computeHistograms;
  std::string m_option_numberOfInternalBinsForHists;
  std::string m_option_computeKDEs;
  std::string m_option_numberOfEvaluationPosForKDEs;

  bool                         m_likelihoodObjComputesMisfits;
  V                            m_paramInitials;
  bool                         m_proposalIsSymmetric;
  po::options_description*     m_optionsDesc;
  std::vector<unsigned int>    m_sizesOfChains;
  unsigned int                 m_maxNumberOfExtraStages;
  std::vector<double>          m_scalesForCovMProposals;
  unsigned int                 m_initialNonAdaptInterval;
  unsigned int                 m_adaptInterval;
  double                       m_eta;
  double                       m_epsilon;
  bool                         m_useChain2;
  bool                         m_generateUniqueChain;
  bool                         m_generateExtraChains;
  std::vector<std::string>     m_namesOfOutputFiles;
  unsigned int                 m_chainDisplayPeriod;
  bool                         m_measureRunTimes;
  bool                         m_generateWhiteNoise;

  std::vector<double>          m_finalPercentsForStats;
  bool                         m_runBMM;
  std::vector<unsigned int>    m_lengthsForBMM;
  bool                         m_printBMM;
  bool                         m_writeBMM;

  bool                         m_computePSDs;
  std::vector<unsigned int>    m_numBlocksForPSD;
  double                       m_hopSizeRatioForPSD;
  bool                         m_printPSD;
  bool                         m_writePSD;

  bool                         m_computeGewekeCoefs;
  double                       m_ratioNaForGeweke;
  double                       m_ratioNbForGeweke;

  bool                         m_computeCorrelations;
  unsigned int                 m_secondLagForCorrs;
  unsigned int                 m_lagSpacingForCorrs;
  unsigned int                 m_numberOfLagsForCorrs;
  bool                         m_printCorrs;
  bool                         m_writeCorrs;

  unsigned int                 m_initialPosForUncorrelation;   // set during run time
  unsigned int                 m_spacingForUncorrelation;      // set during run time
  V*                           m_minPositionsForStatistics;    // set during run time
  V*                           m_maxPositionsForStatistics;    // set during run time
  bool                         m_computeHistograms;
  unsigned int                 m_numberOfInternalBinsForHists;
  std::vector<V*>              m_centersForAllHistogramBins;   // set during run time
  std::vector<V*>              m_histogramBinsForAllParams;    // set during run time
  bool                         m_computeKDEs;
  unsigned int                 m_numberOfEvaluationPosForKDEs;
  std::vector<V*>              m_evaluationPositionsForKDEs;   // set during run time
  V*                           m_scalesForKDEs;                // set during run time 
  std::vector<V*>              m_densityValuesFromGaussianKDE; // set during run time

  std::vector<      M*>        m_lowerCholProposalCovMatrices;
  std::vector<      M*>        m_proposalCovMatrices;
#ifdef UQ_DRAM_MCG_REQUIRES_INVERTED_COV_MATRICES
  std::vector<      M*>        m_upperCholProposalPrecMatrices;
  std::vector<      M*>        m_proposalPrecMatrices;
#endif

  std::vector<const V*>        m_chain1;
  std::vector<const V*>        m_uniqueChain1;
  unsigned int                 m_uniqueChain1Pos;
  uqArrayOfSequencesClass<V>   m_chain2;
  uqArrayOfSequencesClass<V>   m_uniqueChain2;
  unsigned int                 m_uniqueChain2Pos;
  std::vector<const V*>        m_misfitChain;         // Sum of squares of differences between model and experiments: computed by user supplied likelihood obj
  std::vector<const V*>        m_misfitVarianceChain;
  std::vector<const V*>        m_m2lLikelihoodChain;
  std::vector<double>          m_alphaQuotients;
  double                       m_chainRunTime;
  double                       m_candidateRunTime;
  double                       m_priorRunTime;
  double                       m_lhRunTime;
  double                       m_mhAlphaRunTime;
  double                       m_drAlphaRunTime;
  double                       m_drRunTime;
  double                       m_amRunTime;
  unsigned int                 m_numRejections;
  unsigned int                 m_numOutOfBounds;
  double                       m_lastChainSize;
  V*                           m_lastMean;
  M*                           m_lastAdaptedCovMatrix;
};

template <class V, class M>
std::ostream& operator<<(std::ostream& os, const uqDRAM_MarkovChainGeneratorClass<V,M>& obj);

#include <uqDRAM_mcg1.h>
#include <uqDRAM_mcg2.h>

template <class V, class M>
uqDRAM_MarkovChainGeneratorClass<V,M>::uqDRAM_MarkovChainGeneratorClass(
  const uqEnvironmentClass&                   env,
  const char*                                 prefix,
  const uqParamSpaceClass<V,M>&               paramSpace,
  const uqObservableSpaceClass<V,M>&          observableSpace,
  const uq_ProbDensity_BaseClass<V,M>&        m2lPriorProbDensity_Obj,
  const uq_LikelihoodFunction_BaseClass<V,M>& m2lLikelihoodFunction_Obj)
  :
  m_env                          (env),
  m_prefix                       (""),
  m_paramSpace                   (paramSpace),
  m_observableSpace              (observableSpace),
  m_m2lPriorProbDensity_Obj      (m2lPriorProbDensity_Obj),
  m_m2lLikelihoodFunction_Obj    (m2lLikelihoodFunction_Obj),
  m_likelihoodObjComputesMisfits (dynamic_cast<const uq_MisfitLikelihoodFunction_Class<V,M>*>(&m2lLikelihoodFunction_Obj) != NULL),
  m_paramInitials                (m_paramSpace.initialValues()),
  m_proposalIsSymmetric          (true),
  m_optionsDesc                  (new po::options_description("Markov chain Monte Carlo options")),
  m_sizesOfChains                (1,(unsigned int) strtod(UQ_MCMC_MH_CHAIN_SIZE_ODV,NULL)),
  m_maxNumberOfExtraStages       (UQ_MCMC_DR_MAX_NUM_EXTRA_STAGES_ODV),
  m_scalesForCovMProposals       (0),//,0.),
  m_initialNonAdaptInterval      (UQ_MCMC_AM_INIT_NON_ADAPT_INT_ODV),
  m_adaptInterval                (UQ_MCMC_AM_ADAPT_INTERVAL_ODV),
  m_eta                          (UQ_MCMC_AM_ETA_ODV),
  m_epsilon                      (UQ_MCMC_AM_EPSILON_ODV),
  m_useChain2                    (UQ_MCMC_MH_USE_CHAIN2),
  m_generateUniqueChain          (UQ_MCMC_MH_GENERATE_UNIQUE_CHAIN_ODV),
  m_generateExtraChains          (UQ_MCMC_MH_GENERATE_EXTRA_CHAINS_ODV),
  m_namesOfOutputFiles           (1,UQ_MCMC_MH_OUTPUT_FILE_NAME_ODV),
  m_chainDisplayPeriod           (UQ_MCMC_MH_CHAIN_DISPLAY_PERIOD_ODV),
  m_measureRunTimes              (UQ_MCMC_MH_MEASURE_RUN_TIMES_ODV),
  m_generateWhiteNoise           (UQ_MCMC_MH_GENERATE_WHITE_NOISE_ODV),
  m_finalPercentsForStats        (0),//,0.),
  m_runBMM                       (UQ_MCMC_RUN_BMM_ODV),
  m_lengthsForBMM                (0),//,0),
  m_computePSDs                  (UQ_MCMC_COMPUTE_PSDS_ODV),
  m_numBlocksForPSD              (0),//,0),
  m_hopSizeRatioForPSD           (1.),
  m_computeGewekeCoefs           (UQ_MCMC_COMPUTE_GEWEKE_COEFS_ODV),
  m_ratioNaForGeweke             (0.),
  m_ratioNbForGeweke             (0.),
  m_computeCorrelations          (UQ_MCMC_COMPUTE_CORRELATIONS_ODV),
  m_secondLagForCorrs            (0),
  m_lagSpacingForCorrs           (0),
  m_numberOfLagsForCorrs         (0),
  m_printCorrs                   (false),
  m_writeCorrs                   (false),
  m_initialPosForUncorrelation   (0),
  m_spacingForUncorrelation      (0),
  m_minPositionsForStatistics    (m_paramSpace.newVector()),
  m_maxPositionsForStatistics    (m_paramSpace.newVector()),
  m_computeHistograms            (UQ_MCMC_COMPUTE_HISTOGRAMS_ODV),
  m_numberOfInternalBinsForHists (0),
  m_centersForAllHistogramBins   (0),//,NULL),
  m_histogramBinsForAllParams    (0),//,NULL),
  m_computeKDEs                  (UQ_MCMC_COMPUTE_KDES_ODV),
  m_numberOfEvaluationPosForKDEs (0),
  m_evaluationPositionsForKDEs   (0),//,NULL),
  m_scalesForKDEs                (m_paramSpace.newVector()),
  m_densityValuesFromGaussianKDE (0),//,NULL),
  m_lowerCholProposalCovMatrices (1),//,NULL),
  m_proposalCovMatrices          (1),//,NULL),
#ifdef UQ_DRAM_MCG_REQUIRES_INVERTED_COV_MATRICES
  m_upperCholProposalPrecMatrices(1),//,NULL),
  m_proposalPrecMatrices         (1),//,NULL),
#endif
  m_chain1                       (0),//,NULL),
  m_uniqueChain1                 (0),//,NULL),
  m_uniqueChain1Pos              (0),
  m_chain2                       (0,m_paramSpace.zeroVector()),
  m_uniqueChain2                 (0,m_paramSpace.zeroVector()),
  m_uniqueChain2Pos              (0),
  m_misfitChain                  (0),//,NULL),
  m_misfitVarianceChain          (0),//,NULL),
  m_m2lLikelihoodChain           (0),//,NULL),
  m_alphaQuotients               (0),//,0.),
  m_chainRunTime                 (0.),
  m_candidateRunTime             (0.),
  m_priorRunTime                 (0.),
  m_lhRunTime                    (0.),
  m_mhAlphaRunTime               (0.),
  m_drAlphaRunTime               (0.),
  m_drRunTime                    (0.),
  m_amRunTime                    (0.), 
  m_numRejections                (0),
  m_numOutOfBounds               (0),
  m_lastChainSize                (0),
  m_lastMean                     (NULL),
  m_lastAdaptedCovMatrix         (NULL)
{
  if (m_env.rank() == 0) std::cout << "Entering uqDRAM_MarkovChainGeneratorClass<V,M>::constructor()"
                                   << std::endl;

  if ((prefix         != NULL) && 
      (strlen(prefix) != 0   )) {
    std::string tmpString(prefix);
    m_prefix = tmpString + "_";
  }

  m_option_help                         = m_prefix + "MCMC_help";
  m_option_mh_sizesOfChains             = m_prefix + "MCMC_mh_sizesOfChains";
  m_option_dr_maxNumberOfExtraStages    = m_prefix + "MCMC_dr_maxNumberOfExtraStages";
  m_option_dr_scalesForExtraStages      = m_prefix + "MCMC_dr_scalesForExtraStages";
  m_option_am_initialNonAdaptInterval   = m_prefix + "MCMC_am_initialNonAdaptInterval";
  m_option_am_adaptInterval             = m_prefix + "MCMC_am_adaptInterval";
  m_option_am_eta                       = m_prefix + "MCMC_am_eta";
  m_option_am_epsilon                   = m_prefix + "MCMC_am_epsilon";
  m_option_mh_useChain2                 = m_prefix + "MCMC_mh_useChain2";
  m_option_mh_generateUniqueChain       = m_prefix + "MCMC_mh_generateUniqueChain";
  m_option_mh_generateExtraChains       = m_prefix + "MCMC_mh_generateExtraChains";
  m_option_mh_namesOfOutputFiles        = m_prefix + "MCMC_mh_namesOfOutputFiles";
  m_option_mh_chainDisplayPeriod        = m_prefix + "MCMC_mh_chainDisplayPeriod";
  m_option_mh_measureRunTimes           = m_prefix + "MCMC_mh_measureRunTimes";
  m_option_mh_generateWhiteNoise        = m_prefix + "MCMC_mh_generateWhiteNoise";
  m_option_finalPercentsForStats        = m_prefix + "MCMC_finalPercentsForStats";
  m_option_runBMM                       = m_prefix + "MCMC_runBMM";
  m_option_lengthsForBMM                = m_prefix + "MCMC_lengthsForBMM";
  m_option_computePSDs                  = m_prefix + "MCMC_computePSDs";
  m_option_numBlocksForPSD              = m_prefix + "MCMC_numBlocksForPSD";
  m_option_hopSizeRatioForPSD           = m_prefix + "MCMC_hopSizeRatioForPSD";
  m_option_computeGewekeCoefs           = m_prefix + "MCMC_computeGewekeCoefficients";
  m_option_ratioNaForGeweke             = m_prefix + "MCMC_ratioNaForGeweke";
  m_option_ratioNbForGeweke             = m_prefix + "MCMC_ratioNbForGeweke";
  m_option_computeCorrelations          = m_prefix + "MCMC_computeCorrelations";
  m_option_secondLagForCorrs            = m_prefix + "MCMC_secondLagForCorrs";
  m_option_lagSpacingForCorrs           = m_prefix + "MCMC_lagSpacingForCorrs";
  m_option_numberOfLagsForCorrs         = m_prefix + "MCMC_numberOfLagsForCorrs";
  m_option_printCorrs                   = m_prefix + "MCMC_printCorrs";
  m_option_writeCorrs                   = m_prefix + "MCMC_writeCorrs";
  m_option_computeHistograms            = m_prefix + "MCMC_computeHistograms";
  m_option_numberOfInternalBinsForHists = m_prefix + "MCMC_numberOfInternalBinsForHistograms";
  m_option_computeKDEs                  = m_prefix + "MCMC_computeKDEs";
  m_option_numberOfEvaluationPosForKDEs = m_prefix + "MCMC_numberOfEvaluationPositionsForKDEs";

  m_numBlocksForPSD.resize(1);
  m_numBlocksForPSD[0] = 8;

  m_hopSizeRatioForPSD = .50;

  m_ratioNaForGeweke = .1;
  m_ratioNbForGeweke = .5;

  m_initialPosForUncorrelation   = 20000;
  m_spacingForUncorrelation      = 50;

  m_numberOfInternalBinsForHists = 100;

  m_numberOfEvaluationPosForKDEs = 100;

  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.rank() == 0) std::cout << "After getting values of options with prefix '" << m_prefix
                                   << "', state of uqDRAM_MarkovChainGeneratorClass object is:"
                                   << "\n" << *this
                                   << std::endl;

  if (m_env.rank() == 0) std::cout << "Leaving uqDRAM_MarkovChainGeneratorClass<V,M>::constructor()"
                                   << std::endl;
}

template <class V, class M>
uqDRAM_MarkovChainGeneratorClass<V,M>::~uqDRAM_MarkovChainGeneratorClass()
{
  //std::cout << "Entering uqDRAM_MarkovChainGeneratorClass<V,M>::destructor()"
  //          << std::endl;

  resetChainAndRelatedInfo();

  for (unsigned int i = 0; i < m_densityValuesFromGaussianKDE.size(); ++i) {
    if (m_densityValuesFromGaussianKDE[i] != NULL) delete m_densityValuesFromGaussianKDE[i];
  }
  if (m_scalesForKDEs != NULL) delete m_scalesForKDEs;
  for (unsigned int i = 0; i < m_evaluationPositionsForKDEs.size(); ++i) {
    if (m_evaluationPositionsForKDEs[i] != NULL) delete m_evaluationPositionsForKDEs[i];
  }
  for (unsigned int i = 0; i < m_histogramBinsForAllParams.size(); ++i) {
    if (m_histogramBinsForAllParams[i] != NULL) delete m_histogramBinsForAllParams[i];
  }
  for (unsigned int i = 0; i < m_centersForAllHistogramBins.size(); ++i) {
    if (m_centersForAllHistogramBins[i] != NULL) delete m_centersForAllHistogramBins[i];
  }
  if (m_maxPositionsForStatistics != NULL) delete m_maxPositionsForStatistics;
  if (m_minPositionsForStatistics != NULL) delete m_minPositionsForStatistics;

  m_numBlocksForPSD.clear();

  if (m_optionsDesc            ) delete m_optionsDesc;

  //std::cout << "Leaving uqDRAM_MarkovChainGeneratorClass<V,M>::destructor()"
  //          << std::endl;
}

template <class V, class M>
void
uqDRAM_MarkovChainGeneratorClass<V,M>::resetChainAndRelatedInfo()
{
  if (m_lastAdaptedCovMatrix) delete m_lastAdaptedCovMatrix;
  if (m_lastMean)             delete m_lastMean;
  m_lastChainSize     = 0;
  m_numOutOfBounds    = 0;
  m_chainRunTime      = 0.;
  m_candidateRunTime  = 0.;
  m_priorRunTime      = 0.;
  m_lhRunTime         = 0.;
  m_mhAlphaRunTime    = 0.;
  m_drAlphaRunTime    = 0.;
  m_drRunTime         = 0.;
  m_amRunTime         = 0.;
  m_numRejections     = 0;
  m_alphaQuotients.clear();
  for (unsigned int i = 0; i < m_m2lLikelihoodChain.size(); ++i) {
    if (m_m2lLikelihoodChain[i]) delete m_m2lLikelihoodChain[i];
  }
  for (unsigned int i = 0; i < m_misfitVarianceChain.size(); ++i) {
    if (m_misfitVarianceChain[i]) delete m_misfitVarianceChain[i];
  }
  for (unsigned int i = 0; i < m_misfitChain.size(); ++i) {
    if (m_misfitChain[i]) delete m_misfitChain[i];
  }

  m_uniqueChain2Pos = 0;
  //m_uniqueChain2.resetValues();

  m_uniqueChain1Pos = 0;
  for (unsigned int i = 0; i < m_uniqueChain1.size(); ++i) {
    if (m_uniqueChain1[i]) delete m_uniqueChain1[i];
  }
  for (unsigned int i = 0; i < m_chain1.size(); ++i) {
    if (m_chain1[i]) delete m_chain1[i];
  }

#ifdef UQ_DRAM_MCG_REQUIRES_INVERTED_COV_MATRICES
  for (unsigned int i = 0; i < m_proposalPrecMatrices.size(); ++i) {
    if (m_proposalPrecMatrices[i]) delete m_proposalPrecMatrices[i];
  }
  for (unsigned int i = 0; i < m_upperCholProposalPrecMatrices.size(); ++i) {
    if (m_upperCholProposalPrecMatrices[i]) delete m_upperCholProposalPrecMatrices[i];
  }
#endif
  for (unsigned int i = 0; i < m_proposalCovMatrices.size(); ++i) {
    if (m_proposalCovMatrices[i]) delete m_proposalCovMatrices[i];
  }
  m_proposalCovMatrices.clear();
  for (unsigned int i = 0; i < m_lowerCholProposalCovMatrices.size(); ++i) {
    if (m_lowerCholProposalCovMatrices[i]) delete m_lowerCholProposalCovMatrices[i];
  }
  m_lowerCholProposalCovMatrices.clear();

  return;
}

template <class V, class M>
void
uqDRAM_MarkovChainGeneratorClass<V,M>::defineMyOptions(
  po::options_description& optionsDesc) const
{
  optionsDesc.add_options()
    (m_option_help.c_str(),                                                                                                           "produce help message for DRAM Markov chain generator"   )
    (m_option_mh_sizesOfChains.c_str(),             po::value<std::string >()->default_value(UQ_MCMC_MH_CHAIN_SIZE_ODV             ), "'mh' size(s) of chain(s)"                               )
    (m_option_dr_maxNumberOfExtraStages.c_str(),    po::value<unsigned int>()->default_value(UQ_MCMC_DR_MAX_NUM_EXTRA_STAGES_ODV   ), "'dr' maximum number of extra stages"                    )
    (m_option_dr_scalesForExtraStages.c_str(),      po::value<std::string >()->default_value(UQ_MCMC_DR_SCALES_FOR_EXTRA_STAGES_ODV), "'dr' scales for proposal cov matrices from 2nd stage on")
    (m_option_am_initialNonAdaptInterval.c_str(),   po::value<unsigned int>()->default_value(UQ_MCMC_AM_INIT_NON_ADAPT_INT_ODV     ), "'am' initial non adaptation interval"                   )
    (m_option_am_adaptInterval.c_str(),             po::value<unsigned int>()->default_value(UQ_MCMC_AM_ADAPT_INTERVAL_ODV         ), "'am' adaptation interval"                               )
    (m_option_am_eta.c_str(),                       po::value<double      >()->default_value(UQ_MCMC_AM_ETA_ODV                    ), "'am' eta"                                               )
    (m_option_am_epsilon.c_str(),                   po::value<double      >()->default_value(UQ_MCMC_AM_EPSILON_ODV                ), "'am' epsilon"                                           )
    (m_option_mh_useChain2.c_str(),                 po::value<bool        >()->default_value(UQ_MCMC_MH_USE_CHAIN2                 ), "'mh' use chain2"                                        )
    (m_option_mh_generateUniqueChain.c_str(),       po::value<bool        >()->default_value(UQ_MCMC_MH_GENERATE_UNIQUE_CHAIN_ODV  ), "'mh' generate unique chain"                             )
    (m_option_mh_generateExtraChains.c_str(),       po::value<bool        >()->default_value(UQ_MCMC_MH_GENERATE_EXTRA_CHAINS_ODV  ), "'mh' generate extra chains"                             )
    (m_option_mh_namesOfOutputFiles.c_str(),        po::value<std::string >()->default_value(UQ_MCMC_MH_OUTPUT_FILE_NAME_ODV       ), "'mh' name(s) of output file(s)"                         )
    (m_option_mh_chainDisplayPeriod.c_str(),        po::value<unsigned int>()->default_value(UQ_MCMC_MH_CHAIN_DISPLAY_PERIOD_ODV   ), "'mh' period of message display during chain generation" )
    (m_option_mh_measureRunTimes.c_str(),           po::value<bool        >()->default_value(UQ_MCMC_MH_MEASURE_RUN_TIMES_ODV      ), "'mh' measure run times"                                 )
    (m_option_mh_generateWhiteNoise.c_str(),        po::value<bool        >()->default_value(UQ_MCMC_MH_GENERATE_WHITE_NOISE_ODV   ), "'mh' generate white noise"                              )
    (m_option_finalPercentsForStats.c_str(),        po::value<std::string >()->default_value(UQ_MCMC_FINAL_PERCENTS_FOR_STATS_ODV  ), "final percentages for computation of chain statistics"  )
    (m_option_runBMM.c_str(),                       po::value<bool        >()->default_value(UQ_MCMC_RUN_BMM_ODV                   ), "compute variance of sample mean with batch means method")
    (m_option_lengthsForBMM.c_str(),                po::value<std::string >()->default_value(UQ_MCMC_LENGTHS_FOR_BMM_ODV           ), "batch lenghts for BMM")
    (m_option_computePSDs.c_str(),                  po::value<bool        >()->default_value(UQ_MCMC_COMPUTE_PSDS_ODV              ), "compute power spectral densities"                       )
#if 0
    (m_option_numBlocksForPSD.c_str(),              po::value<            >()->default_value(), "")
    (m_option_hopSizeRatioForPSD.c_str(),           po::value<            >()->default_value(), "")
#endif
    (m_option_computeGewekeCoefs.c_str(),           po::value<bool        >()->default_value(UQ_MCMC_COMPUTE_GEWEKE_COEFS_ODV      ), "compute Geweke coefficients"                            )
#if 0
    (m_option_ratioNaForGeweke.c_str(),             po::value<            >()->default_value(), "")
    (m_option_ratioNbForGeweke.c_str(),             po::value<            >()->default_value(), "")
#endif
    (m_option_computeCorrelations.c_str(),          po::value<bool        >()->default_value(UQ_MCMC_COMPUTE_CORRELATIONS_ODV      ), "compute correlations"                                   )
    (m_option_secondLagForCorrs.c_str(),            po::value<unsigned int>()->default_value(UQ_MCMC_SECOND_LAG_FOR_CORRS_ODV      ), "second lag for computation of autocorrelations"         )
    (m_option_lagSpacingForCorrs.c_str(),           po::value<unsigned int>()->default_value(UQ_MCMC_LAG_SPACING_FOR_CORRS_ODV     ), "lag spacing for computation of autocorrelations"        )
    (m_option_numberOfLagsForCorrs.c_str(),         po::value<unsigned int>()->default_value(UQ_MCMC_NUMBER_OF_LAGS_FOR_CORRS_ODV  ), "number of lags for computation of autocorrelations"     )
    (m_option_printCorrs.c_str(),                   po::value<bool        >()->default_value(UQ_MCMC_PRINT_CORRS_ODV               ), "print computed autocorrelations on the screen"          )
    (m_option_writeCorrs.c_str(),                   po::value<bool        >()->default_value(UQ_MCMC_WRITE_CORRS_ODV               ), "write computed autocorrelations to the output file"     )
    (m_option_computeHistograms.c_str(),            po::value<bool        >()->default_value(UQ_MCMC_COMPUTE_HISTOGRAMS_ODV        ), "compute histograms"                                     )
#if 0
    (m_option_numberOfInternalBinsForHists.c_str(), po::value<            >()->default_value(), "")
#endif
    (m_option_computeKDEs.c_str(),                  po::value<bool        >()->default_value(UQ_MCMC_COMPUTE_KDES_ODV              ), "compute kernel density estimators"                      )
#if 0
    (m_option_numberOfEvaluationPosForKDEs.c_str(), po::value<            >()->default_value(), "")
#endif
  ;

  return;
}

template <class V, class M>
void
uqDRAM_MarkovChainGeneratorClass<V,M>::getMyOptionValues(
  po::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count(m_option_help.c_str())) {
    std::cout << optionsDesc
              << std::endl;
  }

  if (m_env.allOptionsMap().count(m_option_mh_sizesOfChains.c_str())) {
    std::string inputString = m_env.allOptionsMap()[m_option_mh_sizesOfChains.c_str()].as<std::string>();
    std::vector<double> inputDoubles(0,0.);
    uqMiscReadDoublesFromString(inputString,inputDoubles);
    m_sizesOfChains.clear();

    m_sizesOfChains.resize(inputDoubles.size(),0);
    for (unsigned int i = 0; i < inputDoubles.size(); ++i) {
      m_sizesOfChains[i] = (unsigned int) inputDoubles[i];
    }
  }

  if (m_env.allOptionsMap().count(m_option_dr_maxNumberOfExtraStages.c_str())) {
    m_maxNumberOfExtraStages = m_env.allOptionsMap()[m_option_dr_maxNumberOfExtraStages.c_str()].as<unsigned int>();
  }

  std::vector<double> tmpScales(0,0.);
  if (m_env.allOptionsMap().count(m_option_dr_scalesForExtraStages.c_str())) {
    std::string inputString = m_env.allOptionsMap()[m_option_dr_scalesForExtraStages.c_str()].as<std::string>();
    uqMiscReadDoublesFromString(inputString,tmpScales);
    //std::cout << "In uqDRAM_MarkovChainGeneratorClass<V,M>::getMyOptionValues(): scales =";
    //for (unsigned int i = 0; i < tmpScales.size(); ++i) {
    //  std::cout << " " << tmpScales[i];
    //}
    //std::cout << std::endl;
  }

  if (m_maxNumberOfExtraStages > 0) {
    double scale = 1.0;
    unsigned int tmpSize = tmpScales.size();

    m_scalesForCovMProposals.resize      (m_maxNumberOfExtraStages+1,1.);
    m_lowerCholProposalCovMatrices.resize(m_maxNumberOfExtraStages+1,NULL);
    m_proposalCovMatrices.resize         (m_maxNumberOfExtraStages+1,NULL);

    for (unsigned int i = 1; i < (m_maxNumberOfExtraStages+1); ++i) {
      if (i <= tmpSize) scale = tmpScales[i-1];
      m_scalesForCovMProposals[i] = scale;
    }
    //updateCovMatrices();
  }

  if (m_env.allOptionsMap().count(m_option_am_initialNonAdaptInterval.c_str())) {
    m_initialNonAdaptInterval = m_env.allOptionsMap()[m_option_am_initialNonAdaptInterval.c_str()].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_am_adaptInterval.c_str())) {
    m_adaptInterval = m_env.allOptionsMap()[m_option_am_adaptInterval.c_str()].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_am_eta.c_str())) {
    m_eta = m_env.allOptionsMap()[m_option_am_eta.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_am_epsilon.c_str())) {
    m_epsilon = m_env.allOptionsMap()[m_option_am_epsilon.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_mh_useChain2.c_str())) {
    m_useChain2 = m_env.allOptionsMap()[m_option_mh_useChain2.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_mh_generateUniqueChain.c_str())) {
    m_generateUniqueChain = m_env.allOptionsMap()[m_option_mh_generateUniqueChain.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_mh_generateExtraChains.c_str())) {
    m_generateExtraChains = m_env.allOptionsMap()[m_option_mh_generateExtraChains.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_mh_namesOfOutputFiles.c_str())) {
    std::string inputString = m_env.allOptionsMap()[m_option_mh_namesOfOutputFiles.c_str()].as<std::string>();
    m_namesOfOutputFiles.clear();
    uqMiscReadWordsFromString(inputString,m_namesOfOutputFiles);

    UQ_FATAL_TEST_MACRO(m_namesOfOutputFiles.size() != m_sizesOfChains.size(),
                        m_env.rank(),
                        "uqDRAM_MarkovChainGeneratorClass<V,M>::readMyOptionsValues()",
                        "size of array for 'outputFileNames' is not equal to size of array for 'chainSizes'");
  }

  if (m_env.allOptionsMap().count(m_option_mh_chainDisplayPeriod.c_str())) {
    m_chainDisplayPeriod = m_env.allOptionsMap()[m_option_mh_chainDisplayPeriod.c_str()].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_mh_measureRunTimes.c_str())) {
    m_measureRunTimes = m_env.allOptionsMap()[m_option_mh_measureRunTimes.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_mh_generateWhiteNoise.c_str())) {
    m_generateWhiteNoise = m_env.allOptionsMap()[m_option_mh_generateWhiteNoise.c_str()].as<bool>();
  }

  std::vector<double> tmpPercents(0,0.);
  if (m_env.allOptionsMap().count(m_option_finalPercentsForStats.c_str())) {
    std::string inputString = m_env.allOptionsMap()[m_option_finalPercentsForStats.c_str()].as<std::string>();
    uqMiscReadDoublesFromString(inputString,tmpPercents);
    //std::cout << "In uqDRAM_MarkovChainGeneratorClass<V,M>::getMyOptionValues(): percents = ";
    //for (unsigned int i = 0; i < tmpPercents.size(); ++i) {
    //  std::cout << " " << tmpPercents[i];
    //}
    //std::cout << std::endl;

    if (tmpPercents.size() > 0) {
      m_finalPercentsForStats.resize(tmpPercents.size(),0.);
      for (unsigned int i = 0; i < m_finalPercentsForStats.size(); ++i) {
        m_finalPercentsForStats[i] = tmpPercents[i];
      }
    }
  }

  if (m_env.allOptionsMap().count(m_option_runBMM.c_str())) {
    m_runBMM = m_env.allOptionsMap()[m_option_runBMM.c_str()].as<bool>();
  }
  std::cout << "Aqui" << std::endl;
  std::vector<double> tmpLengths(0,0.);
  if (m_env.allOptionsMap().count(m_option_lengthsForBMM.c_str())) {
    std::string inputString = m_env.allOptionsMap()[m_option_lengthsForBMM.c_str()].as<std::string>();
    uqMiscReadDoublesFromString(inputString,tmpLengths);
    //std::cout << "In uqDRAM_MarkovChainGeneratorClass<V,M>::getMyOptionValues(): lengths for BMM = ";
    //for (unsigned int i = 0; i < tmpLengths.size(); ++i) {
    //  std::cout << " " << tmpLengths[i];
    //}
    //std::cout << std::endl;

    if (tmpLengths.size() > 0) {
      m_lengthsForBMM.resize(tmpLengths.size(),0);
      for (unsigned int i = 0; i < m_lengthsForBMM.size(); ++i) {
        m_lengthsForBMM[i] = (unsigned int) tmpLengths[i];
      }
    }
  }

  if (m_env.allOptionsMap().count(m_option_computePSDs.c_str())) {
    m_computePSDs = m_env.allOptionsMap()[m_option_computePSDs.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_computeGewekeCoefs.c_str())) {
    m_computeGewekeCoefs = m_env.allOptionsMap()[m_option_computeGewekeCoefs.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_computeCorrelations.c_str())) {
    m_computeCorrelations = m_env.allOptionsMap()[m_option_computeCorrelations.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_secondLagForCorrs.c_str())) {
    m_secondLagForCorrs = m_env.allOptionsMap()[m_option_secondLagForCorrs.c_str()].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_lagSpacingForCorrs.c_str())) {
    m_lagSpacingForCorrs = m_env.allOptionsMap()[m_option_lagSpacingForCorrs.c_str()].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_numberOfLagsForCorrs.c_str())) {
    m_numberOfLagsForCorrs = m_env.allOptionsMap()[m_option_numberOfLagsForCorrs.c_str()].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_printCorrs.c_str())) {
    m_printCorrs = m_env.allOptionsMap()[m_option_printCorrs.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_writeCorrs.c_str())) {
    m_writeCorrs = m_env.allOptionsMap()[m_option_writeCorrs.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_computeHistograms.c_str())) {
    m_computeHistograms = m_env.allOptionsMap()[m_option_computeHistograms.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_computeKDEs.c_str())) {
    m_computeKDEs = m_env.allOptionsMap()[m_option_computeKDEs.c_str()].as<bool>();
  }

  return;
}

template <class V, class M>
void
uqDRAM_MarkovChainGeneratorClass<V,M>::generateChains(
  const M* proposalCovMatrix,
  const M* mahalanobisMatrix,
  bool     applyMahalanobisInvert)
{
  if (m_useChain2) {
    generateChains2(proposalCovMatrix,
                    mahalanobisMatrix,
                    applyMahalanobisInvert);
  }
  else {
    generateChains1(proposalCovMatrix,
                    mahalanobisMatrix,
                    applyMahalanobisInvert);
  }

  return;
}

template <class V, class M>
int
uqDRAM_MarkovChainGeneratorClass<V,M>::prepareForNextChain(
  const M* proposalCovMatrix)
  //const M* proposalPrecMatrix)
{
  //if (m_env.rank() == 0) std::cout << "Entering uqDRAM_MarkovChainGeneratorClass<V,M>::prepareForNextChain()..."
  //                                 << std::endl;

  int iRC = UQ_OK_RC;

  const M* internalProposalCovMatrix = proposalCovMatrix;
  if (proposalCovMatrix == NULL) {
    V tmpVec(m_paramSpace.zeroVector());
    for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
      double sigma = m_paramSpace.parameter(i).priorSigma();
      if ((sigma == INFINITY) ||
          (sigma == NAN     )) {
        tmpVec[i] = pow( fabs(m_paramInitials[i])*0.05,2. );
        if ( tmpVec[i] == 0 ) tmpVec[i] = 1.;
      }
      else if (sigma == 0.) {
        tmpVec[i] = 1.;
      }
      else {
        tmpVec[i] = sigma*sigma;
      }
    }
    internalProposalCovMatrix = m_paramSpace.uqFinDimLinearSpaceClass<V,M>::newDiagMatrix(tmpVec);

    if (m_env.rank() == 0) std::cout << "In uqDRAM_MarkovChainGeneratorClass<V,M>::prepareForNextChain()"
                                     << ", contents of internally generated proposal cov matrix are:"
                                     << std::endl;
    std::cout << *internalProposalCovMatrix;
    if (m_env.rank() == 0) std::cout << std::endl;
  }
  else {
    if (m_env.rank() == 0)  std::cout << "In uqDRAM_MarkovChainGeneratorClass<V,M>::prepareForNextChain()"
                                      << "using suplied proposalCovMatrix, whose contents are:"
                                      << std::endl;
    std::cout << *internalProposalCovMatrix;
    if (m_env.rank() == 0) std::cout << std::endl;
  }

  m_lowerCholProposalCovMatrices[0] = new M(*internalProposalCovMatrix); 
  iRC = m_lowerCholProposalCovMatrices[0]->chol();
  UQ_FATAL_RC_MACRO(iRC,
                    m_env.rank(),
                    "uqDRAM_MarkovChainGeneratorClass<V,M>::prepareForNextChain()",
                    "proposalCovMatrix is not positive definite");
  m_lowerCholProposalCovMatrices[0]->zeroUpper(false);

  m_proposalCovMatrices[0] = new M(*internalProposalCovMatrix);

  if (m_env.rank() == 0) std::cout << "In uqDRAM_MarkovChainGeneratorClass<V,M>::prepareForNextChain()"
                                   << ", m_lowerCholProposalCovMatrices[0] contents are:"
                                   << std::endl;
  std::cout << *(m_lowerCholProposalCovMatrices[0]);
  if (m_env.rank() == 0) std::cout << std::endl;

#ifdef UQ_DRAM_MCG_REQUIRES_INVERTED_COV_MATRICES
  const M* internalProposalPrecMatrix = proposalPrecMatrix;
  if (proposalPrecMatrix == NULL) {
    UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                      m_env.rank(),
                      "uqDRAM_MarkovChainGeneratorClass<V,M>::prepareForNextChain()",
                      "not yet implemented for the case 'proposalPrecMatrix == NULL'");
  }

  m_upperCholProposalPrecMatrices[0] = new M(*internalProposalPrecMatrix); 
  iRC = m_upperCholProposalPrecMatrices[0]->chol();
  UQ_FATAL_RC_MACRO(iRC,
                    m_env.rank(),
                    "uqDRAM_MarkovChainGeneratorClass<V,M>::prepareForNextChain()",
                    "proposalPrecMatrix is not positive definite");
  m_upperCholProposalPrecMatrices[0]->zeroLower(false);

  m_proposalPrecMatrices[0] = new M(*internalProposalPrecMatrix);

  //if (m_env.rank() == 0) std::cout << "In uqDRAM_MarkovChainGeneratorClass<V,M>::prepareForNextChain()"
  //                                 << ", m_upperCholProposalPrecMatrices[0] contents are:"
  //                                 << std::endl;
  //std::cout << *(m_upperCholProposalPrecMatrices[0]);
  //if (m_env.rank() == 0) std::cout << std::endl;
#endif

  if (m_maxNumberOfExtraStages > 0) {
    updateCovMatrices();
  }

  if (internalProposalCovMatrix != NULL) delete internalProposalCovMatrix;

  //if (m_env.rank() == 0) std::cout << "Leaving uqDRAM_MarkovChainGeneratorClass<V,M>::prepareForNextChain()"
  //                                 << std::endl;

  return iRC;
}

template <class V, class M>
void
uqDRAM_MarkovChainGeneratorClass<V,M>::updateCovMatrices()
{
  for (unsigned int i = 1; i < (m_maxNumberOfExtraStages+1); ++i) {
    double scale = m_scalesForCovMProposals[i];
    if (m_lowerCholProposalCovMatrices[i]) delete m_lowerCholProposalCovMatrices[i];
    m_lowerCholProposalCovMatrices [i]   = new M(*(m_lowerCholProposalCovMatrices[i-1]));
  *(m_lowerCholProposalCovMatrices [i]) /= scale;
    if (m_proposalCovMatrices[i]) delete m_proposalCovMatrices[i];
    m_proposalCovMatrices[i]             = new M(*(m_proposalCovMatrices[i-1]));
  *(m_proposalCovMatrices[i])           /= (scale*scale);
#ifdef UQ_DRAM_MCG_REQUIRES_INVERTED_COV_MATRICES
    m_upperCholProposalPrecMatrices[i]   = new M(*(m_upperCholProposalPrecMatrices[i-1]));
  *(m_upperCholProposalPrecMatrices[i]) *= scale;
    m_proposalPrecMatrices[i]            = new M(*(m_proposalPrecMatrices[i-1]));
  *(m_proposalPrecMatrices[i])          *= (scale*scale);
#endif
  }

  return;
}

template <class V, class M>
double
uqDRAM_MarkovChainGeneratorClass<V,M>::logProposal(
  const uqChainPositionClass<V>& x,
  const uqChainPositionClass<V>& y,
  unsigned int                   idOfProposalCovMatrix)
{
  V diffVec(y.paramValues() - x.paramValues());
#ifdef UQ_DRAM_MCG_REQUIRES_INVERTED_COV_MATRICES
  double value = -0.5 * scalarProduct(diffVec, *(m_proposalPrecMatrices[idOfProposalCovMatrix]) * diffVec);
#else
  double value = -0.5 * scalarProduct(diffVec, m_proposalCovMatrices[idOfProposalCovMatrix]->invertMultiply(diffVec));
#endif
  return value;
}

template <class V, class M>
double
uqDRAM_MarkovChainGeneratorClass<V,M>::logProposal(const std::vector<uqChainPositionClass<V>*>& inputPositions)
{
  unsigned int inputSize = inputPositions.size();
  UQ_FATAL_TEST_MACRO((inputSize < 2),
                      m_env.rank(),
                      "uqDRAM_MarkovChainGeneratorClass<V,M>::logProposal()",
                      "inputPositions has size < 2");

  return this->logProposal(*(inputPositions[0            ]),
                           *(inputPositions[inputSize - 1]),
                           inputSize-2);
}

template <class V, class M>
double
uqDRAM_MarkovChainGeneratorClass<V,M>::alpha(
  const uqChainPositionClass<V>& x,
  const uqChainPositionClass<V>& y,
  double*                        alphaQuotientPtr)
{
  double alphaQuotient = 0.;
  bool xOutOfBounds = x.outOfBounds();
  bool yOutOfBounds = y.outOfBounds();
  if ((xOutOfBounds == false) &&
      (yOutOfBounds == false)) {
    double yLogPosteriorToUse = y.logPosterior();
    if (m_likelihoodObjComputesMisfits &&
        m_observableSpace.shouldVariancesBeUpdated()) {
      // Divide the misfitVector of 'y' by the misfitVarianceVector of 'x'
      yLogPosteriorToUse = -0.5 * ( y.m2lPrior() + (y.misfitVector()/x.misfitVarianceVector()).sumOfComponents() );
    }
    if (m_proposalIsSymmetric) {
      alphaQuotient = exp(yLogPosteriorToUse - x.logPosterior());
      if ((m_env.verbosity() >= 10) && (m_env.rank() == 0)) {
        std::cout << "In uqDRAM_MarkovChainGeneratorClass<V,M>::alpha()"
                  << ": symmetric proposal case"
                  << ", yLogPosteriorToUse = " << yLogPosteriorToUse
                  << ", x.logPosterior() = "   << x.logPosterior()
                  << ", alpha = "              << alphaQuotient
                  << std::endl;
      }
    }
    else {
      alphaQuotient = exp(yLogPosteriorToUse + logProposal(y,x,0) - x.logPosterior() - logProposal(x,y,0));
    }
  }
  else {
      if ((m_env.verbosity() >= 10) && (m_env.rank() == 0)) {
        std::cout << "In uqDRAM_MarkovChainGeneratorClass<V,M>::alpha()"
                  << ": xOutOfBounds = " << xOutOfBounds
                  << ", yOutOfBounds = " << yOutOfBounds
                  << std::endl;
      }
  }
  if (alphaQuotientPtr != NULL) *alphaQuotientPtr = alphaQuotient;

  return std::min(1.,alphaQuotient);
}

template <class V, class M>
double
uqDRAM_MarkovChainGeneratorClass<V,M>::alpha(const std::vector<uqChainPositionClass<V>*>& inputPositions)
{
  unsigned int inputSize = inputPositions.size();
  UQ_FATAL_TEST_MACRO((inputSize < 2),
                      m_env.rank(),
                      "uqDRAM_MarkovChainGeneratorClass<V,M>::alpha()",
                      "inputPositions has size < 2");

  // If necessary, return 0. right away
  if (inputPositions[0          ]->outOfBounds()) return 0.;
  if (inputPositions[inputSize-1]->outOfBounds()) return 0.;

  // If inputSize is 2, recursion is not needed
  if (inputSize == 2) return this->alpha(*(inputPositions[0            ]),
                                         *(inputPositions[inputSize - 1]));

  // Prepare two vectors of positions
  std::vector<uqChainPositionClass<V>*>         positions(inputSize,NULL);
  std::vector<uqChainPositionClass<V>*> backwardPositions(inputSize,NULL);
  for (unsigned int i = 0; i < inputSize; ++i) {
            positions[i] = inputPositions[i];
    backwardPositions[i] = inputPositions[inputSize-i-1];
  }

  // Initialize cumulative variables
  double logNumerator      = 0.;
  double logDenominator    = 0.;
  double alphasNumerator   = 1.;
  double alphasDenominator = 1.;

  // Compute cumulative variables
  logNumerator   += logProposal(backwardPositions);
  logDenominator += logProposal(        positions);

  for (unsigned int i = 0; i < (inputSize-2); ++i) { // That is why size must be >= 2
            positions.pop_back();
    backwardPositions.pop_back();

    logNumerator   += logProposal(backwardPositions);
    logDenominator += logProposal(        positions);

    alphasNumerator   *= (1 - this->alpha(backwardPositions));
    alphasDenominator *= (1 - this->alpha(        positions));
  }

  double numeratorLogPosteriorToUse = backwardPositions[0]->logPosterior();
  if (m_likelihoodObjComputesMisfits &&
      m_observableSpace.shouldVariancesBeUpdated()) {
    // Divide the misfitVector of 'back[0]' by the misfitVarianceVector of 'pos[0]'
    numeratorLogPosteriorToUse = -0.5 * ( backwardPositions[0]->m2lPrior() +
      (backwardPositions[0]->misfitVector()/positions[0]->misfitVarianceVector()).sumOfComponents() );
  }
  logNumerator   += numeratorLogPosteriorToUse;
  logDenominator += positions[0]->logPosterior();

  // Return result
  return std::min(1.,(alphasNumerator/alphasDenominator)*exp(logNumerator-logDenominator));
}

template <class V, class M>
bool
uqDRAM_MarkovChainGeneratorClass<V,M>::acceptAlpha(double alpha)
{
  bool result = false;

  if      (alpha <= 0.                          ) result = false;
  else if (alpha >= 1.                          ) result = true;
  else if (alpha >= gsl_rng_uniform(m_env.rng())) result = true;
  else                                            result = false;

  return result;
}

//template <class V, class M>
//void
//uqDRAM_MarkovChainGeneratorClass<V,M>::gammar(
//  double a,
//  double b,
//  M&     mat)
//{
//  for (unsigned int i = 0; i < mat.numRows(); ++i) {
//    for (unsigned int j = 0; j < mat.numCols(); ++j) {
//      mat(i,j) = uqMiscGammar(a,b,m_env.rng());
//    }
//  }
//
//  return mat;
//}

#if 0
template <class V, class M>
const std::vector<const V*>&
uqDRAM_MarkovChainGeneratorClass<V,M>::chain() const
{
  return m_chain1;
}

template <class V, class M>
const std::vector<const V*>&
uqDRAM_MarkovChainGeneratorClass<V,M>::misfitChain() const
{
  return m_misfitChain;
}

template <class V, class M>
const std::vector<const V*>&
uqDRAM_MarkovChainGeneratorClass<V,M>::misfitVarianceChain() const
{
  return m_misfitVarianceChain;
}

template <class V, class M>
const std::string&
uqDRAM_MarkovChainGeneratorClass<V,M>::outputFileName() const
{
  return m_namesOfOutputFiles[m_chain1.size()-1];
}
#endif

template <class V, class M>
void
uqDRAM_MarkovChainGeneratorClass<V,M>::print(std::ostream& os) const
{
  os << m_option_mh_sizesOfChains << " = ";
  for (unsigned int i = 0; i < m_sizesOfChains.size(); ++i) {
    os << m_sizesOfChains[i] << " ";
  }
  os << "\n" << m_option_dr_maxNumberOfExtraStages << " = " << m_maxNumberOfExtraStages
     << "\n" << m_option_dr_scalesForExtraStages << " = ";
  for (unsigned int i = 0; i < m_scalesForCovMProposals.size(); ++i) {
    os << m_scalesForCovMProposals[i] << " ";
  }
  os << "\n" << m_option_am_initialNonAdaptInterval << " = " << m_initialNonAdaptInterval
     << "\n" << m_option_am_adaptInterval           << " = " << m_adaptInterval
     << "\n" << m_option_am_eta                     << " = " << m_eta
     << "\n" << m_option_am_epsilon                 << " = " << m_epsilon
     << "\n" << m_option_mh_useChain2               << " = " << m_useChain2
     << "\n" << m_option_mh_generateUniqueChain     << " = " << m_generateUniqueChain
     << "\n" << m_option_mh_generateExtraChains     << " = " << m_generateExtraChains
     << "\n" << m_option_mh_namesOfOutputFiles << " = ";
  for (unsigned int i = 0; i < m_namesOfOutputFiles.size(); ++i) {
    os << m_namesOfOutputFiles[i] << " ";
  }
  os << "\n" << m_option_mh_chainDisplayPeriod      << " = " << m_chainDisplayPeriod
     << "\n" << m_option_mh_measureRunTimes         << " = " << m_measureRunTimes
     << "\n" << m_option_mh_generateWhiteNoise      << " = " << m_generateWhiteNoise
     << "\n" << m_option_finalPercentsForStats << " = ";
  for (unsigned int i = 0; i < m_finalPercentsForStats.size(); ++i) {
    os << m_finalPercentsForStats[i] << " ";
  }
  os << "\n" << m_option_runBMM                     << " = " << m_runBMM;
  for (unsigned int i = 0; i < m_lengthsForBMM.size(); ++i) {
    os << m_lengthsForBMM[i] << " ";
  }
  os << "\n" << m_option_computePSDs                << " = " << m_computePSDs
     << "\n" << m_option_computeGewekeCoefs         << " = " << m_computeGewekeCoefs
     << "\n" << m_option_computeCorrelations        << " = " << m_computeCorrelations
     << "\n" << m_option_secondLagForCorrs          << " = " << m_secondLagForCorrs
     << "\n" << m_option_lagSpacingForCorrs         << " = " << m_lagSpacingForCorrs
     << "\n" << m_option_numberOfLagsForCorrs       << " = " << m_numberOfLagsForCorrs
     << "\n" << m_option_printCorrs                 << " = " << m_printCorrs
     << "\n" << m_option_writeCorrs                 << " = " << m_writeCorrs
     << "\n" << m_option_computeHistograms          << " = " << m_computeHistograms
     << "\n" << m_option_computeKDEs                << " = " << m_computeKDEs
     << "\n" << "(internal variable) m_likelihoodObjComputesMisfits = " << m_likelihoodObjComputesMisfits
     << std::endl;

  return;
}

template <class V, class M>
std::ostream& operator<<(std::ostream& os, const uqDRAM_MarkovChainGeneratorClass<V,M>& obj)
{
  obj.print(os);

  return os;
}
#endif // __UQ_DRAM_MCG_H__
