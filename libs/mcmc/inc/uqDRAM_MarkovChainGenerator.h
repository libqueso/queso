/* uq/libs/mcmc/inc/uqDRAM_MarkovChainGenerator.h
 *
 * Copyright (C) 2008 The PECOS Team, http://pecos.ices.utexas.edu
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
//   LR = likelihood results
#define UQ_MCMC_MH_CHAIN_SIZE_ODV              "100"
#define UQ_MCMC_MH_LR_UPDATE_SIGMA2_ODV        false
#define UQ_MCMC_MH_LR_SIGMA2_PRIORS_ODV        "1."
#define UQ_MCMC_MH_LR_SIGMA2_ACCURACIES_ODV    "1."
#define UQ_MCMC_MH_LR_NUMBERS_OF_OBS_ODV       "1."
#define UQ_MCMC_DR_MAX_NUM_EXTRA_STAGES_ODV    0
#define UQ_MCMC_DR_SCALES_FOR_EXTRA_STAGES_ODV "1."
#define UQ_MCMC_AM_INIT_NON_ADAPT_INT_ODV      0
#define UQ_MCMC_AM_ADAPT_INTERVAL_ODV          0
#define UQ_MCMC_AM_ETA_ODV                     1.
#define UQ_MCMC_AM_EPSILON_ODV                 1.e-5
#define UQ_MCMC_MH_OUTPUT_FILE_NAME_ODV        "."
#define UQ_MCMC_MH_CHAIN_DISPLAY_PERIOD_ODV    500
#define UQ_MCMC_RUN_BMM_ODV                    0
#define UQ_MCMC_COMPUTE_PSDS_ODV               0
#define UQ_MCMC_COMPUTE_GEWEKE_COEFS_ODV       0
#define UQ_MCMC_COMPUTE_CORRELATIONS_ODV       0
#define UQ_MCMC_COMPUTE_HISTOGRAMS_ODV         0
#define UQ_MCMC_COMPUTE_KDES_ODV               0

#include <uqPriorAndLikelihoodInterfaces.h>
#include <uqParamSpace.h>
#include <uqStateSpace.h>
#include <uqOutputSpace.h>
#include <uqChainPosition.h>
#include <uqMiscellaneous.h>
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
  uqDRAM_MarkovChainGeneratorClass(const uqEnvironmentClass&                  env,                        /*! The QUESO toolkit environment.   */
                                   const char*                                prefix,                     /*! Prefix for the validation phase. */
                                   const uqParamSpaceClass<V,M>&              paramSpace,                 /*! The parameter space.             */
                                   const uqFinDimLinearSpaceClass<V,M>&       outputSpace,                /*! The output space.                */
                                         uq_M2lPriorFunction_Class<V,M>*      m2lPriorFunction_ObjPtr,    /*! -2*ln(prior())                   */
                                         uq_M2lLikelihoodFunction_Class<V,M>& m2lLikelihoodFunction_Obj); /*! -2*ln(likelihood())              */
 ~uqDRAM_MarkovChainGeneratorClass();

  void generateChains             (const M* proposalCovMatrix,
                                   //const M* proposalPrecMatrix,
                                   void*    m2lPriorFunction_DataPtr,
                                   void*    m2lLikelihoodFunction_DataPtr,
                                   const M* mahalanobisMatrix = NULL,
                                   bool     applyMahalanobisInvert = true);
  void print                      (std::ostream& os) const;

  const std::vector<const V*>& chain         () const;
  const std::vector<const V*>& lrSigma2Chain () const;
  const std::string&           outputFileName() const;

private:
  void   resetChainAndRelatedInfo();
  double computeM2lPrior         (const V& paramValues, const void* functionDataPtr);
  void   defineMyOptions         (po::options_description& optionsDesc) const;
  void   getMyOptionValues       (po::options_description& optionsDesc);

  int    prepareForNextChain     (const M* proposalCovMatrix);
                                  //const M* proposalPrecMatrix,
  int    generateChain           (unsigned int chainId,
                                  const V&     valuesOf1stPosition,
                                  const M*     proposalCovMatrix,
                                  void*        m2lPriorFunction_DataPtr,
                                  void*        m2lLikelihoodFunction_DataPtr,
                                  const M*     mahalanobisMatrix = NULL,
                                  bool         applyMahalanobisInvert = true);
  int    writeChainInfoOut       (unsigned int chainId,
                                  const M*     mahalanobisMatrix = NULL,
                                  bool         applyMahalanobisInvert = true) const;

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
  void   updateCovMatrix         (const std::vector<V*>& subChain,
                                  unsigned int           idOfFirstPositionInSubChain,
                                  double&                lastChainSize,
                                  V&                     lastMean,
                                  M&                     lastAdaptedCovMatrix);
  void   gammar                  (double a,
                                  double b,
                                  M&     mat);

  const uqEnvironmentClass&                  m_env;
        std::string                          m_prefix;
  const uqParamSpaceClass<V,M>&              m_paramSpace;
  const uqFinDimLinearSpaceClass<V,M>&       m_outputSpace;
        uq_M2lPriorFunction_Class<V,M>*      m_m2lPriorFunction_ObjPtr;
        uq_M2lLikelihoodFunction_Class<V,M>& m_m2lLikelihoodFunction_Obj;

  std::string m_option_help;
  std::string m_option_mh_sizesOfChains;
  std::string m_option_mh_lrUpdateSigma2;
  std::string m_option_mh_lrSigma2Priors;
  std::string m_option_mh_lrSigma2Accuracies;
  std::string m_option_mh_lrNumbersOfObs;
  std::string m_option_dr_maxNumberOfExtraStages;
  std::string m_option_dr_scalesForExtraStages;
  std::string m_option_am_initialNonAdaptInterval;
  std::string m_option_am_adaptInterval;
  std::string m_option_am_eta;
  std::string m_option_am_epsilon;
  std::string m_option_mh_namesOfOutputFiles;
  std::string m_option_mh_chainDisplayPeriod;
  std::string m_option_runBMM;
  std::string m_option_initialPosForBMM;
  std::string m_option_lengthsForBMM;
  std::string m_option_computePSDs;
  std::string m_option_initialPosForPSD;
  std::string m_option_numBlocksForPSD;
  std::string m_option_hopSizeRatioForPSD;
  std::string m_option_computeGewekeCoefs;
  std::string m_option_initialPosForGeweke;
  std::string m_option_ratioNaForGeweke;
  std::string m_option_ratioNbForGeweke;
  std::string m_option_computeCorrelations;
  std::string m_option_initialPosForCorrs;
  std::string m_option_lagsForCorrs;
  std::string m_option_computeHistograms;
  std::string m_option_numberOfInternalBinsForHists;
  std::string m_option_computeKDEs;
  std::string m_option_numberOfEvaluationPosForKDEs;

  bool                      m_useDefaultM2lPriorFunction;
  uqDefault_M2lPriorFunction_DataType<V,M> m_defaultM2lPriorFunction_Data;
  V                         m_paramInitials;
  bool                      m_proposalIsSymmetric;
  po::options_description*  m_optionsDesc;
  std::vector<unsigned int> m_sizesOfChains;
  bool                      m_lrUpdateSigma2; // lr = Likelihood Results
  V*                        m_lrSigma2Priors;
  V*                        m_lrSigma2Accuracies;
  V*                        m_lrNumbersOfObservations;
  unsigned int              m_maxNumberOfExtraStages;
  std::vector<double>       m_scalesForCovMProposals;
  unsigned int              m_initialNonAdaptInterval;
  unsigned int              m_adaptInterval;
  double                    m_eta;
  double                    m_epsilon;
  std::vector<std::string>  m_namesOfOutputFiles;
  unsigned int              m_chainDisplayPeriod;
  bool                      m_runBMM;
  std::vector<unsigned int> m_initialPosForBMM;
  std::vector<unsigned int> m_lengthsForBMM;
  bool                      m_computePSDs;
  std::vector<unsigned int> m_initialPosForPSD;
  std::vector<unsigned int> m_numBlocksForPSD;
  double                    m_hopSizeRatioForPSD;
  bool                      m_computeGewekeCoefs;
  std::vector<unsigned int> m_initialPosForGeweke;
  double                    m_ratioNaForGeweke;
  double                    m_ratioNbForGeweke;
  bool                      m_computeCorrelations;
  std::vector<unsigned int> m_initialPosForCorrs;
  std::vector<unsigned int> m_lagsForCorrs;
  unsigned int              m_initialPosForUncorrelation;   // set during run time
  unsigned int              m_spacingForUncorrelation;      // set during run time
  V*                        m_minPositionsForStatistics;    // set during run time
  V*                        m_maxPositionsForStatistics;    // set during run time
  bool                      m_computeHistograms;
  unsigned int              m_numberOfInternalBinsForHists;
  std::vector<V*>           m_centersForAllHistogramBins;   // set during run time
  std::vector<V*>           m_histogramBinsForAllParams;    // set during run time
  bool                      m_computeKDEs;
  unsigned int              m_numberOfEvaluationPosForKDEs;
  std::vector<V*>           m_evaluationPositionsForKDEs;   // set during run time
  V*                        m_scalesForKDEs;                // set during run time 
  std::vector<V*>           m_densityValuesFromGaussianKDE; // set during run time

  std::vector<      M*>     m_lowerCholProposalCovMatrices;
  std::vector<      M*>     m_proposalCovMatrices;
#ifdef UQ_DRAM_MCG_REQUIRES_INVERTED_COV_MATRICES
  std::vector<      M*>     m_upperCholProposalPrecMatrices;
  std::vector<      M*>     m_proposalPrecMatrices;
#endif

  std::vector<const V*>     m_chain;
  std::vector<const V*>     m_lrChain;       // lr = likelihood results
  std::vector<const V*>     m_lrSigma2Chain; // lr = likelihood results
  std::vector<double>       m_alphaQuotients;
  double                    m_chainRunTime;
  unsigned int              m_numRejections;
  unsigned int              m_numOutOfBounds;
  double                    m_lastChainSize;
  V*                        m_lastMean;
  M*                        m_lastAdaptedCovMatrix;
};

template <class V, class M>
std::ostream& operator<<(std::ostream& os, const uqDRAM_MarkovChainGeneratorClass<V,M>& obj);

template <class V, class M>
uqDRAM_MarkovChainGeneratorClass<V,M>::uqDRAM_MarkovChainGeneratorClass(
  const uqEnvironmentClass&                  env,
  const char*                                prefix,
  const uqParamSpaceClass<V,M>&              paramSpace,
  const uqFinDimLinearSpaceClass<V,M>&       outputSpace,
        uq_M2lPriorFunction_Class<V,M>*      m2lPriorFunction_ObjPtr,
        uq_M2lLikelihoodFunction_Class<V,M>& m2lLikelihoodFunction_Obj)
  :
  m_env                          (env),
  m_prefix                       (""),
  m_paramSpace                   (paramSpace),
  m_outputSpace                  (outputSpace),
  m_m2lPriorFunction_ObjPtr      (m2lPriorFunction_ObjPtr),
  m_m2lLikelihoodFunction_Obj    (m2lLikelihoodFunction_Obj),
  m_useDefaultM2lPriorFunction   (m2lPriorFunction_ObjPtr == NULL),
  m_paramInitials                (m_paramSpace.initialValues()),
  m_proposalIsSymmetric          (true),
  m_optionsDesc                  (new po::options_description("Markov chain Monte Carlo options")),
  m_sizesOfChains                (1,(unsigned int) strtod(UQ_MCMC_MH_CHAIN_SIZE_ODV,NULL)),
  m_lrUpdateSigma2               (UQ_MCMC_MH_LR_UPDATE_SIGMA2_ODV),
  m_lrSigma2Priors               (m_outputSpace.newVector()),
  m_lrSigma2Accuracies           (m_outputSpace.newVector()),
  m_lrNumbersOfObservations      (m_outputSpace.newVector()),
  m_maxNumberOfExtraStages       (UQ_MCMC_DR_MAX_NUM_EXTRA_STAGES_ODV),
  m_scalesForCovMProposals       (0,0.),
  m_initialNonAdaptInterval      (UQ_MCMC_AM_INIT_NON_ADAPT_INT_ODV),
  m_adaptInterval                (UQ_MCMC_AM_ADAPT_INTERVAL_ODV),
  m_eta                          (UQ_MCMC_AM_ETA_ODV),
  m_epsilon                      (UQ_MCMC_AM_EPSILON_ODV),
  m_namesOfOutputFiles           (1,UQ_MCMC_MH_OUTPUT_FILE_NAME_ODV),
  m_chainDisplayPeriod           (UQ_MCMC_MH_CHAIN_DISPLAY_PERIOD_ODV),
  m_runBMM                       (false),
  m_initialPosForBMM             (0,0),
  m_lengthsForBMM                (0,0),
  m_computePSDs                  (false),
  m_initialPosForPSD             (0,0),
  m_numBlocksForPSD              (0,0),
  m_hopSizeRatioForPSD           (1.),
  m_computeGewekeCoefs           (false),
  m_initialPosForGeweke          (0,0),
  m_ratioNaForGeweke             (0.),
  m_ratioNbForGeweke             (0.),
  m_computeCorrelations          (false),
  m_initialPosForCorrs           (0,0),
  m_lagsForCorrs                 (0,0),
  m_initialPosForUncorrelation   (0),
  m_spacingForUncorrelation      (0),
  m_minPositionsForStatistics    (m_paramSpace.newVector()),
  m_maxPositionsForStatistics    (m_paramSpace.newVector()),
  m_computeHistograms            (false),
  m_numberOfInternalBinsForHists (0),
  m_centersForAllHistogramBins   (0,NULL),
  m_histogramBinsForAllParams    (0,NULL),
  m_computeKDEs                  (false),
  m_numberOfEvaluationPosForKDEs (0),
  m_evaluationPositionsForKDEs   (0,NULL),
  m_scalesForKDEs                (m_paramSpace.newVector()),
  m_densityValuesFromGaussianKDE (0,NULL),
  m_lowerCholProposalCovMatrices (1,NULL),
  m_proposalCovMatrices          (1,NULL),
#ifdef UQ_DRAM_MCG_REQUIRES_INVERTED_COV_MATRICES
  m_upperCholProposalPrecMatrices(1,NULL),
  m_proposalPrecMatrices         (1,NULL),
#endif
  m_chain                        (0,NULL),
  m_lrChain                      (0,NULL),
  m_lrSigma2Chain                (0,NULL),
  m_alphaQuotients               (0,0.),
  m_chainRunTime                 (0.),
  m_numRejections                (0),
  m_numOutOfBounds               (0),
  m_lastChainSize                (0),
  m_lastMean                     (NULL),
  m_lastAdaptedCovMatrix         (NULL)
{
  if ((prefix         != NULL) && 
      (strlen(prefix) != 0   )) {
    std::string tmpString(prefix);
    m_prefix = tmpString + "_";
  }

  m_option_help                         = m_prefix + "MCMC_help";
  m_option_mh_sizesOfChains             = m_prefix + "MCMC_mh_sizesOfChains";
  m_option_mh_lrUpdateSigma2            = m_prefix + "MCMC_mh_lrUpdateSigma2";
  m_option_mh_lrSigma2Priors            = m_prefix + "MCMC_mh_lrSigma2Priors";
  m_option_mh_lrSigma2Accuracies        = m_prefix + "MCMC_mh_lrSigma2Accuracies";
  m_option_mh_lrNumbersOfObs            = m_prefix + "MCMC_mh_lrNumbersOfObs";
  m_option_dr_maxNumberOfExtraStages    = m_prefix + "MCMC_dr_maxNumberOfExtraStages";
  m_option_dr_scalesForExtraStages      = m_prefix + "MCMC_dr_scalesForExtraStages";
  m_option_am_initialNonAdaptInterval   = m_prefix + "MCMC_am_initialNonAdaptInterval";
  m_option_am_adaptInterval             = m_prefix + "MCMC_am_adaptInterval";
  m_option_am_eta                       = m_prefix + "MCMC_am_eta";
  m_option_am_epsilon                   = m_prefix + "MCMC_am_epsilon";
  m_option_mh_namesOfOutputFiles        = m_prefix + "MCMC_mh_namesOfOutputFiles";
  m_option_mh_chainDisplayPeriod        = m_prefix + "MCMC_mh_chainDisplayPeriod";
  m_option_runBMM                       = m_prefix + "MCMC_runBMM";
  m_option_initialPosForBMM             = m_prefix + "MCMC_initialPositionsForBMM";
  m_option_lengthsForBMM                = m_prefix + "MCMC_lengthsForBMM";
  m_option_computePSDs                  = m_prefix + "MCMC_computePSDs";
  m_option_initialPosForPSD             = m_prefix + "MCMC_initialPositionsForPSD";
  m_option_numBlocksForPSD              = m_prefix + "MCMC_numBlocksForPSD";
  m_option_hopSizeRatioForPSD           = m_prefix + "MCMC_hopSizeRatioForPSD";
  m_option_computeGewekeCoefs           = m_prefix + "MCMC_computeGewekeCoefficients";
  m_option_initialPosForGeweke          = m_prefix + "MCMC_initialPosForGeweke";
  m_option_ratioNaForGeweke             = m_prefix + "MCMC_ratioNaForGeweke";
  m_option_ratioNbForGeweke             = m_prefix + "MCMC_ratioNbForGeweke";
  m_option_computeCorrelations          = m_prefix + "MCMC_computeCorrelations";
  m_option_initialPosForCorrs           = m_prefix + "MCMC_initialPositionsForCorrs";
  m_option_lagsForCorrs                 = m_prefix + "MCMC_lagsForCorrs";
  m_option_computeHistograms            = m_prefix + "MCMC_computeHistograms";
  m_option_numberOfInternalBinsForHists = m_prefix + "MCMC_numberOfInternalBinsForHistograms";
  m_option_computeKDEs                  = m_prefix + "MCMC_computeKDEs";
  m_option_numberOfEvaluationPosForKDEs = m_prefix + "MCMC_numberOfEvaluationPositionsForKDEs";

  if (m_useDefaultM2lPriorFunction) {
    m_m2lPriorFunction_ObjPtr = new uq_M2lPriorFunction_Class<V,M>();
    m_defaultM2lPriorFunction_Data.paramPriorMus    = new V(m_paramSpace.priorMuValues   ());
    m_defaultM2lPriorFunction_Data.paramPriorSigmas = new V(m_paramSpace.priorSigmaValues());
  }
  m_lrSigma2Priors->cwSet         (strtod(UQ_MCMC_MH_LR_SIGMA2_PRIORS_ODV    ,NULL));
  m_lrSigma2Accuracies->cwSet     (strtod(UQ_MCMC_MH_LR_SIGMA2_ACCURACIES_ODV,NULL));
  m_lrNumbersOfObservations->cwSet(strtod(UQ_MCMC_MH_LR_NUMBERS_OF_OBS_ODV   ,NULL));

  m_initialPosForBMM.resize(5);
  m_initialPosForBMM[0] =  1000;
  m_initialPosForBMM[1] =  5000;
  m_initialPosForBMM[2] = 10000;
  m_initialPosForBMM[3] = 15000;
  m_initialPosForBMM[4] = 20000;

  m_lengthsForBMM.resize(5);
  m_lengthsForBMM[0] =  100;
  m_lengthsForBMM[1] =  200;
  m_lengthsForBMM[2] =  500;
  m_lengthsForBMM[3] = 1000;
  m_lengthsForBMM[4] = 5000;

  m_initialPosForPSD.resize(1);
  m_initialPosForPSD[0] = 20000;

  m_numBlocksForPSD.resize(1);
  m_numBlocksForPSD[0] = 8;

  m_hopSizeRatioForPSD = .50;

  m_initialPosForGeweke.resize(4);
  m_initialPosForGeweke[0] =  5000;
  m_initialPosForGeweke[1] = 10000;
  m_initialPosForGeweke[2] = 15000;
  m_initialPosForGeweke[3] = 20000;

  m_ratioNaForGeweke = .1;
  m_ratioNbForGeweke = .5;

  m_initialPosForCorrs.resize(5);
  m_initialPosForCorrs[0] =  1000;
  m_initialPosForCorrs[1] =  5000;
  m_initialPosForCorrs[2] = 10000;
  m_initialPosForCorrs[3] = 15000;
  m_initialPosForCorrs[4] = 20000;

  m_lagsForCorrs.resize(12);
  m_lagsForCorrs[ 0] =   1;
  m_lagsForCorrs[ 1] =   5;
  m_lagsForCorrs[ 2] =  10;
  m_lagsForCorrs[ 3] =  20;
  m_lagsForCorrs[ 4] =  30;
  m_lagsForCorrs[ 5] =  40;
  m_lagsForCorrs[ 6] =  50;
  m_lagsForCorrs[ 7] =  60;
  m_lagsForCorrs[ 8] =  70;
  m_lagsForCorrs[ 9] =  80;
  m_lagsForCorrs[10] =  90;
  m_lagsForCorrs[11] = 100;

  m_initialPosForUncorrelation   = 20000;
  m_spacingForUncorrelation      = 50;

  m_numberOfInternalBinsForHists = 100;

  m_numberOfEvaluationPosForKDEs = 100;

  std::cout << "Entering uqDRAM_MarkovChainGeneratorClass<V,M>::constructor()"
            << std::endl;

  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.rank() == 0) std::cout << "After getting option values, state of uqDRAM_MarkovChainGeneratorClass object is:"
                                   << "\n" << *this
                                   << std::endl;

  std::cout << "Leaving uqDRAM_MarkovChainGeneratorClass<V,M>::constructor()"
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

  m_lagsForCorrs.clear();
  m_initialPosForCorrs.clear();
  m_initialPosForGeweke.clear();
  m_numBlocksForPSD.clear();
  m_initialPosForPSD.clear();
  m_lengthsForBMM.clear();
  m_initialPosForBMM.clear();

  if (m_lrNumbersOfObservations) delete m_lrNumbersOfObservations;
  if (m_lrSigma2Accuracies     ) delete m_lrSigma2Accuracies;
  if (m_lrSigma2Priors         ) delete m_lrSigma2Priors;
  if (m_optionsDesc            ) delete m_optionsDesc;

  if (m_useDefaultM2lPriorFunction) {
    delete m_m2lPriorFunction_ObjPtr;
    delete m_defaultM2lPriorFunction_Data.paramPriorMus;
    delete m_defaultM2lPriorFunction_Data.paramPriorSigmas;
  }

  //std::cout << "Leaving uqDRAM_MarkovChainGeneratorClass<V,M>::destructor()"
  //          << std::endl;
}

template <class V, class M>
void
uqDRAM_MarkovChainGeneratorClass<V,M>::resetChainAndRelatedInfo()
{
  if (m_lastAdaptedCovMatrix) delete m_lastAdaptedCovMatrix;
  if (m_lastMean)             delete m_lastMean;
  m_lastChainSize = 0;
  m_numOutOfBounds   = 0;
  m_chainRunTime     = 0.;
  m_numRejections    = 0;
  m_alphaQuotients.clear();
  for (unsigned int i = 0; i < m_lrSigma2Chain.size(); ++i) {
    if (m_lrSigma2Chain[i]) delete m_lrSigma2Chain[i];
  }
  for (unsigned int i = 0; i < m_lrChain.size(); ++i) {
    if (m_lrChain[i]) delete m_lrChain[i];
  }
  for (unsigned int i = 0; i < m_chain.size(); ++i) {
    if (m_chain[i]) delete m_chain[i];
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
  for (unsigned int i = 0; i < m_lowerCholProposalCovMatrices.size(); ++i) {
    if (m_lowerCholProposalCovMatrices[i]) delete m_lowerCholProposalCovMatrices[i];
  }

  return;
}

template <class V, class M>
double
uqDRAM_MarkovChainGeneratorClass<V,M>::computeM2lPrior(
  const V&    paramValues,
  const void* functionDataPtr)
{
  double value = 0.;
  if (m_useDefaultM2lPriorFunction) {
    value = (*m_m2lPriorFunction_ObjPtr)(paramValues, (void *) &m_defaultM2lPriorFunction_Data, true);
  }
  else {
    value = (*m_m2lPriorFunction_ObjPtr)(paramValues, functionDataPtr, false);
  }

  return value;
}

template <class V, class M>
void
uqDRAM_MarkovChainGeneratorClass<V,M>::defineMyOptions(
  po::options_description& optionsDesc) const
{
  optionsDesc.add_options()
    (m_option_help.c_str(),                                                                                                           "produce help message for DRAM Markov chain generator"   )
    (m_option_mh_sizesOfChains.c_str(),             po::value<std::string >()->default_value(UQ_MCMC_MH_CHAIN_SIZE_ODV             ), "'mh' size(s) of chain(s)"                               )
    (m_option_mh_lrUpdateSigma2.c_str(),            po::value<bool        >()->default_value(UQ_MCMC_MH_LR_UPDATE_SIGMA2_ODV       ), "'mh' update the sigma2 of likelihood results"           )
    (m_option_mh_lrSigma2Priors.c_str(),            po::value<std::string >()->default_value(UQ_MCMC_MH_LR_SIGMA2_PRIORS_ODV       ), "'mh' priors of sigma2 on likelihood results"            )
    (m_option_mh_lrSigma2Accuracies.c_str(),        po::value<std::string >()->default_value(UQ_MCMC_MH_LR_SIGMA2_ACCURACIES_ODV   ), "'mh' accuracies of sigma2 on likelihood results"        )
    (m_option_mh_lrNumbersOfObs.c_str(),            po::value<std::string >()->default_value(UQ_MCMC_MH_LR_NUMBERS_OF_OBS_ODV      ), "'mh' numbers of observations for likelihood results"    )
    (m_option_dr_maxNumberOfExtraStages.c_str(),    po::value<unsigned int>()->default_value(UQ_MCMC_DR_MAX_NUM_EXTRA_STAGES_ODV   ), "'dr' maximum number of extra stages"                    )
    (m_option_dr_scalesForExtraStages.c_str(),      po::value<std::string >()->default_value(UQ_MCMC_DR_SCALES_FOR_EXTRA_STAGES_ODV), "'dr' scales for proposal cov matrices from 2nd stage on")
    (m_option_am_initialNonAdaptInterval.c_str(),   po::value<unsigned int>()->default_value(UQ_MCMC_AM_INIT_NON_ADAPT_INT_ODV     ), "'am' initial non adaptation interval"                   )
    (m_option_am_adaptInterval.c_str(),             po::value<unsigned int>()->default_value(UQ_MCMC_AM_ADAPT_INTERVAL_ODV         ), "'am' adaptation interval"                               )
    (m_option_am_eta.c_str(),                       po::value<double      >()->default_value(UQ_MCMC_AM_ETA_ODV                    ), "'am' eta"                                               )
    (m_option_am_epsilon.c_str(),                   po::value<double      >()->default_value(UQ_MCMC_AM_EPSILON_ODV                ), "'am' epsilon"                                           )
    (m_option_mh_namesOfOutputFiles.c_str(),        po::value<std::string >()->default_value(UQ_MCMC_MH_OUTPUT_FILE_NAME_ODV       ), "'mh' name(s) of output file(s)"                         )
    (m_option_mh_chainDisplayPeriod.c_str(),        po::value<unsigned int>()->default_value(UQ_MCMC_MH_CHAIN_DISPLAY_PERIOD_ODV   ), "'mh' period of message display during chain generation" )
    (m_option_runBMM.c_str(),                       po::value<bool        >()->default_value(UQ_MCMC_RUN_BMM_ODV                   ), "compute sample variance with batch means method"        )
#if 0
    (m_option_initialPosForBMM.c_str(),             po::value<            >()->default_value(), "")
    (m_option_lengthsForBMM.c_str(),                po::value<            >()->default_value(), "")
#endif
    (m_option_computePSDs.c_str(),                  po::value<bool        >()->default_value(UQ_MCMC_COMPUTE_PSDS_ODV              ), "compute power spectral densities"                       )
#if 0
    (m_option_initialPosForPSD.c_str(),             po::value<            >()->default_value(), "")
    (m_option_numBlocksForPSD.c_str(),              po::value<            >()->default_value(), "")
    (m_option_hopSizeRatioForPSD.c_str(),           po::value<            >()->default_value(), "")
#endif
    (m_option_computeGewekeCoefs.c_str(),           po::value<bool        >()->default_value(UQ_MCMC_COMPUTE_GEWEKE_COEFS_ODV      ), "compute Geweke coefficients"                            )
#if 0
    (m_option_initialPosForGeweke.c_str(),          po::value<            >()->default_value(), "")
    (m_option_ratioNaForGeweke.c_str(),             po::value<            >()->default_value(), "")
    (m_option_ratioNbForGeweke.c_str(),             po::value<            >()->default_value(), "")
#endif
    (m_option_computeCorrelations.c_str(),          po::value<bool        >()->default_value(UQ_MCMC_COMPUTE_CORRELATIONS_ODV      ), "compute correlations"                                   )
#if 0
    (m_option_initialPosForCorrs.c_str(),           po::value<            >()->default_value(), "")
    (m_option_lagsForCorrs.c_str(),                 po::value<            >()->default_value(), "")
#endif
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

  if (m_env.allOptionsMap().count(m_option_mh_lrUpdateSigma2.c_str())) {
    m_lrUpdateSigma2 = m_env.allOptionsMap()[m_option_mh_lrUpdateSigma2.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_mh_lrSigma2Priors.c_str())) {
    std::string inputString = m_env.allOptionsMap()[m_option_mh_lrSigma2Priors.c_str()].as<std::string>();
    std::vector<double> inputDoubles(0,0.);
    uqMiscReadDoublesFromString(inputString,inputDoubles);
    delete m_lrSigma2Priors;

    UQ_FATAL_TEST_MACRO(inputDoubles.size() != m_outputSpace.dim(),
                        m_env.rank(),
                        "uqDRAM_MarkovChainGeneratorClass<V,M>::readMyOptionsValues()",
                        "size of array for 'lrSigma2Priors' is not equal to the dimension of the output space");

    m_lrSigma2Priors = m_outputSpace.newVector();
    for (unsigned int i = 0; i < inputDoubles.size(); ++i) {
      (*m_lrSigma2Priors)[i] = inputDoubles[i];
    }
  }

  if (m_env.allOptionsMap().count(m_option_mh_lrSigma2Accuracies.c_str())) {
    std::string inputString = m_env.allOptionsMap()[m_option_mh_lrSigma2Accuracies.c_str()].as<std::string>();
    std::vector<double> inputDoubles(0,0.);
    uqMiscReadDoublesFromString(inputString,inputDoubles);
    delete m_lrSigma2Accuracies;

    UQ_FATAL_TEST_MACRO(inputDoubles.size() != m_outputSpace.dim(),
                        m_env.rank(),
                        "uqDRAM_MarkovChainGeneratorClass<V,M>::readMyOptionsValues()",
                        "size of array for 'lrSigma2Accuracies' is not equal to the dimension of the output space");

    m_lrSigma2Accuracies = m_outputSpace.newVector();
    for (unsigned int i = 0; i < inputDoubles.size(); ++i) {
      (*m_lrSigma2Accuracies)[i] = inputDoubles[i];
    }
  }

  if (m_env.allOptionsMap().count(m_option_mh_lrNumbersOfObs.c_str())) {
    std::string inputString = m_env.allOptionsMap()[m_option_mh_lrNumbersOfObs.c_str()].as<std::string>();
    std::vector<double> inputDoubles(0,0.);
    uqMiscReadDoublesFromString(inputString,inputDoubles);
    delete m_lrNumbersOfObservations;

    UQ_FATAL_TEST_MACRO(inputDoubles.size() != m_outputSpace.dim(),
                        m_env.rank(),
                        "uqDRAM_MarkovChainGeneratorClass<V,M>::readMyOptionsValues()",
                        "size of array for 'lrNumbersOfObs' is not equal to the dimension of the output space");

    m_lrNumbersOfObservations = m_outputSpace.newVector();
    for (unsigned int i = 0; i < inputDoubles.size(); ++i) {
      (*m_lrNumbersOfObservations)[i] = inputDoubles[i];
    }
  }

  bool bRC = ((m_lrSigma2Priors->size() == m_lrSigma2Accuracies->size()     ) &&
              (m_lrSigma2Priors->size() == m_lrNumbersOfObservations->size()));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqDRAM_MarkovChainGeneratorClass<V,M>::readMyOptionsValues()",
                      "sizes of value arrays for 'lr' are not all the same");

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
    m_scalesForCovMProposals.resize(m_maxNumberOfExtraStages+1,1.);
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

  return;
}

template <class V, class M>
void
uqDRAM_MarkovChainGeneratorClass<V,M>::generateChains(
  const M* proposalCovMatrix,
  void*    m2lPriorFunction_DataPtr,
  void*    m2lLikelihoodFunction_DataPtr,
  const M* mahalanobisMatrix,
  bool     applyMahalanobisInvert)
{
  //if (m_env.rank() == 0) std::cout << "Entering uqDRAM_MarkovChainGeneratorClass<V,M>::generateChains()..."
  //                                 << std::endl;

  V valuesOf1stPosition(m_paramInitials);
  for (unsigned int chainId = 0; chainId < m_sizesOfChains.size(); ++chainId) {
    int iRC = UQ_OK_RC;

    //****************************************************
    // Initialize variables before chain loop
    //****************************************************
    if (chainId > 0) {
      valuesOf1stPosition = *(m_chain[m_chain.size()-1]);
      resetChainAndRelatedInfo();
    }

    //****************************************************
    // Initialize m_lowerCholProposalCovMatrices[0]
    // Initialize m_proposalCovMatrices[0]
    //****************************************************
    iRC = prepareForNextChain(proposalCovMatrix);
    UQ_FATAL_RC_MACRO(iRC,
                      m_env.rank(),
                      "uqDRAM_MarkovChainGeneratorClass<V,M>::generateChains()",
                      "improper prepareForNextChain() return");

    //****************************************************
    // Generate chain
    //****************************************************
    iRC = generateChain(chainId,
                        valuesOf1stPosition,
                        proposalCovMatrix,
                        m2lPriorFunction_DataPtr,
                        m2lLikelihoodFunction_DataPtr,
                        mahalanobisMatrix,
                        applyMahalanobisInvert);
    UQ_FATAL_RC_MACRO(iRC,
                      m_env.rank(),
                      "uqDRAM_MarkovChainGeneratorClass<V,M>::generateChains()",
                      "improper generateChain() return");

    //****************************************************
    // Write chain out
    //****************************************************
    iRC = writeChainInfoOut(chainId,
                            mahalanobisMatrix,
                            applyMahalanobisInvert);
    UQ_FATAL_RC_MACRO(iRC,
                      m_env.rank(),
                      "uqDRAM_MarkovChainGeneratorClass<V,M>::writeChainInfoOut()",
                      "improper writeChainInfoOut() return");

    if (m_env.rank() == 0) {
      std::cout << std::endl;
    }
  }

  //if (m_env.rank() == 0) std::cout << "Leaving uqDRAM_MarkovChainGeneratorClass<V,M>::generateChains()"
  //                                 << std::endl;

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
    V* tmpVec = m_paramSpace.newVector();
    for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
      double sigma = m_paramSpace.parameter(i).priorSigma();
      if ((sigma == INFINITY) ||
          (sigma == NAN     )) {
        (*tmpVec)[i] = pow( fabs(m_paramInitials[i])*0.05,2. );
        if ( (*tmpVec)[i] == 0 ) (*tmpVec)[i] = 1.;
      }
      else if (sigma == 0.) {
        (*tmpVec)[i] = 1.;
      }
      else {
        (*tmpVec)[i] = sigma*sigma;
      }
    }
    internalProposalCovMatrix = m_paramSpace.uqFinDimLinearSpaceClass<V,M>::newDiagMatrix(*tmpVec);
    delete tmpVec;

    //if (m_env.rank() == 0) std::cout << "In uqDRAM_MarkovChainGeneratorClass<V,M>::prepareForNextChain()"
    //                                 << ", contents of internally generated proposal cov matrix are:"
    //                                 << std::endl;
    //std::cout << *internalProposalCovMatrix;
    //if (m_env.rank() == 0) std::cout << std::endl;
  }

  m_lowerCholProposalCovMatrices[0] = new M(*internalProposalCovMatrix); 
  iRC = m_lowerCholProposalCovMatrices[0]->chol();
  UQ_FATAL_RC_MACRO(iRC,
                    m_env.rank(),
                    "uqDRAM_MarkovChainGeneratorClass<V,M>::prepareForNextChain()",
                    "proposalCovMatrix is not positive definite");
  m_lowerCholProposalCovMatrices[0]->zeroUpper(false);

  m_proposalCovMatrices[0] = new M(*internalProposalCovMatrix);

  //if (m_env.rank() == 0) std::cout << "In uqDRAM_MarkovChainGeneratorClass<V,M>::prepareForNextChain()"
  //                                 << ", m_lowerCholProposalCovMatrices[0] contents are:"
  //                                 << std::endl;
  //std::cout << *(m_lowerCholProposalCovMatrices[0]);
  //if (m_env.rank() == 0) std::cout << std::endl;

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
    m_lowerCholProposalCovMatrices [i]   = new M(*(m_lowerCholProposalCovMatrices[i-1]));
  *(m_lowerCholProposalCovMatrices [i]) /= scale;
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
int
uqDRAM_MarkovChainGeneratorClass<V,M>::generateChain(
  unsigned int chainId,
  const V&     valuesOf1stPosition,
  const M*     proposalCovMatrix,
  void*        m2lPriorFunction_DataPtr,
  void*        m2lLikelihoodFunction_DataPtr,
  const M*     mahalanobisMatrix,
  bool         applyMahalanobisInvert)
{
  if (m_env.rank() == 0) {
    std::cout << "Generating chain of id " << chainId
              << ", with "                 << m_sizesOfChains[chainId]
              << " positions..."
              << std::endl;
  }

  int iRC = UQ_OK_RC;

  struct timeval timeval0;
  iRC = gettimeofday(&timeval0, NULL);

  bool   outOfBounds = m_paramSpace.outOfBounds(valuesOf1stPosition);
  UQ_FATAL_TEST_MACRO(outOfBounds,
                      m_env.rank(),
                      "uqDRAM_MarkovChainGeneratorClass<V,M>::generateChains()",
                      "paramInitials should not be out of bound");
  double m2lPrior             = computeM2lPrior(valuesOf1stPosition, m2lPriorFunction_DataPtr);
  V*     m2lLikelihoodResults = m_outputSpace.newVector();
  m_m2lLikelihoodFunction_Obj(valuesOf1stPosition, m2lLikelihoodFunction_DataPtr, *m2lLikelihoodResults);
  V      lrSigma2(*m_lrSigma2Priors);
  double logPosterior         = -0.5 * ( m2lPrior + (*m2lLikelihoodResults/lrSigma2).sumOfComponents() );
  uqChainPositionClass<V> currentPosition(m_env,
                                          valuesOf1stPosition,
                                          outOfBounds,
                                          m2lPrior,
                                          *m2lLikelihoodResults,
                                          lrSigma2,
                                          logPosterior);

  V* gaussianVector = m_paramSpace.newVector();
  V* tmpParamValues = m_paramSpace.newVector();
  uqChainPositionClass<V> currentCandidate(m_env);

  //****************************************************
  // Begin chain loop from positionId = 1
  //****************************************************
  m_chain.resize         (m_sizesOfChains[chainId],NULL); 
  m_lrChain.resize       (m_sizesOfChains[chainId],NULL); 
  m_lrSigma2Chain.resize (m_sizesOfChains[chainId],NULL); 
  m_alphaQuotients.resize(m_sizesOfChains[chainId],0.);

  m_chain         [0] = m_paramSpace.uqFinDimLinearSpaceClass<V,M>::newVector(currentPosition.paramValues());
  m_lrSigma2Chain [0] = m_outputSpace.newVector(lrSigma2);
  m_lrChain       [0] = m_outputSpace.newVector(currentPosition.m2lLikelihood());
  m_alphaQuotients[0] = 1.;

  for (unsigned int positionId = 1; positionId < m_sizesOfChains[chainId]; ++positionId) {
    //if (m_env.rank() == 0) std::cout << "In uqDRAM_MarkovChainGeneratorClass<V,M>::generateChain()"
    //                                 << ": beginning chain position of id = " << positionId
    //                                 << ", m_maxNumberOfExtraStages = "       << m_maxNumberOfExtraStages
    //                                 << std::endl;
    unsigned int stageId = 0;

    //****************************************************
    // Loop: generate new parameters
    //****************************************************
    gaussianVector->cwSetGaussian(m_env.rng(),0.,1.);

    *tmpParamValues = currentPosition.paramValues() + *(m_lowerCholProposalCovMatrices[stageId]) * *gaussianVector;
    outOfBounds     = m_paramSpace.outOfBounds(*tmpParamValues);
    if (outOfBounds) {
      m_numOutOfBounds++;
      m2lPrior      = 0.;
      m2lLikelihoodResults->cwSet(INFINITY);
      logPosterior  = -INFINITY;
    }
    else {
      m2lPrior      = computeM2lPrior(*tmpParamValues, m2lPriorFunction_DataPtr);
      m_m2lLikelihoodFunction_Obj(*tmpParamValues, m2lLikelihoodFunction_DataPtr, *m2lLikelihoodResults);
      logPosterior  = -0.5 * ( m2lPrior + (*m2lLikelihoodResults/lrSigma2).sumOfComponents() );
    }
    currentCandidate.set(*tmpParamValues,
                         outOfBounds,
                         m2lPrior,
                         *m2lLikelihoodResults,
                         lrSigma2,
                         logPosterior);

    bool accept = false;
    if (outOfBounds) {
      m_alphaQuotients[positionId] = 0.;
    }
    else {
      accept = acceptAlpha(this->alpha(currentPosition,currentCandidate,&m_alphaQuotients[positionId]));
    }
#if 0 // For debug only
    if (m_env.rank() == 0) std::cout << "In uqDRAM_MarkovChainGeneratorClass<V,M>::generateChain()"
                                     << ": for chain position of id = " << positionId
                                     << " contents of currentCandidate.paramValues() are:"
                                     << std::endl;
    std::cout << currentCandidate.paramValues();
    if (m_env.rank() == 0) std::cout << std::endl;

    if (m_env.rank() == 0) std::cout << "In uqDRAM_MarkovChainGeneratorClass<V,M>::generateChain()"
                                     << ": for chain position of id = " << positionId
                                     << ", outOfBounds = "              << outOfBounds
                                     << ", curM2lPrior = "              << currentPosition.m2lPrior()
                                     << ", curM2lLikelihood = "         << currentPosition.m2lLikelihood()
                                     << ", canM2lPrior = "              << currentCandidate.m2lPrior()
                                     << ", canM2lLikelihoodResults = "  << currentCandidate.m2lLikelihood()
                                     << ", canLrSigma2 = "              << currentCandidate.lrSigma2()
                                     << ", canLogPosterior = "          << currentCandidate.logPosterior()
                                     << ", accept = "                   << accept
                                     << std::endl;
#endif

    //****************************************************
    // Loop: delayed rejection
    //****************************************************
    std::vector<uqChainPositionClass<V>*> drPositions(stageId+2,NULL);
    if ((accept == false) && (outOfBounds == false) && (m_maxNumberOfExtraStages > 0)) {
      drPositions[0] = new uqChainPositionClass<V>(currentPosition);
      drPositions[1] = new uqChainPositionClass<V>(currentCandidate);

      while ((accept == false) && (stageId < m_maxNumberOfExtraStages)) {
        stageId++;

        gaussianVector->cwSetGaussian(m_env.rng(),0.,1.);

        *tmpParamValues = currentPosition.paramValues() + *(m_lowerCholProposalCovMatrices[stageId]) * *gaussianVector;
        outOfBounds   = m_paramSpace.outOfBounds(*tmpParamValues);
        if (outOfBounds) {
          m2lPrior      = 0.;
          m2lLikelihoodResults->cwSet(INFINITY);
          logPosterior  = -INFINITY;
        }
        else {
          m2lPrior      = computeM2lPrior(*tmpParamValues, m2lPriorFunction_DataPtr);
          m_m2lLikelihoodFunction_Obj(*tmpParamValues, m2lLikelihoodFunction_DataPtr, *m2lLikelihoodResults);
          logPosterior  = -0.5 * ( m2lPrior + (*m2lLikelihoodResults/lrSigma2).sumOfComponents() );
        }
        currentCandidate.set(*tmpParamValues,
                             outOfBounds,
                             m2lPrior,
                             *m2lLikelihoodResults,
                             lrSigma2,
                             logPosterior);

        drPositions.push_back(new uqChainPositionClass<V>(currentCandidate));
        if (outOfBounds == false) {
          accept = acceptAlpha(this->alpha(drPositions));
        }
      } // while
    } // end of 'delayed rejection' logic

    for (unsigned int i = 0; i < drPositions.size(); ++i) {
      if (drPositions[i]) delete drPositions[i];
    }

    //****************************************************
    // Loop: update chain
    //****************************************************
    if (accept) {
      m_chain[positionId]   = m_paramSpace.uqFinDimLinearSpaceClass<V,M>::newVector(currentCandidate.paramValues());
      m_lrChain[positionId] = m_outputSpace.newVector(currentCandidate.m2lLikelihood());
      currentPosition = currentCandidate;
    }
    else {
      m_chain[positionId]   = m_paramSpace.uqFinDimLinearSpaceClass<V,M>::newVector(currentPosition.paramValues());
      m_lrChain[positionId] = m_outputSpace.newVector(currentPosition.m2lLikelihood());
      m_numRejections++;
    }

    if (m_lrUpdateSigma2) {
      for (unsigned int i = 0; i < lrSigma2.size(); ++i) {
        double term1 = 0.5*(double)( (*m_lrSigma2Accuracies)[i] + (*m_lrNumbersOfObservations)[i] );
        double term2 = 2./ (double)( (*m_lrSigma2Accuracies)[i] * (*m_lrSigma2Priors)[i] + (*m_lrChain[positionId])[i] );
        lrSigma2[i] = 1./uqMiscGammar(term1,term2,m_env.rng());
        //if (m_env.rank() == 0) {
        //  std::cout << "In uqDRAM_MarkovChainGeneratorClass<V,M>::generateChain()"
        //            << ": for chain position of id = "     << positionId
        //            << ", (*m_lrSigma2Accuracies) = "      << (*m_lrSigma2Accuracies)
        //            << ", (*m_lrNumbersOfObservations) = " << (*m_lrNumbersOfObservations)
        //            << ", (*m_lrSigma2Priors) = "          << (*m_lrSigma2Priors)
        //            << ", (*m_lrChain[positionId]) = "     << (*m_lrChain[positionId])
        //            << ", term1 = "                        << term1
        //            << ", term2 = "                        << term2
        //            << std::endl;
        //}
      }
      //if (m_env.rank() == 0) {
      //  std::cout << "In uqDRAM_MarkovChainGeneratorClass<V,M>::generateChain()"
      //            << ": for chain position of id = " << positionId
      //            << ", lrSigma2 changed from "      << *(m_lrSigma2Chain[positionId])
      //            << " to "                          << lrSigma2
      //            << std::endl;
      //}
    }
    m_lrSigma2Chain[positionId] = m_outputSpace.newVector(lrSigma2);

    //****************************************************
    // Loop: adaptive Metropolis (adaptation of covariance matrix)
    //****************************************************
    if ((m_initialNonAdaptInterval > 0) &&
        (m_adaptInterval           > 0)) {
      // Now might be the moment to adapt
      unsigned int idOfFirstPositionInSubChain = 0;
      std::vector<V*> subChain(0,NULL);

      // Check if now is indeed the moment to adapt
      if (positionId < m_initialNonAdaptInterval) {
        // Do nothing
      }
      else if (positionId == m_initialNonAdaptInterval) {
        idOfFirstPositionInSubChain = 0;
        subChain.resize(m_initialNonAdaptInterval+1,NULL);
        m_lastMean             = m_paramSpace.newVector();
        m_lastAdaptedCovMatrix = m_paramSpace.newMatrix();
      }
      else {
        unsigned int interval = positionId - m_initialNonAdaptInterval;
        if ((interval % m_adaptInterval) == 0) {
          idOfFirstPositionInSubChain = positionId - m_adaptInterval;
          subChain.resize(m_adaptInterval,NULL);
        }
      }

      // If now is indeed the moment to adapt, then do it!
      if (subChain.size() > 0) {
        for (unsigned int i = 0; i < subChain.size(); ++i) {
          subChain[i] = new V(*(m_chain[idOfFirstPositionInSubChain+i]));
        }
        updateCovMatrix(subChain,
			idOfFirstPositionInSubChain,
                        m_lastChainSize,
                        *m_lastMean,
                        *m_lastAdaptedCovMatrix);

        bool tmpCholIsPositiveDefinite = false;
        M tmpChol(*m_lastAdaptedCovMatrix);
        //if (m_env.rank() == 0) {
        //  std::cout << "DRAM: chainId = " << chainId
        //            << ", positionId = "  << positionId
        //            << ": 'am' calling first tmpChol.chol()"
        //            << std::endl;
        //}
        iRC = tmpChol.chol();
        //if (m_env.rank() == 0) {
        //  std::cout << "DRAM: chainId = " << chainId
        //            << ", positionId = "  << positionId
        //            << ": 'am' got first tmpChol.chol() with iRC = " << iRC
        //            << std::endl;
        //}
        if (iRC) {
          UQ_FATAL_TEST_MACRO(iRC != UQ_MATRIX_IS_NOT_POS_DEFINITE_RC,
                              m_env.rank(),
                              "uqDRAM_MarkovChainGeneratorClass<V,M>::generateChain()",
                              "invalid iRC returned from first chol()");
          // Matrix is not positive definite
          M* tmpDiag = m_paramSpace.newDiagMatrix(m_epsilon);
          tmpChol = *m_lastAdaptedCovMatrix + *tmpDiag;
          delete tmpDiag;
          //if (m_env.rank() == 0) {
          //  std::cout << "DRAM: chainId = " << chainId
          //            << ", positionId = "  << positionId
          //            << ": 'am' calling second tmpChol.chol()"
          //            << std::endl;
          //}
          iRC = tmpChol.chol();
          //if (m_env.rank() == 0) {
          //  std::cout << "DRAM: chainId = " << chainId
          //            << ", positionId = "  << positionId
          //            << ": 'am' got second tmpChol.chol() with iRC = " << iRC
          //            << std::endl;
          //}
          if (iRC) {
            UQ_FATAL_TEST_MACRO(iRC != UQ_MATRIX_IS_NOT_POS_DEFINITE_RC,
                                m_env.rank(),
                                "uqDRAM_MarkovChainGeneratorClass<V,M>::generateChain()",
                                "invalid iRC returned from second chol()");
            // Do nothing
          }
          else {
            tmpCholIsPositiveDefinite = true;
          }
        }
        else {
          tmpCholIsPositiveDefinite = true;
        }
        if (tmpCholIsPositiveDefinite) {
          *(m_lowerCholProposalCovMatrices[0]) = tmpChol;
          *(m_lowerCholProposalCovMatrices[0]) *= sqrt(m_eta);
          m_lowerCholProposalCovMatrices[0]->zeroUpper(false);

#ifdef UQ_DRAM_MCG_REQUIRES_INVERTED_COV_MATRICES
          UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                            m_env.rank(),
                            "uqDRAM_MarkovChainGeneratorClass<V,M>::generateChain()",
                            "need to code the update of m_upperCholProposalPrecMatrices");
#endif

          if (m_maxNumberOfExtraStages > 0) updateCovMatrices();
        }

        for (unsigned int i = 0; i < subChain.size(); ++i) {
          if (subChain[i]) delete subChain[i];
        }
      }
    } // End of 'adaptive Metropolis' logic

    if ((m_chainDisplayPeriod                     > 0) && 
        (((positionId+1) % m_chainDisplayPeriod) == 0)) {
      if (m_env.rank() == 0) {
        std::cout << "Finished generating " << positionId+1
                  << " positions"
                  << std::endl;
      }
    }
  } // end chain loop

  //****************************************************
  // Print basic statistics
  //****************************************************
  m_chainRunTime = uqMiscGetEllapsedSeconds(&timeval0);
  if (m_env.rank() == 0) {
    std::cout << "Finished generating the chain of id " << chainId
              << ". Chain statistics are:";
    std::cout << "\n  Run time             = " << m_chainRunTime
              << " seconds";
    std::cout << "\n  Rejection percentage = " << 100. * (double) m_numRejections/(double) m_chain.size()
              << " %";
    std::cout << "\n   Outbound percentage = " << 100. * (double) m_numOutOfBounds/(double) m_chain.size()
              << " %";
    std::cout << std::endl;
  }

  //****************************************************
  // Compute mean, sample std, population std
  //****************************************************
  V* chainMean = m_paramSpace.newVector();
  uqVectorSequenceMean(m_chain,
                       0,
                       m_chain.size(),
                       *chainMean);

  V* chainSampleVariance = m_paramSpace.newVector();
  uqVectorSequenceSampleVariance(m_chain,
                                 0,
                                 m_chain.size(),
                                 *chainMean,
                                 *chainSampleVariance);

  V* chainPopulationVariance = m_paramSpace.newVector();
  uqVectorSequencePopulationVariance(m_chain,
                                     0,
                                     m_chain.size(),
                                     *chainMean,
                                     *chainPopulationVariance);

  if (m_env.rank() == 0) {
    std::cout << "\nMean, sample std, population std"
              << std::endl;
    char line[512];
    sprintf(line,"%s%4s%s%9s%s%9s%s",
	    "Parameter",
            " ",
            "Mean",
            " ",
            "SampleStd",
            " ",
            "Popul.Std");
    std::cout << line;

    for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
      sprintf(line,"\n%8.8s%2s%11.4e%2s%11.4e%2s%11.4e",
              m_paramSpace.parameter(i).name().c_str(),
              " ",
	      (*chainMean)[i],
              " ",
              sqrt((*chainSampleVariance)[i]),
              " ",
              sqrt((*chainPopulationVariance)[i]));
      std::cout << line;
    }
    std::cout << std::endl;
  }

  //****************************************************
  // Compute variance of sample mean through the 'batch means method' (BMM)
  //****************************************************
  if ((m_runBMM                     ) &&
      (m_initialPosForBMM.size() > 0) &&
      (m_lengthsForBMM.size()    > 0)) { 
    if (m_env.rank() == 0) {
      std::cout << "\nComputing variance of sample mean through BMM..."
                << std::endl;
    }
    uq2dArrayOfStuff<V> _2dArrayOfBMM(m_initialPosForBMM.size(),m_lengthsForBMM.size());
    for (unsigned int i = 0; i < _2dArrayOfBMM.numRows(); ++i) {
      for (unsigned int j = 0; j < _2dArrayOfBMM.numCols(); ++j) {
        _2dArrayOfBMM.setLocation(i,j,m_paramSpace.newVector());
      }
    }
    uqVectorSequenceBMM(m_chain,
                        m_initialPosForBMM,
                        m_lengthsForBMM,
                        _2dArrayOfBMM);

    if (m_env.rank() == 0) {
      for (unsigned int initialPosId = 0; initialPosId < m_initialPosForBMM.size(); initialPosId++) {
        std::cout << "\nEstimated covariances of sample mean, through batch means method, for subchain beggining at position " << m_initialPosForBMM[initialPosId]
                  << " (each column corresponds to a batch length)"
                  << std::endl;

        char line[512];
        sprintf(line,"%s",
                "Parameter");
        std::cout << line;
        for (unsigned int batchLengthId = 0; batchLengthId < m_lengthsForBMM.size(); batchLengthId++) {
          sprintf(line,"%9s%3d",
                  " ",
                  m_lengthsForBMM[batchLengthId]);
          std::cout << line;
        }

        for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
          sprintf(line,"\n%8.8s",
                  m_paramSpace.parameter(i).name().c_str());
          std::cout << line;
          for (unsigned int batchLengthId = 0; batchLengthId < m_lengthsForBMM.size(); batchLengthId++) {
            sprintf(line,"%2s%11.4e",
                    " ",
                    _2dArrayOfBMM(initialPosId,batchLengthId)[i]);
            std::cout << line;
          }
        }
        std::cout << std::endl;
      }
    }
  }

  //****************************************************
  // Compute power spectral density (PSD) of chain at zero frequency
  //****************************************************
  if ((m_computePSDs                ) &&
      (m_initialPosForPSD.size() > 0) &&
      (m_numBlocksForPSD.size()  > 0)) { 
    if (m_env.rank() == 0) {
      std::cout << "\nComputing PSD at frequency zero..."
                << std::endl;
    }
    uq2dArrayOfStuff<V> _2dArrayOfPSDAtZero(m_initialPosForPSD.size(),m_numBlocksForPSD.size());
    for (unsigned int i = 0; i < _2dArrayOfPSDAtZero.numRows(); ++i) {
      for (unsigned int j = 0; j < _2dArrayOfPSDAtZero.numCols(); ++j) {
        _2dArrayOfPSDAtZero.setLocation(i,j,m_paramSpace.newVector());
      }
    }
    uqVectorSequencePSD(m_chain,
                        m_initialPosForPSD,
                        m_numBlocksForPSD,
                        m_hopSizeRatioForPSD,
                        _2dArrayOfPSDAtZero);

    if (m_env.rank() == 0) {
      for (unsigned int initialPosId = 0; initialPosId < m_initialPosForPSD.size(); initialPosId++) {
        double sizeForPSD = m_chain.size() - m_initialPosForPSD[initialPosId];
        std::cout << "\nEstimated covariances of sample mean, through psd (fft), for subchain beggining at position " << m_initialPosForPSD[initialPosId]
                  << ", so effective data size = " << sizeForPSD
                  << " (each column corresponds to a number of blocks)"
                  << std::endl;

        char line[512];
        sprintf(line,"%s",
                "Parameter");
        std::cout << line;
        for (unsigned int numBlocksId = 0; numBlocksId < m_numBlocksForPSD.size(); numBlocksId++) {
          sprintf(line,"%9s%3d",
                  " ",
                  m_numBlocksForPSD[numBlocksId]);
          std::cout << line;
        }

        for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
          sprintf(line,"\n%8.8s",
                  m_paramSpace.parameter(i).name().c_str());
          std::cout << line;
          for (unsigned int numBlocksId = 0; numBlocksId < m_numBlocksForPSD.size(); numBlocksId++) {
            sprintf(line,"%2s%11.4e",
                    " ",
                    M_PI*_2dArrayOfPSDAtZero(initialPosId,numBlocksId)[i]/sizeForPSD); // CHECK
            std::cout << line;
          }
        }
        std::cout << std::endl;
      }
    }
  }

  //****************************************************
  // Compute Geweke
  //****************************************************
  if ((m_computeGewekeCoefs            ) &&
      (m_initialPosForGeweke.size() > 0)) {
    if (m_env.rank() == 0) {
      std::cout << "\nComputing Geweke coefficients..."
                << std::endl;
    }
    std::vector<V*> vectorOfGeweke(m_initialPosForGeweke.size(),NULL);
    uqVectorSequenceGeweke(m_chain,
                           m_initialPosForGeweke,
                           m_ratioNaForGeweke,
                           m_ratioNbForGeweke,
                           vectorOfGeweke);

    if (m_env.rank() == 0) {
      std::cout << "\nComputed Geweke coefficients with 10% and 50% percentages"
                  << " (each column corresponds to a different initial position on the full chain)"
                  << std::endl;

      char line[512];
      sprintf(line,"%s",
              "Parameter");
      std::cout << line;
      for (unsigned int initialPosId = 0; initialPosId < m_initialPosForGeweke.size(); initialPosId++) {
        sprintf(line,"%9s%3d",
                " ",
                m_initialPosForGeweke[initialPosId]);
        std::cout << line;
      }

      for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
        sprintf(line,"\n%8.8s",
                m_paramSpace.parameter(i).name().c_str());
        std::cout << line;
        for (unsigned int initialPosId = 0; initialPosId < m_initialPosForGeweke.size(); initialPosId++) {
          sprintf(line,"%2s%11.4e",
                  " ",
                  (*(vectorOfGeweke[initialPosId]))[i]);
          std::cout << line;
        }
      }
      std::cout << std::endl;
    }
  }

  //****************************************************
  // Compute autocorrelation coefficients
  //****************************************************
  if ((m_computeCorrelations          ) &&
      (m_initialPosForCorrs.size() > 0) &&
      (m_lagsForCorrs.size()       > 0)) { 
    if (m_env.rank() == 0) {
      std::cout << "\nComputing autocorrelation coefficients..."
                << std::endl;
    }
    uq2dArrayOfStuff<V> _2dArrayOfAutoCorrs(m_initialPosForCorrs.size(),m_lagsForCorrs.size());
    for (unsigned int i = 0; i < _2dArrayOfAutoCorrs.numRows(); ++i) {
      for (unsigned int j = 0; j < _2dArrayOfAutoCorrs.numCols(); ++j) {
        _2dArrayOfAutoCorrs.setLocation(i,j,m_paramSpace.newVector());
      }
    }
    uqVectorSequenceAutoCorrelations(m_chain,
                                     m_initialPosForCorrs,
                                     m_lagsForCorrs,
                                     _2dArrayOfAutoCorrs);

    if (m_env.rank() == 0) {
      for (unsigned int initialPosId = 0; initialPosId < m_initialPosForCorrs.size(); initialPosId++) {
        std::cout << "\nEstimated autocorrelation coefficients, for subchain beggining at position " << m_initialPosForCorrs[initialPosId]
                  << " (each column corresponds to a different lag)"
                  << std::endl;

        char line[512];
        sprintf(line,"%s",
	        "Parameter");
        std::cout << line;
        for (unsigned int lagId = 0; lagId < m_lagsForCorrs.size(); lagId++) {
          sprintf(line,"%9s%3d",
                  " ",
                  m_lagsForCorrs[lagId]);
          std::cout << line;
        }

        for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
          sprintf(line,"\n%8.8s",
                  m_paramSpace.parameter(i).name().c_str());
          std::cout << line;
          for (unsigned int lagId = 0; lagId < m_lagsForCorrs.size(); lagId++) {
            sprintf(line,"%2s%11.4e",
                    " ",
	            _2dArrayOfAutoCorrs(initialPosId,lagId)[i]);
            std::cout << line;
          }
        }
        std::cout << std::endl;
      }
    }
  }

  //****************************************************
  // Compute MIN and MAX: for histograms and KDE
  //****************************************************
  uqVectorSequenceMinMax(m_chain,
                         m_initialPosForUncorrelation,
                         *m_minPositionsForStatistics,
                         *m_maxPositionsForStatistics);

  //****************************************************
  // Compute histograms
  //****************************************************
  if ((m_computeHistograms               ) &&
      (m_numberOfInternalBinsForHists > 0)) {
    if (m_env.rank() == 0) {
      std::cout << "\nComputing histograms..."
                << std::endl;
    }
    for (unsigned int i = 0; i < m_maxPositionsForStatistics->size(); ++i) {
      (*m_maxPositionsForStatistics)[i] *= (1. + 1.e-15);
    }

    m_centersForAllHistogramBins.resize(m_numberOfInternalBinsForHists+2,NULL);
    m_histogramBinsForAllParams.resize (m_numberOfInternalBinsForHists+2,NULL);
    uqVectorSequenceHistogram(m_chain,
                              m_initialPosForUncorrelation,
                              m_spacingForUncorrelation,
                              *m_minPositionsForStatistics,
                              *m_maxPositionsForStatistics,
                              m_centersForAllHistogramBins,
                              m_histogramBinsForAllParams);
  }

  //****************************************************
  // Compute estimations of probability densities
  //****************************************************
  if ((m_computeKDEs                     ) &&
      (m_numberOfEvaluationPosForKDEs > 0)) {
    if (m_env.rank() == 0) {
      std::cout << "\nComputing KDE..."
                << std::endl;
    }

    m_evaluationPositionsForKDEs.resize(m_numberOfEvaluationPosForKDEs,NULL);
    uqMiscComputePositionsBetweenMinMax(*m_minPositionsForStatistics,
                                        *m_maxPositionsForStatistics,
                                        m_evaluationPositionsForKDEs);

    V iqrs(*(m_chain[0]));
    uqVectorSequenceInterQuantileRange(m_chain,
                                       m_initialPosForUncorrelation,
                                       m_spacingForUncorrelation,
                                       iqrs);
    if (m_env.rank() == 0) {
      std::cout << "\nComputed min values, max values and inter quantile ranges, for subchain beggining at position " << m_initialPosForUncorrelation
                  << std::endl;

      char line[512];
      sprintf(line,"%s",
              "Parameter");
      std::cout << line;

      sprintf(line,"%9s%s%9s%s%9s%s",
              " ",
              "min",
              " ",
              "max",
              " ",
              "iqr");
      std::cout << line;

      for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
        sprintf(line,"\n%8.8s",
                m_paramSpace.parameter(i).name().c_str());
        std::cout << line;

        sprintf(line,"%2s%11.4e%2s%11.4e%2s%11.4e",
                " ",
                (*m_minPositionsForStatistics)[i],
                " ",
                (*m_maxPositionsForStatistics)[i],
                " ",
                iqrs[i]);
        std::cout << line;
      }
      std::cout << std::endl;
    }

    uqVectorSequenceScaleForKDE(m_chain,
                                m_initialPosForUncorrelation,
                                m_spacingForUncorrelation,
                                iqrs,
                                *m_scalesForKDEs);

    m_densityValuesFromGaussianKDE.resize(m_numberOfEvaluationPosForKDEs,NULL);
    uqVectorSequenceGaussianKDE(m_chain,
                                m_initialPosForUncorrelation,
                                m_spacingForUncorrelation,
                                m_evaluationPositionsForKDEs,
                                *m_scalesForKDEs,
                                m_densityValuesFromGaussianKDE);
  }

  //****************************************************
  // Release memory before leaving routine
  //****************************************************
  if (chainPopulationVariance) delete chainPopulationVariance;
  if (chainSampleVariance    ) delete chainSampleVariance;
  if (chainMean              ) delete chainMean;
  if (gaussianVector         ) delete gaussianVector;
  if (tmpParamValues         ) delete tmpParamValues;

  return iRC;
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
  if ((x.outOfBounds() == false) &&
      (y.outOfBounds() == false)) {
    double yLogPosteriorToUse = y.logPosterior();
    if (m_lrUpdateSigma2) {
      yLogPosteriorToUse = -0.5 * ( y.m2lPrior() + (y.m2lLikelihood()/x.lrSigma2()).sumOfComponents() );
    }
    if (m_proposalIsSymmetric) {
      alphaQuotient = exp(yLogPosteriorToUse - x.logPosterior());
    }
    else {
      alphaQuotient = exp(yLogPosteriorToUse + logProposal(y,x,0) - x.logPosterior() - logProposal(x,y,0));
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
  if (m_lrUpdateSigma2) {
    numeratorLogPosteriorToUse = -0.5 * ( backwardPositions[0]->m2lPrior() + (backwardPositions[0]->m2lLikelihood()/positions[0]->lrSigma2()).sumOfComponents() );
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

template <class V, class M>
void
uqDRAM_MarkovChainGeneratorClass<V,M>::updateCovMatrix(
  const std::vector<V*>& subChain,
  unsigned int           idOfFirstPositionInSubChain,
  double&                lastChainSize,
  V&                     lastMean,
  M&                     lastAdaptedCovMatrix)
{
  double doubleSubChainSize = (double) subChain.size();
  if (lastChainSize == 0) {
    UQ_FATAL_TEST_MACRO(subChain.size() < 2,
                        m_env.rank(),
                        "uqDRAM_MarkovChainGeneratorClass<V,M>::updateCovMatrix()",
                        "'subChain.size()' should be >= 2");

    lastMean.cwSet(0.);
    double ratio = 1./doubleSubChainSize;
    for (unsigned int i = 0; i < subChain.size(); ++i) {
      lastMean += ratio * *(subChain[i]);
    }

    lastAdaptedCovMatrix = -doubleSubChainSize * matrixProduct(lastMean,lastMean);
    for (unsigned int i = 0; i < subChain.size(); ++i) {
      lastAdaptedCovMatrix += matrixProduct(*(subChain[i]),*(subChain[i]));
    }
    lastAdaptedCovMatrix /= (doubleSubChainSize - 1.); // That is why subChain size must be >= 2
  }
  else {
    UQ_FATAL_TEST_MACRO(subChain.size() < 1,
                        m_env.rank(),
                        "uqDRAM_MarkovChainGeneratorClass<V,M>::updateCovMatrix()",
                        "'subChain.size()' should be >= 1");

    UQ_FATAL_TEST_MACRO(idOfFirstPositionInSubChain < 1,
                        m_env.rank(),
                        "uqDRAM_MarkovChainGeneratorClass<V,M>::updateCovMatrix()",
                        "'idOfFirstPositionInSubChain' should be >= 1");

    for (unsigned int i = 0; i < subChain.size(); ++i) {
      double doubleCurrentId  = (double) (idOfFirstPositionInSubChain+i);
      V diffVec(*(subChain[i]) - lastMean);

      double ratio1         = (1. - 1./doubleCurrentId); // That is why idOfFirstPositionInSubChain must be >= 1
      double ratio2         = (1./(1.+doubleCurrentId));
      lastAdaptedCovMatrix  = ratio1 * lastAdaptedCovMatrix + ratio2 * matrixProduct(diffVec,diffVec);
      lastMean             += ratio2 * diffVec;
    } 
  }
  lastChainSize += doubleSubChainSize;

  return;
}

template <class V, class M>
void
uqDRAM_MarkovChainGeneratorClass<V,M>::gammar(
  double a,
  double b,
  M&     mat)
{
  for (unsigned int i = 0; i < mat.numRows(); ++i) {
    for (unsigned int j = 0; j < mat.numCols(); ++j) {
      mat(i,j) = uqMiscGammar(a,b,m_env.rng());
    }
  }

  return mat;
}

template <class V, class M>
int
uqDRAM_MarkovChainGeneratorClass<V,M>::writeChainInfoOut(
  unsigned int chainId,
  const M*     mahalanobisMatrix,
  bool         applyMahalanobisInvert) const
{
  int iRC = UQ_OK_RC;

  if (m_namesOfOutputFiles[chainId] == ".") {
    if (m_env.rank() == 0) {
      std::cout << "No info written out for chain of id " << chainId
                << std::endl;
    }
    return iRC;
  }

  if (m_env.rank() == 0) {
    std::cout << "Writing out info about the chain of id " << chainId
              << " into file '"                            << m_namesOfOutputFiles[chainId]
              << "'..."
              << std::endl;
  }
  
  // Open file
  std::ofstream ofs(m_namesOfOutputFiles[chainId].c_str(), std::ofstream::out | std::ofstream::trunc);
  UQ_TEST_MACRO((ofs && ofs.is_open()) == false,
                m_env.rank(),
                "uqDRAM_MarkovChainGeneratorClass<V,M>::writeChainInfoOut()",
                "failed to open file'",
                UQ_FAILED_TO_OPEN_FILE_RC);

  // Write m_chain
  ofs << "queso_" << m_prefix << "chain = [";
  for (unsigned int i = 0; i < m_chain.size(); ++i) {
    ofs << *(m_chain[i])
        << std::endl;
  }
  ofs << "];\n";

  // Write m_lrChain
  ofs << "queso_" << m_prefix << "sschain = [";
  for (unsigned int i = 0; i < m_lrChain.size(); ++i) {
    ofs << *(m_lrChain[i])
        << std::endl;
  }
  ofs << "];\n";

  // Write m_lrSigma2Chain
  ofs << "queso_" << m_prefix << "s2chain = [";
  for (unsigned int i = 0; i < m_lrSigma2Chain.size(); ++i) {
    ofs << *(m_lrSigma2Chain[i])
        << std::endl;
  }
  ofs << "];\n";

  // Write m_alphaQuotients
  ofs << "queso_" << m_prefix << "alphaQuotients = [";
  for (unsigned int i = 0; i < m_alphaQuotients.size(); ++i) {
    ofs << m_alphaQuotients[i]
        << std::endl;
  }
  ofs << "];\n";

  // Write names of parameters
  ofs << "queso_" << m_prefix << "results.names = {";
  m_paramSpace.printParameterNames(ofs,false);
  ofs << "};\n";

  if (mahalanobisMatrix != NULL) {
    // Write mahalanobis distances
    V* diffVec = m_paramSpace.newVector();
    ofs << "queso_" << m_prefix << "d = [";
    if (applyMahalanobisInvert) {
      for (unsigned int i = 0; i < m_chain.size(); ++i) {
        *diffVec = *(m_chain[i]) - *(m_chain[0]);
        ofs << scalarProduct(*diffVec, mahalanobisMatrix->invertMultiply(*diffVec))
            << std::endl;
      }
    }
    else {
      for (unsigned int i = 0; i < m_chain.size(); ++i) {
        *diffVec = *(m_chain[i]) - *(m_chain[0]);
        ofs << scalarProduct(*diffVec, *mahalanobisMatrix * *diffVec)
            << std::endl;
      }
    }
    ofs << "];\n";
    delete diffVec;
  }

  // Write prior mean values
  ofs << "queso_" << m_prefix << "priorMeanValues = ["
      << m_paramSpace.priorMuValues()
      << "];\n";

  // Write prior sigma values
  ofs << "queso_" << m_prefix << "priorSigValues = ["
      << m_paramSpace.priorSigmaValues()
      << "];\n";

  ofs << "queso_" << m_prefix << "results.prior = [queso_priorMeanValues',queso_priorSigValues'];\n";

  // Write param lower bounds
  ofs << "queso_" << m_prefix << "low = ["
      << m_paramSpace.minValues()
      << "];\n";

  // Write param upper bounds
  ofs << "queso_" << m_prefix << "upp = ["
      << m_paramSpace.maxValues()
      << "];\n";

  ofs << "queso_" << m_prefix << "results.limits = [queso_low',queso_upp'];\n";

  // Write out data for mcmcpred.m
  ofs << "queso_" << m_prefix << "results.parind = ["; // FIXME
  for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
    ofs << i+1
        << std::endl;
  }
  ofs << "];\n";

  ofs << "queso_" << m_prefix << "results.local = [\n"; // FIXME
  for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
    ofs << " 0";
    //<< std::endl;
  }
  ofs << "];\n";

  bool savedVectorPrintState = m_chain[m_chain.size()-1]->getPrintHorizontally();
  m_chain[m_chain.size()-1]->setPrintHorizontally(false);
  ofs << "queso_" << m_prefix << "results.theta = ["
      << *(m_chain[m_chain.size()-1])
      << "];\n";
  m_chain[m_chain.size()-1]->setPrintHorizontally(savedVectorPrintState);
  
  ofs << "queso_" << m_prefix << "results.nbatch = 1.;\n"; // FIXME

  if (mahalanobisMatrix != NULL) {
    // Write covMatrix
    ofs << "queso_" << m_prefix << "mahalanobisMatrix = ["
        << *mahalanobisMatrix
        << "];\n";
  }

  // Write number of rejections
  ofs << "queso_" << m_prefix << "results.rejected = "  << (double) m_numRejections/(double) m_chain.size()
      << ";\n"
      << std::endl;

  // Write number of rejections
  ofs << "queso_" << m_prefix << "results.simutime = "  << m_chainRunTime
      << ";\n"
      << std::endl;

  // Write number of outbounds
  ofs << "queso_" << m_prefix << "results.ulrejected = " << (double) m_numOutOfBounds/(double) m_chain.size()
      << ";\n"
      << std::endl;

  // Write histograms
  if ((m_computeHistograms               ) &&
      (m_numberOfInternalBinsForHists > 0)) {
    // plot(queso_centersOfHistBins(1,:)',queso_histBins(1,:)','r-');
    ofs << "queso_" << m_prefix << "centersOfHistBins = zeros(" << m_paramSpace.dim()
        << ","                                                  << m_centersForAllHistogramBins.size()
        << ");"
        << std::endl;
    for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
      for (unsigned int j = 0; j < m_centersForAllHistogramBins.size(); ++j) {
         ofs << "queso_" << m_prefix << "centersOfHistBins(" << i+1
             << ","                                          << j+1
             << ") = "                                       << (*(m_centersForAllHistogramBins[j]))[i]
             << ";"
             << std::endl;
      }
    }

    ofs << "queso_" << m_prefix << "histBins = zeros(" << m_paramSpace.dim()
        << ","                       << m_histogramBinsForAllParams.size()
        << ");"
        << std::endl;
    for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
      for (unsigned int j = 0; j < m_histogramBinsForAllParams.size(); ++j) {
         ofs << "queso_" << m_prefix << "histBins(" << i+1
             << ","                                 << j+1
             << ") = "                              << (*(m_histogramBinsForAllParams[j]))[i]
             << ";"
             << std::endl;
      }
    }
  }

  // Write estimations of probability densities
  if ((m_computeKDEs                     ) &&
      (m_numberOfEvaluationPosForKDEs > 0)) {
    // hold
    // plot(queso_evalPosForKDE(1,:)',7*queso_densFromGaussianKDE(1,:)','r-');
    ofs << "queso_" << m_prefix << "evalPosForKDE = zeros(" << m_paramSpace.dim()
        << ","                                              << m_evaluationPositionsForKDEs.size()
        << ");"
        << std::endl;
    for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
      for (unsigned int j = 0; j < m_evaluationPositionsForKDEs.size(); ++j) {
        ofs << "queso_" << m_prefix << "evalPosForKDE(" << i+1
            << ","                                      << j+1
            << ") = "                                   << (*(m_evaluationPositionsForKDEs[j]))[i]
            << ";"
            << std::endl;
      }
    }

    ofs << "queso_" << m_prefix << "scalesForKDE = zeros(" << m_paramSpace.dim()
        << ");"
        << std::endl;
    for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
      ofs << "queso_" << m_prefix << "scalesForKDE(" << i+1
          << ") = "                                  << (*m_scalesForKDEs)[i]
          << ";"
          << std::endl;
    }

    ofs << "queso_" << m_prefix << "densFromGaussianKDE = zeros(" << m_paramSpace.dim()
        << ","                                                    << m_densityValuesFromGaussianKDE.size()
        << ");"
        << std::endl;
    for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
      for (unsigned int j = 0; j < m_densityValuesFromGaussianKDE.size(); ++j) {
        ofs << "queso_" << m_prefix << "densFromGaussianKDE(" << i+1
            << ","                                            << j+1
            << ") = "                                         << (*(m_densityValuesFromGaussianKDE[j]))[i]
            << ";"
            << std::endl;
      }
    }
  }

  // Close file
  ofs.close();

  if (m_env.rank() == 0) {
    std::cout << "Finished writing out info about the chain of id " << chainId
              << std::endl;
  }

  return iRC;
}

template <class V, class M>
const std::vector<const V*>&
uqDRAM_MarkovChainGeneratorClass<V,M>::chain() const
{
  return m_chain;
}

template <class V, class M>
const std::vector<const V*>&
uqDRAM_MarkovChainGeneratorClass<V,M>::lrSigma2Chain() const
{
  return m_lrSigma2Chain;
}

template <class V, class M>
const std::string&
uqDRAM_MarkovChainGeneratorClass<V,M>::outputFileName() const
{
  return m_namesOfOutputFiles[m_chain.size()-1];
}

template <class V, class M>
void
uqDRAM_MarkovChainGeneratorClass<V,M>::print(std::ostream& os) const
{
  os << m_prefix << "sizesOfChains = ";
  for (unsigned int i = 0; i < m_sizesOfChains.size(); ++i) {
    os << m_sizesOfChains[i] << " ";
  }
  os << "\n" << m_prefix << "lrUpdateSigma2 = " << m_lrUpdateSigma2
     << "\n" << m_prefix << "lrSigma2Priors = ";
  for (unsigned int i = 0; i < m_lrSigma2Priors->size(); ++i) {
    os << (*m_lrSigma2Priors)[i] << " ";
  }
  os << "\n" << m_prefix << "lrSigma2Accuracies = ";
  for (unsigned int i = 0; i < m_lrSigma2Accuracies->size(); ++i) {
    os << (*m_lrSigma2Accuracies)[i] << " ";
  }
  os << "\n" << m_prefix << "lrNumbersOfObservations = ";
  for (unsigned int i = 0; i < m_lrNumbersOfObservations->size(); ++i) {
    os << (*m_lrNumbersOfObservations)[i] << " ";
  }
  os << "\n" << m_prefix << "maxNumberOfExtraStages = " << m_maxNumberOfExtraStages
     << "\n" << m_prefix << "scalesForCovMProposals = ";
  for (unsigned int i = 0; i < m_scalesForCovMProposals.size(); ++i) {
    os << m_scalesForCovMProposals[i] << " ";
  }
  os << "\n" << m_prefix << "initialNonAdaptInterval = " << m_initialNonAdaptInterval
     << "\n" << m_prefix << "adaptInterval = "           << m_adaptInterval
     << "\n" << m_prefix << "eta = "                     << m_eta
     << "\n" << m_prefix << "epsilon = "                 << m_epsilon;
  os << "\n" << m_prefix << "namesOfOutputFiles = ";
  for (unsigned int i = 0; i < m_namesOfOutputFiles.size(); ++i) {
    os << m_namesOfOutputFiles[i] << " ";
  }
  os << "\n" << m_prefix << "chainDisplayPeriod = " << m_chainDisplayPeriod;
  os << std::endl;

  return;
}

template <class V, class M>
std::ostream& operator<<(std::ostream& os, const uqDRAM_MarkovChainGeneratorClass<V,M>& obj)
{
  obj.print(os);

  return os;
}
// Possible comparisons between Mcmc++ and Mcmc matlab toolbox

// semilogy([1:1000],alphaQuotientsCpp,'*');         
// hold                                      
// semilogy([1:1000],results.alphaQuotients,'r.');

// mean(alphaQuotientsCpp)
// mean(results.alphaQuotients)

// semilogy([1:1000],s2chainCpp,'.');    
// hold                              
// semilogy([1:1000],s2chain,'r.');    

// mean(s2chainCpp)
// mean(s2chain)   

// [zCpp,pCpp] = geweke(chainCpp)
// [z,p] = geweke(chain)

// tauCpp = iact(chainCpp)
// tau = iact(chain)

#endif // __UQ_DRAM_MCG_H__
