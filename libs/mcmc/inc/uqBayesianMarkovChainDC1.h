/* uq/libs/mcmc/inc/uqBayesianMarkovChainDC1.h
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

#ifndef __UQ_BMC_DC1_H__
#define __UQ_BMC_DC1_H__

#undef  UQ_BMC_DC_REQUIRES_INVERTED_COV_MATRICES

#define UQ_BMC_DC_MARKOV_CHAIN_TYPE           1
#define UQ_BMC_DC_WHITE_NOISE_CHAIN_TYPE      2
#define UQ_BMC_DC_UNIFORM_CHAIN_TYPE          3
#define UQ_BMC_DC_FILENAME_FOR_NO_OUTPUT_FILE "."

// _ODV = option default value
#define UQ_BMC_DC_CHAIN_TYPE_ODV                       UQ_BMC_DC_MARKOV_CHAIN_TYPE
#define UQ_BMC_DC_CHAIN_NUMBER_ODV                     1
#define UQ_BMC_DC_CHAIN_SIZES_ODV                      "100"
#define UQ_BMC_DC_CHAIN_USE2_ODV                       0
#define UQ_BMC_DC_CHAIN_GENERATE_EXTRA_ODV             0
#define UQ_BMC_DC_CHAIN_DISPLAY_PERIOD_ODV             500
#define UQ_BMC_DC_CHAIN_MEASURE_RUN_TIMES_ODV          0
#define UQ_BMC_DC_CHAIN_WRITE_ODV                      0
#define UQ_BMC_DC_CHAIN_COMPUTE_STATS_ODV              0
#define UQ_BMC_DC_CHAIN_OUTPUT_FILE_NAMES_ODV          UQ_BMC_DC_FILENAME_FOR_NO_OUTPUT_FILE
#define UQ_BMC_DC_UNIQUE_CHAIN_GENERATE_ODV            0
#define UQ_BMC_DC_UNIQUE_CHAIN_WRITE_ODV               0
#define UQ_BMC_DC_UNIQUE_CHAIN_COMPUTE_STATS_ODV       0
#define UQ_BMC_DC_FILTERED_CHAIN_GENERATE_ODV          0
#define UQ_BMC_DC_FILTERED_CHAIN_DISCARDED_PORTION_ODV 0.
#define UQ_BMC_DC_FILTERED_CHAIN_LAG_ODV               1
#define UQ_BMC_DC_FILTERED_CHAIN_WRITE_ODV             0
#define UQ_BMC_DC_FILTERED_CHAIN_COMPUTE_STATS_ODV     0
#define UQ_BMC_DC_AVG_CHAIN_COMPUTE_ODV                "0"
#define UQ_BMC_DC_AVG_CHAIN_WRITE_ODV                  0
#define UQ_BMC_DC_AVG_CHAIN_COMPUTE_STATS_ODV          0
#define UQ_BMC_DC_DR_MAX_NUM_EXTRA_STAGES_ODV          0
#define UQ_BMC_DC_DR_SCALES_FOR_EXTRA_STAGES_ODV       "1."
#define UQ_BMC_DC_AM_INIT_NON_ADAPT_INT_ODV            0
#define UQ_BMC_DC_AM_ADAPT_INTERVAL_ODV                0
#define UQ_BMC_DC_AM_ETA_ODV                           1.
#define UQ_BMC_DC_AM_EPSILON_ODV                       1.e-5

#include <uqChainStatisticalOptions.h>
#include <uqProbDensity.h>
#include <uqLikelihoodFunction.h>
#include <uqParamSpace.h>
#include <uqObservableSpace.h>
#include <uqChainPosition.h>
#include <uqMiscellaneous.h>
#include <uqSequenceOfVectors.h>
#include <uqArrayOfSequences.h>
#include <uq2dArrayOfStuff.h>
#include <sys/time.h>
#include <fstream>

/*! A templated class that implements a Bayesian Markov Chain Distribution Calculator
 */
template <class P_V,class P_M,class L_V,class L_M>
class uqBayesianMarkovChainDCClass
{
public:
  uqBayesianMarkovChainDCClass(const uqEnvironmentClass&                              env,                        /*! The QUESO toolkit environment.   */
                               const char*                                            prefix,                     /*! Prefix for the validation phase. */
                               const uqParamSpaceClass<P_V,P_M>&                      paramSpace,                 /*! The parameter space.             */
                               const uqObservableSpaceClass<L_V,L_M>&                 observableSpace,            /*! The observable space.            */
                               const uqProbDensity_BaseClass<P_V,P_M>&                m2lPriorProbDensity_Obj,    /*! -2*ln(prior())                   */
                               const uqLikelihoodFunction_BaseClass<P_V,P_M,L_V,L_M>& m2lLikelihoodFunction_Obj); /*! -2*ln(likelihood())              */
 ~uqBayesianMarkovChainDCClass();

  void calculatePosterior      (const P_M* proposalCovMatrix,
                              //const P_M* proposalPrecMatrix,
                                const P_M* mahalanobisMatrix = NULL,
                                bool       applyMahalanobisInvert = true);

  void print                   (std::ostream& os) const;

//const uqSequenceOfVectorsClass<P_V>& chain              () const;
//const std::vector<const L_V*>&       misfitChain        () const;
//const std::vector<const L_V*>&       misfitVarianceChain() const;
//const std::string&                   outputFileName     () const;

private:
  void   resetChainAndRelatedInfo();
  void   defineMyOptions         (po::options_description&     optionsDesc);
  void   getMyOptionValues       (po::options_description&     optionsDesc);

  int    prepareForNextChain     (const P_M*                   proposalCovMatrix);
                                //const P_M*                   proposalPrecMatrix,

  void   calculatePosterior      (const P_M*                   proposalCovMatrix,
                                //const P_M*                   proposalPrecMatrix,
                                  const P_M*                   mahalanobisMatrix,
                                  bool                         applyMahalanobisInvert,
                                  uqChainBaseClass<P_V>&       workingChain);
  int    generateChain           (unsigned int                 chainSize,
                                  const P_V&                   valuesOf1stPosition,
                                  const P_M*                   proposalCovMatrix,
                                  const P_M*                   mahalanobisMatrix,
                                  bool                         applyMahalanobisInvert,
                                  uqChainBaseClass<P_V>&       workingChain,
                                  const std::string&           chainName,
                                  const std::string&           prefixName,
                                  std::ofstream*               ofs);
  int    generateWhiteNoise      (unsigned int                 chainSize,
                                  uqChainBaseClass<P_V>&       workingChain,
                                  const std::string&           chainName);
  int    generateUniform         (unsigned int                 chainSize,
                                  uqChainBaseClass<P_V>&       workingChain,
                                  const std::string&           chainName);
  void   updateCovMatrix         (const uqChainBaseClass<P_V>& subChain,
                                  unsigned int                 idOfFirstPositionInSubChain,
                                  double&                      lastChainSize,
                                  P_V&                         lastMean,
                                  P_M&                         lastAdaptedCovMatrix);

  double logProposal             (const uqChainPositionClass<P_V>&               x,
                                  const uqChainPositionClass<P_V>&               y,
                                  unsigned int                                   idOfProposalCovMatrix);
  double logProposal             (const std::vector<uqChainPositionClass<P_V>*>& inputPositions);
  double alpha                   (const uqChainPositionClass<P_V>&               x,
                                  const uqChainPositionClass<P_V>&               y,
                                  double*                                        alphaQuotientPtr = NULL);
  double alpha                   (const std::vector<uqChainPositionClass<P_V>*>& inputPositions);
  bool   acceptAlpha             (double                                         alpha);
  void   updateCovMatrices       ();
  void   computeStatistics       (const uqChainBaseClass<P_V>&          workingChain,
                                  const uqChainStatisticalOptionsClass& statisticalOptions,
                                  const std::string&                    chainName,
                                  std::ofstream*                        passedOfs);
  void   computeMeanVars         (const uqChainBaseClass<P_V>&          workingChain,
                                  const uqChainStatisticalOptionsClass& statisticalOptions,
                                  const std::string&                    chainName,
                                  std::ofstream*                        passedOfs,
                                  P_V*                                  meanPtr,
                                  P_V*                                  sampleVarPtr,
                                  P_V*                                  populVarPtr);
  void   computeBMM              (const uqChainBaseClass<P_V>&          workingChain,
                                  const uqChainStatisticalOptionsClass& statisticalOptions,
                                  const std::vector<unsigned int>&      initialPosForStatistics,
                                  const std::string&                    chainName,
                                  std::ofstream*                        passedOfs);
  void   computeFFT              (const uqChainBaseClass<P_V>&          workingChain,
                                  const uqChainStatisticalOptionsClass& statisticalOptions,
                                  const std::vector<unsigned int>&      initialPosForStatistics,
                                  const std::string&                    chainName,
                                  std::ofstream*                        passedOfs);
  void   computePSD              (const uqChainBaseClass<P_V>&          workingChain,
                                  const uqChainStatisticalOptionsClass& statisticalOptions,
                                  const std::vector<unsigned int>&      initialPosForStatistics,
                                  const std::string&                    chainName,
                                  std::ofstream*                        passedOfs);
  void   computePSDAtZero        (const uqChainBaseClass<P_V>&          workingChain,
                                  const uqChainStatisticalOptionsClass& statisticalOptions,
                                  const std::vector<unsigned int>&      initialPosForStatistics,
                                  const std::string&                    chainName,
                                  std::ofstream*                        passedOfs);
  void   computeGeweke           (const uqChainBaseClass<P_V>&          workingChain,
                                  const uqChainStatisticalOptionsClass& statisticalOptions,
                                  const std::vector<unsigned int>&      initialPosForStatistics,
                                  const std::string&                    chainName,
                                  std::ofstream*                        passedOfs);
  void   computeCorrViaDef       (const uqChainBaseClass<P_V>&          workingChain,
                                  const uqChainStatisticalOptionsClass& statisticalOptions,
                                  const std::vector<unsigned int>&      initialPosForStatistics,
                                  const std::vector<unsigned int>&      lagsForCorrs,
                                  const std::string&                    chainName,
                                  std::ofstream*                        passedOfs);
  void   computeCorrViaFFT       (const uqChainBaseClass<P_V>&          workingChain,
                                  const uqChainStatisticalOptionsClass& statisticalOptions,
                                  const std::vector<unsigned int>&      initialPosForStatistics,
                                  const std::vector<unsigned int>&      lagsForCorrs,
                                  const std::string&                    chainName,
                                  std::ofstream*                        passedOfs);
  void   computeFilterParameters (const uqChainBaseClass<P_V>&          workingChain,
                                  const uqChainStatisticalOptionsClass& statisticalOptions,
                                  const std::string&                    chainName,
                                  std::ofstream*                        passedOfs,
                                  unsigned int&                         initialPos,
                                  unsigned int&                         spacing);
  void   computeHistKde          (const uqChainBaseClass<P_V>&          workingChain,
                                  const uqChainStatisticalOptionsClass& statisticalOptions,
                                  const std::string&                    chainName,
                                  std::ofstream*                        passedOfs);

  int    writeInfo               (const uqChainBaseClass<P_V>&          workingChain,
                                  const std::string&                    chainName,
                                  const std::string&                    prefixName,
                                  std::ofstream&                        ofs,
                                  const P_M*                            mahalanobisMatrix = NULL,
                                  bool                                  applyMahalanobisInvert = true) const;

  const uqEnvironmentClass&                              m_env;
        std::string                                      m_prefix;
  const uqParamSpaceClass<P_V,P_M>&                      m_paramSpace;
  const uqObservableSpaceClass<L_V,L_M>&                 m_observableSpace;
  const uqProbDensity_BaseClass<P_V,P_M>&                m_m2lPriorProbDensity_Obj;
  const uqLikelihoodFunction_BaseClass<P_V,P_M,L_V,L_M>& m_m2lLikelihoodFunction_Obj;

  std::string m_option_help;
  std::string m_option_chain_type;
  std::string m_option_chain_number;
  std::string m_option_chain_sizes;
  std::string m_option_chain_use2;
  std::string m_option_chain_generateExtra;
  std::string m_option_chain_displayPeriod;
  std::string m_option_chain_measureRunTimes;
  std::string m_option_chain_write;
  std::string m_option_chain_computeStats;
  std::string m_option_uniqueChain_generate;
  std::string m_option_uniqueChain_write;
  std::string m_option_uniqueChain_computeStats;
  std::string m_option_filteredChain_generate;
  std::string m_option_filteredChain_discardedPortion;
  std::string m_option_filteredChain_lag;
  std::string m_option_filteredChain_write;
  std::string m_option_filteredChain_computeStats;
  std::string m_option_avgChain_compute;
  std::string m_option_avgChain_write;
  std::string m_option_avgChain_computeStats;
  std::string m_option_dr_maxNumExtraStages;
  std::string m_option_dr_scalesForExtraStages;
  std::string m_option_am_initialNonAdaptInterval;
  std::string m_option_am_adaptInterval;
  std::string m_option_am_eta;
  std::string m_option_am_epsilon;
  std::string m_option_chain_outputFileNames;

  bool                            m_likelihoodObjComputesMisfits;
  P_V                             m_paramInitials;
  bool                            m_proposalIsSymmetric;
  po::options_description*        m_optionsDesc;

  unsigned int                    m_chainType;
  unsigned int                    m_chainNumber;
  std::vector<unsigned int>       m_chainSizes;
  bool                            m_chainUse2;
  bool                            m_chainGenerateExtra;
  unsigned int                    m_chainDisplayPeriod;
  bool                            m_chainMeasureRunTimes;
  bool                            m_chainWrite;
  bool                            m_chainComputeStats;
  uqChainStatisticalOptionsClass* m_chainStatisticalOptions;
  std::vector<std::string>        m_chainOutputFileNames;

  bool                            m_uniqueChainGenerate;
  bool                            m_uniqueChainWrite;
  bool                            m_uniqueChainComputeStats;
  uqChainStatisticalOptionsClass* m_uniqueChainStatisticalOptions;

  bool                            m_filteredChainGenerate;
  double                          m_filteredChainDiscardedPortion; // input or set during run time
  unsigned int                    m_filteredChainLag;              // input or set during run time
  bool                            m_filteredChainWrite;
  bool                            m_filteredChainComputeStats;
  uqChainStatisticalOptionsClass* m_filteredChainStatisticalOptions;

  std::vector<unsigned int>       m_avgChainCompute;
  bool                            m_avgChainWrite;
  bool                            m_avgChainComputeStats;

  unsigned int                    m_maxNumExtraStages;
  std::vector<double>             m_scalesForCovMProposals;
  unsigned int                    m_initialNonAdaptInterval;
  unsigned int                    m_adaptInterval;
  double                          m_eta;
  double                          m_epsilon;

  std::vector<      P_M*>         m_lowerCholProposalCovMatrices;
  std::vector<      P_M*>         m_proposalCovMatrices;
#ifdef UQ_BMC_DC_REQUIRES_INVERTED_COV_MATRICES
  std::vector<      P_M*>         m_upperCholProposalPrecMatrices;
  std::vector<      P_M*>         m_proposalPrecMatrices;
#endif

  uqSequenceOfVectorsClass<P_V>   m_chain1;
  uqArrayOfSequencesClass<P_V>    m_chain2;
  std::vector<unsigned int>       m_idsOfUniquePositions;
  std::vector<const L_V*>         m_misfitChain; // Sum of squares of differences between model and experiments: computed by user supplied likelihood obj
  std::vector<const L_V*>         m_misfitVarianceChain;
  std::vector<const L_V*>         m_m2lLikelihoodChain;
  std::vector<double>             m_alphaQuotients;
  double                          m_chainRunTime;
  double                          m_candidateRunTime;
  double                          m_priorRunTime;
  double                          m_lhRunTime;
  double                          m_mhAlphaRunTime;
  double                          m_drAlphaRunTime;
  double                          m_drRunTime;
  double                          m_amRunTime;
  unsigned int                    m_numRejections;
  unsigned int                    m_numOutOfBounds;
  double                          m_lastChainSize;
  P_V*                            m_lastMean;
  P_M*                            m_lastAdaptedCovMatrix;
};

template<class P_V,class P_M,class L_V,class L_M>
std::ostream& operator<<(std::ostream& os, const uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>& obj);

#include <uqBayesianMarkovChainDC2.h>

template<class P_V,class P_M,class L_V,class L_M>
uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::uqBayesianMarkovChainDCClass(
  const uqEnvironmentClass&                              env,
  const char*                                            prefix,
  const uqParamSpaceClass<P_V,P_M>&                      paramSpace,
  const uqObservableSpaceClass<L_V,L_M>&                 observableSpace,
  const uqProbDensity_BaseClass<P_V,P_M>&                m2lPriorProbDensity_Obj,
  const uqLikelihoodFunction_BaseClass<P_V,P_M,L_V,L_M>& m2lLikelihoodFunction_Obj)
  :
  m_env                            (env),
  m_prefix                         (prefix),
  m_paramSpace                     (paramSpace),
  m_observableSpace                (observableSpace),
  m_m2lPriorProbDensity_Obj        (m2lPriorProbDensity_Obj),
  m_m2lLikelihoodFunction_Obj      (m2lLikelihoodFunction_Obj),
  m_likelihoodObjComputesMisfits   (dynamic_cast<const uqMisfitLikelihoodFunction_Class<P_V,P_M,L_V,L_M>*>(&m2lLikelihoodFunction_Obj) != NULL),
  m_paramInitials                  (m_paramSpace.initialValues()),
  m_proposalIsSymmetric            (true),
  m_optionsDesc                    (new po::options_description("Markov chain Monte Carlo options")),
  m_chainType                      (UQ_BMC_DC_CHAIN_TYPE_ODV),
  m_chainNumber                    (UQ_BMC_DC_CHAIN_NUMBER_ODV),
  m_chainSizes                     (1,(unsigned int) strtod(UQ_BMC_DC_CHAIN_SIZES_ODV,NULL)),
  m_chainUse2                      (UQ_BMC_DC_CHAIN_USE2_ODV),
  m_chainGenerateExtra             (UQ_BMC_DC_CHAIN_GENERATE_EXTRA_ODV),
  m_chainDisplayPeriod             (UQ_BMC_DC_CHAIN_DISPLAY_PERIOD_ODV),
  m_chainMeasureRunTimes           (UQ_BMC_DC_CHAIN_MEASURE_RUN_TIMES_ODV),
  m_chainWrite                     (UQ_BMC_DC_CHAIN_WRITE_ODV),
  m_chainComputeStats              (UQ_BMC_DC_CHAIN_COMPUTE_STATS_ODV),
  m_chainStatisticalOptions        (NULL),
  m_chainOutputFileNames           (1,UQ_BMC_DC_CHAIN_OUTPUT_FILE_NAMES_ODV),
  m_uniqueChainGenerate            (UQ_BMC_DC_UNIQUE_CHAIN_GENERATE_ODV),
  m_uniqueChainWrite               (UQ_BMC_DC_UNIQUE_CHAIN_WRITE_ODV),
  m_uniqueChainComputeStats        (UQ_BMC_DC_UNIQUE_CHAIN_COMPUTE_STATS_ODV),
  m_uniqueChainStatisticalOptions  (NULL),
  m_filteredChainGenerate          (UQ_BMC_DC_FILTERED_CHAIN_GENERATE_ODV),
  m_filteredChainDiscardedPortion  (UQ_BMC_DC_FILTERED_CHAIN_DISCARDED_PORTION_ODV),
  m_filteredChainLag               (UQ_BMC_DC_FILTERED_CHAIN_LAG_ODV),
  m_filteredChainWrite             (UQ_BMC_DC_FILTERED_CHAIN_WRITE_ODV),
  m_filteredChainComputeStats      (UQ_BMC_DC_FILTERED_CHAIN_COMPUTE_STATS_ODV),
  m_filteredChainStatisticalOptions(NULL),
  m_avgChainCompute                (0),//,0.),
  m_avgChainWrite                  (UQ_BMC_DC_AVG_CHAIN_WRITE_ODV),
  m_avgChainComputeStats           (UQ_BMC_DC_AVG_CHAIN_COMPUTE_STATS_ODV),
  m_maxNumExtraStages              (UQ_BMC_DC_DR_MAX_NUM_EXTRA_STAGES_ODV),
  m_scalesForCovMProposals         (0),//,0.),
  m_initialNonAdaptInterval        (UQ_BMC_DC_AM_INIT_NON_ADAPT_INT_ODV),
  m_adaptInterval                  (UQ_BMC_DC_AM_ADAPT_INTERVAL_ODV),
  m_eta                            (UQ_BMC_DC_AM_ETA_ODV),
  m_epsilon                        (UQ_BMC_DC_AM_EPSILON_ODV),
  m_lowerCholProposalCovMatrices   (1),//,NULL),
  m_proposalCovMatrices            (1),//,NULL),
#ifdef UQ_BMC_DC_REQUIRES_INVERTED_COV_MATRICES
  m_upperCholProposalPrecMatrices  (1),//,NULL),
  m_proposalPrecMatrices           (1),//,NULL),
#endif
  m_chain1                         (0,m_paramSpace.zeroVector()),
  m_chain2                         (0,m_paramSpace.zeroVector()),
  m_idsOfUniquePositions           (0),//0),
  m_misfitChain                    (0),//,NULL),
  m_misfitVarianceChain            (0),//,NULL),
  m_m2lLikelihoodChain             (0),//,NULL),
  m_alphaQuotients                 (0),//,0.),
  m_chainRunTime                   (0.),
  m_candidateRunTime               (0.),
  m_priorRunTime                   (0.),
  m_lhRunTime                      (0.),
  m_mhAlphaRunTime                 (0.),
  m_drAlphaRunTime                 (0.),
  m_drRunTime                      (0.),
  m_amRunTime                      (0.), 
  m_numRejections                  (0),
  m_numOutOfBounds                 (0),
  m_lastChainSize                  (0),
  m_lastMean                       (NULL),
  m_lastAdaptedCovMatrix           (NULL)
{
  if (m_env.rank() == 0) std::cout << "Entering uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::constructor()"
                                   << std::endl;

  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.rank() == 0) std::cout << "In uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::constructor()"
                                   << ": after getting values of options with prefix '" << m_prefix
                                   << "', state of  object is:"
                                   << "\n" << *this
                                   << std::endl;

  if (m_chainComputeStats        ) m_chainStatisticalOptions         = new uqChainStatisticalOptionsClass(m_env,m_prefix+"BMC_DC_chain_"        );
  if (m_uniqueChainComputeStats  ) m_uniqueChainStatisticalOptions   = new uqChainStatisticalOptionsClass(m_env,m_prefix+"BMC_DC_uniqueChain_"  );
  if (m_filteredChainComputeStats) m_filteredChainStatisticalOptions = new uqChainStatisticalOptionsClass(m_env,m_prefix+"BMC_DC_filteredChain_");

  if (m_env.rank() == 0) std::cout << "Leaving uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::constructor()"
                                   << std::endl;
}

template<class P_V,class P_M,class L_V,class L_M>
uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::~uqBayesianMarkovChainDCClass()
{
  //std::cout << "Entering uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::destructor()"
  //          << std::endl;

  resetChainAndRelatedInfo();
  //if (m_lowerCholProposalCovMatrices[0]) delete m_lowerCholProposalCovMatrices[0]; // Loop inside 'resetChainAndRelatedInfo()' deletes just from position '1' on

  if (m_filteredChainStatisticalOptions) delete m_filteredChainStatisticalOptions;
  if (m_uniqueChainStatisticalOptions  ) delete m_uniqueChainStatisticalOptions;
  if (m_chainStatisticalOptions        ) delete m_chainStatisticalOptions;
  if (m_optionsDesc                    ) delete m_optionsDesc;

  //std::cout << "Leaving uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::destructor()"
  //          << std::endl;
}

template<class P_V,class P_M,class L_V,class L_M>
void
uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::resetChainAndRelatedInfo()
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

  m_idsOfUniquePositions.clear();

#ifdef UQ_BMC_DC_REQUIRES_INVERTED_COV_MATRICES
  for (unsigned int i = 0; i < m_proposalPrecMatrices.size(); ++i) {
    if (m_proposalPrecMatrices[i]) delete m_proposalPrecMatrices[i];
  }
  for (unsigned int i = 0; i < m_upperCholProposalPrecMatrices.size(); ++i) {
    if (m_upperCholProposalPrecMatrices[i]) delete m_upperCholProposalPrecMatrices[i];
  }
#endif
  for (unsigned int i = 0; i < m_proposalCovMatrices.size(); ++i) {
    if (m_proposalCovMatrices[i]) {
      delete m_proposalCovMatrices[i];
      m_proposalCovMatrices[i] = NULL;
    }
  }
//m_proposalCovMatrices.clear(); // Do not clear
  for (unsigned int i = 0; i < m_lowerCholProposalCovMatrices.size(); ++i) { // Yes, from '1' on
    if (m_lowerCholProposalCovMatrices[i]) {
      delete m_lowerCholProposalCovMatrices[i];
      m_lowerCholProposalCovMatrices[i] = NULL;
    }
  }
//m_lowerCholProposalCovMatrices.clear(); // Do not clear

  return;
}

template<class P_V,class P_M,class L_V,class L_M>
void
uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::defineMyOptions(
  po::options_description& optionsDesc)
{
  m_option_help                           = m_prefix + "BMC_DC_help";

  m_option_chain_type                     = m_prefix + "BMC_DC_chain_type";
  m_option_chain_number                   = m_prefix + "BMC_DC_chain_number";
  m_option_chain_sizes                    = m_prefix + "BMC_DC_chain_sizes";
  m_option_chain_use2                     = m_prefix + "BMC_DC_chain_use2";
  m_option_chain_generateExtra            = m_prefix + "BMC_DC_chain_generateExtra";
  m_option_chain_displayPeriod            = m_prefix + "BMC_DC_chain_displayPeriod";
  m_option_chain_measureRunTimes          = m_prefix + "BMC_DC_chain_measureRunTimes";
  m_option_chain_write                    = m_prefix + "BMC_DC_chain_write";
  m_option_chain_computeStats             = m_prefix + "BMC_DC_chain_computeStats";
  m_option_chain_outputFileNames          = m_prefix + "BMC_DC_chain_outputFileNames";

  m_option_uniqueChain_generate           = m_prefix + "BMC_DC_uniqueChain_generate";
  m_option_uniqueChain_write              = m_prefix + "BMC_DC_uniqueChain_write";
  m_option_uniqueChain_computeStats       = m_prefix + "BMC_DC_uniqueChain_computeStats";

  m_option_filteredChain_generate         = m_prefix + "BMC_DC_filteredChain_generate";
  m_option_filteredChain_discardedPortion = m_prefix + "BMC_DC_filteredChain_discardedPortion";
  m_option_filteredChain_lag              = m_prefix + "BMC_DC_filteredChain_lag";
  m_option_filteredChain_write            = m_prefix + "BMC_DC_filteredChain_write";
  m_option_filteredChain_computeStats     = m_prefix + "BMC_DC_filteredChain_computeStats";

  m_option_avgChain_compute               = m_prefix + "BMC_DC_avgChain_compute";
  m_option_avgChain_write                 = m_prefix + "BMC_DC_avgChain_write";
  m_option_avgChain_computeStats          = m_prefix + "BMC_DC_avgChain_computeStats";

  m_option_dr_maxNumExtraStages           = m_prefix + "BMC_DC_dr_maxNumExtraStages";
  m_option_dr_scalesForExtraStages        = m_prefix + "BMC_DC_dr_scalesForExtraStages";

  m_option_am_initialNonAdaptInterval     = m_prefix + "BMC_DC_am_initialNonAdaptInterval";
  m_option_am_adaptInterval               = m_prefix + "BMC_DC_am_adaptInterval";
  m_option_am_eta                         = m_prefix + "BMC_DC_am_eta";
  m_option_am_epsilon                     = m_prefix + "BMC_DC_am_epsilon";

  optionsDesc.add_options()
    (m_option_help.c_str(),                                                                                                                     "produce help message for Bayesian Markov chain distr. calculator")
    (m_option_chain_type.c_str(),                     po::value<unsigned int>()->default_value(UQ_BMC_DC_CHAIN_TYPE_ODV                      ), "type of chain (1=Markov, 2=White noise)"                         )
    (m_option_chain_number.c_str(),                   po::value<unsigned int>()->default_value(UQ_BMC_DC_CHAIN_NUMBER_ODV                    ), "number of chain(s)"                                              )
    (m_option_chain_sizes.c_str(),                    po::value<std::string >()->default_value(UQ_BMC_DC_CHAIN_SIZES_ODV                     ), "list of size(s) of chain(s)"                                     )
    (m_option_chain_use2.c_str(),                     po::value<bool        >()->default_value(UQ_BMC_DC_CHAIN_USE2_ODV                      ), "use chain2"                                                      )
    (m_option_chain_generateExtra.c_str(),            po::value<bool        >()->default_value(UQ_BMC_DC_CHAIN_GENERATE_EXTRA_ODV            ), "generate extra chains"                                           )
    (m_option_chain_displayPeriod.c_str(),            po::value<unsigned int>()->default_value(UQ_BMC_DC_CHAIN_DISPLAY_PERIOD_ODV            ), "period of message display during chain generation"               )
    (m_option_chain_measureRunTimes.c_str(),          po::value<bool        >()->default_value(UQ_BMC_DC_CHAIN_MEASURE_RUN_TIMES_ODV         ), "measure run times"                                               )
    (m_option_chain_write.c_str(),                    po::value<bool        >()->default_value(UQ_BMC_DC_CHAIN_WRITE_ODV                     ), "write chain values to the output file"                           )
    (m_option_chain_computeStats.c_str(),             po::value<bool        >()->default_value(UQ_BMC_DC_CHAIN_COMPUTE_STATS_ODV             ), "compute statistics on chain"                                     )
    (m_option_chain_outputFileNames.c_str(),          po::value<std::string >()->default_value(UQ_BMC_DC_CHAIN_OUTPUT_FILE_NAMES_ODV         ), "list of name(s) of output file(s)"                               )
  //(m_option_uniqueChain_generate.c_str(),           po::value<bool        >()->default_value(UQ_BMC_DC_UNIQUE_CHAIN_GENERATE_ODV           ), "generate unique chain"                                           )
  //(m_option_uniqueChain_write.c_str(),              po::value<bool        >()->default_value(UQ_BMC_DC_UNIQUE_CHAIN_WRITE_ODV              ), "write unique chain"                                              )
  //(m_option_uniqueChain_computeStats.c_str(),       po::value<bool        >()->default_value(UQ_BMC_DC_UNIQUE_CHAIN_COMPUTE_STATS_ODV      ), "compute statistics on unique chain"                              )
    (m_option_filteredChain_generate.c_str(),         po::value<bool        >()->default_value(UQ_BMC_DC_FILTERED_CHAIN_GENERATE_ODV         ), "generate filtered chain"                                         )
    (m_option_filteredChain_discardedPortion.c_str(), po::value<double      >()->default_value(UQ_BMC_DC_FILTERED_CHAIN_DISCARDED_PORTION_ODV), "initial discarded portion for chain filtering"                   )
    (m_option_filteredChain_lag.c_str(),              po::value<unsigned int>()->default_value(UQ_BMC_DC_FILTERED_CHAIN_LAG_ODV              ), "spacing for chain filtering"                                     )
    (m_option_filteredChain_write.c_str(),            po::value<bool        >()->default_value(UQ_BMC_DC_FILTERED_CHAIN_WRITE_ODV            ), "write filtered chain"                                            )
    (m_option_filteredChain_computeStats.c_str(),     po::value<bool        >()->default_value(UQ_BMC_DC_FILTERED_CHAIN_COMPUTE_STATS_ODV    ), "compute statistics on filtered chain"                            )
  //(m_option_avgChain_compute.c_str(),               po::value<std::string >()->default_value(UQ_BMC_DC_AVG_CHAIN_COMPUTE_ODV               ), "list of amounts of chains involved in chain averages"            )
  //(m_option_avgChain_write.c_str(),                 po::value<bool        >()->default_value(UQ_BMC_DC_AVG_CHAIN_WRITE_ODV                 ), "write averages of chains"                                        )
  //(m_option_avgChain_computeStats.c_str(),          po::value<bool        >()->default_value(UQ_BMC_DC_AVG_CHAIN_COMPUTE_STATS_ODV         ), "compute statistics on the averages of chains"                    )
    (m_option_dr_maxNumExtraStages.c_str(),           po::value<unsigned int>()->default_value(UQ_BMC_DC_DR_MAX_NUM_EXTRA_STAGES_ODV         ), "'dr' maximum number of extra stages"                             )
    (m_option_dr_scalesForExtraStages.c_str(),        po::value<std::string >()->default_value(UQ_BMC_DC_DR_SCALES_FOR_EXTRA_STAGES_ODV      ), "'dr' list of scales for proposal cov matrices from 2nd stage on" )
    (m_option_am_initialNonAdaptInterval.c_str(),     po::value<unsigned int>()->default_value(UQ_BMC_DC_AM_INIT_NON_ADAPT_INT_ODV           ), "'am' initial non adaptation interval"                            )
    (m_option_am_adaptInterval.c_str(),               po::value<unsigned int>()->default_value(UQ_BMC_DC_AM_ADAPT_INTERVAL_ODV               ), "'am' adaptation interval"                                        )
    (m_option_am_eta.c_str(),                         po::value<double      >()->default_value(UQ_BMC_DC_AM_ETA_ODV                          ), "'am' eta"                                                        )
    (m_option_am_epsilon.c_str(),                     po::value<double      >()->default_value(UQ_BMC_DC_AM_EPSILON_ODV                      ), "'am' epsilon"                                                    )
  ;

  return;
}

template<class P_V,class P_M,class L_V,class L_M>
void
uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::getMyOptionValues(
  po::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count(m_option_help.c_str())) {
    std::cout << optionsDesc
              << std::endl;
  }

  if (m_env.allOptionsMap().count(m_option_chain_type.c_str())) {
    m_chainType = m_env.allOptionsMap()[m_option_chain_type.c_str()].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_chain_number.c_str())) {
    m_chainNumber = m_env.allOptionsMap()[m_option_chain_number.c_str()].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_chain_sizes.c_str())) {
    m_chainSizes.clear();
    std::string inputString = m_env.allOptionsMap()[m_option_chain_sizes.c_str()].as<std::string>();
    std::vector<double> inputDoubles(0,0.);
    uqMiscReadDoublesFromString(inputString,inputDoubles);

    if ((inputDoubles.size() == 1               ) && 
        (inputDoubles.size() != m_chainNumber)) {
      inputDoubles.resize(m_chainNumber,inputDoubles[0]);
    }

    m_chainSizes.resize(inputDoubles.size(),0);
    for (unsigned int i = 0; i < inputDoubles.size(); ++i) {
      m_chainSizes[i] = (unsigned int) inputDoubles[i];
    }
  }

  if (m_env.allOptionsMap().count(m_option_chain_use2.c_str())) {
    m_chainUse2 = m_env.allOptionsMap()[m_option_chain_use2.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_chain_displayPeriod.c_str())) {
    m_chainDisplayPeriod = m_env.allOptionsMap()[m_option_chain_displayPeriod.c_str()].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_chain_measureRunTimes.c_str())) {
    m_chainMeasureRunTimes = m_env.allOptionsMap()[m_option_chain_measureRunTimes.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_chain_write.c_str())) {
    m_chainWrite = m_env.allOptionsMap()[m_option_chain_write.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_chain_computeStats.c_str())) {
    m_chainComputeStats = m_env.allOptionsMap()[m_option_chain_computeStats.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_chain_generateExtra.c_str())) {
    m_chainGenerateExtra = m_env.allOptionsMap()[m_option_chain_generateExtra.c_str()].as<bool>();
  }

  //if (m_env.allOptionsMap().count(m_option_uniqueChain_generate.c_str())) {
  //  m_uniqueChainGenerate = m_env.allOptionsMap()[m_option_uniqueChain_generate.c_str()].as<bool>();
  //}

  //if (m_env.allOptionsMap().count(m_option_uniqueChain_write.c_str())) {
  //  m_uniqueChainWrite = m_env.allOptionsMap()[m_option_uniqueChain_write.c_str()].as<bool>();
  //}

  //if (m_env.allOptionsMap().count(m_option_uniqueChain_computeStats.c_str())) {
  //  m_uniqueChainComputeStats = m_env.allOptionsMap()[m_option_uniqueChain_computeStats.c_str()].as<bool>();
  //}

  if (m_env.allOptionsMap().count(m_option_filteredChain_generate.c_str())) {
    m_filteredChainGenerate = m_env.allOptionsMap()[m_option_filteredChain_generate.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_filteredChain_discardedPortion.c_str())) {
    m_filteredChainDiscardedPortion = m_env.allOptionsMap()[m_option_filteredChain_discardedPortion.c_str()].as<double>();
  }

  if (m_env.allOptionsMap().count(m_option_filteredChain_lag.c_str())) {
    m_filteredChainLag = m_env.allOptionsMap()[m_option_filteredChain_lag.c_str()].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_filteredChain_write.c_str())) {
    m_filteredChainWrite = m_env.allOptionsMap()[m_option_filteredChain_write.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_filteredChain_computeStats.c_str())) {
    m_filteredChainComputeStats = m_env.allOptionsMap()[m_option_filteredChain_computeStats.c_str()].as<bool>();
  }
#if 0
  if (m_env.allOptionsMap().count(m_option_avgChain_compute.c_str())) {
    m_avgChainCompute.clear();
    std::vector<double> tmpDenominators(0,0.);
    std::string inputString = m_env.allOptionsMap()[m_option_avgChain_compute.c_str()].as<std::string>();
    uqMiscReadDoublesFromString(inputString,tmpDenominators);

    if (tmpDenominators.size() > 0) {
      m_avgChainCompute.resize(tmpDenominators.size(),0);
      for (unsigned int i = 0; i < m_avgChainCompute.size(); ++i) {
        m_avgChainCompute[i] = (unsigned int) tmpDenominators[i];
      }
    }
  }

  if (m_env.allOptionsMap().count(m_option_avgChain_write.c_str())) {
    m_avgChainWrite = m_env.allOptionsMap()[m_option_avgChain_write.c_str()].as<bool>();
  }

  if (m_env.allOptionsMap().count(m_option_avgChain_computeStats.c_str())) {
    m_avgChainComputeStats = m_env.allOptionsMap()[m_option_avgChain_computeStats.c_str()].as<bool>();
  }
#endif
  if (m_env.allOptionsMap().count(m_option_chain_outputFileNames.c_str())) {
    m_chainOutputFileNames.clear();
    std::string inputString = m_env.allOptionsMap()[m_option_chain_outputFileNames.c_str()].as<std::string>();
    uqMiscReadWordsFromString(inputString,m_chainOutputFileNames);

    if ((m_chainOutputFileNames.size() == 1                     ) && 
        (m_chainOutputFileNames.size() != m_chainSizes.size())) {
      m_chainOutputFileNames.resize(m_chainSizes.size(),m_chainOutputFileNames[0]);
    }

    UQ_FATAL_TEST_MACRO(m_chainOutputFileNames.size() != m_chainSizes.size(),
                        m_env.rank(),
                        "uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::readMyOptionsValues()",
                        "size of array for 'outputFileNames' is not equal to size of array for 'chainSizes'");
  }

  if (m_env.allOptionsMap().count(m_option_dr_maxNumExtraStages.c_str())) {
    m_maxNumExtraStages = m_env.allOptionsMap()[m_option_dr_maxNumExtraStages.c_str()].as<unsigned int>();
  }

  std::vector<double> tmpScales(0,0.);
  if (m_env.allOptionsMap().count(m_option_dr_scalesForExtraStages.c_str())) {
    std::string inputString = m_env.allOptionsMap()[m_option_dr_scalesForExtraStages.c_str()].as<std::string>();
    uqMiscReadDoublesFromString(inputString,tmpScales);
    //std::cout << "In uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::getMyOptionValues(): scales =";
    //for (unsigned int i = 0; i < tmpScales.size(); ++i) {
    //  std::cout << " " << tmpScales[i];
    //}
    //std::cout << std::endl;
  }

  if (m_maxNumExtraStages > 0) {
    m_scalesForCovMProposals.clear();
    m_lowerCholProposalCovMatrices.clear();
    m_proposalCovMatrices.clear();

    double scale = 1.0;
    unsigned int tmpSize = tmpScales.size();

    m_scalesForCovMProposals.resize      (m_maxNumExtraStages+1,1.);
    m_lowerCholProposalCovMatrices.resize(m_maxNumExtraStages+1,NULL);
    m_proposalCovMatrices.resize         (m_maxNumExtraStages+1,NULL);

    for (unsigned int i = 1; i < (m_maxNumExtraStages+1); ++i) {
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

  return;
}

template<class P_V,class P_M,class L_V,class L_M>
void
uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::calculatePosterior(
  const P_M* proposalCovMatrix,
  const P_M* mahalanobisMatrix,
  bool       applyMahalanobisInvert)
{
  if (m_chainUse2) {
    calculatePosterior(proposalCovMatrix,
                       mahalanobisMatrix,
                       applyMahalanobisInvert,
                       m_chain2);
  }
  else {
    calculatePosterior(proposalCovMatrix,
                       mahalanobisMatrix,
                       applyMahalanobisInvert,
                       m_chain1);
  }

  return;
}

template<class P_V,class P_M,class L_V,class L_M>
int
uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::prepareForNextChain(
  const P_M* proposalCovMatrix)
//const P_M* proposalPrecMatrix)
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::prepareForNextChain()..."
              << std::endl;
  }

  int iRC = UQ_OK_RC;

  const P_M* internalProposalCovMatrix = proposalCovMatrix;
  if (proposalCovMatrix == NULL) {
    P_V tmpVec(m_paramSpace.zeroVector());
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
    internalProposalCovMatrix = m_paramSpace.uqFinDimLinearSpaceClass<P_V,P_M>::newDiagMatrix(tmpVec);

    if (m_env.rank() == 0) std::cout << "In uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::prepareForNextChain()"
                                     << ", contents of internally generated proposal cov matrix are:"
                                     << std::endl;
    std::cout << *internalProposalCovMatrix;
    if (m_env.rank() == 0) std::cout << std::endl;
  }
  else {
    if (m_env.rank() == 0)  std::cout << "In uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::prepareForNextChain()"
                                      << "using suplied proposalCovMatrix, whose contents are:"
                                      << std::endl;
    std::cout << *internalProposalCovMatrix;
    if (m_env.rank() == 0) std::cout << std::endl;
  }

  m_lowerCholProposalCovMatrices[0] = new P_M(*internalProposalCovMatrix); 
  iRC = m_lowerCholProposalCovMatrices[0]->chol();
  UQ_FATAL_RC_MACRO(iRC,
                    m_env.rank(),
                    "uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::prepareForNextChain()",
                    "proposalCovMatrix is not positive definite");
  m_lowerCholProposalCovMatrices[0]->zeroUpper(false);

  m_proposalCovMatrices[0] = new P_M(*internalProposalCovMatrix);

  if (m_env.rank() == 0) std::cout << "In uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::prepareForNextChain()"
                                   << ", m_lowerCholProposalCovMatrices[0] contents are:"
                                   << std::endl;
  std::cout << *(m_lowerCholProposalCovMatrices[0]);
  if (m_env.rank() == 0) std::cout << std::endl;

#ifdef UQ_BMC_DC_REQUIRES_INVERTED_COV_MATRICES
  const P_M* internalProposalPrecMatrix = proposalPrecMatrix;
  if (proposalPrecMatrix == NULL) {
    UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                      m_env.rank(),
                      "uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::prepareForNextChain()",
                      "not yet implemented for the case 'proposalPrecMatrix == NULL'");
  }

  m_upperCholProposalPrecMatrices[0] = new P_M(*internalProposalPrecMatrix); 
  iRC = m_upperCholProposalPrecMatrices[0]->chol();
  UQ_FATAL_RC_MACRO(iRC,
                    m_env.rank(),
                    "uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::prepareForNextChain()",
                    "proposalPrecMatrix is not positive definite");
  m_upperCholProposalPrecMatrices[0]->zeroLower(false);

  m_proposalPrecMatrices[0] = new P_M(*internalProposalPrecMatrix);

  //if (m_env.rank() == 0) std::cout << "In uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::prepareForNextChain()"
  //                                 << ", m_upperCholProposalPrecMatrices[0] contents are:"
  //                                 << std::endl;
  //std::cout << *(m_upperCholProposalPrecMatrices[0]);
  //if (m_env.rank() == 0) std::cout << std::endl;
#endif

  if (m_maxNumExtraStages > 0) {
    updateCovMatrices();
  }

  if (proposalCovMatrix == NULL) delete internalProposalCovMatrix;

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::prepareForNextChain()"
              << std::endl;
  }

  return iRC;
}

template<class P_V,class P_M,class L_V,class L_M>
void
uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::updateCovMatrices()
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::updateCovMatrices()"
              << std::endl;
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "In uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::updateCovMatrices()"
              << ": m_maxNumExtraStages = "                   << m_maxNumExtraStages
              << ", m_scalesForCovMProposals.size() = "       << m_scalesForCovMProposals.size()
              << ", m_lowerCholProposalCovMatrices.size() = " << m_lowerCholProposalCovMatrices.size()
              << std::endl;
  }

  for (unsigned int i = 1; i < (m_maxNumExtraStages+1); ++i) {
    double scale = m_scalesForCovMProposals[i];
    if (m_lowerCholProposalCovMatrices[i]) delete m_lowerCholProposalCovMatrices[i];
    m_lowerCholProposalCovMatrices [i]   = new P_M(*(m_lowerCholProposalCovMatrices[i-1]));
  *(m_lowerCholProposalCovMatrices [i]) /= scale;
    if (m_proposalCovMatrices[i]) delete m_proposalCovMatrices[i];
    m_proposalCovMatrices[i]             = new P_M(*(m_proposalCovMatrices[i-1]));
  *(m_proposalCovMatrices[i])           /= (scale*scale);
#ifdef UQ_BMC_DC_REQUIRES_INVERTED_COV_MATRICES
    m_upperCholProposalPrecMatrices[i]   = new P_M(*(m_upperCholProposalPrecMatrices[i-1]));
  *(m_upperCholProposalPrecMatrices[i]) *= scale;
    m_proposalPrecMatrices[i]            = new P_M(*(m_proposalPrecMatrices[i-1]));
  *(m_proposalPrecMatrices[i])          *= (scale*scale);
#endif
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::updateCovMatrices()"
              << std::endl;
  }

  return;
}

template<class P_V,class P_M,class L_V,class L_M>
double
uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::logProposal(
  const uqChainPositionClass<P_V>& x,
  const uqChainPositionClass<P_V>& y,
  unsigned int                     idOfProposalCovMatrix)
{
  P_V diffVec(y.paramValues() - x.paramValues());
#ifdef UQ_BMC_DC_REQUIRES_INVERTED_COV_MATRICES
  double value = -0.5 * scalarProduct(diffVec, *(m_proposalPrecMatrices[idOfProposalCovMatrix]) * diffVec);
#else
  double value = -0.5 * scalarProduct(diffVec, m_proposalCovMatrices[idOfProposalCovMatrix]->invertMultiply(diffVec));
#endif
  return value;
}

template<class P_V,class P_M,class L_V,class L_M>
double
uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::logProposal(const std::vector<uqChainPositionClass<P_V>*>& inputPositions)
{
  unsigned int inputSize = inputPositions.size();
  UQ_FATAL_TEST_MACRO((inputSize < 2),
                      m_env.rank(),
                      "uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::logProposal()",
                      "inputPositions has size < 2");

  return this->logProposal(*(inputPositions[0            ]),
                           *(inputPositions[inputSize - 1]),
                           inputSize-2);
}

template<class P_V,class P_M,class L_V,class L_M>
double
uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::alpha(
  const uqChainPositionClass<P_V>& x,
  const uqChainPositionClass<P_V>& y,
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
        std::cout << "In uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::alpha()"
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
        std::cout << "In uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::alpha()"
                  << ": xOutOfBounds = " << xOutOfBounds
                  << ", yOutOfBounds = " << yOutOfBounds
                  << std::endl;
      }
  }
  if (alphaQuotientPtr != NULL) *alphaQuotientPtr = alphaQuotient;

  return std::min(1.,alphaQuotient);
}

template<class P_V,class P_M,class L_V,class L_M>
double
uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::alpha(const std::vector<uqChainPositionClass<P_V>*>& inputPositions)
{
  unsigned int inputSize = inputPositions.size();
  UQ_FATAL_TEST_MACRO((inputSize < 2),
                      m_env.rank(),
                      "uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::alpha()",
                      "inputPositions has size < 2");

  // If necessary, return 0. right away
  if (inputPositions[0          ]->outOfBounds()) return 0.;
  if (inputPositions[inputSize-1]->outOfBounds()) return 0.;

  // If inputSize is 2, recursion is not needed
  if (inputSize == 2) return this->alpha(*(inputPositions[0            ]),
                                         *(inputPositions[inputSize - 1]));

  // Prepare two vectors of positions
  std::vector<uqChainPositionClass<P_V>*>         positions(inputSize,NULL);
  std::vector<uqChainPositionClass<P_V>*> backwardPositions(inputSize,NULL);
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

template<class P_V,class P_M,class L_V,class L_M>
bool
uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::acceptAlpha(double alpha)
{
  bool result = false;

  if      (alpha <= 0.                          ) result = false;
  else if (alpha >= 1.                          ) result = true;
  else if (alpha >= gsl_rng_uniform(m_env.rng())) result = true;
  else                                            result = false;

  return result;
}

template<class P_V,class P_M,class L_V,class L_M>
void
uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::computeStatistics(
  const uqChainBaseClass<P_V>&          workingChain,
  const uqChainStatisticalOptionsClass& statisticalOptions,
  const std::string&                    chainName,
  std::ofstream*                        passedOfs)
{
  if (m_env.rank() == 0) {
    std::cout << "\n"
              << "\n-----------------------------------------------------"
              << "\n Computing statistics for chain " << chainName << " ..."
              << "\n-----------------------------------------------------"
              << "\n"
              << std::endl;
  }

  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  iRC = gettimeofday(&timevalTmp, NULL);
  double tmpRunTime = 0.;

  // Set initial positions for the computation of chain statistics
  std::vector<unsigned int> initialPosForStatistics(statisticalOptions.initialDiscardedPortions().size(),0);
  for (unsigned int i = 0; i < initialPosForStatistics.size(); ++i) {
    initialPosForStatistics[i] = (unsigned int) (statisticalOptions.initialDiscardedPortions()[i] * (double) workingChain.sequenceSize());
  }
  std::cout << "In uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::computeStatistics(): initial positions for statistics =";
  for (unsigned int i = 0; i < initialPosForStatistics.size(); ++i) {
    std::cout << " " << initialPosForStatistics[i];
  }
  std::cout << std::endl;

  //****************************************************
  // Compute mean, sample std, population std
  //****************************************************
  computeMeanVars(workingChain,
                  statisticalOptions,
                  chainName,
                  passedOfs,
                  NULL,
                  NULL,
                  NULL);

  //****************************************************
  // Compute variance of sample mean through the 'batch means method' (BMM)
  //****************************************************
  if ((statisticalOptions.bmmRun()               ) &&
      (initialPosForStatistics.size()         > 0) &&
      (statisticalOptions.bmmLengths().size() > 0)) { 
    computeBMM(workingChain,
               statisticalOptions,
               initialPosForStatistics,
               chainName,
               passedOfs);
  }

  //****************************************************
  // Compute FFT of chain, for one parameter only
  //****************************************************
  if ((statisticalOptions.fftCompute()   ) &&
      (initialPosForStatistics.size() > 0)) {
    computeFFT(workingChain,
               statisticalOptions,
               initialPosForStatistics,
               chainName,
               passedOfs);
  }

  //****************************************************
  // Compute power spectral density (PSD) of chain, for one parameter only
  //****************************************************
  if ((statisticalOptions.psdCompute()   ) &&
      (initialPosForStatistics.size() > 0)) {
    computePSD(workingChain,
               statisticalOptions,
               initialPosForStatistics,
               chainName,
               passedOfs);
  }

  //****************************************************
  // Compute power spectral density (PSD) of chain at zero frequency
  //****************************************************
  if ((statisticalOptions.psdAtZeroCompute()             ) &&
      (initialPosForStatistics.size()                 > 0) &&
      (statisticalOptions.psdAtZeroNumBlocks().size() > 0)) { 
    computePSDAtZero(workingChain,
                     statisticalOptions,
                     initialPosForStatistics,
                     chainName,
                     passedOfs);
  }

  //****************************************************
  // Compute Geweke
  //****************************************************
  if ((statisticalOptions.gewekeCompute()) &&
      (initialPosForStatistics.size() > 0)) {
    computeGeweke(workingChain,
                  statisticalOptions,
                  initialPosForStatistics,
                  chainName,
                  passedOfs);
  }

  // Set lags for the computation of chain autocorrelations
  std::vector<unsigned int> lagsForCorrs(statisticalOptions.corrNumLags(),1);
  for (unsigned int i = 1; i < lagsForCorrs.size(); ++i) {
    lagsForCorrs[i] = statisticalOptions.corrSecondLag() + (i-1)*statisticalOptions.corrLagSpacing();
  }

  //****************************************************
  // Compute autocorrelation coefficients via definition
  //****************************************************
  if ((statisticalOptions.corrComputeViaDef()) &&
      (initialPosForStatistics.size() > 0    ) &&
      (lagsForCorrs.size()            > 0    )) { 
    computeCorrViaDef(workingChain,
                      statisticalOptions,
                      initialPosForStatistics,
                      lagsForCorrs,
                      chainName,
                      passedOfs);
  }

  //****************************************************
  // Compute autocorrelation coefficients via FFT
  //****************************************************
  if ((statisticalOptions.corrComputeViaFft()) &&
      (initialPosForStatistics.size() > 0    ) &&
      (lagsForCorrs.size()            > 0    )) { 
    computeCorrViaFFT(workingChain,
                      statisticalOptions,
                      initialPosForStatistics,
                      lagsForCorrs,
                      chainName,
                      passedOfs);
  }

  //****************************************************
  // Compute histogram and/or Kde
  //****************************************************
  if ((statisticalOptions.histCompute()) ||
      (statisticalOptions.kdeCompute() )) {
    computeHistKde(workingChain,
                   statisticalOptions,
                   chainName,
                   passedOfs);
  }

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.rank() == 0) {
    std::cout << "All statistics took " << tmpRunTime
              << " seconds"
              << std::endl;
  }

  if (m_env.rank() == 0) {
    std::cout << "\n-----------------------------------------------------"
              << "\n Finished computing statistics for chain " << chainName
              << "\n-----------------------------------------------------"
              << "\n"
              << std::endl;
  }

  return;
}

template<class P_V,class P_M,class L_V,class L_M>
void
uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::computeMeanVars(
  const uqChainBaseClass<P_V>&          workingChain,
  const uqChainStatisticalOptionsClass& statisticalOptions,
  const std::string&                    chainName,
  std::ofstream*                        passedOfs,
  P_V*                                  meanPtr,
  P_V*                                  sampleVarPtr,
  P_V*                                  populVarPtr)
{
  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  iRC = gettimeofday(&timevalTmp, NULL);
  double tmpRunTime = 0.;

  if (m_env.rank() == 0) {
    std::cout << "\n-----------------------------------------------------"
              << "\nComputing mean, sample variance and population variance"
              << std::endl;
  }

  P_V chainMean(m_paramSpace.zeroVector());
  workingChain.mean(0,
                    workingChain.sequenceSize(),
                    chainMean);

  P_V chainSampleVariance(m_paramSpace.zeroVector());
  workingChain.sampleVariance(0,
                              workingChain.sequenceSize(),
                              chainMean,
                              chainSampleVariance);

  if (m_env.rank() == 0) {
    std::cout << "\nEstimated variance of sample mean for the whole chain " << chainName
              << ", under independence assumption:"
              << std::endl;
  }
  P_V estimatedVarianceOfSampleMean(chainSampleVariance);
  estimatedVarianceOfSampleMean /= (double) workingChain.sequenceSize();
  bool savedVectorPrintState = estimatedVarianceOfSampleMean.getPrintHorizontally();
  estimatedVarianceOfSampleMean.setPrintHorizontally(false);
  std::cout << estimatedVarianceOfSampleMean;
  estimatedVarianceOfSampleMean.setPrintHorizontally(savedVectorPrintState);
  if (m_env.rank() == 0) {
    std::cout << std::endl;
  }

  P_V chainPopulationVariance(m_paramSpace.zeroVector());
  workingChain.populationVariance(0,
                                  workingChain.sequenceSize(),
                                  chainMean,
                                  chainPopulationVariance);

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.rank() == 0) {
    std::cout << "Mean and variances took " << tmpRunTime
              << " seconds"
              << std::endl;
  }

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
	      chainMean[i],
              " ",
              sqrt(chainSampleVariance[i]),
              " ",
              sqrt(chainPopulationVariance[i]));
      std::cout << line;
    }
    std::cout << std::endl;
  }

  if (meanPtr     ) *meanPtr      = chainMean;
  if (sampleVarPtr) *sampleVarPtr = chainSampleVariance;
  if (populVarPtr ) *populVarPtr  = chainPopulationVariance;

  return;
}

template<class P_V,class P_M,class L_V,class L_M>
void
uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::computeBMM(
  const uqChainBaseClass<P_V>&          workingChain,
  const uqChainStatisticalOptionsClass& statisticalOptions,
  const std::vector<unsigned int>&      initialPosForStatistics,
  const std::string&                    chainName,
  std::ofstream*                        passedOfs)
{
  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  iRC = gettimeofday(&timevalTmp, NULL);
  double tmpRunTime = 0.;

  if (m_env.rank() == 0) {
    std::cout << "\n-----------------------------------------------------"
              << "\nComputing variance of sample mean through BMM"
              << std::endl;
  }

  std::cout << "In uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::computeBMM(): lengths for batchs in BMM =";
  for (unsigned int i = 0; i < statisticalOptions.bmmLengths().size(); ++i) {
    std::cout << " " << statisticalOptions.bmmLengths()[i];
  }
  std::cout << std::endl;

  uq2dArrayOfStuff<P_V> _2dArrayOfBMM(initialPosForStatistics.size(),statisticalOptions.bmmLengths().size());
  for (unsigned int i = 0; i < _2dArrayOfBMM.numRows(); ++i) {
    for (unsigned int j = 0; j < _2dArrayOfBMM.numCols(); ++j) {
      _2dArrayOfBMM.setLocation(i,j,m_paramSpace.newVector());
    }
  }
  P_V bmmVec(m_paramSpace.zeroVector());
  for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
    unsigned int initialPos = initialPosForStatistics[initialPosId];
    for (unsigned int batchLengthId = 0; batchLengthId < statisticalOptions.bmmLengths().size(); batchLengthId++) {
      unsigned int batchLength = statisticalOptions.bmmLengths()[batchLengthId];
      workingChain.bmm(initialPos,
                       batchLength,
                       bmmVec);
      _2dArrayOfBMM(initialPosId,batchLengthId) = bmmVec;
    }
  }

  if (m_env.rank() == 0) {
    for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
      std::cout << "\nEstimated variance of sample mean, through batch means method, for subchain beggining at position " << initialPosForStatistics[initialPosId]
                << " (each column corresponds to a batch length)"
                << std::endl;

      char line[512];
      sprintf(line,"%s",
              "Parameter");
      std::cout << line;
      for (unsigned int batchLengthId = 0; batchLengthId < statisticalOptions.bmmLengths().size(); batchLengthId++) {
        sprintf(line,"%10s%3d",
                " ",
                statisticalOptions.bmmLengths()[batchLengthId]);
        std::cout << line;
      }

      for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
        sprintf(line,"\n%9.9s",
                m_paramSpace.parameter(i).name().c_str());
        std::cout << line;
        for (unsigned int batchLengthId = 0; batchLengthId < statisticalOptions.bmmLengths().size(); batchLengthId++) {
          sprintf(line,"%2s%11.4e",
                  " ",
                  _2dArrayOfBMM(initialPosId,batchLengthId)[i]);
          std::cout << line;
        }
      }
      std::cout << std::endl;
    }
  }

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.rank() == 0) {
    std::cout << "Chain BMM took " << tmpRunTime
              << " seconds"
              << std::endl;
  }

  return;
}

template<class P_V,class P_M,class L_V,class L_M>
void
uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::computeFFT(
  const uqChainBaseClass<P_V>&          workingChain,
  const uqChainStatisticalOptionsClass& statisticalOptions,
  const std::vector<unsigned int>&      initialPosForStatistics,
  const std::string&                    chainName,
  std::ofstream*                        passedOfs)
{
  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  iRC = gettimeofday(&timevalTmp, NULL);
  double tmpRunTime = 0.;

  if (m_env.rank() == 0) {
    std::cout << "\n-----------------------------------------------------"
              << "\nComputing FFT of chain on parameter of id = " << statisticalOptions.fftParamId()
              << std::endl;
  }

  std::vector<std::complex<double> > forwardResult(0,std::complex<double>(0.,0.));
  std::vector<std::complex<double> > inverseResult(0,std::complex<double>(0.,0.));
  uqFftClass<std::complex<double> > fftObj(m_env);
  for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
    unsigned int initialPosition = initialPosForStatistics[initialPosId];
    workingChain.fftForward(initialPosition,
                            statisticalOptions.fftSize(),
                            statisticalOptions.fftParamId(),
                            forwardResult);

    if (statisticalOptions.fftWrite() && passedOfs) {
      std::ofstream& ofs = *passedOfs;
      ofs << chainName << "_fft_initPos" << initialPosForStatistics[initialPosId] << " = zeros(" << 1
          << ","                                                                                 << forwardResult.size()
          << ");"
          << std::endl;
      for (unsigned int j = 0; j < forwardResult.size(); ++j) {
        ofs << chainName << "_fft_initPos" << initialPosForStatistics[initialPosId] << "(" << 1
            << ","                                                                         << j+1
            << ") = "                                                                      << forwardResult[j].real()
            << " + i*"                                                                     << forwardResult[j].imag()
            << ";"
            << std::endl;
      }
    } // if write

    if (statisticalOptions.fftTestInversion()) {
      fftObj.inverse(forwardResult,
                     statisticalOptions.fftSize(),
                     inverseResult);
      if (statisticalOptions.fftWrite() && passedOfs) {
        std::ofstream& ofs = *passedOfs;
        ofs << chainName << "_inv_initPos" << initialPosForStatistics[initialPosId] << " = zeros(" << 1
            << ","                                                                                 << inverseResult.size()
            << ");"
            << std::endl;
        for (unsigned int j = 0; j < inverseResult.size(); ++j) {
          ofs << chainName << "_inv_initPos" << initialPosForStatistics[initialPosId] << "(" << 1
              << ","                                                                         << j+1
              << ") = "                                                                      << inverseResult[j].real()
              << " + i*"                                                                     << inverseResult[j].imag()
              << ";"
              << std::endl;
        }
      } // if write
    }
  } // for initialPosId

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.rank() == 0) {
    std::cout << "Chain FFT took " << tmpRunTime
              << " seconds"
              << std::endl;
  }

  return;
}

template<class P_V,class P_M,class L_V,class L_M>
void
uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::computePSD(
  const uqChainBaseClass<P_V>&          workingChain,
  const uqChainStatisticalOptionsClass& statisticalOptions,
  const std::vector<unsigned int>&      initialPosForStatistics,
  const std::string&                    chainName,
  std::ofstream*                        passedOfs)
{
  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  iRC = gettimeofday(&timevalTmp, NULL);
  double tmpRunTime = 0.;

  if (m_env.rank() == 0) {
    std::cout << "\n-----------------------------------------------------"
              << "\nComputing PSD of chain on parameter of id = " << statisticalOptions.psdParamId()
              << std::endl;
  }

  std::vector<double> psdResult(0,0.);
  for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
    unsigned int initialPosition = initialPosForStatistics[initialPosId];
    workingChain.psd(initialPosition,
                     statisticalOptions.psdNumBlocks(),
                     statisticalOptions.psdHopSizeRatio(),
                     statisticalOptions.psdParamId(),
                     psdResult);

    if (statisticalOptions.psdWrite() && passedOfs) {
      std::ofstream& ofs = *passedOfs;
      ofs << chainName << "_psd_initPos" << initialPosForStatistics[initialPosId] << " = zeros(" << 1
          << ","                                                                                 << psdResult.size()
          << ");"
          << std::endl;
      for (unsigned int j = 0; j < psdResult.size(); ++j) {
        ofs << chainName << "_psd_initPos" << initialPosForStatistics[initialPosId] << "(" << 1
            << ","                                                                         << j+1
            << ") = "                                                                      << psdResult[j]
            << ";"
            << std::endl;
      }
    } // if write
  }

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.rank() == 0) {
    std::cout << "Chain PSD took " << tmpRunTime
              << " seconds"
              << std::endl;
  }

  return;
}

template<class P_V,class P_M,class L_V,class L_M>
void
uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::computePSDAtZero(
  const uqChainBaseClass<P_V>&          workingChain,
  const uqChainStatisticalOptionsClass& statisticalOptions,
  const std::vector<unsigned int>&      initialPosForStatistics,
  const std::string&                    chainName,
  std::ofstream*                        passedOfs)
{
  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  iRC = gettimeofday(&timevalTmp, NULL);
  double tmpRunTime = 0.;

  if (m_env.rank() == 0) {
    std::cout << "\n-----------------------------------------------------"
              << "\nComputing PSD at frequency zero for all parameters"
              << std::endl;
  }

  uq2dArrayOfStuff<P_V> _2dArrayOfPSDAtZero(initialPosForStatistics.size(),statisticalOptions.psdAtZeroNumBlocks().size());
  for (unsigned int i = 0; i < _2dArrayOfPSDAtZero.numRows(); ++i) {
    for (unsigned int j = 0; j < _2dArrayOfPSDAtZero.numCols(); ++j) {
      _2dArrayOfPSDAtZero.setLocation(i,j,m_paramSpace.newVector());
    }
  }
  P_V psdVec(m_paramSpace.zeroVector());
  for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
    unsigned int initialPosition = initialPosForStatistics[initialPosId];
    for (unsigned int numBlocksId = 0; numBlocksId < statisticalOptions.psdAtZeroNumBlocks().size(); numBlocksId++) {
      unsigned int numBlocks = statisticalOptions.psdAtZeroNumBlocks()[numBlocksId];
      workingChain.psdAtZero(initialPosition,
                             numBlocks,
                             statisticalOptions.psdAtZeroHopSizeRatio(),
                             psdVec);
      _2dArrayOfPSDAtZero(initialPosId,numBlocksId) = psdVec;
    }
  }

  // Display PSD at frequency zero
  if ((statisticalOptions.psdAtZeroDisplay()) && (m_env.rank() == 0)) {
    for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
      unsigned int initialPos = initialPosForStatistics[initialPosId];
      std::cout << "\nComputed PSD at frequency zero for subchain beggining at position " << initialPos
                << ", so effective data size = " << workingChain.sequenceSize() - initialPos
                << " (each column corresponds to a number of blocks)"
                << std::endl;

      char line[512];
      sprintf(line,"%s",
              "Parameter");
      std::cout << line;
      for (unsigned int numBlocksId = 0; numBlocksId < statisticalOptions.psdAtZeroNumBlocks().size(); numBlocksId++) {
        sprintf(line,"%10s%3d",
                " ",
                statisticalOptions.psdAtZeroNumBlocks()[numBlocksId]);
        std::cout << line;
      }

      for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
        sprintf(line,"\n%9.9s",
                m_paramSpace.parameter(i).name().c_str());
        std::cout << line;
        for (unsigned int numBlocksId = 0; numBlocksId < statisticalOptions.psdAtZeroNumBlocks().size(); numBlocksId++) {
          sprintf(line,"%2s%11.4e",
                  " ",
                  _2dArrayOfPSDAtZero(initialPosId,numBlocksId)[i]);
          std::cout << line;
        }
      }
      std::cout << std::endl;
    }
  }

  // Display estimated variance of sample mean through PSD
  if (/*(statisticalOptions.psdAtZeroDisplay()) &&*/ (m_env.rank() == 0)) {
    for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
      unsigned int initialPos = initialPosForStatistics[initialPosId];
      std::cout << "\nEstimated variance of sample mean, through psd, for subchain beggining at position " << initialPos
                << ", so effective data size = " << workingChain.sequenceSize() - initialPos
                << " (each column corresponds to a number of blocks)"
                << std::endl;

      char line[512];
      sprintf(line,"%s",
              "Parameter");
      std::cout << line;
      for (unsigned int numBlocksId = 0; numBlocksId < statisticalOptions.psdAtZeroNumBlocks().size(); numBlocksId++) {
        sprintf(line,"%10s%3d",
                " ",
                statisticalOptions.psdAtZeroNumBlocks()[numBlocksId]);
        std::cout << line;
      }

      for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
        sprintf(line,"\n%9.9s",
                m_paramSpace.parameter(i).name().c_str());
        std::cout << line;
        for (unsigned int numBlocksId = 0; numBlocksId < statisticalOptions.psdAtZeroNumBlocks().size(); numBlocksId++) {
          sprintf(line,"%2s%11.4e",
                  " ",
                  2.*M_PI*_2dArrayOfPSDAtZero(initialPosId,numBlocksId)[i]/(double) (workingChain.sequenceSize() - initialPos));
          std::cout << line;
        }
      }
      std::cout << std::endl;
    }
  }

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.rank() == 0) {
    std::cout << "Chain PSD at frequency zero took " << tmpRunTime
              << " seconds"
              << std::endl;
  }

  // Write PSD at frequency zero
  if (statisticalOptions.psdAtZeroWrite() && passedOfs) {
    std::ofstream& ofs = *passedOfs;
    ofs << chainName << "_psdAtZero_numBlocks = zeros(" << 1
        << ","                                          << statisticalOptions.psdAtZeroNumBlocks().size()
        << ");"
        << std::endl;
    for (unsigned int numBlocksId = 0; numBlocksId < statisticalOptions.psdAtZeroNumBlocks().size(); numBlocksId++) {
      ofs << chainName << "_psdAtZero_numBlocks(" << 1
          << ","                                  << numBlocksId+1
          << ") = "                               << statisticalOptions.psdAtZeroNumBlocks()[numBlocksId]
          << ";"
          << std::endl;
    }

    for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
      ofs << chainName << "_psdAtZero_initPos" << initialPosForStatistics[initialPosId] << " = zeros(" << m_paramSpace.dim()
          << ","                                                                                       << statisticalOptions.psdAtZeroNumBlocks().size()
          << ");"
          << std::endl;
      for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
        for (unsigned int numBlocksId = 0; numBlocksId < statisticalOptions.psdAtZeroNumBlocks().size(); numBlocksId++) {
          ofs << chainName << "_psdAtZero_initPos" << initialPosForStatistics[initialPosId] << "(" << i+1
              << ","                                                                               << numBlocksId+1
              << ") = "                                                                            << _2dArrayOfPSDAtZero(initialPosId,numBlocksId)[i]
              << ";"
              << std::endl;
        }
      }
    }
  } 

  return;
}

template<class P_V,class P_M,class L_V,class L_M>
void
uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::computeGeweke(
  const uqChainBaseClass<P_V>&          workingChain,
  const uqChainStatisticalOptionsClass& statisticalOptions,
  const std::vector<unsigned int>&      initialPosForStatistics,
  const std::string&                    chainName,
  std::ofstream*                        passedOfs)
{
  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  iRC = gettimeofday(&timevalTmp, NULL);
  double tmpRunTime = 0.;

  if (m_env.rank() == 0) {
    std::cout << "\n-----------------------------------------------------"
              << "\nComputing Geweke coefficients"
              << std::endl;
  }

  std::vector<P_V*> vectorOfGeweke(initialPosForStatistics.size(),NULL);
  P_V gewVec(m_paramSpace.zeroVector());
  for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
    unsigned int initialPosition = initialPosForStatistics[initialPosId];
    workingChain.geweke(initialPosition,
                        statisticalOptions.gewekeNaRatio(),
                        statisticalOptions.gewekeNbRatio(),
                        gewVec);
    vectorOfGeweke[initialPosId] = new P_V(gewVec);
  }

  if (m_env.rank() == 0) {
    std::cout << "\nComputed Geweke coefficients with 10% and 50% percentages"
              << " (each column corresponds to a different initial position on the full chain)"
              << std::endl;

    char line[512];
    sprintf(line,"%s",
            "Parameter");
    std::cout << line;
    for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
      sprintf(line,"%10s%3d",
              " ",
              initialPosForStatistics[initialPosId]);
      std::cout << line;
    }

    for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
      sprintf(line,"\n%9.9s",
              m_paramSpace.parameter(i).name().c_str());
      std::cout << line;
      for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
        sprintf(line,"%2s%11.4e",
                " ",
                (*(vectorOfGeweke[initialPosId]))[i]);
        std::cout << line;
      }
    }
    std::cout << std::endl;
  }

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.rank() == 0) {
    std::cout << "Chain Geweke took " << tmpRunTime
              << " seconds"
              << std::endl;
  }

  return;
}

template<class P_V,class P_M,class L_V,class L_M>
void
uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::computeCorrViaDef(
  const uqChainBaseClass<P_V>&          workingChain,
  const uqChainStatisticalOptionsClass& statisticalOptions,
  const std::vector<unsigned int>&      initialPosForStatistics,
  const std::vector<unsigned int>&      lagsForCorrs,
  const std::string&                    chainName,
  std::ofstream*                        passedOfs)
{
  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  iRC = gettimeofday(&timevalTmp, NULL);
  double tmpRunTime = 0.;

  if (m_env.rank() == 0) {
    std::cout << "\n-----------------------------------------------------"
              << "\nComputing autocorrelation coefficients (via def)"
              << std::endl;
  }

  if (statisticalOptions.corrDisplay() && (m_env.rank() == 0)) {
    std::cout << "In uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::computeCorrViaDef(): lags for autocorrelation (via def) = ";
    for (unsigned int i = 0; i < lagsForCorrs.size(); ++i) {
      std::cout << " " << lagsForCorrs[i];
    }
    std::cout << std::endl;
  }

  uq2dArrayOfStuff<P_V> _2dArrayOfAutoCorrs(initialPosForStatistics.size(),lagsForCorrs.size());
  for (unsigned int i = 0; i < _2dArrayOfAutoCorrs.numRows(); ++i) {
    for (unsigned int j = 0; j < _2dArrayOfAutoCorrs.numCols(); ++j) {
      _2dArrayOfAutoCorrs.setLocation(i,j,m_paramSpace.newVector());
    }
  }
  //V corrVec(m_paramSpace.zeroVector());
  for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
    unsigned int initialPos = initialPosForStatistics[initialPosId];
    for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
      unsigned int lag = lagsForCorrs[lagId];
      workingChain.autoCorrViaDef(initialPos,
                                  workingChain.sequenceSize()-initialPos,
                                  lag,
                                  _2dArrayOfAutoCorrs(initialPosId,lagId));
      //_2dArrayOfAutoCorrs(initialPosId,lagId) = corrVec;
    }
  }

  // It is not practical to compute the variance of sample mean by computing the autocorrelations via definition for each lag
  // The code computes the variance of sample mean by computing the autocorrelations via fft, below, in another routine

  if ((statisticalOptions.corrDisplay()) && (m_env.rank() == 0)) {
    for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
      std::cout << "\nComputed autocorrelation coefficients (via def), for subchain beggining at position " << initialPosForStatistics[initialPosId]
                << " (each column corresponds to a different lag)"
                << std::endl;
      char line[512];
      sprintf(line,"%s",
              "Parameter");
      std::cout << line;
      for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
        sprintf(line,"%10s%3d",
                " ",
                lagsForCorrs[lagId]);
        std::cout << line;
      }

      for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
        sprintf(line,"\n%9.9s",
                m_paramSpace.parameter(i).name().c_str());
        std::cout << line;
        for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
          sprintf(line,"%2s%11.4e",
                  " ",
                  _2dArrayOfAutoCorrs(initialPosId,lagId)[i]);
          std::cout << line;
        }
      }
      std::cout << std::endl;
    }
  }

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.rank() == 0) {
    std::cout << "Chain autocorrelation (via def) took " << tmpRunTime
              << " seconds"
              << std::endl;
  }

  // Write autocorrelations
  if (statisticalOptions.corrWrite() && passedOfs) {
    std::ofstream& ofs = *passedOfs;
    ofs << chainName << "_corrViaDef_lags = zeros(" << 1
        << ","                                      << lagsForCorrs.size()
        << ");"
        << std::endl;
    for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
      ofs << chainName << "_corrViaDef_lags(" << 1
          << ","                              << lagId+1
          << ") = "                           << lagsForCorrs[lagId]
          << ";"
          << std::endl;
    }

    for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
      ofs << chainName << "_corrViaDef_initPos" << initialPosForStatistics[initialPosId] << " = zeros(" << m_paramSpace.dim()
          << ","                                                                                        << lagsForCorrs.size()
          << ");"
          << std::endl;
      for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
        for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
          ofs << chainName << "_corrViaDef_initPos" << initialPosForStatistics[initialPosId] << "(" << i+1
              << ","                                                                                << lagId+1
              << ") = "                                                                             << _2dArrayOfAutoCorrs(initialPosId,lagId)[i]
              << ";"
              << std::endl;
        }
      }
    }
  } 

  return;
}

template<class P_V,class P_M,class L_V,class L_M>
void
uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::computeCorrViaFFT(
  const uqChainBaseClass<P_V>&          workingChain,
  const uqChainStatisticalOptionsClass& statisticalOptions,
  const std::vector<unsigned int>&      initialPosForStatistics,
  const std::vector<unsigned int>&      lagsForCorrs,
  const std::string&                    chainName,
  std::ofstream*                        passedOfs)
{
  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  iRC = gettimeofday(&timevalTmp, NULL);
  double tmpRunTime = 0.;

  if (m_env.rank() == 0) {
    std::cout << "\n-----------------------------------------------------"
              << "\nComputing autocorrelation coefficients (via fft)"
              << std::endl;
  }

  if (statisticalOptions.corrDisplay() && (m_env.rank() == 0)) {
    std::cout << "In uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::computeCorrViaFFT(): lags for autocorrelation (via fft) = ";
    for (unsigned int i = 0; i < lagsForCorrs.size(); ++i) {
      std::cout << " " << lagsForCorrs[i];
     }
     std::cout << std::endl;
  }

  uq2dArrayOfStuff<P_V> _2dArrayOfAutoCorrs(initialPosForStatistics.size(),lagsForCorrs.size());
  for (unsigned int i = 0; i < _2dArrayOfAutoCorrs.numRows(); ++i) {
    for (unsigned int j = 0; j < _2dArrayOfAutoCorrs.numCols(); ++j) {
      _2dArrayOfAutoCorrs.setLocation(i,j,m_paramSpace.newVector());
    }
  }
  std::vector<P_V*> corrVecs(lagsForCorrs.size(),NULL);
  std::vector<P_V*> corrSumVecs(initialPosForStatistics.size(),NULL);
  for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
    corrSumVecs[initialPosId] = m_paramSpace.newVector();
    unsigned int initialPos = initialPosForStatistics[initialPosId];
    for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
      corrVecs[lagId] = m_paramSpace.newVector();
    }
    if (m_env.rank() == 0) {
      std::cout << "In uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::computeCorrViaFFT()"
                << ": about to call chain.autoCorrViaFft()"
                << " with initialPos = "      << initialPos
                << ", numPos = "              << workingChain.sequenceSize()-initialPos
                << ", lagsForCorrs.size() = " << lagsForCorrs.size()
                << ", corrVecs.size() = "     << corrVecs.size()
                << std::endl;
    }
    workingChain.autoCorrViaFft(initialPos,
                                workingChain.sequenceSize()-initialPos, // Use all possible data positions
                                lagsForCorrs,
                                corrVecs);
    workingChain.autoCorrViaFft(initialPos,
                                workingChain.sequenceSize()-initialPos, // Use all possible data positions
                                (unsigned int) (1.0 * (double) (workingChain.sequenceSize()-initialPos)), // CHECK
                                *corrSumVecs[initialPosId]); // Sum of all possibly computable autocorrelations, not only the asked ones in lagsForCorrs
    for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
      _2dArrayOfAutoCorrs(initialPosId,lagId) = *(corrVecs[lagId]);
    }
  }
  for (unsigned int j = 0; j < corrVecs.size(); ++j) {
    if (corrVecs[j] != NULL) delete corrVecs[j];
  }

  if ((statisticalOptions.corrDisplay()) && (m_env.rank() == 0)) {
    P_V chainMean                    (m_paramSpace.zeroVector());
    P_V chainSampleVariance          (m_paramSpace.zeroVector());
    P_V estimatedVarianceOfSampleMean(m_paramSpace.zeroVector());
    for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
      unsigned int initialPos = initialPosForStatistics[initialPosId];

      workingChain.mean(initialPos,
                        workingChain.sequenceSize()-initialPos,
                        chainMean);

      workingChain.sampleVariance(initialPos,
                                  workingChain.sequenceSize()-initialPos,
                                  chainMean,
                                  chainSampleVariance);

      std::cout << "\nEstimated variance of sample mean, through autocorrelation (via fft), for subchain beggining at position " << initialPosForStatistics[initialPosId]
                << std::endl;
      estimatedVarianceOfSampleMean.cwSet(-1.); // Yes, '-1' because the autocorrelation at lag 0, which values '+1', is already counted in the sum
      estimatedVarianceOfSampleMean += 2.* (*corrSumVecs[initialPosId]);
      estimatedVarianceOfSampleMean *= chainSampleVariance;
      estimatedVarianceOfSampleMean /= (double) (workingChain.sequenceSize() - initialPos);
      bool savedVectorPrintState = estimatedVarianceOfSampleMean.getPrintHorizontally();
      estimatedVarianceOfSampleMean.setPrintHorizontally(false);
      std::cout << estimatedVarianceOfSampleMean;
      estimatedVarianceOfSampleMean.setPrintHorizontally(savedVectorPrintState);
      if (m_env.rank() == 0) {
        std::cout << std::endl;
      }

      std::cout << "\nComputed autocorrelation coefficients (via fft), for subchain beggining at position " << initialPosForStatistics[initialPosId]
                << " (each column corresponds to a different lag)"
                << std::endl;

      char line[512];
      sprintf(line,"%s",
              "Parameter");
      std::cout << line;
      for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
        sprintf(line,"%10s%3d",
                " ",
                lagsForCorrs[lagId]);
        std::cout << line;
      }

      for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
        sprintf(line,"\n%9.9s",
                m_paramSpace.parameter(i).name().c_str());
        std::cout << line;
        for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
          sprintf(line,"%2s%11.4e",
                  " ",
                  _2dArrayOfAutoCorrs(initialPosId,lagId)[i]);
          std::cout << line;
        }
      }
      std::cout << std::endl;
    }
  }

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.rank() == 0) {
    std::cout << "Chain autocorrelation (via fft) took " << tmpRunTime
              << " seconds"
              << std::endl;
  }

  // Write autocorrelations
  if (statisticalOptions.corrWrite() && passedOfs) {
    std::ofstream& ofs = *passedOfs;
    ofs << chainName << "_corrViaFft_lags = zeros(" << 1
        << ","                                      << lagsForCorrs.size()
        << ");"
        << std::endl;
    for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
      ofs << chainName << "_corrViaFft_lags(" << 1
          << ","                              << lagId+1
          << ") = "                           << lagsForCorrs[lagId]
          << ";"
          << std::endl;
    }

    for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
      ofs << chainName << "_corrViaFft_initPos" << initialPosForStatistics[initialPosId] << " = zeros(" << m_paramSpace.dim()
          << ","                                                                                        << lagsForCorrs.size()
          << ");"
          << std::endl;
      for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
        for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
          ofs << chainName << "_corrViaFft_initPos" << initialPosForStatistics[initialPosId] << "(" << i+1
              << ","                                                                                << lagId+1
              << ") = "                                                                             << _2dArrayOfAutoCorrs(initialPosId,lagId)[i]
              << ";"
              << std::endl;
        }
      }
    }
  } 

  return;
}

template<class P_V,class P_M,class L_V,class L_M>
void
uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::computeFilterParameters(
  const uqChainBaseClass<P_V>&          workingChain,
  const uqChainStatisticalOptionsClass& statisticalOptions,
  const std::string&                    chainName,
  std::ofstream*                        passedOfs,
  unsigned int&                         initialPos,
  unsigned int&                         spacing)
{
  if (m_env.rank() == 0) {
    std::cout << "\n"
              << "\n-----------------------------------------------------"
              << "\n Computing filter parameters for chain " << chainName << " ..."
              << "\n-----------------------------------------------------"
              << "\n"
              << std::endl;
  }

  initialPos = 0;
  spacing    = 1;

  if (m_env.rank() == 0) {
    std::cout << "\n-----------------------------------------------------"
              << "\n Finished computing filter parameters for chain " << chainName
              << ": initialPos = " << initialPos
              << ", spacing = "    << spacing
              << "\n-----------------------------------------------------"
              << "\n"
              << std::endl;
  }

  return;
}

template<class P_V,class P_M,class L_V,class L_M>
void
uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::computeHistKde(
  const uqChainBaseClass<P_V>&          workingChain, // Use the whole chain
  const uqChainStatisticalOptionsClass& statisticalOptions,
  const std::string&                    chainName,
  std::ofstream*                        passedOfs)
{
  if (m_env.rank() == 0) {
    std::cout << "\n"
              << "\n-----------------------------------------------------"
              << "\n Computing histogram and/or KDE for chain " << chainName << " ..."
              << "\n-----------------------------------------------------"
              << "\n"
              << std::endl;
  }

  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;

  //****************************************************
  // Compute MIN and MAX: for histograms and KDE
  //****************************************************
  double tmpRunTime = 0.;
  iRC = gettimeofday(&timevalTmp, NULL);
  if (m_env.rank() == 0) {
    std::cout << "\n-----------------------------------------------------"
              << "\nComputing min and max for histograms and KDE"
              << std::endl;
  }

  P_V statsMinPositions(m_paramSpace.zeroVector());
  P_V statsMaxPositions(m_paramSpace.zeroVector());
  workingChain.minMax(0, // Use the whole chain
                      statsMinPositions,
                      statsMaxPositions);

  if (m_env.rank() == 0) {
    std::cout << "\nComputed min values and max values for chain " << chainName
              << std::endl;

    char line[512];
    sprintf(line,"%s",
            "Parameter");
    std::cout << line;

    sprintf(line,"%9s%s%9s%s",
            " ",
            "min",
            " ",
            "max");
    std::cout << line;

    for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
      sprintf(line,"\n%8.8s",
              m_paramSpace.parameter(i).name().c_str());
      std::cout << line;

      sprintf(line,"%2s%11.4e%2s%11.4e",
              " ",
              statsMinPositions[i],
              " ",
              statsMaxPositions[i]);
      std::cout << line;
    }
    std::cout << std::endl;
  }

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.rank() == 0) {
    std::cout << "Chain min and max took " << tmpRunTime
              << " seconds"
              << std::endl;
  }

  //****************************************************
  // Compute histograms
  //****************************************************
  if ((statisticalOptions.histCompute()            ) &&
      (statisticalOptions.histNumInternalBins() > 0)) {
    tmpRunTime = 0.;
    iRC = gettimeofday(&timevalTmp, NULL);
    if (m_env.rank() == 0) {
      std::cout << "\n-----------------------------------------------------"
                << "\nComputing histograms"
                << std::endl;
    }

    std::vector<P_V*> histCentersForAllBins(0);
    std::vector<P_V*> histBinsForAllParams(0);

    for (unsigned int i = 0; i < statsMaxPositions.size(); ++i) {
      statsMaxPositions[i] *= (1. + 1.e-15);
    }

    histCentersForAllBins.resize(statisticalOptions.histNumInternalBins()+2,NULL);
    histBinsForAllParams.resize (statisticalOptions.histNumInternalBins()+2,NULL);
    workingChain.histogram(0, // Use the whole chain
                           statsMinPositions,
                           statsMaxPositions,
                           histCentersForAllBins,
                           histBinsForAllParams);

    tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
    if (m_env.rank() == 0) {
      std::cout << "Chain histograms took " << tmpRunTime
                << " seconds"
                << std::endl;
    }

    // Write histograms
    // plot(queso_centersOfHistBins(1,:)',queso_histBins(1,:)','r-');
    if (passedOfs) {
      std::ofstream& ofs = *passedOfs;
      ofs << chainName << "_centersOfHistBins = zeros(" << m_paramSpace.dim()
          << ","                                        << histCentersForAllBins.size()
          << ");"
          << std::endl;
      for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
        for (unsigned int j = 0; j < histCentersForAllBins.size(); ++j) {
           ofs << chainName << "_centersOfHistBins(" << i+1
               << ","                                << j+1
               << ") = "                             << (*(histCentersForAllBins[j]))[i]
               << ";"
               << std::endl;
        }
      }

      ofs << chainName << "_histBins = zeros(" << m_paramSpace.dim()
          << ","                               << histBinsForAllParams.size()
          << ");"
          << std::endl;
      for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
        for (unsigned int j = 0; j < histBinsForAllParams.size(); ++j) {
           ofs << chainName << "_histBins(" << i+1
               << ","                       << j+1
               << ") = "                    << (*(histBinsForAllParams[j]))[i]
               << ";"
               << std::endl;
        }
      }
    }

    for (unsigned int i = 0; i < histBinsForAllParams.size(); ++i) {
      if (histBinsForAllParams[i] != NULL) delete histBinsForAllParams[i];
    }
    for (unsigned int i = 0; i < histCentersForAllBins.size(); ++i) {
      if (histCentersForAllBins[i] != NULL) delete histCentersForAllBins[i];
    }
  }

  //****************************************************
  // Compute estimations of probability densities
  //****************************************************
  if ((statisticalOptions.kdeCompute()             ) &&
      (statisticalOptions.kdeNumEvalPositions() > 0)) {
    tmpRunTime = 0.;
    iRC = gettimeofday(&timevalTmp, NULL);
    if (m_env.rank() == 0) {
      std::cout << "\n-----------------------------------------------------"
                << "\nComputing KDE"
                << std::endl;
    }

    std::vector<P_V*> kdeEvalPositions(0);
    P_V               gaussianKdeScaleVec(m_paramSpace.zeroVector());
    std::vector<P_V*> gaussianKdeDensities(0);

    kdeEvalPositions.resize(statisticalOptions.kdeNumEvalPositions(),NULL);
    uqMiscComputePositionsBetweenMinMax(statsMinPositions,
                                        statsMaxPositions,
                                        kdeEvalPositions);

    P_V iqrVec(m_paramSpace.zeroVector());
    workingChain.interQuantileRange(0, // Use the whole chain
                                    iqrVec);

    if (m_env.rank() == 0) {
      std::cout << "\nComputed inter quantile ranges for chain " << chainName
                  << std::endl;

      char line[512];
      sprintf(line,"%s",
              "Parameter");
      std::cout << line;

      sprintf(line,"%9s%s",
              " ",
              "iqr");
      std::cout << line;

      for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
        sprintf(line,"\n%8.8s",
                m_paramSpace.parameter(i).name().c_str());
        std::cout << line;

        sprintf(line,"%2s%11.4e",
                " ",
                iqrVec[i]);
        std::cout << line;
      }
      std::cout << std::endl;
    }

    workingChain.scalesForKDE(0, // Use the whole chain
                              iqrVec,
                              gaussianKdeScaleVec);

    gaussianKdeDensities.resize(statisticalOptions.kdeNumEvalPositions(),NULL);
    workingChain.gaussianKDE(0, // Use the whole chain
                             gaussianKdeScaleVec,
                             kdeEvalPositions,
                             gaussianKdeDensities);

    tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
    if (m_env.rank() == 0) {
      std::cout << "Chain KDE took " << tmpRunTime
                << " seconds"
                << std::endl;
    }

    // Write estimations of probability densities
    // hold
    // plot(queso_kdeEvalPositions(1,:)',7*queso_gaussianKdeDensities(1,:)','r-');
    if (passedOfs) {
      std::ofstream& ofs = *passedOfs;
      ofs << chainName << "_kdeEvalPositions = zeros(" << m_paramSpace.dim()
          << ","                                       << kdeEvalPositions.size()
          << ");"
          << std::endl;
      for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
        for (unsigned int j = 0; j < kdeEvalPositions.size(); ++j) {
          ofs << chainName << "_kdeEvalPositions(" << i+1
              << ","                               << j+1
              << ") = "                            << (*(kdeEvalPositions[j]))[i]
              << ";"
              << std::endl;
        }
      }

      ofs << chainName << "_gaussianKdeScaleVec = zeros(" << m_paramSpace.dim()
          << ");"
          << std::endl;
      for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
        ofs << chainName << "_gaussianKdeScaleVec(" << i+1
            << ") = "                               << gaussianKdeScaleVec[i]
            << ";"
            << std::endl;
      }

      ofs << chainName << "_gaussianKdeDensities = zeros(" << m_paramSpace.dim()
          << ","                                           << gaussianKdeDensities.size()
          << ");"
          << std::endl;
      for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
        for (unsigned int j = 0; j < gaussianKdeDensities.size(); ++j) {
          ofs << chainName << "_gaussianKdeDensities(" << i+1
              << ","                                   << j+1
              << ") = "                                << (*(gaussianKdeDensities[j]))[i]
              << ";"
              << std::endl;
        }
      }
    }

    for (unsigned int i = 0; i < gaussianKdeDensities.size(); ++i) {
      if (gaussianKdeDensities[i] != NULL) delete gaussianKdeDensities[i];
    }
    for (unsigned int i = 0; i < kdeEvalPositions.size(); ++i) {
      if (kdeEvalPositions[i] != NULL) delete kdeEvalPositions[i];
    }
  }

  if (m_env.rank() == 0) {
    std::cout << "\n-----------------------------------------------------"
              << "\n Finished computing histogram and/or KDE for chain " << chainName
              << "\n-----------------------------------------------------"
              << "\n"
              << std::endl;
  }

  return;
}

template<class P_V,class P_M,class L_V,class L_M>
int
uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::writeInfo(
  const uqChainBaseClass<P_V>&     workingChain,
  const std::string&               chainName,
  const std::string&               prefixName,
  std::ofstream&                   ofs,
  const P_M*                       mahalanobisMatrix,
  bool                             applyMahalanobisInvert) const
{
  if (m_env.rank() == 0) {
    std::cout << "\n"
              << "\n-----------------------------------------------------"
              << "\n Writing extra information about the Markov chain " << chainName << " to output file ..."
              << "\n-----------------------------------------------------"
              << "\n"
              << std::endl;
  }

  int iRC = UQ_OK_RC;

  if (m_chainGenerateExtra) {
    if (m_likelihoodObjComputesMisfits) {
      // Write m_misfitChain
      ofs << prefixName << "misfitChain = zeros(" << m_misfitChain.size()
          << ","                                   << m_misfitChain[0]->size()
          << ");"
          << std::endl;
      ofs << prefixName << "misfitChain = [";
      for (unsigned int i = 0; i < m_misfitChain.size(); ++i) {
        ofs << *(m_misfitChain[i])
            << std::endl;
      }
      ofs << "];\n";

      // Write m_misfitVarianceChain
      ofs << prefixName << "misfitVarianceChain = zeros(" << m_misfitVarianceChain.size()
          << ","                                           << m_misfitVarianceChain[0]->size()
          << ");"
          << std::endl;
      ofs << prefixName << "misfitVarianceChain = [";
      for (unsigned int i = 0; i < m_misfitVarianceChain.size(); ++i) {
        ofs << *(m_misfitVarianceChain[i])
            << std::endl;
      }
      ofs << "];\n";
    }

    // Write m_m2lLikelihoodChain
    ofs << prefixName << "m2lLikelihoodChain = zeros(" << m_m2lLikelihoodChain.size()
        << ","                                          << m_m2lLikelihoodChain[0]->size()
        << ");"
        << std::endl;
    ofs << prefixName << "m2lLikelihoodChain = [";
    for (unsigned int i = 0; i < m_m2lLikelihoodChain.size(); ++i) {
      ofs << *(m_m2lLikelihoodChain[i])
          << std::endl;
    }
    ofs << "];\n";

    // Write m_alphaQuotients
    ofs << prefixName << "alphaQuotients = zeros(" << m_alphaQuotients.size()
        << ","                                      << 1
        << ");"
        << std::endl;
    ofs << prefixName << "alphaQuotients = [";
    for (unsigned int i = 0; i < m_alphaQuotients.size(); ++i) {
      ofs << m_alphaQuotients[i]
          << std::endl;
    }
    ofs << "];\n";
  }

  // Write names of parameters
  ofs << prefixName << "paramNames = {";
  m_paramSpace.printParameterNames(ofs,false);
  ofs << "};\n";

  // Write mahalanobis distances
  if (mahalanobisMatrix != NULL) {
    P_V diffVec(m_paramSpace.zeroVector());
    ofs << prefixName << "d = [";
    if (applyMahalanobisInvert) {
      P_V tmpVec(m_paramSpace.zeroVector());
      P_V vec0(m_paramSpace.zeroVector());
      workingChain.getPositionValues(0,vec0);
      for (unsigned int i = 0; i < workingChain.sequenceSize(); ++i) {
        workingChain.getPositionValues(i,tmpVec);
        diffVec = tmpVec - vec0;
        //diffVec = *(workingChain[i]) - *(workingChain[0]);
        ofs << scalarProduct(diffVec, mahalanobisMatrix->invertMultiply(diffVec))
            << std::endl;
      }
    }
    else {
      P_V tmpVec(m_paramSpace.zeroVector());
      P_V vec0(m_paramSpace.zeroVector());
      workingChain.getPositionValues(0,vec0);
      for (unsigned int i = 0; i < workingChain.sequenceSize(); ++i) {
        workingChain.getPositionValues(i,tmpVec);
        diffVec = tmpVec - vec0;
        //diffVec = *(workingChain[i]) - *(workingChain[0]);
        ofs << scalarProduct(diffVec, *mahalanobisMatrix * diffVec)
            << std::endl;
      }
    }
    ofs << "];\n";
  }

  // Write prior mean values
  ofs << prefixName << "priorMeanValues = ["
      << m_paramSpace.priorMuValues()
      << "];\n";

  // Write prior sigma values
  ofs << prefixName << "priorSigmaValues = ["
      << m_paramSpace.priorSigmaValues()
      << "];\n";

#if 0
  ofs << prefixName << "results.prior = [queso_priorMeanValues',queso_priorSigmaValues'];\n";
#endif

  // Write param lower bounds
  ofs << prefixName << "minValues = ["
      << m_paramSpace.minValues()
      << "];\n";

  // Write param upper bounds
  ofs << prefixName << "maxValues = ["
      << m_paramSpace.maxValues()
      << "];\n";

#if 0
  ofs << prefixName << "results.limits = [queso_low',queso_upp'];\n";

  // Write out data for mcmcpred.m
  ofs << prefixName << "results.parind = ["; // FIXME
  for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
    ofs << i+1
        << std::endl;
  }
  ofs << "];\n";

  ofs << prefixName << "results.local = [\n"; // FIXME
  for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
    ofs << " 0";
    //<< std::endl;
  }
  ofs << "];\n";

  if (m_chainUse2) {
  }
  else {
    bool savedVectorPrintState = workingChain[workingChain.sequenceSize()-1]->getPrintHorizontally();
    workingChain[workingChain.sequenceSize()-1]->setPrintHorizontally(false);
    ofs << prefixName << "results.theta = ["
        << *(workingChain[workingChain.sequenceSize()-1])
        << "];\n";
    workingChain[workingChain.sequenceSize()-1]->setPrintHorizontally(savedVectorPrintState);
  }
  
  ofs << prefixName << "results.nbatch = 1.;\n"; // FIXME

  if (mahalanobisMatrix != NULL) {
    // Write covMatrix
    ofs << prefixName << "mahalanobisMatrix = ["
        << *mahalanobisMatrix
        << "];\n";
  }
#endif

  // Write number of rejections
  ofs << prefixName << "rejected = " << (double) m_numRejections/(double) (workingChain.sequenceSize()-1)
      << ";\n"
      << std::endl;

  // Write number of outbounds
  ofs << prefixName << "outbounds = " << (double) m_numOutOfBounds/(double) (workingChain.sequenceSize()-1)
      << ";\n"
      << std::endl;

  // Write chain run time
  ofs << prefixName << "runTime = " << m_chainRunTime
      << ";\n"
      << std::endl;

  if (m_env.rank() == 0) {
    std::cout << "\n-----------------------------------------------------"
              << "\n Finished writing extra information about the Markov chain " << chainName
              << "\n-----------------------------------------------------"
              << "\n"
              << std::endl;
  }

  return iRC;
}

#if 0
template<class P_V,class P_M,class L_V,class L_M>
const uqSequenceOfVectorsClass<P_V>&
uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::chain() const
{
  return m_chain1;
}

template<class P_V,class P_M,class L_V,class L_M>
const std::vector<const L_V*>&
uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::misfitChain() const
{
  return m_misfitChain;
}

template<class P_V,class P_M,class L_V,class L_M>
const std::vector<const L_V*>&
uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::misfitVarianceChain() const
{
  return m_misfitVarianceChain;
}

template<class P_V,class P_M,class L_V,class L_M>
const std::string&
uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::outputFileName() const
{
  return m_chainOutputFileNames[nothing yet];
}
#endif

template<class P_V,class P_M,class L_V,class L_M>
void
uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::print(std::ostream& os) const
{
  os <<         m_option_chain_type   << " = " << m_chainType
     << "\n" << m_option_chain_number << " = " << m_chainNumber
     << "\n" << m_option_chain_sizes  << " = ";
  for (unsigned int i = 0; i < m_chainSizes.size(); ++i) {
    os << m_chainSizes[i] << " ";
  }
  os << "\n" << m_option_chain_use2            << " = " << m_chainUse2
     << "\n" << m_option_chain_generateExtra   << " = " << m_chainGenerateExtra
     << "\n" << m_option_chain_displayPeriod   << " = " << m_chainDisplayPeriod
     << "\n" << m_option_chain_measureRunTimes << " = " << m_chainMeasureRunTimes
     << "\n" << m_option_chain_write           << " = " << m_chainWrite
     << "\n" << m_option_chain_computeStats    << " = " << m_chainComputeStats
     << "\n" << m_option_chain_outputFileNames << " = ";
  for (unsigned int i = 0; i < m_chainOutputFileNames.size(); ++i) {
    os << m_chainOutputFileNames[i] << " ";
  }
//os << "\n" << m_option_uniqueChain_generate           << " = " << m_uniqueChainGenerate
//   << "\n" << m_option_uniqueChain_write              << " = " << m_uniqueChainWrite
//   << "\n" << m_option_uniqueChain_computeStats       << " = " << m_uniqueChainComputeStats
  os << "\n" << m_option_filteredChain_generate         << " = " << m_filteredChainGenerate
     << "\n" << m_option_filteredChain_discardedPortion << " = " << m_filteredChainDiscardedPortion
     << "\n" << m_option_filteredChain_lag              << " = " << m_filteredChainLag
     << "\n" << m_option_filteredChain_write            << " = " << m_filteredChainWrite
     << "\n" << m_option_filteredChain_computeStats     << " = " << m_filteredChainComputeStats
#if 0
     << "\n" << m_option_avgChain_compute               << " = ";
  for (unsigned int i = 0; i < m_avgChainCompute.size(); ++i) {
    os << m_avgChainCompute[i] << " ";
  }
  os << "\n" << m_option_avgChain_write          << " = " << m_avgChainWrite
     << "\n" << m_option_avgChain_computeStats   << " = " << m_avgChainComputeStats
#endif
     << "\n" << m_option_dr_maxNumExtraStages    << " = " << m_maxNumExtraStages
     << "\n" << m_option_dr_scalesForExtraStages << " = ";
  for (unsigned int i = 0; i < m_scalesForCovMProposals.size(); ++i) {
    os << m_scalesForCovMProposals[i] << " ";
  }
  os << "\n" << m_option_am_initialNonAdaptInterval << " = " << m_initialNonAdaptInterval
     << "\n" << m_option_am_adaptInterval           << " = " << m_adaptInterval
     << "\n" << m_option_am_eta                     << " = " << m_eta
     << "\n" << m_option_am_epsilon                 << " = " << m_epsilon
     << "\n" << "(internal variable) m_likelihoodObjComputesMisfits = " << m_likelihoodObjComputesMisfits
     << std::endl;

  return;
}

template<class P_V,class P_M,class L_V,class L_M>
std::ostream& operator<<(std::ostream& os, const uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>& obj)
{
  obj.print(os);

  return os;
}
#endif // __UQ_BMC_DC1_H__
