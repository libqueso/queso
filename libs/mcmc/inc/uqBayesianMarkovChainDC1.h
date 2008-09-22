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

#ifndef __UQ_BMCDC1_H__
#define __UQ_BMCDC1_H__

#undef  UQ_BMCDC_REQUIRES_INVERTED_COV_MATRICES
#define UQ_BMCDC_REQUIRES_TARGET_DISTRIBUTION_ONLY

#define UQ_BMCDC_MARKOV_CHAIN_TYPE           1
#define UQ_BMCDC_WHITE_NOISE_CHAIN_TYPE      2
#define UQ_BMCDC_UNIFORM_CHAIN_TYPE          3
#define UQ_BMCDC_FILENAME_FOR_NO_OUTPUT_FILE "."

// _ODV = option default value
#define UQ_BMCDC_CHAIN_TYPE_ODV                       UQ_BMCDC_MARKOV_CHAIN_TYPE
#define UQ_BMCDC_CHAIN_NUMBER_ODV                     1
#define UQ_BMCDC_CHAIN_SIZES_ODV                      "100"
#define UQ_BMCDC_CHAIN_OUTPUT_FILE_NAMES_ODV          UQ_BMCDC_FILENAME_FOR_NO_OUTPUT_FILE
#define UQ_BMCDC_CHAIN_USE2_ODV                       0
#define UQ_BMCDC_CHAIN_GENERATE_EXTRA_ODV             0
#define UQ_BMCDC_CHAIN_DISPLAY_PERIOD_ODV             500
#define UQ_BMCDC_CHAIN_MEASURE_RUN_TIMES_ODV          0
#define UQ_BMCDC_CHAIN_WRITE_ODV                      0
#define UQ_BMCDC_CHAIN_COMPUTE_STATS_ODV              0
#define UQ_BMCDC_UNIQUE_CHAIN_GENERATE_ODV            0
#define UQ_BMCDC_UNIQUE_CHAIN_WRITE_ODV               0
#define UQ_BMCDC_UNIQUE_CHAIN_COMPUTE_STATS_ODV       0
#define UQ_BMCDC_FILTERED_CHAIN_GENERATE_ODV          0
#define UQ_BMCDC_FILTERED_CHAIN_DISCARDED_PORTION_ODV 0.
#define UQ_BMCDC_FILTERED_CHAIN_LAG_ODV               1
#define UQ_BMCDC_FILTERED_CHAIN_WRITE_ODV             0
#define UQ_BMCDC_FILTERED_CHAIN_COMPUTE_STATS_ODV     0
#define UQ_BMCDC_AVG_CHAIN_COMPUTE_ODV                "0"
#define UQ_BMCDC_AVG_CHAIN_WRITE_ODV                  0
#define UQ_BMCDC_AVG_CHAIN_COMPUTE_STATS_ODV          0
#define UQ_BMCDC_DR_MAX_NUM_EXTRA_STAGES_ODV          0
#define UQ_BMCDC_DR_SCALES_FOR_EXTRA_STAGES_ODV       "1."
#define UQ_BMCDC_AM_INIT_NON_ADAPT_INT_ODV            0
#define UQ_BMCDC_AM_ADAPT_INTERVAL_ODV                0
#define UQ_BMCDC_AM_ETA_ODV                           1.
#define UQ_BMCDC_AM_EPSILON_ODV                       1.e-5

#include <uqChainStatisticalOptions.h>
#include <uqProbDensity.h>
#include <uqParamSpace.h>
#ifdef UQ_BMCDC_REQUIRES_TARGET_DISTRIBUTION_ONLY
#else
#include <uqVectorLhFunction.h>
#include <uqObservableSpace.h>
#endif
#include <uqChainPosition.h>
#include <uqMiscellaneous.h>
#include <uqSequenceOfVectors.h>
#include <uqArrayOfSequences.h>
#include <sys/time.h>
#include <fstream>

/*! A templated class that implements a Bayesian Markov Chain Distribution Calculator
 */
template <class P_V,class P_M,class L_V,class L_M>
class uqBayesianMarkovChainDCClass
{
public:
  uqBayesianMarkovChainDCClass(const uqEnvironmentClass&                             env,                      /*! The QUESO toolkit environment. */
                               const char*                                           prefix,                   /*! Prefix.                        */
                               const uqParamSpaceClass            <P_V,P_M>&         paramSpace,               /*! The parameter space.           */
#ifdef UQ_BMCDC_REQUIRES_TARGET_DISTRIBUTION_ONLY
                               const uqBayesProbDensity_BaseClass <P_V,P_M>&         targetParamDensityObj,
#else
                               const uqObservableSpaceClass       <L_V,L_M>&         observableSpace,          /*! The observable space.          */
                               const uqProbDensity_BaseClass      <P_V,P_M>&         m2lPriorParamDensityObj,  /*! -2*ln(prior()).                */
                               const uqVectorLhFunction_BaseClass <P_V,P_M,L_V,L_M>& m2lVectorLhFunctionObj,   /*! -2*ln(likelihood()).           */
#endif
                                     P_M*                                            proposalCovMatrix,        /*! */
                               const uqProposalDensity_BaseClass  <P_V,P_M>*         proposalDensityObj,       /*! */
                               const uqProposalGenerator_BaseClass<P_V,P_M>*         proposalGeneratorObj);    /*! */
 ~uqBayesianMarkovChainDCClass();

  void calculateDistributions ();
                            //(const P_M* proposalCovMatrix;
                            // const P_M* proposalPrecMatrix,
                            // const P_M* mahalanobisMatrix = NULL,
                            // bool       applyMahalanobisInvert = true);
  void calculateDistributions (const uqProbDensity_BaseClass<P_V,P_M>& priorParamDensityObj);

  void print                  (std::ostream& os) const;

  const uqChainBaseClass<P_V>&         chain              () const;
//const std::vector<const L_V*>&       misfitChain        () const;
//const std::vector<const L_V*>&       misfitVarianceChain() const;
//const std::string&                   outputFileName     () const;

private:
  void   resetChainAndRelatedInfo();
  void   defineMyOptions         (po::options_description&                       optionsDesc);
  void   getMyOptionValues       (po::options_description&                       optionsDesc);

  int    prepareForNextChain     (const P_M*                                     proposalCovMatrix);
                                //const P_M*                                     proposalPrecMatrix,

  void   calculateDistributions  (const P_M*                                     proposalCovMatrix,
                                //const P_M*                                     proposalPrecMatrix,
                                //const P_M*                                     mahalanobisMatrix,
                                //bool                                           applyMahalanobisInvert,
                                  uqChainBaseClass<P_V>&                         workingChain);
  void   generateMarkovChain     (unsigned int                                   chainSize,
                                  const P_V&                                     valuesOf1stPosition,
                                  const P_M*                                     proposalCovMatrix,
                                  uqChainBaseClass<P_V>&                         workingChain,
                                  const std::string&                             chainName);
  void   generateWhiteNoiseChain (unsigned int                                   chainSize,
                                  uqChainBaseClass<P_V>&                         workingChain,
                                  const std::string&                             chainName);
  void   generateUniformChain    (unsigned int                                   chainSize,
                                  uqChainBaseClass<P_V>&                         workingChain,
                                  const std::string&                             chainName);
  void   updateCovMatrix         (const uqChainBaseClass<P_V>&                   subChain,
                                  unsigned int                                   idOfFirstPositionInSubChain,
                                  double&                                        lastChainSize,
                                  P_V&                                           lastMean,
                                  P_M&                                           lastAdaptedCovMatrix);

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

  int    writeInfo               (const uqChainBaseClass<P_V>&                   workingChain,
                                  const std::string&                             chainName,
                                  const std::string&                             prefixName,
                                  std::ofstream&                                 ofs) const;
                                //const P_M*                                     mahalanobisMatrix = NULL,
                                //bool                                           applyMahalanobisInvert = true) const;

  const uqEnvironmentClass&                             m_env;
        std::string                                     m_prefix;
  const uqParamSpaceClass            <P_V,P_M>&         m_paramSpace;
#ifdef UQ_BMCDC_REQUIRES_TARGET_DISTRIBUTION_ONLY
  const uqBayesProbDensity_BaseClass <P_V,P_M>&         m_targetParamDensityObj;
#else
  const uqObservableSpaceClass       <L_V,L_M>&         m_observableSpace;
  const uqProbDensity_BaseClass      <P_V,P_M>&         m_m2lPriorParamDensityObj;
  const uqVectorLhFunction_BaseClass <P_V,P_M,L_V,L_M>& m_m2lVectorLhFunctionObj;
#endif
        P_M*                                            m_proposalCovMatrix;
  const uqProposalDensity_BaseClass  <P_V,P_M>*         m_proposalDensityObj;
  const uqProposalGenerator_BaseClass<P_V,P_M>*         m_proposalGeneratorObj;

  po::options_description*        m_optionsDesc;
  std::string                     m_option_help;
  std::string                     m_option_chain_type;
  std::string                     m_option_chain_number;
  std::string                     m_option_chain_sizes;
  std::string                     m_option_chain_outputFileNames;
  std::string                     m_option_chain_use2;
  std::string                     m_option_chain_generateExtra;
  std::string                     m_option_chain_displayPeriod;
  std::string                     m_option_chain_measureRunTimes;
  std::string                     m_option_chain_write;
  std::string                     m_option_chain_computeStats;
  std::string                     m_option_uniqueChain_generate;
  std::string                     m_option_uniqueChain_write;
  std::string                     m_option_uniqueChain_computeStats;
  std::string                     m_option_filteredChain_generate;
  std::string                     m_option_filteredChain_discardedPortion;
  std::string                     m_option_filteredChain_lag;
  std::string                     m_option_filteredChain_write;
  std::string                     m_option_filteredChain_computeStats;
  std::string                     m_option_avgChain_compute;
  std::string                     m_option_avgChain_write;
  std::string                     m_option_avgChain_computeStats;
  std::string                     m_option_dr_maxNumExtraStages;
  std::string                     m_option_dr_scalesForExtraStages;
  std::string                     m_option_am_initialNonAdaptInterval;
  std::string                     m_option_am_adaptInterval;
  std::string                     m_option_am_eta;
  std::string                     m_option_am_epsilon;

#ifdef UQ_BMCDC_REQUIRES_TARGET_DISTRIBUTION_ONLY
#else
  bool                            m_likelihoodObjComputesMisfits;
#endif
  P_V                             m_paramInitials;
  bool                            m_proposalIsSymmetric;

  unsigned int                    m_chainType;
  unsigned int                    m_chainNumber;
  std::vector<unsigned int>       m_chainSizes;
  std::vector<std::string>        m_chainOutputFileNames;
  bool                            m_chainUse2;
  bool                            m_chainGenerateExtra;
  unsigned int                    m_chainDisplayPeriod;
  bool                            m_chainMeasureRunTimes;
  bool                            m_chainWrite;
  bool                            m_chainComputeStats;
  uqChainStatisticalOptionsClass* m_chainStatisticalOptions;

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
#ifdef UQ_BMCDC_REQUIRES_INVERTED_COV_MATRICES
  std::vector<      P_M*>         m_upperCholProposalPrecMatrices;
  std::vector<      P_M*>         m_proposalPrecMatrices;
#endif

  uqSequenceOfVectorsClass<P_V>   m_chain1;
  uqArrayOfSequencesClass<P_V>    m_chain2;
  std::vector<unsigned int>       m_idsOfUniquePositions;
#ifdef UQ_BMCDC_REQUIRES_TARGET_DISTRIBUTION_ONLY
#else
  std::vector<const L_V*>         m_misfitChain; // Sum of squares of differences between model and experiments: computed by user supplied likelihood obj
  std::vector<const L_V*>         m_misfitVarianceChain;
  std::vector<const L_V*>         m_m2lLikelihoodChain;
#endif
  std::vector<double>             m_alphaQuotients;
  double                          m_chainRunTime;
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
  const uqEnvironmentClass&                             env,
  const char*                                           prefix,
  const uqParamSpaceClass            <P_V,P_M>&         paramSpace,
#ifdef UQ_BMCDC_REQUIRES_TARGET_DISTRIBUTION_ONLY
  const uqBayesProbDensity_BaseClass <P_V,P_M>&         targetParamDensityObj,
#else
  const uqObservableSpaceClass       <L_V,L_M>&         observableSpace,
  const uqProbDensity_BaseClass      <P_V,P_M>&         m2lPriorParamDensityObj,
  const uqVectorLhFunction_BaseClass <P_V,P_M,L_V,L_M>& m2lVectorLhFunctionObj,
#endif
        P_M*                                            proposalCovMatrix,
  const uqProposalDensity_BaseClass  <P_V,P_M>*         proposalDensityObj,
  const uqProposalGenerator_BaseClass<P_V,P_M>*         proposalGeneratorObj)
  :
  m_env                                  (env),
  m_prefix                               ((std::string)(prefix) + "bmcdc_"),
  m_paramSpace                           (paramSpace),
#ifdef UQ_BMCDC_REQUIRES_TARGET_DISTRIBUTION_ONLY
  m_targetParamDensityObj                (targetParamDensityObj),
#else
  m_observableSpace                      (observableSpace),
  m_m2lPriorParamDensityObj              (m2lPriorParamDensityObj),
  m_m2lVectorLhFunctionObj               (m2lVectorLhFunctionObj),
#endif
  m_proposalCovMatrix                    (proposalCovMatrix),
  m_proposalDensityObj                   (proposalDensityObj),
  m_proposalGeneratorObj                 (proposalGeneratorObj),
  m_optionsDesc                          (new po::options_description("Bayesian Markov chain options")),
  m_option_help                          (m_prefix + "help"                          ),
  m_option_chain_type                    (m_prefix + "chain_type"                    ),
  m_option_chain_number                  (m_prefix + "chain_number"                  ),
  m_option_chain_sizes                   (m_prefix + "chain_sizes"                   ),
  m_option_chain_outputFileNames         (m_prefix + "chain_outputFileNames"         ),
  m_option_chain_use2                    (m_prefix + "chain_use2"                    ),
  m_option_chain_generateExtra           (m_prefix + "chain_generateExtra"           ),
  m_option_chain_displayPeriod           (m_prefix + "chain_displayPeriod"           ),
  m_option_chain_measureRunTimes         (m_prefix + "chain_measureRunTimes"         ),
  m_option_chain_write                   (m_prefix + "chain_write"                   ),
  m_option_chain_computeStats            (m_prefix + "chain_computeStats"            ),
  m_option_uniqueChain_generate          (m_prefix + "uniqueChain_generate"          ),
  m_option_uniqueChain_write             (m_prefix + "uniqueChain_write"             ),
  m_option_uniqueChain_computeStats      (m_prefix + "uniqueChain_computeStats"      ),
  m_option_filteredChain_generate        (m_prefix + "filteredChain_generate"        ),
  m_option_filteredChain_discardedPortion(m_prefix + "filteredChain_discardedPortion"),
  m_option_filteredChain_lag             (m_prefix + "filteredChain_lag"             ),
  m_option_filteredChain_write           (m_prefix + "filteredChain_write"           ),
  m_option_filteredChain_computeStats    (m_prefix + "filteredChain_computeStats"    ),
  m_option_avgChain_compute              (m_prefix + "avgChain_compute"              ),
  m_option_avgChain_write                (m_prefix + "avgChain_write"                ),
  m_option_avgChain_computeStats         (m_prefix + "avgChain_computeStats"         ),
  m_option_dr_maxNumExtraStages          (m_prefix + "dr_maxNumExtraStages"          ),
  m_option_dr_scalesForExtraStages       (m_prefix + "dr_scalesForExtraStages"       ),
  m_option_am_initialNonAdaptInterval    (m_prefix + "am_initialNonAdaptInterval"    ),
  m_option_am_adaptInterval              (m_prefix + "am_adaptInterval"              ),
  m_option_am_eta                        (m_prefix + "am_eta"                        ),
  m_option_am_epsilon                    (m_prefix + "am_epsilon"                    ),
#ifdef UQ_BMCDC_REQUIRES_TARGET_DISTRIBUTION_ONLY
#else
  m_likelihoodObjComputesMisfits         (dynamic_cast<const uqMisfitVectorLhFunction_Class<P_V,P_M,L_V,L_M>*>(&m2lVectorLhFunctionObj) != NULL),
#endif
  m_paramInitials                        (m_paramSpace.initialValues()),
  m_proposalIsSymmetric                  (true),
  m_chainType                            (UQ_BMCDC_CHAIN_TYPE_ODV),
  m_chainNumber                          (UQ_BMCDC_CHAIN_NUMBER_ODV),
  m_chainSizes                           (1,(unsigned int) strtod(UQ_BMCDC_CHAIN_SIZES_ODV,NULL)),
  m_chainOutputFileNames                 (1,UQ_BMCDC_CHAIN_OUTPUT_FILE_NAMES_ODV),
  m_chainUse2                            (UQ_BMCDC_CHAIN_USE2_ODV),
  m_chainGenerateExtra                   (UQ_BMCDC_CHAIN_GENERATE_EXTRA_ODV),
  m_chainDisplayPeriod                   (UQ_BMCDC_CHAIN_DISPLAY_PERIOD_ODV),
  m_chainMeasureRunTimes                 (UQ_BMCDC_CHAIN_MEASURE_RUN_TIMES_ODV),
  m_chainWrite                           (UQ_BMCDC_CHAIN_WRITE_ODV),
  m_chainComputeStats                    (UQ_BMCDC_CHAIN_COMPUTE_STATS_ODV),
  m_chainStatisticalOptions              (NULL),
  m_uniqueChainGenerate                  (UQ_BMCDC_UNIQUE_CHAIN_GENERATE_ODV),
  m_uniqueChainWrite                     (UQ_BMCDC_UNIQUE_CHAIN_WRITE_ODV),
  m_uniqueChainComputeStats              (UQ_BMCDC_UNIQUE_CHAIN_COMPUTE_STATS_ODV),
  m_uniqueChainStatisticalOptions        (NULL),
  m_filteredChainGenerate                (UQ_BMCDC_FILTERED_CHAIN_GENERATE_ODV),
  m_filteredChainDiscardedPortion        (UQ_BMCDC_FILTERED_CHAIN_DISCARDED_PORTION_ODV),
  m_filteredChainLag                     (UQ_BMCDC_FILTERED_CHAIN_LAG_ODV),
  m_filteredChainWrite                   (UQ_BMCDC_FILTERED_CHAIN_WRITE_ODV),
  m_filteredChainComputeStats            (UQ_BMCDC_FILTERED_CHAIN_COMPUTE_STATS_ODV),
  m_filteredChainStatisticalOptions      (NULL),
  m_avgChainCompute                      (0),//,0.),
  m_avgChainWrite                        (UQ_BMCDC_AVG_CHAIN_WRITE_ODV),
  m_avgChainComputeStats                 (UQ_BMCDC_AVG_CHAIN_COMPUTE_STATS_ODV),
  m_maxNumExtraStages                    (UQ_BMCDC_DR_MAX_NUM_EXTRA_STAGES_ODV),
  m_scalesForCovMProposals               (0),//,0.),
  m_initialNonAdaptInterval              (UQ_BMCDC_AM_INIT_NON_ADAPT_INT_ODV),
  m_adaptInterval                        (UQ_BMCDC_AM_ADAPT_INTERVAL_ODV),
  m_eta                                  (UQ_BMCDC_AM_ETA_ODV),
  m_epsilon                              (UQ_BMCDC_AM_EPSILON_ODV),
  m_lowerCholProposalCovMatrices         (1),//,NULL),
  m_proposalCovMatrices                  (1),//,NULL),
#ifdef UQ_BMCDC_REQUIRES_INVERTED_COV_MATRICES
  m_upperCholProposalPrecMatrices        (1),//,NULL),
  m_proposalPrecMatrices                 (1),//,NULL),
#endif
  m_chain1                               (0,m_paramSpace.zeroVector()),
  m_chain2                               (0,m_paramSpace.zeroVector()),
  m_idsOfUniquePositions                 (0),//0),
#ifdef UQ_BMCDC_REQUIRES_TARGET_DISTRIBUTION_ONLY
#else
  m_misfitChain                          (0),//,NULL),
  m_misfitVarianceChain                  (0),//,NULL),
  m_m2lLikelihoodChain                   (0),//,NULL),
#endif
  m_alphaQuotients                       (0),//,0.),
  m_chainRunTime                         (0.),
  m_numRejections                        (0),
  m_numOutOfBounds                       (0),
  m_lastChainSize                        (0),
  m_lastMean                             (NULL),
  m_lastAdaptedCovMatrix                 (NULL)
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

  if (m_chainComputeStats        ) m_chainStatisticalOptions         = new uqChainStatisticalOptionsClass(m_env,m_prefix + "chain_"        );
  if (m_uniqueChainComputeStats  ) m_uniqueChainStatisticalOptions   = new uqChainStatisticalOptionsClass(m_env,m_prefix + "uniqueChain_"  );
  if (m_filteredChainComputeStats) m_filteredChainStatisticalOptions = new uqChainStatisticalOptionsClass(m_env,m_prefix + "filteredChain_");

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
  m_lastChainSize  = 0;
  m_numOutOfBounds = 0;
  m_chainRunTime   = 0.;
  m_numRejections  = 0;
  m_alphaQuotients.clear();
#ifdef UQ_BMCDC_REQUIRES_TARGET_DISTRIBUTION_ONLY
#else
  for (unsigned int i = 0; i < m_m2lLikelihoodChain.size(); ++i) {
    if (m_m2lLikelihoodChain[i]) delete m_m2lLikelihoodChain[i];
  }
  for (unsigned int i = 0; i < m_misfitVarianceChain.size(); ++i) {
    if (m_misfitVarianceChain[i]) delete m_misfitVarianceChain[i];
  }
  for (unsigned int i = 0; i < m_misfitChain.size(); ++i) {
    if (m_misfitChain[i]) delete m_misfitChain[i];
  }
#endif
  m_idsOfUniquePositions.clear();
  m_chain1.clear();
  m_chain2.clear();

#ifdef UQ_BMCDC_REQUIRES_INVERTED_COV_MATRICES
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
  optionsDesc.add_options()
    (m_option_help.c_str(),                                                                                                                     "produce help message for Bayesian Markov chain distr. calculator")
    (m_option_chain_type.c_str(),                     po::value<unsigned int>()->default_value(UQ_BMCDC_CHAIN_TYPE_ODV                      ), "type of chain (1=Markov, 2=White noise)"                         )
    (m_option_chain_number.c_str(),                   po::value<unsigned int>()->default_value(UQ_BMCDC_CHAIN_NUMBER_ODV                    ), "number of chain(s)"                                              )
    (m_option_chain_sizes.c_str(),                    po::value<std::string >()->default_value(UQ_BMCDC_CHAIN_SIZES_ODV                     ), "list of size(s) of chain(s)"                                     )
    (m_option_chain_outputFileNames.c_str(),          po::value<std::string >()->default_value(UQ_BMCDC_CHAIN_OUTPUT_FILE_NAMES_ODV         ), "list of name(s) of output file(s)"                               )
    (m_option_chain_use2.c_str(),                     po::value<bool        >()->default_value(UQ_BMCDC_CHAIN_USE2_ODV                      ), "use chain2"                                                      )
    (m_option_chain_generateExtra.c_str(),            po::value<bool        >()->default_value(UQ_BMCDC_CHAIN_GENERATE_EXTRA_ODV            ), "generate extra chains"                                           )
    (m_option_chain_displayPeriod.c_str(),            po::value<unsigned int>()->default_value(UQ_BMCDC_CHAIN_DISPLAY_PERIOD_ODV            ), "period of message display during chain generation"               )
    (m_option_chain_measureRunTimes.c_str(),          po::value<bool        >()->default_value(UQ_BMCDC_CHAIN_MEASURE_RUN_TIMES_ODV         ), "measure run times"                                               )
    (m_option_chain_write.c_str(),                    po::value<bool        >()->default_value(UQ_BMCDC_CHAIN_WRITE_ODV                     ), "write chain values to the output file"                           )
    (m_option_chain_computeStats.c_str(),             po::value<bool        >()->default_value(UQ_BMCDC_CHAIN_COMPUTE_STATS_ODV             ), "compute statistics on chain"                                     )
  //(m_option_uniqueChain_generate.c_str(),           po::value<bool        >()->default_value(UQ_BMCDC_UNIQUE_CHAIN_GENERATE_ODV           ), "generate unique chain"                                           )
  //(m_option_uniqueChain_write.c_str(),              po::value<bool        >()->default_value(UQ_BMCDC_UNIQUE_CHAIN_WRITE_ODV              ), "write unique chain"                                              )
  //(m_option_uniqueChain_computeStats.c_str(),       po::value<bool        >()->default_value(UQ_BMCDC_UNIQUE_CHAIN_COMPUTE_STATS_ODV      ), "compute statistics on unique chain"                              )
    (m_option_filteredChain_generate.c_str(),         po::value<bool        >()->default_value(UQ_BMCDC_FILTERED_CHAIN_GENERATE_ODV         ), "generate filtered chain"                                         )
    (m_option_filteredChain_discardedPortion.c_str(), po::value<double      >()->default_value(UQ_BMCDC_FILTERED_CHAIN_DISCARDED_PORTION_ODV), "initial discarded portion for chain filtering"                   )
    (m_option_filteredChain_lag.c_str(),              po::value<unsigned int>()->default_value(UQ_BMCDC_FILTERED_CHAIN_LAG_ODV              ), "spacing for chain filtering"                                     )
    (m_option_filteredChain_write.c_str(),            po::value<bool        >()->default_value(UQ_BMCDC_FILTERED_CHAIN_WRITE_ODV            ), "write filtered chain"                                            )
    (m_option_filteredChain_computeStats.c_str(),     po::value<bool        >()->default_value(UQ_BMCDC_FILTERED_CHAIN_COMPUTE_STATS_ODV    ), "compute statistics on filtered chain"                            )
  //(m_option_avgChain_compute.c_str(),               po::value<std::string >()->default_value(UQ_BMCDC_AVG_CHAIN_COMPUTE_ODV               ), "list of amounts of chains involved in chain averages"            )
  //(m_option_avgChain_write.c_str(),                 po::value<bool        >()->default_value(UQ_BMCDC_AVG_CHAIN_WRITE_ODV                 ), "write averages of chains"                                        )
  //(m_option_avgChain_computeStats.c_str(),          po::value<bool        >()->default_value(UQ_BMCDC_AVG_CHAIN_COMPUTE_STATS_ODV         ), "compute statistics on the averages of chains"                    )
    (m_option_dr_maxNumExtraStages.c_str(),           po::value<unsigned int>()->default_value(UQ_BMCDC_DR_MAX_NUM_EXTRA_STAGES_ODV         ), "'dr' maximum number of extra stages"                             )
    (m_option_dr_scalesForExtraStages.c_str(),        po::value<std::string >()->default_value(UQ_BMCDC_DR_SCALES_FOR_EXTRA_STAGES_ODV      ), "'dr' list of scales for proposal cov matrices from 2nd stage on" )
    (m_option_am_initialNonAdaptInterval.c_str(),     po::value<unsigned int>()->default_value(UQ_BMCDC_AM_INIT_NON_ADAPT_INT_ODV           ), "'am' initial non adaptation interval"                            )
    (m_option_am_adaptInterval.c_str(),               po::value<unsigned int>()->default_value(UQ_BMCDC_AM_ADAPT_INTERVAL_ODV               ), "'am' adaptation interval"                                        )
    (m_option_am_eta.c_str(),                         po::value<double      >()->default_value(UQ_BMCDC_AM_ETA_ODV                          ), "'am' eta"                                                        )
    (m_option_am_epsilon.c_str(),                     po::value<double      >()->default_value(UQ_BMCDC_AM_EPSILON_ODV                      ), "'am' epsilon"                                                    )
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
uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::calculateDistributions()
//const P_M* proposalCovMatrix,
//const P_M* mahalanobisMatrix,
//bool       applyMahalanobisInvert)
{
  // FIX ME: complement code logic

  if (m_chainUse2) {
    calculateDistributions(m_proposalCovMatrix,
                         //mahalanobisMatrix,
                         //applyMahalanobisInvert,
                           m_chain2);
  }
  else {
    calculateDistributions(m_proposalCovMatrix,
                         //mahalanobisMatrix,
                         //applyMahalanobisInvert,
                           m_chain1);
  }

  return;
}

template<class P_V,class P_M,class L_V,class L_M>
void
uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::calculateDistributions(const uqProbDensity_BaseClass<P_V,P_M>& priorParamDensityObj)
{
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

#ifdef UQ_BMCDC_REQUIRES_INVERTED_COV_MATRICES
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
#ifdef UQ_BMCDC_REQUIRES_INVERTED_COV_MATRICES
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
#ifdef UQ_BMCDC_REQUIRES_INVERTED_COV_MATRICES
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
  double*                          alphaQuotientPtr)
{
  double alphaQuotient = 0.;
  bool xOutOfBounds = x.outOfBounds();
  bool yOutOfBounds = y.outOfBounds();
  if ((xOutOfBounds == false) &&
      (yOutOfBounds == false)) {
    double yLogPosteriorToUse = y.logPosterior();
#ifdef UQ_BMCDC_REQUIRES_TARGET_DISTRIBUTION_ONLY
#else
    if (m_likelihoodObjComputesMisfits &&
        m_observableSpace.shouldVariancesBeUpdated()) {
      // Divide the misfitVector of 'y' by the misfitVarianceVector of 'x'
      yLogPosteriorToUse = -0.5 * ( y.m2lPrior() + (y.misfitVector()/x.misfitVarianceVector()).sumOfComponents() );
    }
#endif
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
#ifdef UQ_BMCDC_REQUIRES_TARGET_DISTRIBUTION_ONLY
#else
  if (m_likelihoodObjComputesMisfits &&
      m_observableSpace.shouldVariancesBeUpdated()) {
    // Divide the misfitVector of 'back[0]' by the misfitVarianceVector of 'pos[0]'
    numeratorLogPosteriorToUse = -0.5 * ( backwardPositions[0]->m2lPrior() +
      (backwardPositions[0]->misfitVector()/positions[0]->misfitVarianceVector()).sumOfComponents() );
  }
#endif
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
int
uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::writeInfo(
  const uqChainBaseClass<P_V>&     workingChain,
  const std::string&               chainName,
  const std::string&               prefixName,
  std::ofstream&                   ofs) const
//const P_M*                       mahalanobisMatrix,
//bool                             applyMahalanobisInvert) const
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
#ifdef UQ_BMCDC_REQUIRES_TARGET_DISTRIBUTION_ONLY
#else
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
        << ","                                         << m_m2lLikelihoodChain[0]->size()
        << ");"
        << std::endl;
    ofs << prefixName << "m2lLikelihoodChain = [";
    for (unsigned int i = 0; i < m_m2lLikelihoodChain.size(); ++i) {
      ofs << *(m_m2lLikelihoodChain[i])
          << std::endl;
    }
    ofs << "];\n";
#endif
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
#if 0
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
#endif
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

template<class P_V,class P_M,class L_V,class L_M>
const uqChainBaseClass<P_V>&
uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>::chain() const
{
  if (m_chainUse2) return m_chain2;
  return m_chain1;
}

#if 0
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
  os << "\n" << m_option_chain_outputFileNames << " = ";
  for (unsigned int i = 0; i < m_chainOutputFileNames.size(); ++i) {
    os << m_chainOutputFileNames[i] << " ";
  }
  os << "\n" << m_option_chain_use2            << " = " << m_chainUse2
     << "\n" << m_option_chain_generateExtra   << " = " << m_chainGenerateExtra
     << "\n" << m_option_chain_displayPeriod   << " = " << m_chainDisplayPeriod
     << "\n" << m_option_chain_measureRunTimes << " = " << m_chainMeasureRunTimes
     << "\n" << m_option_chain_write           << " = " << m_chainWrite
     << "\n" << m_option_chain_computeStats    << " = " << m_chainComputeStats;
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
#ifdef UQ_BMCDC_REQUIRES_TARGET_DISTRIBUTION_ONLY
#else
     << "\n" << "(internal variable) m_likelihoodObjComputesMisfits = " << m_likelihoodObjComputesMisfits
#endif
     << std::endl;

  return;
}

template<class P_V,class P_M,class L_V,class L_M>
std::ostream& operator<<(std::ostream& os, const uqBayesianMarkovChainDCClass<P_V,P_M,L_V,L_M>& obj)
{
  obj.print(os);

  return os;
}
#endif // __UQ_BMCDC1_H__
