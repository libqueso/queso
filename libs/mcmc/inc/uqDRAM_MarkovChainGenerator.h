#ifndef __UQ_DRAM_MCG_H__
#define __UQ_DRAM_MCG_H__

#undef UQ_DRAM_MCG_REQUIRES_INVERTED_COV_MATRICES

// _ODV = option default value
//   LR = likelihood results
#define UQ_DRAM_MCG_MH_CHAIN_SIZE_ODV              "100"
#define UQ_DRAM_MCG_MH_LR_UPDATE_SIGMA2_ODV        false
#define UQ_DRAM_MCG_MH_LR_SIGMA2_PRIORS_ODV        "1."
#define UQ_DRAM_MCG_MH_LR_SIGMA2_ACCURACIES_ODV    "1."
#define UQ_DRAM_MCG_MH_LR_NUMBERS_OF_OBS_ODV       "1."
#define UQ_DRAM_MCG_DR_MAX_NUM_EXTRA_STAGES_ODV    0
#define UQ_DRAM_MCG_DR_SCALES_FOR_EXTRA_STAGES_ODV "1."
#define UQ_DRAM_MCG_AM_INIT_NON_ADAPT_INT_ODV      0
#define UQ_DRAM_MCG_AM_ADAPT_INTERVAL_ODV          0
#define UQ_DRAM_MCG_AM_SD_ODV                      1.
#define UQ_DRAM_MCG_AM_EPSILON_ODV                 1.e-5
#define UQ_DRAM_MCG_MH_OUTPUT_FILE_NAME_ODV        "."
#define UQ_DRAM_MCG_MH_CHAIN_DISPLAY_PERIOD_ODV    500

#include <uqPriorAndLikelihoodInterfaces.h>
#include <uqParamSpace.h>
#include <uqStateSpace.h>
#include <uqOutputSpace.h>
#include <uqChainPosition.h>
#include <uqMiscellaneous.h>
#include <uqChainStats.h>
#include <fstream>

template <class V, class M>
class uqDRAM_MarkovChainGeneratorClass
{
public:
  uqDRAM_MarkovChainGeneratorClass(const uqEnvironmentClass&                  env,
                                   const uqParamSpaceClass<V,M>&              paramSpace,
                                   const uqFinDimLinearSpaceClass<V,M>&       outputSpace,
                                         uq_M2lPriorFunction_Class<V,M>*      m2lPriorFunction_ObjPtr,
                                         uq_M2lLikelihoodFunction_Class<V,M>& m2lLikelihoodFunction_Obj);
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
  const uqParamSpaceClass<V,M>&              m_paramSpace;
  const uqFinDimLinearSpaceClass<V,M>&       m_outputSpace;
        uq_M2lPriorFunction_Class<V,M>*      m_m2lPriorFunction_ObjPtr;
        uq_M2lLikelihoodFunction_Class<V,M>& m_m2lLikelihoodFunction_Obj;

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
  double                    m_sd;
  double                    m_epsilon;
  std::vector<std::string>  m_namesOfOutputFiles;
  unsigned int              m_chainDisplayPeriod;

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
  const uqParamSpaceClass<V,M>&              paramSpace,
  const uqFinDimLinearSpaceClass<V,M>&       outputSpace,
        uq_M2lPriorFunction_Class<V,M>*      m2lPriorFunction_ObjPtr,
        uq_M2lLikelihoodFunction_Class<V,M>& m2lLikelihoodFunction_Obj)
  :
  m_env                          (env),
  m_paramSpace                   (paramSpace),
  m_outputSpace                  (outputSpace),
  m_m2lPriorFunction_ObjPtr      (m2lPriorFunction_ObjPtr),
  m_m2lLikelihoodFunction_Obj    (m2lLikelihoodFunction_Obj),
  m_useDefaultM2lPriorFunction   (m2lPriorFunction_ObjPtr == NULL),
  m_paramInitials                (m_paramSpace.initialValues()),
  m_proposalIsSymmetric          (true),
  m_optionsDesc                  (new po::options_description("Markov chain Monte Carlo options")),
  m_sizesOfChains                (1,(unsigned int) strtod(UQ_DRAM_MCG_MH_CHAIN_SIZE_ODV,NULL)),
  m_lrUpdateSigma2               (UQ_DRAM_MCG_MH_LR_UPDATE_SIGMA2_ODV),
  m_lrSigma2Priors               (m_outputSpace.newVector()),
  m_lrSigma2Accuracies           (m_outputSpace.newVector()),
  m_lrNumbersOfObservations      (m_outputSpace.newVector()),
  m_maxNumberOfExtraStages       (UQ_DRAM_MCG_DR_MAX_NUM_EXTRA_STAGES_ODV),
  m_scalesForCovMProposals       (0,0.),
  m_initialNonAdaptInterval      (UQ_DRAM_MCG_AM_INIT_NON_ADAPT_INT_ODV),
  m_adaptInterval                (UQ_DRAM_MCG_AM_ADAPT_INTERVAL_ODV),
  m_sd                           (UQ_DRAM_MCG_AM_SD_ODV),
  m_epsilon                      (UQ_DRAM_MCG_AM_EPSILON_ODV),
  m_namesOfOutputFiles           (1,UQ_DRAM_MCG_MH_OUTPUT_FILE_NAME_ODV),
  m_chainDisplayPeriod           (UQ_DRAM_MCG_MH_CHAIN_DISPLAY_PERIOD_ODV),
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
  m_numRejections                (0),
  m_numOutOfBounds               (0),
  m_lastChainSize                (0),
  m_lastMean                     (NULL),
  m_lastAdaptedCovMatrix         (NULL)
{
  if (m_useDefaultM2lPriorFunction) {
    m_m2lPriorFunction_ObjPtr = new uq_M2lPriorFunction_Class<V,M>();
    m_defaultM2lPriorFunction_Data.paramPriorMus    = new V(m_paramSpace.priorMuValues   ());
    m_defaultM2lPriorFunction_Data.paramPriorSigmas = new V(m_paramSpace.priorSigmaValues());
  }
  m_lrSigma2Priors->cwSet         (strtod(UQ_DRAM_MCG_MH_LR_SIGMA2_PRIORS_ODV    ,NULL));
  m_lrSigma2Accuracies->cwSet     (strtod(UQ_DRAM_MCG_MH_LR_SIGMA2_ACCURACIES_ODV,NULL));
  m_lrNumbersOfObservations->cwSet(strtod(UQ_DRAM_MCG_MH_LR_NUMBERS_OF_OBS_ODV   ,NULL));

  //std::cout << "Entering uqDRAM_MarkovChainGeneratorClass<V,M>::constructor()"
  //          << std::endl;

  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.rank() == 0) std::cout << "After getting option values, state of uqDRAM_MarkovChainGeneratorClass object is:"
                                   << "\n" << *this
                                   << std::endl;

  //std::cout << "Leaving uqDRAM_MarkovChainGeneratorClass<V,M>::constructor()"
  //          << std::endl;
}

template <class V, class M>
uqDRAM_MarkovChainGeneratorClass<V,M>::~uqDRAM_MarkovChainGeneratorClass()
{
  //std::cout << "Entering uqDRAM_MarkovChainGeneratorClass<V,M>::destructor()"
  //          << std::endl;

  resetChainAndRelatedInfo();

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
    ("uqDRAM_help",                                                                                                             "produce help message for DRAM Markov chain generator"   )
    ("uqDRAM_mh_sizesOfChains",           po::value<std::string >()->default_value(UQ_DRAM_MCG_MH_CHAIN_SIZE_ODV             ), "'mh' size(s) of chain(s)"                               )
    ("uqDRAM_mh_lrUpdateSigma2",          po::value<bool        >()->default_value(UQ_DRAM_MCG_MH_LR_UPDATE_SIGMA2_ODV       ), "'mh' update the sigma2 of likelihood results"           )
    ("uqDRAM_mh_lrSigma2Priors",          po::value<std::string >()->default_value(UQ_DRAM_MCG_MH_LR_SIGMA2_PRIORS_ODV       ), "'mh' priors of sigma2 on likelihood results"            )
    ("uqDRAM_mh_lrSigma2Accuracies",      po::value<std::string >()->default_value(UQ_DRAM_MCG_MH_LR_SIGMA2_ACCURACIES_ODV   ), "'mh' accuracies of sigma2 on likelihood results"        )
    ("uqDRAM_mh_lrNumbersOfObs",          po::value<std::string >()->default_value(UQ_DRAM_MCG_MH_LR_NUMBERS_OF_OBS_ODV      ), "'mh' numbers of observations for likelihood results"    )
    ("uqDRAM_dr_maxNumberOfExtraStages",  po::value<unsigned int>()->default_value(UQ_DRAM_MCG_DR_MAX_NUM_EXTRA_STAGES_ODV   ), "'dr' maximum number of extra stages"                    )
    ("uqDRAM_dr_scalesForExtraStages",    po::value<std::string >()->default_value(UQ_DRAM_MCG_DR_SCALES_FOR_EXTRA_STAGES_ODV), "'dr' scales for proposal cov matrices from 2nd stage on")
    ("uqDRAM_am_initialNonAdaptInterval", po::value<unsigned int>()->default_value(UQ_DRAM_MCG_AM_INIT_NON_ADAPT_INT_ODV     ), "'am' initial non adaptation interval"                   )
    ("uqDRAM_am_adaptInterval",           po::value<unsigned int>()->default_value(UQ_DRAM_MCG_AM_ADAPT_INTERVAL_ODV         ), "'am' adaptation interval"                               )
    ("uqDRAM_am_sd",                      po::value<double      >()->default_value(UQ_DRAM_MCG_AM_SD_ODV                     ), "'am' sd"                                                )
    ("uqDRAM_am_epsilon",                 po::value<double      >()->default_value(UQ_DRAM_MCG_AM_EPSILON_ODV                ), "'am' epsilon"                                           )
    ("uqDRAM_mh_namesOfOutputFiles",      po::value<std::string >()->default_value(UQ_DRAM_MCG_MH_OUTPUT_FILE_NAME_ODV       ), "'mh' name(s) of output file(s)"                         )
    ("uqDRAM_mh_chainDisplayPeriod",      po::value<unsigned int>()->default_value(UQ_DRAM_MCG_MH_CHAIN_DISPLAY_PERIOD_ODV   ), "'mh' period of message display during chain generation" )
  ;

  return;
}

template <class V, class M>
void
uqDRAM_MarkovChainGeneratorClass<V,M>::getMyOptionValues(
  po::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count("uqDRAM_help")) {
    std::cout << optionsDesc
              << std::endl;
  }

  if (m_env.allOptionsMap().count("uqDRAM_mh_sizesOfChains")) {
    std::string inputString = m_env.allOptionsMap()["uqDRAM_mh_sizesOfChains"].as<std::string>();
    std::vector<double> inputDoubles(0,0.);
    uqMiscReadDoublesFromString(inputString,inputDoubles);
    m_sizesOfChains.clear();

    m_sizesOfChains.resize(inputDoubles.size(),0);
    for (unsigned int i = 0; i < inputDoubles.size(); ++i) {
      m_sizesOfChains[i] = (unsigned int) inputDoubles[i];
    }
  }

  if (m_env.allOptionsMap().count("uqDRAM_mh_lrUpdateSigma2")) {
    m_lrUpdateSigma2 = m_env.allOptionsMap()["uqDRAM_mh_lrUpdateSigma2"].as<bool>();
  }

  if (m_env.allOptionsMap().count("uqDRAM_mh_lrSigma2Priors")) {
    std::string inputString = m_env.allOptionsMap()["uqDRAM_mh_lrSigma2Priors"].as<std::string>();
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

  if (m_env.allOptionsMap().count("uqDRAM_mh_lrSigma2Accuracies")) {
    std::string inputString = m_env.allOptionsMap()["uqDRAM_mh_lrSigma2Accuracies"].as<std::string>();
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

  if (m_env.allOptionsMap().count("uqDRAM_mh_lrNumbersOfObs")) {
    std::string inputString = m_env.allOptionsMap()["uqDRAM_mh_lrNumbersOfObs"].as<std::string>();
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

  if (m_env.allOptionsMap().count("uqDRAM_dr_maxNumberOfExtraStages")) {
    m_maxNumberOfExtraStages = m_env.allOptionsMap()["uqDRAM_dr_maxNumberOfExtraStages"].as<unsigned int>();
  }

  std::vector<double> tmpScales(0,0.);
  if (m_env.allOptionsMap().count("uqDRAM_dr_scalesForExtraStages")) {
    std::string inputString = m_env.allOptionsMap()["uqDRAM_dr_scalesForExtraStages"].as<std::string>();
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

  if (m_env.allOptionsMap().count("uqDRAM_am_initialNonAdaptInterval")) {
    m_initialNonAdaptInterval = m_env.allOptionsMap()["uqDRAM_am_initialNonAdaptInterval"].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count("uqDRAM_am_adaptInterval")) {
    m_adaptInterval = m_env.allOptionsMap()["uqDRAM_am_adaptInterval"].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count("uqDRAM_am_sd")) {
    m_sd = m_env.allOptionsMap()["uqDRAM_am_sd"].as<double>();
  }

  if (m_env.allOptionsMap().count("uqDRAM_am_epsilon")) {
    m_epsilon = m_env.allOptionsMap()["uqDRAM_am_epsilon"].as<double>();
  }

  if (m_env.allOptionsMap().count("uqDRAM_mh_namesOfOutputFiles")) {
    std::string inputString = m_env.allOptionsMap()["uqDRAM_mh_namesOfOutputFiles"].as<std::string>();
    m_namesOfOutputFiles.clear();
    uqMiscReadWordsFromString(inputString,m_namesOfOutputFiles);

    UQ_FATAL_TEST_MACRO(m_namesOfOutputFiles.size() != m_sizesOfChains.size(),
                        m_env.rank(),
                        "uqDRAM_MarkovChainGeneratorClass<V,M>::readMyOptionsValues()",
                        "size of array for 'outputFileNames' is not equal to size of array for 'chainSizes'");
  }

  if (m_env.allOptionsMap().count("uqDRAM_mh_chainDisplayPeriod")) {
    m_chainDisplayPeriod = m_env.allOptionsMap()["uqDRAM_mh_chainDisplayPeriod"].as<unsigned int>();
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

    bool accept = acceptAlpha(this->alpha(currentPosition,currentCandidate,&m_alphaQuotients[positionId]));

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
        accept = acceptAlpha(this->alpha(drPositions));
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
          M* tmpI = m_paramSpace.newDiagMatrix(m_epsilon);
          tmpChol = *m_lastAdaptedCovMatrix + *tmpI;
          delete tmpI;
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
          *(m_lowerCholProposalCovMatrices[0]) *= sqrt(m_sd);
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

  V* chainMean = m_paramSpace.newVector();
  V* chainStd  = m_paramSpace.newVector();
  uqChainStats(m_chain,
               *chainMean,
               *chainStd);

  if (m_env.rank() == 0) {
    std::cout << "Finished generating the chain of id " << chainId
              << ". Chain statistics are:";
    std::cout << "\n  Rejection percentage = " << 100. * (double) m_numRejections/(double) m_chain.size()
              << " %";
    std::cout << "\n   Outbound percentage = " << 100. * (double) m_numOutOfBounds/(double) m_chain.size()
              << " %";
    char line[512];
    sprintf(line,"\n%s%4s%s%9s%s",
	    "Parameter",
            " ",
            "Mean",
            " ",
            "Std");
    std::cout << line;

    for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
      sprintf(line,"\n%8.8s%2s%11.4e%2s%11.4e",
              m_paramSpace.parameter(i).name().c_str(),
              " ",
	      (*chainMean)[i],
              " ",
              (*chainStd)[i]);
      std::cout << line;
    }
    std::cout << std::endl;
  }

  //****************************************************
  // Release memory before leaving routine
  //****************************************************
  if (chainStd      ) delete chainStd;
  if (chainMean     ) delete chainMean;
  if (gaussianVector) delete gaussianVector;
  if (tmpParamValues) delete tmpParamValues;

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
  ofs << "chainCpp = [";
  for (unsigned int i = 0; i < m_chain.size(); ++i) {
    ofs << *(m_chain[i])
        << std::endl;
  }
  ofs << "];\n";

  // Write m_lrChain
  ofs << "sschainCpp = [";
  for (unsigned int i = 0; i < m_lrChain.size(); ++i) {
    ofs << *(m_lrChain[i])
        << std::endl;
  }
  ofs << "];\n";

  // Write m_lrSigma2Chain
  ofs << "s2chainCpp = [";
  for (unsigned int i = 0; i < m_lrSigma2Chain.size(); ++i) {
    ofs << *(m_lrSigma2Chain[i])
        << std::endl;
  }
  ofs << "];\n";

  // Write m_alphaQuotients
  ofs << "alphaQuotientsCpp = [";
  for (unsigned int i = 0; i < m_alphaQuotients.size(); ++i) {
    ofs << m_alphaQuotients[i]
        << std::endl;
  }
  ofs << "];\n";

  // Write names of parameters
  ofs << "resultsCpp.names = {";
  m_paramSpace.printParameterNames(ofs,false);
  ofs << "};\n";

  if (mahalanobisMatrix != NULL) {
    // Write mahalanobis distances
    V* diffVec = m_paramSpace.newVector();
    ofs << "dCpp = [";
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
  ofs << "priorMeanValuesCpp = ["
      << m_paramSpace.priorMuValues()
      << "];\n";

  // Write prior sigma values
  ofs << "priorSigValuesCpp = ["
      << m_paramSpace.priorSigmaValues()
      << "];\n";

  ofs << "resultsCpp.prior = [priorMeanValuesCpp',priorSigValuesCpp'];\n";

  // Write param lower bounds
  ofs << "lowCpp = ["
      << m_paramSpace.minValues()
      << "];\n";

  // Write param upper bounds
  ofs << "uppCpp = ["
      << m_paramSpace.maxValues()
      << "];\n";

  ofs << "resultsCpp.limits = [lowCpp',uppCpp'];\n";

  // Write out data for mcmcpred.m
  ofs << "resultsCpp.parind = ["; // FIXME
  for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
    ofs << i+1
        << std::endl;
  }
  ofs << "];\n";

  ofs << "resultsCpp.local = [\n"; // FIXME
  for (unsigned int i = 0; i < m_paramSpace.dim(); ++i) {
    ofs << " 0";
    //<< std::endl;
  }
  ofs << "];\n";

  bool savedVectorPrintState = m_chain[m_chain.size()-1]->getPrintHorizontally();
  m_chain[m_chain.size()-1]->setPrintHorizontally(false);
  ofs << "resultsCpp.theta = ["
      << *(m_chain[m_chain.size()-1])
      << "];\n";
  m_chain[m_chain.size()-1]->setPrintHorizontally(savedVectorPrintState);
  
  ofs << "resultsCpp.nbatch = 1.;\n"; // FIXME

  if (mahalanobisMatrix != NULL) {
    // Write covMatrix
    ofs << "mahalanobisMatrixCpp = ["
        << *mahalanobisMatrix
        << "];\n";
  }

  // Write number of rejections
  ofs << "resultsCpp.rejected = "  << (double) m_numRejections/(double) m_chain.size()
      << ";\n"
      << std::endl;

  // Write number of outbounds
  ofs << "resultsCpp.ulrejected = " << (double) m_numOutOfBounds/(double) m_chain.size()
      << ";\n"
      << std::endl;

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
  os << "m_sizesOfChains = ";
  for (unsigned int i = 0; i < m_sizesOfChains.size(); ++i) {
    os << m_sizesOfChains[i] << " ";
  }
  os << "\nm_lrUpdateSigma2 = " << m_lrUpdateSigma2
     << "\nm_lrSigma2Priors = ";
  for (unsigned int i = 0; i < m_lrSigma2Priors->size(); ++i) {
    os << (*m_lrSigma2Priors)[i] << " ";
  }
  os << "\nm_lrSigma2Accuracies = ";
  for (unsigned int i = 0; i < m_lrSigma2Accuracies->size(); ++i) {
    os << (*m_lrSigma2Accuracies)[i] << " ";
  }
  os << "\nm_lrNumbersOfObservations = ";
  for (unsigned int i = 0; i < m_lrNumbersOfObservations->size(); ++i) {
    os << (*m_lrNumbersOfObservations)[i] << " ";
  }
  os << "\nm_maxNumberOfExtraStages = " << m_maxNumberOfExtraStages
     << "\nm_scalesForCovMProposals = ";
  for (unsigned int i = 0; i < m_scalesForCovMProposals.size(); ++i) {
    os << m_scalesForCovMProposals[i] << " ";
  }
  os << "\nm_initialNonAdaptInterval = " << m_initialNonAdaptInterval
     << "\nm_adaptInterval = "           << m_adaptInterval
     << "\nm_sd = "                      << m_sd
     << "\nm_epsilon = "                 << m_epsilon;
  os << "\nm_namesOfOutputFiles = ";
  for (unsigned int i = 0; i < m_namesOfOutputFiles.size(); ++i) {
    os << m_namesOfOutputFiles[i] << " ";
  }
  os << "\nm_chainDisplayPeriod = " << m_chainDisplayPeriod;
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
