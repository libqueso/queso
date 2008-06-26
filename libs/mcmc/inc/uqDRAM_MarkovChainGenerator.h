#ifndef __UQ_DRAM_MCG_H__
#define __UQ_DRAM_MCG_H__

#undef UQ_DRAM_MCG_REQUIRES_INVERTED_COV_MATRICES

// _ODV = option default value
//   LR = likelihood results
#define UQ_DRAM_MCG_CHAIN_SIZE_ODV           "100"
#define UQ_DRAM_MCG_LR_UPDATE_SIGMA2_ODV     false
#define UQ_DRAM_MCG_LR_SIGMA2_PRIORS_ODV     "1."
#define UQ_DRAM_MCG_LR_SIGMA2_ACCURACIES_ODV "1."
#define UQ_DRAM_MCG_LR_NUMBERS_OF_OBS_ODV    "1."
#define UQ_DRAM_MCG_NUM_EXTRA_STAGES_ODV     0
#define UQ_DRAM_MCG_COV_MATRIX_SCALES_ODV    "1."
#define UQ_DRAM_MCG_INIT_NON_ADAPT_INT_ODV   0
#define UQ_DRAM_MCG_ADAPT_INTERVAL_ODV       0
#define UQ_DRAM_MCG_SD_ODV                   1.
#define UQ_DRAM_MCG_EPSILON_ODV              1.e-5
#define UQ_DRAM_MCG_OUTPUT_FILE_NAME_ODV     "."

#include <uqParamSpace.h>
#include <uqStateSpace.h>
#include <uqOutputSpace.h>
#include <uqCandidate.h>
#include <uqMiscellaneous.h>
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
  const std::vector<const V*>& chain() const;

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

  double logProposal             (const uqCandidateClass<V>& x,
                                  const uqCandidateClass<V>& y,
                                  unsigned int               idOfMatrix);
  double logProposal             (const std::vector<uqCandidateClass<V>*>& inputPositions);
  double alpha                   (const uqCandidateClass<V>& x,
                                  const uqCandidateClass<V>& y,
                                  double* alphaQuotientPtr = NULL);
  double alpha                   (const std::vector<uqCandidateClass<V>*>& inputPositions);
  bool   acceptAlpha             (double alpha);
  void   updateCovMatrices       ();
  void   updateCovMatrix         (const std::vector<V*>& subChain,
                                  double&                lastSubChainSize,
                                  V*&                    lastSubMean,
                                  M*&                    lastAdaptedCovMatrix);
  void   gammar                  (double       a,
                                  double       b,
                                  M&           mat);

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
  unsigned int              m_numExtraStages;
  std::vector<double>       m_scalesForCovMProposals;
  unsigned int              m_initialNonAdaptInterval;
  unsigned int              m_adaptInterval;
  double                    m_sd;
  double                    m_epsilon;
  std::vector<std::string>  m_namesOfOutputFiles;
  gsl_rng*                  m_rng;

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
  double                    m_lastSubChainSize;
  V*                        m_lastSubMean;
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
  m_sizesOfChains                (1,(unsigned int) strtod(UQ_DRAM_MCG_CHAIN_SIZE_ODV,NULL)),
  m_lrUpdateSigma2               (UQ_DRAM_MCG_LR_UPDATE_SIGMA2_ODV),
  m_lrSigma2Priors               (m_outputSpace.newVector()),
  m_lrSigma2Accuracies           (m_outputSpace.newVector()),
  m_lrNumbersOfObservations      (m_outputSpace.newVector()),
  m_numExtraStages               (UQ_DRAM_MCG_NUM_EXTRA_STAGES_ODV),
  m_scalesForCovMProposals       (0,0.),
  m_initialNonAdaptInterval      (UQ_DRAM_MCG_INIT_NON_ADAPT_INT_ODV),
  m_adaptInterval                (UQ_DRAM_MCG_ADAPT_INTERVAL_ODV),
  m_sd                           (UQ_DRAM_MCG_SD_ODV),
  m_epsilon                      (UQ_DRAM_MCG_EPSILON_ODV),
  m_namesOfOutputFiles           (1,UQ_DRAM_MCG_OUTPUT_FILE_NAME_ODV),
  m_rng                          (NULL),
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
  m_lastSubChainSize             (0),
  m_lastSubMean                  (NULL),
  m_lastAdaptedCovMatrix         (NULL)
{
  if (m_useDefaultM2lPriorFunction) {
    m_m2lPriorFunction_ObjPtr = new uq_M2lPriorFunction_Class<V,M>();
    m_defaultM2lPriorFunction_Data.paramPriorMus    = new V(m_paramSpace.priorMuValues   ());
    m_defaultM2lPriorFunction_Data.paramPriorSigmas = new V(m_paramSpace.priorSigmaValues());
  }
  m_lrSigma2Priors->cwSet         (strtod(UQ_DRAM_MCG_LR_SIGMA2_PRIORS_ODV    ,NULL));
  m_lrSigma2Accuracies->cwSet     (strtod(UQ_DRAM_MCG_LR_SIGMA2_ACCURACIES_ODV,NULL));
  m_lrNumbersOfObservations->cwSet(strtod(UQ_DRAM_MCG_LR_NUMBERS_OF_OBS_ODV   ,NULL));

  //std::cout << "Entering uqDRAM_MarkovChainGeneratorClass<V,M>::constructor()"
  //          << std::endl;

  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.rank() == 0) std::cout << "After getting option values, state of uqDRAM_MarkovChainGeneratorClass object is:"
                                   << "\n" << *this
                                   << std::endl;

  m_rng = gsl_rng_alloc(gsl_rng_ranlxd2);
  UQ_FATAL_TEST_MACRO((m_rng == NULL),
                      m_env.rank(),
                      "uqDRAM_MarkovChainGeneratorClass<V,M>::constructor()",
                      "null m_rng");

  //std::cout << "Leaving uqDRAM_MarkovChainGeneratorClass<V,M>::constructor()"
  //          << std::endl;
}

template <class V, class M>
uqDRAM_MarkovChainGeneratorClass<V,M>::~uqDRAM_MarkovChainGeneratorClass()
{
  //std::cout << "Entering uqDRAM_MarkovChainGeneratorClass<V,M>::destructor()"
  //          << std::endl;

  resetChainAndRelatedInfo();

  if (m_rng                    ) gsl_rng_free(m_rng);
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
  if (m_lastSubMean)          delete m_lastSubMean;
  m_lastSubChainSize = 0;
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
    ("uqDRAM_MCG_help",                                                                                               "produce help message for DRAM Markov chain generator"   )
    ("uqDRAM_MCG_sizesOfChains",      po::value<std::string >()->default_value(UQ_DRAM_MCG_CHAIN_SIZE_ODV          ), "'mh' size(s) of chain(s)"                               )
    ("uqDRAM_MCG_lrUpdateSigma2",     po::value<bool        >()->default_value(UQ_DRAM_MCG_LR_UPDATE_SIGMA2_ODV    ), "'mh' update the sigma2 of likelihood results"           )
    ("uqDRAM_MCG_lrSigma2Priors",     po::value<std::string >()->default_value(UQ_DRAM_MCG_LR_SIGMA2_PRIORS_ODV    ), "'mh' priors of sigma2 on likelihood results"            )
    ("uqDRAM_MCG_lrSigma2Accuracies", po::value<std::string >()->default_value(UQ_DRAM_MCG_LR_SIGMA2_ACCURACIES_ODV), "'mh' accuracies of sigma2 on likelihood results"        )
    ("uqDRAM_MCG_lrNumbersOfObs",     po::value<std::string >()->default_value(UQ_DRAM_MCG_LR_NUMBERS_OF_OBS_ODV   ), "'mh' numbers of observations for likelihood results"    )
    ("uqDRAM_MCG_numExtraStages",     po::value<unsigned int>()->default_value(UQ_DRAM_MCG_NUM_EXTRA_STAGES_ODV    ), "'dr' maximum number of extra stages"                    )
    ("uqDRAM_MCG_stageScales",        po::value<std::string >()->default_value(UQ_DRAM_MCG_COV_MATRIX_SCALES_ODV   ), "'dr' scales for proposal cov matrices from 2nd stage on")
    ("uqDRAM_MCG_initialNonAdaptInt", po::value<unsigned int>()->default_value(UQ_DRAM_MCG_INIT_NON_ADAPT_INT_ODV  ), "'am' initial non adaptation interval"                   )
    ("uqDRAM_MCG_adaptInterval",      po::value<unsigned int>()->default_value(UQ_DRAM_MCG_ADAPT_INTERVAL_ODV      ), "'am' adaptation interval"                               )
    ("uqDRAM_MCG_sd",                 po::value<double      >()->default_value(UQ_DRAM_MCG_SD_ODV                  ), "'am' sd"                                                )
    ("uqDRAM_MCG_epsilon",            po::value<double      >()->default_value(UQ_DRAM_MCG_EPSILON_ODV             ), "'am' epsilon"                                           )
    ("uqDRAM_MCG_namesOfOutputFiles", po::value<std::string >()->default_value(UQ_DRAM_MCG_OUTPUT_FILE_NAME_ODV    ), "'mh' name(s) of output file(s)"                         )
  ;

  return;
}

template <class V, class M>
void
uqDRAM_MarkovChainGeneratorClass<V,M>::getMyOptionValues(
  po::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count("uqDRAM_MCG_help")) {
    std::cout << optionsDesc
              << std::endl;
  }

  if (m_env.allOptionsMap().count("uqDRAM_MCG_sizesOfChains")) {
    std::string inputString = m_env.allOptionsMap()["uqDRAM_MCG_sizesOfChains"].as<std::string>();
    std::vector<double> inputDoubles(0,0.);
    uqMiscReadDoublesFromString(inputString,inputDoubles);
    m_sizesOfChains.clear();

    m_sizesOfChains.resize(inputDoubles.size(),0);
    for (unsigned int i = 0; i < inputDoubles.size(); ++i) {
      m_sizesOfChains[i] = (unsigned int) inputDoubles[i];
    }
  }

  if (m_env.allOptionsMap().count("uqDRAM_MCG_lrUpdateSigma2")) {
    m_lrUpdateSigma2 = m_env.allOptionsMap()["uqDRAM_MCG_lrUpdateSigma2"].as<bool>();
  }

  if (m_env.allOptionsMap().count("uqDRAM_MCG_lrSigma2Priors")) {
    std::string inputString = m_env.allOptionsMap()["uqDRAM_MCG_lrSigma2Priors"].as<std::string>();
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

  if (m_env.allOptionsMap().count("uqDRAM_MCG_lrSigma2Accuracies")) {
    std::string inputString = m_env.allOptionsMap()["uqDRAM_MCG_lrSigma2Accuracies"].as<std::string>();
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

  if (m_env.allOptionsMap().count("uqDRAM_MCG_lrNumbersOfObs")) {
    std::string inputString = m_env.allOptionsMap()["uqDRAM_MCG_lrNumbersOfObs"].as<std::string>();
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

  if (m_env.allOptionsMap().count("uqDRAM_MCG_numExtraStages")) {
    m_numExtraStages = m_env.allOptionsMap()["uqDRAM_MCG_numExtraStages"].as<unsigned int>();
  }

  std::vector<double> tmpScales(0,0.);
  if (m_env.allOptionsMap().count("uqDRAM_MCG_scales")) {
    std::string inputString = m_env.allOptionsMap()["uqDRAM_MCG_scales"].as<std::string>();
    uqMiscReadDoublesFromString(inputString,tmpScales);
    //std::cout << "In uqDRAM_MarkovChainGeneratorClass<V,M>::getMyOptionValues(): scales =";
    //for (unsigned int i = 0; i < tmpScales.size(); ++i) {
    //  std::cout << " " << tmpScales[i];
    //}
    //std::cout << std::endl;
  }

  if (m_numExtraStages > 0) {
    double scale = 1.0;
    unsigned int tmpSize = tmpScales.size();
    m_scalesForCovMProposals.resize(m_numExtraStages+1,1.);
    for (unsigned int i = 1; i < (m_numExtraStages+1); ++i) {
      if (i <= tmpSize) scale = tmpScales[i-1];
      m_scalesForCovMProposals[i] = scale;
    }
    //updateCovMatrices();
  }

  if (m_env.allOptionsMap().count("uqDRAM_MCG_initialNonAdaptInt")) {
    m_initialNonAdaptInterval = m_env.allOptionsMap()["uqDRAM_MCG_initialNonAdaptInt"].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count("uqDRAM_MCG_adaptInterval")) {
    m_adaptInterval = m_env.allOptionsMap()["uqDRAM_MCG_adaptInterval"].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count("uqDRAM_MCG_sd")) {
    m_sd = m_env.allOptionsMap()["uqDRAM_MCG_sd"].as<double>();
  }

  if (m_env.allOptionsMap().count("uqDRAM_MCG_epsilon")) {
    m_epsilon = m_env.allOptionsMap()["uqDRAM_MCG_epsilon"].as<double>();
  }

  if (m_env.allOptionsMap().count("uqDRAM_MCG_namesOfOutputFiles")) {
    std::string inputString = m_env.allOptionsMap()["uqDRAM_MCG_namesOfOutputFiles"].as<std::string>();
    m_namesOfOutputFiles.clear();
    uqMiscReadWordsFromString(inputString,m_namesOfOutputFiles);

    UQ_FATAL_TEST_MACRO(m_namesOfOutputFiles.size() != m_sizesOfChains.size(),
                        m_env.rank(),
                        "uqDRAM_MarkovChainGeneratorClass<V,M>::readMyOptionsValues()",
                        "size of array for 'outputFileNames' is not equal to size of array for 'chainSizes'");
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

  if (m_numExtraStages > 0) {
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
  for (unsigned int i = 1; i < (m_numExtraStages+1); ++i) {
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
  uqCandidateClass<V> currentPosition(m_env,
                                      valuesOf1stPosition,
                                      outOfBounds,
                                      m2lPrior,
                                      *m2lLikelihoodResults,
                                      lrSigma2,
                                      logPosterior);

  V* gaussianVector = m_paramSpace.newVector();
  V* tmpParamValues = m_paramSpace.newVector();
  uqCandidateClass<V> currentCandidate(m_env);

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
    //                                 << ", m_numExtraStages = "               << m_numExtraStages
    //                                 << std::endl;
    unsigned int stageId = 0;

    //****************************************************
    // Loop: generate new parameters
    //****************************************************
    gaussianVector->cwSetGaussian(m_rng,0.,1.);

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
    if ((accept == false) && (outOfBounds == false) && (m_numExtraStages > 0)) {
      std::vector<uqCandidateClass<V>*> positions(stageId+2,NULL);
      positions[0] = new uqCandidateClass<V>(currentPosition);
      positions[1] = new uqCandidateClass<V>(currentCandidate);

      while ((accept == false) && (stageId < m_numExtraStages)) {
        stageId++;

        gaussianVector->cwSetGaussian(m_rng,0.,1.);

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

        positions.push_back(new uqCandidateClass<V>(currentCandidate));
        accept = acceptAlpha(this->alpha(positions));
      } // while
    }

    //****************************************************
    // Loop: update chain
    //****************************************************
    if (accept) {
      m_chain[positionId]   = m_paramSpace.uqFinDimLinearSpaceClass<V,M>::newVector(currentCandidate.paramValues());
      m_lrChain[positionId] = m_outputSpace.newVector(currentCandidate.m2lLikelihood());
      currentPosition  = currentCandidate;
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
        lrSigma2[i] = 1./uqMiscGammar(term1,term2,m_rng);
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
    // Loop: adaptation of covariance matrix
    //****************************************************
    if (m_adaptInterval > 0) {
      // Now might be the moment to adapt
      std::vector<V*> subChain(0,NULL);

      // Check if now is indeed the moment to adapt
      if (positionId < m_initialNonAdaptInterval) {
        // Do nothing
      }
      else if (positionId == m_initialNonAdaptInterval) {
        subChain.resize(m_initialNonAdaptInterval,NULL);
        for (unsigned int i = 0; i < subChain.size(); ++i) {
          subChain[i] = new V(*(m_chain[i]));
        }
      }
      else {
        unsigned int deltaInterval = positionId - m_initialNonAdaptInterval;
        if ((deltaInterval % m_adaptInterval) == 0) {
          subChain.resize(m_adaptInterval,NULL);
          for (unsigned int i = 0; i < subChain.size(); ++i) {
            subChain[i] = new V(*(m_chain[m_chain.size()-subChain.size()+i]));
          }
        }
      }

      // If now is indeed the moment to adapt, then do it!
      if (subChain.size() > 0) {
        updateCovMatrix(subChain,
                        m_lastSubChainSize,
                        m_lastSubMean,
                        m_lastAdaptedCovMatrix);

        bool tmpCholIsPositiveDefinite = false;
        M tmpChol(*m_lastAdaptedCovMatrix);
        iRC = tmpChol.chol();
        if (iRC) {
          // Matrix is not positive definite
          tmpChol  = *m_lastAdaptedCovMatrix;
          M* tmpI = m_paramSpace.newDiagMatrix(m_epsilon);
          tmpChol += *tmpI;
          delete tmpI;
          iRC = tmpChol.chol();
          if (iRC) {
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

          if (m_numExtraStages > 0) updateCovMatrices();
        }

        for (unsigned int i = 0; i < subChain.size(); ++i) {
          if (subChain[i]) delete subChain[i];
        }
      }
    }
  } // end chain loop

  //****************************************************
  // Release memory before leaving routine
  //****************************************************
  if (gaussianVector) delete gaussianVector;
  if (tmpParamValues) delete tmpParamValues;

  return iRC;
}

template <class V, class M>
double
uqDRAM_MarkovChainGeneratorClass<V,M>::logProposal(
  const uqCandidateClass<V>& x,
  const uqCandidateClass<V>& y,
  unsigned int               idOfMatrix)
{
  V diffVec(y.paramValues() - x.paramValues());
#ifdef UQ_DRAM_MCG_REQUIRES_INVERTED_COV_MATRICES
  double value = -0.5 * scalarProduct(diffVec, *(m_proposalPrecMatrices[idOfMatrix]) * diffVec);
#else
  double value = -0.5 * scalarProduct(diffVec, m_proposalCovMatrices[idOfMatrix]->invertMultiply(diffVec));
#endif
  return value;
}

template <class V, class M>
double
uqDRAM_MarkovChainGeneratorClass<V,M>::logProposal(const std::vector<uqCandidateClass<V>*>& inputPositions)
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
  const uqCandidateClass<V>& x,
  const uqCandidateClass<V>& y,
  double* alphaQuotientPtr)
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
uqDRAM_MarkovChainGeneratorClass<V,M>::alpha(const std::vector<uqCandidateClass<V>*>& inputPositions)
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
  std::vector<uqCandidateClass<V>*>         positions(inputSize,NULL);
  std::vector<uqCandidateClass<V>*> backwardPositions(inputSize,NULL);
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

  for (unsigned int i = 0; i < inputSize; ++i) {
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

  if      (alpha <= 0.                    ) result = false;
  else if (alpha >= 1.                    ) result = true;
  else if (alpha >= gsl_rng_uniform(m_rng)) result = true;
  else                                      result = false;

  return result;
}

template <class V, class M>
void
uqDRAM_MarkovChainGeneratorClass<V,M>::updateCovMatrix(
  const std::vector<V*>& subChain,
  double&                lastSubChainSize,
  V*&                    lastSubMean,
  M*&                    lastAdaptedCovMatrix)
{
  if ((lastSubChainSize     == 0   ) &&
      (lastSubMean          == NULL) &&
      (lastAdaptedCovMatrix == NULL)) {
    if (subChain.size() < 2) {
      UQ_FATAL_RC_MACRO(UQ_INVALID_INPUT_PARAMETER_RC,
                        m_env.rank(),
                        "uqDRAM_MarkovChainGeneratorClass<V,M>::updateCovMatrix()",
                        "'subChain.size()' should be >= 2");
    }
  
    double doubleSubChainSize = (double) subChain.size();
    lastSubMean = m_paramSpace.newVector();
    for (unsigned int i = 0; i < subChain.size(); ++i) {
      *lastSubMean += *(subChain[i]);
    }
    *lastSubMean /= doubleSubChainSize;

    lastAdaptedCovMatrix = m_paramSpace.newMatrix();

    M* tmpMatrix = new M(matrixProduct(*lastSubMean,*lastSubMean));
    *lastAdaptedCovMatrix += *tmpMatrix;
    delete tmpMatrix;
    *lastAdaptedCovMatrix *= doubleSubChainSize;

    for (unsigned int i = 0; i < subChain.size(); ++i) {
      tmpMatrix = new M(matrixProduct(*(subChain[i]),*(subChain[i])));
      *lastAdaptedCovMatrix += *tmpMatrix;
      delete tmpMatrix;
    }
    *lastAdaptedCovMatrix /= (doubleSubChainSize - 1.); // That is why size must be >= 2
  }
  else if ((lastSubChainSize     >  0   ) &&
           (lastSubMean          != NULL) &&
           (lastAdaptedCovMatrix != NULL)) {
    if (subChain.size() < 1) {
      UQ_FATAL_RC_MACRO(UQ_INVALID_INPUT_PARAMETER_RC,
                        m_env.rank(),
                        "uqDRAM_MarkovChainGeneratorClass<V,M>::updateCovMatrix()",
                        "'subChain.size()' should be >= 1");
    }
  }
  else {
    UQ_FATAL_RC_MACRO(UQ_INVALID_INPUT_PARAMETER_RC,
                      m_env.rank(),
                      "uqDRAM_MarkovChainGeneratorClass<V,M>::updateCovMatrix()",
                      "'lastSubChainSize', 'lastSubMean' and 'lastAdaptedCovMatrix' are incompatible with each other");
  }

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
      mat(i,j) = uqMiscGammar(a,b,m_rng);
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

  if (m_namesOfOutputFiles[chainId] == ".") return iRC;

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

  return iRC;
}

template <class V, class M>
const std::vector<const V*>&
uqDRAM_MarkovChainGeneratorClass<V,M>::chain() const
{
  return m_chain;
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
  os << "\nm_numExtraStages = " << m_numExtraStages
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
