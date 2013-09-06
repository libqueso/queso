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

#ifndef __UQ_MH_SG1_H__
#define __UQ_MH_SG1_H__

#include <uqMetropolisHastingsSGOptions.h>
#include <uqTKGroup.h>
#include <uqVectorRV.h>
#include <uqVectorSpace.h>
#include <uqMarkovChainPositionData.h>
#include <uqScalarFunctionSynchronizer.h>
#include <uqSequenceOfVectors.h>
#include <uqArrayOfSequences.h>
#include <sys/time.h>
#include <fstream>
#include <boost/math/special_functions.hpp> // for Boost isnan. Note parentheses are important in function call.

namespace QUESO {

//--------------------------------------------------
// uqMHRawChainInfoStruct --------------------------
//--------------------------------------------------
 /*! \file uqMetropolisHastingsSG1.h
 * \brief A templated class that represents a Metropolis-Hastings generator of samples and a struct which stores its info.
 * 
 * \struct  uqMHRawChainInfoStruct
 * \brief A struct that represents a Metropolis-Hastings sample.
 * 
 * Some of the information about the  Metropolis-Hastings sample generator includes the allowed number
 * of delayed rejections, number of rejections, number of positions in or out of target support, and 
 * so on. This struct is responsible for the storage of such info. */

struct uqMHRawChainInfoStruct
{
 //! @name Constructor/Destructor methods
 //@{
 //! Constructor. 
  uqMHRawChainInfoStruct();
  
  //! Copy constructor. 
  uqMHRawChainInfoStruct(const uqMHRawChainInfoStruct& rhs);
  
  //! Destructor
  ~uqMHRawChainInfoStruct();
  //@}

  //! @name Set methods
  //@{
  //! Assignment operator.
  uqMHRawChainInfoStruct& operator= (const uqMHRawChainInfoStruct& rhs);
  
  //! Addition assignment operator.
  uqMHRawChainInfoStruct& operator+=(const uqMHRawChainInfoStruct& rhs);
  //@}
  
   //! @name Misc methods
  //@{
  //! Copies Metropolis-Hastings chain info from \c src to \c this.
  void copy  (const uqMHRawChainInfoStruct& src);
  
  //! Resets Metropolis-Hastings chain info.
  void reset ();
  
  //! Calculates the MPI sum of \c this.
  void mpiSum(const uqMpiCommClass& comm, uqMHRawChainInfoStruct& sumInfo) const;
  //@}
  
  double       runTime;
  double       candidateRunTime;
  double       targetRunTime;
  double       mhAlphaRunTime;
  double       drAlphaRunTime;
  double       drRunTime;
  double       amRunTime;

  unsigned int numTargetCalls;
  unsigned int numDRs;
  unsigned int numOutOfTargetSupport;
  unsigned int numOutOfTargetSupportInDR;
  unsigned int numRejections;

};

//--------------------------------------------------
// uqMetropolisHastingsSGClass----------------------
//--------------------------------------------------

/*!\class uqMetropolisHastingsSGClass
 * \brief A templated class that represents a Metropolis-Hastings generator of samples.
 *
 * This class implements a Metropolis-Hastings generator of samples. 'SG' stands for 'Sequence Generator'.
 * Options reading is handled by class 'uqMetropolisHastingsOptionsClass'. If options request data to be 
 * written in the output file (MATLAB .m format only, for now), the user can check which MATLAB variables 
 * are defined and set by running 'grep zeros <OUTPUT FILE NAME>' after the solution procedures ends. 
 * The names of the variables are self explanatory. */

template <class P_V,class P_M>
class uqMetropolisHastingsSGClass
{
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Constructor.
  /*! This method reads the options from the options input file. It calls commonConstructor().
   * Requirements: 1) the image set of the vector random variable 'sourceRv' should belong to a 
   * vector space of dimension equal to the size of the vector 'initialPosition' and 2) if 
   * 'inputProposalCovMatrix' is not NULL, is should be square and its size should be equal to the
   * size of 'initialPosition'. If the requirements are satisfied, the constructor then reads input 
   * options that begin with the string '\<prefix\>_mh_'. For instance, if 'prefix' is 
   * 'pROblem_775_ip_', then the constructor will read all options that begin with 'pROblem_775_ip_mh_'.
    Options reading is handled by class 'uqMetropolisHastingsOptionsClass'.*/
  uqMetropolisHastingsSGClass(const char*                         prefix,
                              const uqMhOptionsValuesClass*       alternativeOptionsValues, // dakota
                              const uqBaseVectorRVClass<P_V,P_M>& sourceRv,
                              const P_V&                          initialPosition,
                              const P_M*                          inputProposalCovMatrix);
  
  //! Constructor.
  uqMetropolisHastingsSGClass(const uqMLSamplingLevelOptionsClass& mlOptions,
                              const uqBaseVectorRVClass<P_V,P_M>&  sourceRv,
                              const P_V&                           initialPosition,
                              const P_M*                           inputProposalCovMatrix);
 
  //! Destructor
  ~uqMetropolisHastingsSGClass();
  //@}
 
  //! @name Statistical methods
  //@{
 //! Method to generate the chain. 
 /*! Requirement: the vector space 'm_vectorSpace' should have dimension equal to the size of a
  * vector in 'workingChain'. If the requirement is satisfied, this operation sets the size and 
  * the contents of 'workingChain' using the algorithm options set in the constructor. If not NULL,
  * 'workingLogLikelihoodValues' and 'workingLogTargetValues' are set accordingly. This operation 
  * currently implements the DRAM algorithm (Heikki Haario, Marko Laine, Antonietta Mira and
  * Eero Saksman, "DRAM: Efficient Adaptive MCMC", Statistics and Computing (2006), 16:339-354), 
  * as a translation of the core routine at the MCMC toolbox for MATLAB, available at 
  *  \htmlonly helios.fmi.fi/~lainema/mcmc/‎ \endhtmlonly (accessed in July 3rd, 2013). Indeed, 
  * the example available in examples/statisticalInverseProblem is closely related to the 
  * 'normal example' in the toolbox. Over time, though:
  <list type=number>
  <item> the whole set of QUESO classes took shape, focusing not only on Markov Chains, but on 
  statistical forward problems and model validation as well;
  <item> the interfaces to this Metropolis-Hastings class changed;
  <item> QUESO has parallel capabilities;
  <item> TK (transition kernel) class has been added in order to have both DRAM with Stochastic 
  Newton capabilities.
  </list>
  */
  void         generateSequence   (uqBaseVectorSequenceClass<P_V,P_M>& workingChain,
                                   uqScalarSequenceClass<double>*      workingLogLikelihoodValues,
                                   uqScalarSequenceClass<double>*      workingLogTargetValues);
  
  //! Gets information from the raw chain.
  void         getRawChainInfo    (uqMHRawChainInfoStruct& info) const;

   //@}
  
  //! @name I/O methods
  //@{
  //! TODO: Prints the sequence.
  /*! \todo: implement me!*/
  void   print                    (std::ostream& os) const;
  //@}

private:
  //! Reads the options values from the options input file.
  /*!  This method \b actually reads the options input file, such as the value for the Delayed Rejection
   * scales, the presence of Hessian covariance matrices and reads the user-provided initial proposal 
   * covariance matrix (alternative to local Hessians).*/
  void   commonConstructor        ();

  //!  This method generates the chain.
  /*! Given the size of the chain, the values of the first position in the chain. It checks if the position
   * is out of target pdf support, otherwise it calculates the value of both its likelihood and prior and 
   * adds the position to the chain. For the next positions, once they are generated, some tests are 
   * performed (such as unicity and the value of alpha) and the steps for the first position are repeated,
   * including the optional Delayed Rejection and the adaptive Metropolis (adaptation of covariance matrix) 
   * steps.*/
  void   generateFullChain        (const P_V&                          valuesOf1stPosition,
                                   unsigned int                        chainSize,
                                   uqBaseVectorSequenceClass<P_V,P_M>& workingChain,
                                   uqScalarSequenceClass<double>*      workingLogLikelihoodValues,
                                   uqScalarSequenceClass<double>*      workingLogTargetValues);
  
  //! This method reads the chain contents.
  void   readFullChain            (const std::string&                  inputFileName,
                                   const std::string&                  inputFileType,
                                   unsigned int                        chainSize,
                                   uqBaseVectorSequenceClass<P_V,P_M>& workingChain);
  
  //! This method updates the adapted covariance matrix
  /*! This function is called is the option to used adaptive Metropolis was chosen by the user 
   * (via options input file). It performs an adaptation of covariance matrix. */
  void   updateAdaptedCovMatrix   (const uqBaseVectorSequenceClass<P_V,P_M>&  subChain,
                                   unsigned int                               idOfFirstPositionInSubChain,
                                   double&                                    lastChainSize,
                                   P_V&                                       lastMean,
                                   P_M&                                       lastAdaptedCovMatrix);

  //! Calculates acceptance ration.
  /*! It is called by alpha(const std::vector<uqMarkovChainPositionDataClass<P_V>*>& inputPositions,
      const std::vector<unsigned int>& inputTKStageIds); */
  double alpha                    (const uqMarkovChainPositionDataClass<P_V>& x,
                                   const uqMarkovChainPositionDataClass<P_V>& y,
                                   unsigned int                               xStageId,
                                   unsigned int                               yStageId,
                                   double*                                    alphaQuotientPtr = NULL);
  
  //! Calculates acceptance ration.
  /*! The acceptance ratio is used to decide whether to accept or reject a candidate. */
  double alpha                    (const std::vector<uqMarkovChainPositionDataClass<P_V>*>& inputPositions,
                                   const std::vector<unsigned int                        >& inputTKStageIds);
  
  //! Decides whether or not to accept alpha.
  /*! If either alpha is negative or greater than one, its value will not be accepted.*/
  bool   acceptAlpha              (double                                     alpha);

  //! Writes information about the Markov chain in a file.
  /*! It writes down the alpha quotients, the number of rejected positions, number of positions out of
   * target support, the name of the components and the chain runtime.*/
  int    writeInfo                (const uqBaseVectorSequenceClass<P_V,P_M>&  workingChain,
                                   std::ofstream&                             ofsvar) const;

  const uqBaseEnvironmentClass&                     m_env;
  const uqVectorSpaceClass <P_V,P_M>&               m_vectorSpace;
  const uqBaseJointPdfClass<P_V,P_M>&               m_targetPdf;
        P_V                                         m_initialPosition;
        P_M                                         m_initialProposalCovMatrix;
        bool                                        m_nullInputProposalCovMatrix;
  const uqScalarFunctionSynchronizerClass<P_V,P_M>* m_targetPdfSynchronizer;

        uqBaseTKGroupClass<P_V,P_M>*                m_tk;
        unsigned int                                m_positionIdForDebugging;
        unsigned int                                m_stageIdForDebugging;
        std::vector<unsigned int>                   m_idsOfUniquePositions;
        std::vector<double>                         m_logTargets;
        std::vector<double>                         m_alphaQuotients;
        double                                      m_lastChainSize;
        P_V*                                        m_lastMean;
        P_M*                                        m_lastAdaptedCovMatrix;
        unsigned int                                m_numPositionsNotSubWritten;

        uqMHRawChainInfoStruct                      m_rawChainInfo;

        uqMhOptionsValuesClass                      m_alternativeOptionsValues;
        uqMetropolisHastingsSGOptionsClass*         m_optionsObj;
};
//! Prints the object \c obj, overloading an operator.
template<class P_V,class P_M>
std::ostream& operator<<(std::ostream& os, const uqMetropolisHastingsSGClass<P_V,P_M>& obj);

#include <uqMetropolisHastingsSG2.h>

// Default constructor -----------------------------
template<class P_V,class P_M>
uqMetropolisHastingsSGClass<P_V,P_M>::uqMetropolisHastingsSGClass(
  /*! Prefix                     */ const char*                         prefix,
  /*! Options (if no input file) */ const uqMhOptionsValuesClass*       alternativeOptionsValues, // dakota
  /*! The source RV              */ const uqBaseVectorRVClass<P_V,P_M>& sourceRv,
  /*! Initial chain position     */ const P_V&                          initialPosition,
  /*! Proposal cov. matrix       */ const P_M*                          inputProposalCovMatrix)
  :
  m_env                       (sourceRv.env()),
  m_vectorSpace               (sourceRv.imageSet().vectorSpace()),
  m_targetPdf                 (sourceRv.pdf()),
  m_initialPosition           (initialPosition),
  m_initialProposalCovMatrix  (m_vectorSpace.zeroVector()),
  m_nullInputProposalCovMatrix(inputProposalCovMatrix == NULL),
  m_targetPdfSynchronizer     (new uqScalarFunctionSynchronizerClass<P_V,P_M>(m_targetPdf,m_initialPosition)),
  m_tk                        (NULL),
  m_positionIdForDebugging    (0),
  m_stageIdForDebugging       (0),
  m_idsOfUniquePositions      (0),//0.),
  m_logTargets                (0),//0.),
  m_alphaQuotients            (0),//0.),
  m_lastChainSize             (0),
  m_lastMean                  (NULL),
  m_lastAdaptedCovMatrix      (NULL),
  m_numPositionsNotSubWritten (0),
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_alternativeOptionsValues  (NULL,NULL),
#else
  m_alternativeOptionsValues  (),
#endif
  m_optionsObj                (NULL)
{
  if (inputProposalCovMatrix != NULL) {
    m_initialProposalCovMatrix = *inputProposalCovMatrix;
  }
  if (alternativeOptionsValues) m_alternativeOptionsValues = *alternativeOptionsValues;
  if (m_env.optionsInputFileName() == "") {
    m_optionsObj = new uqMetropolisHastingsSGOptionsClass(m_env,prefix,m_alternativeOptionsValues);
  }
  else {
    m_optionsObj = new uqMetropolisHastingsSGOptionsClass(m_env,prefix);
    m_optionsObj->scanOptionsValues();
  }

  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_ov.m_totallyMute == false)) {
    *m_env.subDisplayFile() << "Entering uqMetropolisHastingsSGClass<P_V,P_M>::constructor(1)"
                            << ": prefix = " << prefix
                            << ", alternativeOptionsValues = " << alternativeOptionsValues
                            << ", m_env.optionsInputFileName() = " << m_env.optionsInputFileName()
                            << ", m_initialProposalCovMatrix = " << m_initialProposalCovMatrix
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO(sourceRv.imageSet().vectorSpace().dimLocal() != initialPosition.sizeLocal(),
                      m_env.worldRank(),
                      "uqMetropolisHastingsSGClass<P_V,P_M>::constructor(1)",
                      "'sourceRv' and 'initialPosition' should have equal dimensions");

  if (inputProposalCovMatrix) {
    UQ_FATAL_TEST_MACRO(sourceRv.imageSet().vectorSpace().dimLocal() != inputProposalCovMatrix->numRowsLocal(),
                        m_env.worldRank(),
                        "uqMetropolisHastingsSGClass<P_V,P_M>::constructor(1)",
                        "'sourceRv' and 'inputProposalCovMatrix' should have equal dimensions");
    UQ_FATAL_TEST_MACRO(inputProposalCovMatrix->numCols() != inputProposalCovMatrix->numRowsGlobal(),
                        m_env.worldRank(),
                        "uqMetropolisHastingsSGClass<P_V,P_M>::constructor(1)",
                        "'inputProposalCovMatrix' should be a square matrix");
  }

  commonConstructor();

  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_ov.m_totallyMute == false)) {
    *m_env.subDisplayFile() << "Leaving uqMetropolisHastingsSGClass<P_V,P_M>::constructor(1)"
                            << std::endl;
  }
}
// Constructor -------------------------------------
template<class P_V,class P_M>
uqMetropolisHastingsSGClass<P_V,P_M>::uqMetropolisHastingsSGClass(
  const uqMLSamplingLevelOptionsClass& mlOptions,
  const uqBaseVectorRVClass<P_V,P_M>&  sourceRv,
  const P_V&                           initialPosition, // KEY
  const P_M*                           inputProposalCovMatrix)
  :
  m_env                       (sourceRv.env()),
  m_vectorSpace               (sourceRv.imageSet().vectorSpace()),
  m_targetPdf                 (sourceRv.pdf()),
  m_initialPosition           (initialPosition),
  m_initialProposalCovMatrix  (m_vectorSpace.zeroVector()),
  m_nullInputProposalCovMatrix(inputProposalCovMatrix == NULL),
  m_targetPdfSynchronizer     (new uqScalarFunctionSynchronizerClass<P_V,P_M>(m_targetPdf,m_initialPosition)),
  m_tk                        (NULL),
  m_positionIdForDebugging    (0),
  m_stageIdForDebugging       (0),
  m_idsOfUniquePositions      (0),//0.),
  m_logTargets                (0),//0.),
  m_alphaQuotients            (0),//0.),
  m_lastChainSize             (0),
  m_lastMean                  (NULL),
  m_lastAdaptedCovMatrix      (NULL),
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_alternativeOptionsValues  (NULL,NULL),
#else
  m_alternativeOptionsValues  (),
#endif
  m_optionsObj                (new uqMetropolisHastingsSGOptionsClass(mlOptions))
{
  if (inputProposalCovMatrix != NULL) {
    m_initialProposalCovMatrix = *inputProposalCovMatrix;
    if ((m_env.subDisplayFile()                   ) &&
        (m_optionsObj->m_ov.m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::constructor(2)"
                              << ": just set m_initialProposalCovMatrix = " << m_initialProposalCovMatrix
                              << std::endl;
    }
  }
  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_ov.m_totallyMute == false)) {
    *m_env.subDisplayFile() << "Entering uqMetropolisHastingsSGClass<P_V,P_M>::constructor(2)"
                            << std::endl;
  }

  commonConstructor();

  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_ov.m_totallyMute == false)) {
    *m_env.subDisplayFile() << "Leaving uqMetropolisHastingsSGClass<P_V,P_M>::constructor(2)"
                            << std::endl;
  }
}
// Destructor ---------------------------------------
template<class P_V,class P_M>
uqMetropolisHastingsSGClass<P_V,P_M>::~uqMetropolisHastingsSGClass()
{
  //if (m_env.subDisplayFile()) {
  //  *m_env.subDisplayFile() << "Entering uqMetropolisHastingsSGClass<P_V,P_M>::destructor()"
  //                          << std::endl;
  //}

  if (m_lastAdaptedCovMatrix) delete m_lastAdaptedCovMatrix;
  if (m_lastMean)             delete m_lastMean;
  m_lastChainSize             = 0;
  m_rawChainInfo.reset();
  m_alphaQuotients.clear();
  m_logTargets.clear();
  m_positionIdForDebugging = 0;
  m_stageIdForDebugging    = 0;
  m_idsOfUniquePositions.clear();

  if (m_tk                   ) delete m_tk;
  if (m_targetPdfSynchronizer) delete m_targetPdfSynchronizer;
  if (m_optionsObj           ) delete m_optionsObj;

  //if (m_env.subDisplayFile()) {
  //  *m_env.subDisplayFile() << "Leaving uqMetropolisHastingsSGClass<P_V,P_M>::destructor()"
  //                          << std::endl;
  //}
}

// Private methods----------------------------------
template<class P_V,class P_M>
void
uqMetropolisHastingsSGClass<P_V,P_M>::commonConstructor()
{
  /////////////////////////////////////////////////////////////////
  // Instantiate the appropriate TK (transition kernel)
  /////////////////////////////////////////////////////////////////
  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_ov.m_totallyMute == false)) {
    *m_env.subDisplayFile() << "Entering uqMetropolisHastingsSGClass<P_V,P_M>::commonConstructor()"
                            << std::endl;
  }

  if (m_optionsObj->m_ov.m_initialPositionDataInputFileName != ".") { // palms
    std::set<unsigned int> tmpSet;
    tmpSet.insert(m_env.subId());
    m_initialPosition.subReadContents((m_optionsObj->m_ov.m_initialPositionDataInputFileName+"_sub"+m_env.subIdString()),
                                      m_optionsObj->m_ov.m_initialPositionDataInputFileType,
                                      tmpSet);
    if ((m_env.subDisplayFile()                   ) &&
        (m_optionsObj->m_ov.m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::commonConstructor()"
                              << ": just read initial position contents = " << m_initialPosition
                              << std::endl;
    }
  }

  std::vector<double> drScalesAll(m_optionsObj->m_ov.m_drScalesForExtraStages.size()+1,1.);
  for (unsigned int i = 1; i < (m_optionsObj->m_ov.m_drScalesForExtraStages.size()+1); ++i) {
    drScalesAll[i] = m_optionsObj->m_ov.m_drScalesForExtraStages[i-1];
  }
  if (m_optionsObj->m_ov.m_tkUseLocalHessian) { // sep2011
    m_tk = new uqHessianCovMatricesTKGroupClass<P_V,P_M>(m_optionsObj->m_prefix.c_str(),
                                                         m_vectorSpace,
                                                         drScalesAll,
                                                         *m_targetPdfSynchronizer);
    if ((m_env.subDisplayFile()                   ) &&
        (m_optionsObj->m_ov.m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::commonConstructor()"
                              << ": just instantiated a 'HessianCovMatrices' TK class"
                              << std::endl;
    }
  }
  else {
    if (m_optionsObj->m_ov.m_initialProposalCovMatrixDataInputFileName != ".") { // palms
      std::set<unsigned int> tmpSet;
      tmpSet.insert(m_env.subId());
      m_initialProposalCovMatrix.subReadContents((m_optionsObj->m_ov.m_initialProposalCovMatrixDataInputFileName+"_sub"+m_env.subIdString()),
                                                 m_optionsObj->m_ov.m_initialProposalCovMatrixDataInputFileType,
                                                 tmpSet);
      if ((m_env.subDisplayFile()                   ) &&
          (m_optionsObj->m_ov.m_totallyMute == false)) {
        *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::commonConstructor()"
                                << ": just read initial proposal cov matrix contents = " << m_initialProposalCovMatrix
                                << std::endl;
      }
    }
    else {
      UQ_FATAL_TEST_MACRO(m_nullInputProposalCovMatrix,
                          m_env.worldRank(),
                          "uqMetropolisHastingsSGClass<P_V,P_M>::commonConstructor()",
                          "proposal cov matrix should have been passed by user, since, according to the input algorithm options, local Hessians will not be used in the proposal");
    }

    m_tk = new uqScaledCovMatrixTKGroupClass<P_V,P_M>(m_optionsObj->m_prefix.c_str(),
                                                      m_vectorSpace,
                                                      drScalesAll,
                                                      m_initialProposalCovMatrix);
    if ((m_env.subDisplayFile()                   ) &&
        (m_optionsObj->m_ov.m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::commonConstructor()"
                              << ": just instantiated a 'ScaledCovMatrix' TK class"
                              << std::endl;
    }
  }

  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_ov.m_totallyMute == false)) {
    *m_env.subDisplayFile() << "Leaving uqMetropolisHastingsSGClass<P_V,P_M>::commonConstructor()"
                            << std::endl;
  }
  return;
}
//--------------------------------------------------
template<class P_V,class P_M>
double
uqMetropolisHastingsSGClass<P_V,P_M>::alpha(
  const uqMarkovChainPositionDataClass<P_V>& x,
  const uqMarkovChainPositionDataClass<P_V>& y,
  unsigned int                               xStageId,
  unsigned int                               yStageId,
  double*                                    alphaQuotientPtr)
{
  double alphaQuotient = 0.;
  if ((x.outOfTargetSupport() == false) &&
      (y.outOfTargetSupport() == false)) {
    if ((x.logTarget() == -INFINITY) ||
        (x.logTarget() ==  INFINITY) ||
        ( (boost::math::isnan)(x.logTarget())      )) {
      std::cerr << "WARNING In uqMetropolisHastingsSGClass<P_V,P_M>::alpha(x,y)"
                << ", worldRank "       << m_env.worldRank()
                << ", fullRank "        << m_env.fullRank()
                << ", subEnvironment "  << m_env.subId()
                << ", subRank "         << m_env.subRank()
                << ", inter0Rank "      << m_env.inter0Rank()
                << ", positionId = "    << m_positionIdForDebugging
                << ", stageId = "       << m_stageIdForDebugging
                << ": x.logTarget() = " << x.logTarget()
                << ", x.values() = "    << x.vecValues()
                << ", y.values() = "    << y.vecValues()
                << std::endl;
    }
    else if ((y.logTarget() == -INFINITY           ) ||
             (y.logTarget() ==  INFINITY           ) ||
             ( (boost::math::isnan)(y.logTarget()) )) {
      std::cerr << "WARNING In uqMetropolisHastingsSGClass<P_V,P_M>::alpha(x,y)"
                << ", worldRank "       << m_env.worldRank()
                << ", fullRank "        << m_env.fullRank()
                << ", subEnvironment "  << m_env.subId()
                << ", subRank "         << m_env.subRank()
                << ", inter0Rank "      << m_env.inter0Rank()
                << ", positionId = "    << m_positionIdForDebugging
                << ", stageId = "       << m_stageIdForDebugging
                << ": y.logTarget() = " << y.logTarget()
                << ", x.values() = "    << x.vecValues()
                << ", y.values() = "    << y.vecValues()
                << std::endl;
    }
    else {
      double yLogTargetToUse = y.logTarget();
      if (m_tk->symmetric()) {
        alphaQuotient = std::exp(yLogTargetToUse - x.logTarget());
        if ((m_env.subDisplayFile()                   ) &&
            (m_env.displayVerbosity() >= 3            ) &&
            (m_optionsObj->m_ov.m_totallyMute == false)) {
          *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::alpha(x,y)"
                                 << ": symmetric proposal case"
                                 << ", x = "               << x.vecValues()
                                 << ", y = "               << y.vecValues()
                                 << ", yLogTargetToUse = " << yLogTargetToUse
                                 << ", x.logTarget() = "   << x.logTarget()
                                 << ", alpha = "           << alphaQuotient
                                 << std::endl;
        }
      }
      else {
#ifdef QUESO_EXPECTS_LN_LIKELIHOOD_INSTEAD_OF_MINUS_2_LN
        double qyx = m_tk->rv(yStageId).pdf().lnValue(x.vecValues(),NULL,NULL,NULL,NULL);
#else
        double qyx = -.5 * m_tk->rv(yStageId).pdf().lnValue(x.vecValues(),NULL,NULL,NULL,NULL);
#endif
        if ((m_env.subDisplayFile()                   ) &&
            (m_env.displayVerbosity() >= 10           ) &&
            (m_optionsObj->m_ov.m_totallyMute == false)) {
          const uqGaussianJointPdfClass<P_V,P_M>* pdfYX = dynamic_cast< const uqGaussianJointPdfClass<P_V,P_M>* >(&(m_tk->rv(yStageId).pdf()));
          *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::alpha(x,y)"
                                 << ", rvYX.lawExpVector = " << pdfYX->lawExpVector()
                                 << ", rvYX.lawVarVector = " << pdfYX->lawVarVector()
                                 << ", rvYX.lawCovMatrix = " << pdfYX->lawCovMatrix()
                                 << std::endl;
        }
#ifdef QUESO_EXPECTS_LN_LIKELIHOOD_INSTEAD_OF_MINUS_2_LN
        double qxy = m_tk->rv(xStageId).pdf().lnValue(y.vecValues(),NULL,NULL,NULL,NULL);
#else
        double qxy = -.5 * m_tk->rv(xStageId).pdf().lnValue(y.vecValues(),NULL,NULL,NULL,NULL);
#endif
        if ((m_env.subDisplayFile()                   ) &&
            (m_env.displayVerbosity() >= 10           ) &&
            (m_optionsObj->m_ov.m_totallyMute == false)) {
          const uqGaussianJointPdfClass<P_V,P_M>* pdfXY = dynamic_cast< const uqGaussianJointPdfClass<P_V,P_M>* >(&(m_tk->rv(xStageId).pdf()));
          *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::alpha(x,y)"
                                 << ", rvXY.lawExpVector = " << pdfXY->lawExpVector()
                                 << ", rvXY.lawVarVector = " << pdfXY->lawVarVector()
                                 << ", rvXY.lawCovMatrix = " << pdfXY->lawCovMatrix()
                                 << std::endl;
        }
        alphaQuotient = std::exp(yLogTargetToUse +
                                 qyx -
                                 x.logTarget() -
                                 qxy);
        if ((m_env.subDisplayFile()                   ) &&
            (m_env.displayVerbosity() >= 3            ) &&
            (m_optionsObj->m_ov.m_totallyMute == false)) {
          *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::alpha(x,y)"
                                 << ": asymmetric proposal case"
                                 << ", xStageId = "        << xStageId
                                 << ", yStageId = "        << yStageId
                                 << ", x = "               << x.vecValues()
                                 << ", y = "               << y.vecValues()
                                 << ", yLogTargetToUse = " << yLogTargetToUse
                                 << ", q(y,x) = "          << qyx
                                 << ", x.logTarget() = "   << x.logTarget()
                                 << ", q(x,y) = "          << qxy
                                 << ", alpha = "           << alphaQuotient
                                 << std::endl;
        }
      }
    } // protection logic against logTarget values
  }
  else {
    if ((m_env.subDisplayFile()                   ) &&
        (m_env.displayVerbosity() >= 10           ) &&
        (m_optionsObj->m_ov.m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::alpha(x,y)"
                             << ": x.outOfTargetSupport = " << x.outOfTargetSupport()
                             << ", y.outOfTargetSupport = " << y.outOfTargetSupport()
                             << std::endl;
    }
  }
  if (alphaQuotientPtr != NULL) *alphaQuotientPtr = alphaQuotient;

  return std::min(1.,alphaQuotient);
}
//--------------------------------------------------
template<class P_V,class P_M>
double
uqMetropolisHastingsSGClass<P_V,P_M>::alpha(
  const std::vector<uqMarkovChainPositionDataClass<P_V>*>& inputPositionsData,
  const std::vector<unsigned int                        >& inputTKStageIds)
{
  unsigned int inputSize = inputPositionsData.size();
  if ((m_env.subDisplayFile()                   ) &&
      (m_env.displayVerbosity() >= 10           ) &&
      (m_optionsObj->m_ov.m_totallyMute == false)) {
    *m_env.subDisplayFile() << "Entering uqMetropolisHastingsSGClass<P_V,P_M>::alpha(vec)"
                           << ", inputSize = " << inputSize
                           << std::endl;
  }
  UQ_FATAL_TEST_MACRO((inputSize < 2),
                      m_env.worldRank(),
                      "uqMetropolisHastingsSGClass<P_V,P_M>::alpha(vec)",
                      "inputPositionsData has size < 2");

  // If necessary, return 0. right away
  if (inputPositionsData[0          ]->outOfTargetSupport()) return 0.;
  if (inputPositionsData[inputSize-1]->outOfTargetSupport()) return 0.;

  if ((inputPositionsData[0]->logTarget() == -INFINITY           ) ||
      (inputPositionsData[0]->logTarget() ==  INFINITY           ) ||
      ( (boost::math::isnan)(inputPositionsData[0]->logTarget()) )) {
    std::cerr << "WARNING In uqMetropolisHastingsSGClass<P_V,P_M>::alpha(vec)"
              << ", worldRank "      << m_env.worldRank()
              << ", fullRank "       << m_env.fullRank()
              << ", subEnvironment " << m_env.subId()
              << ", subRank "        << m_env.subRank()
              << ", inter0Rank "     << m_env.inter0Rank()
              << ", positionId = "   << m_positionIdForDebugging
              << ", stageId = "      << m_stageIdForDebugging
              << ": inputSize = "    << inputSize
              << ", inputPositionsData[0]->logTarget() = " << inputPositionsData[0]->logTarget()
              << ", [0]->values() = "                      << inputPositionsData[0]->vecValues()
              << ", [inputSize - 1]->values() = "          << inputPositionsData[inputSize-1]->vecValues()
              << std::endl;
    return 0.;
  }
  else if ((inputPositionsData[inputSize - 1]->logTarget() == -INFINITY           ) ||
           (inputPositionsData[inputSize - 1]->logTarget() ==  INFINITY           ) ||
           ( (boost::math::isnan)(inputPositionsData[inputSize - 1]->logTarget()) )) {
    std::cerr << "WARNING In uqMetropolisHastingsSGClass<P_V,P_M>::alpha(vec)"
              << ", worldRank "      << m_env.worldRank()
              << ", fullRank "       << m_env.fullRank()
              << ", subEnvironment " << m_env.subId()
              << ", subRank "        << m_env.subRank()
              << ", inter0Rank "     << m_env.inter0Rank()
              << ", positionId = "   << m_positionIdForDebugging
              << ", stageId = "      << m_stageIdForDebugging
              << ": inputSize = "    << inputSize
              << ", inputPositionsData[inputSize - 1]->logTarget() = " << inputPositionsData[inputSize-1]->logTarget()
              << ", [0]->values() = "                                  << inputPositionsData[0]->vecValues()
              << ", [inputSize - 1]->values() = "                      << inputPositionsData[inputSize-1]->vecValues()
              << std::endl;
    return 0.;
  }

  // If inputSize is 2, recursion is not needed
  if (inputSize == 2) return this->alpha(*(inputPositionsData[0            ]),
                                         *(inputPositionsData[inputSize - 1]),
                                         inputTKStageIds[0],
                                         inputTKStageIds[inputSize-1]);

  // Prepare two vectors of positions
  std::vector<uqMarkovChainPositionDataClass<P_V>*>         positionsData  (inputSize,NULL);
  std::vector<uqMarkovChainPositionDataClass<P_V>*> backwardPositionsData  (inputSize,NULL);

  std::vector<unsigned int                        >         tkStageIds     (inputSize,0);
  std::vector<unsigned int                        > backwardTKStageIds     (inputSize,0);

  std::vector<unsigned int                        >         tkStageIdsLess1(inputSize,0);
  std::vector<unsigned int                        > backwardTKStageIdsLess1(inputSize,0);

  for (unsigned int i = 0; i < inputSize; ++i) {
            positionsData  [i] = inputPositionsData[i];
    backwardPositionsData  [i] = inputPositionsData[inputSize-i-1];

            tkStageIds     [i] = inputTKStageIds   [i];
    backwardTKStageIds     [i] = inputTKStageIds   [inputSize-i-1];

            tkStageIdsLess1[i] = inputTKStageIds   [i];
    backwardTKStageIdsLess1[i] = inputTKStageIds   [inputSize-i-1];
  }

          tkStageIdsLess1.pop_back();
  backwardTKStageIdsLess1.pop_back();

  // Initialize cumulative variables
  double logNumerator      = 0.;
  double logDenominator    = 0.;
  double alphasNumerator   = 1.;
  double alphasDenominator = 1.;

  // Compute cumulative variables
  const P_V& _lastTKPosition         = m_tk->preComputingPosition(        tkStageIds[inputSize-1]);
  const P_V& _lastBackwardTKPosition = m_tk->preComputingPosition(backwardTKStageIds[inputSize-1]);

#ifdef QUESO_EXPECTS_LN_LIKELIHOOD_INSTEAD_OF_MINUS_2_LN
  double numContrib = m_tk->rv(backwardTKStageIdsLess1).pdf().lnValue(_lastBackwardTKPosition,NULL,NULL,NULL,NULL);
  double denContrib = m_tk->rv(        tkStageIdsLess1).pdf().lnValue(_lastTKPosition        ,NULL,NULL,NULL,NULL);
#else
  double numContrib = -.5 * m_tk->rv(backwardTKStageIdsLess1).pdf().lnValue(_lastBackwardTKPosition,NULL,NULL,NULL,NULL);
  double denContrib = -.5 * m_tk->rv(        tkStageIdsLess1).pdf().lnValue(_lastTKPosition        ,NULL,NULL,NULL,NULL);
#endif
  if ((m_env.subDisplayFile()                   ) &&
      (m_env.displayVerbosity() >= 10           ) &&
      (m_optionsObj->m_ov.m_totallyMute == false)) {
    *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::alpha(vec)"
                           << ", inputSize = "  << inputSize
                           << ", before loop"
                           << ": numContrib = " << numContrib
                           << ", denContrib = " << denContrib
                           << std::endl;
  }
  logNumerator   += numContrib;
  logDenominator += denContrib;

  for (unsigned int i = 0; i < (inputSize-2); ++i) { // That is why size must be >= 2
            positionsData.pop_back();
    backwardPositionsData.pop_back();

    const P_V& lastTKPosition         = m_tk->preComputingPosition(        tkStageIds[inputSize-2-i]);
    const P_V& lastBackwardTKPosition = m_tk->preComputingPosition(backwardTKStageIds[inputSize-2-i]);

            tkStageIds.pop_back();
    backwardTKStageIds.pop_back();

            tkStageIdsLess1.pop_back();
    backwardTKStageIdsLess1.pop_back();

#ifdef QUESO_EXPECTS_LN_LIKELIHOOD_INSTEAD_OF_MINUS_2_LN
    numContrib = m_tk->rv(backwardTKStageIdsLess1).pdf().lnValue(lastBackwardTKPosition,NULL,NULL,NULL,NULL);
    denContrib = m_tk->rv(        tkStageIdsLess1).pdf().lnValue(lastTKPosition        ,NULL,NULL,NULL,NULL);
#else
    numContrib = -.5 * m_tk->rv(backwardTKStageIdsLess1).pdf().lnValue(lastBackwardTKPosition,NULL,NULL,NULL,NULL);
    denContrib = -.5 * m_tk->rv(        tkStageIdsLess1).pdf().lnValue(lastTKPosition        ,NULL,NULL,NULL,NULL);
#endif
    if ((m_env.subDisplayFile()                   ) &&
        (m_env.displayVerbosity() >= 10           ) &&
        (m_optionsObj->m_ov.m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::alpha(vec)"
                             << ", inputSize = "  << inputSize
                             << ", in loop, i = " << i
                             << ": numContrib = " << numContrib
                             << ", denContrib = " << denContrib
                             << std::endl;
    }
    logNumerator   += numContrib;
    logDenominator += denContrib;

    alphasNumerator   *= (1. - this->alpha(backwardPositionsData,backwardTKStageIds));
    alphasDenominator *= (1. - this->alpha(        positionsData,        tkStageIds));
  }

  double numeratorLogTargetToUse = backwardPositionsData[0]->logTarget();
  numContrib = numeratorLogTargetToUse;
  denContrib = positionsData[0]->logTarget();
  if ((m_env.subDisplayFile()                   ) &&
      (m_env.displayVerbosity() >= 10           ) &&
      (m_optionsObj->m_ov.m_totallyMute == false)) {
    *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::alpha(vec)"
                           << ", inputSize = "  << inputSize
                           << ", after loop"
                           << ": numContrib = " << numContrib
                           << ", denContrib = " << denContrib
                           << std::endl;
  }
  logNumerator   += numContrib;
  logDenominator += denContrib;

  if ((m_env.subDisplayFile()                   ) &&
      (m_env.displayVerbosity() >= 10           ) &&
      (m_optionsObj->m_ov.m_totallyMute == false)) {
    *m_env.subDisplayFile() << "Leaving uqMetropolisHastingsSGClass<P_V,P_M>::alpha(vec)"
                           << ", inputSize = "         << inputSize
                           << ": alphasNumerator = "   << alphasNumerator
                           << ", alphasDenominator = " << alphasDenominator
                           << ", logNumerator = "      << logNumerator
                           << ", logDenominator = "    << logDenominator
                           << std::endl;
  }

  // Return result
  return std::min(1.,(alphasNumerator/alphasDenominator)*std::exp(logNumerator-logDenominator));
}
//--------------------------------------------------
template<class P_V,class P_M>
bool
uqMetropolisHastingsSGClass<P_V,P_M>::acceptAlpha(double alpha)
{
  bool result = false;

  if      (alpha <= 0.                                ) result = false;
  else if (alpha >= 1.                                ) result = true;
  else if (alpha >= m_env.rngObject()->uniformSample()) result = true;
  else                                                  result = false;

  return result;
}
//--------------------------------------------------
template<class P_V,class P_M>
int
uqMetropolisHastingsSGClass<P_V,P_M>::writeInfo(
  const uqBaseVectorSequenceClass<P_V,P_M>& workingChain,
  std::ofstream&                            ofsvar) const
{
  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_ov.m_totallyMute == false)) {
    *m_env.subDisplayFile() << "\n"
                            << "\n-----------------------------------------------------"
                            << "\n Writing more information about the Markov chain " << workingChain.name() << " to output file ..."
                            << "\n-----------------------------------------------------"
                            << "\n"
                            << std::endl;
  }

  int iRC = UQ_OK_RC;

  if (m_optionsObj->m_ov.m_rawChainGenerateExtra) {
    // Write m_logTargets
    ofsvar << m_optionsObj->m_prefix << "logTargets_sub" << m_env.subIdString() << " = zeros(" << m_logTargets.size()
           << ","                    << 1
           << ");"
           << std::endl;
    ofsvar << m_optionsObj->m_prefix << "logTargets_sub" << m_env.subIdString() << " = [";
    for (unsigned int i = 0; i < m_logTargets.size(); ++i) {
      ofsvar << m_logTargets[i]
             << std::endl;
    }
    ofsvar << "];\n";

    // Write m_alphaQuotients
    ofsvar << m_optionsObj->m_prefix << "alphaQuotients_sub" << m_env.subIdString() << " = zeros(" << m_alphaQuotients.size()
           << ","                    << 1
           << ");"
           << std::endl;
    ofsvar << m_optionsObj->m_prefix << "alphaQuotients_sub" << m_env.subIdString() << " = [";
    for (unsigned int i = 0; i < m_alphaQuotients.size(); ++i) {
      ofsvar << m_alphaQuotients[i]
             << std::endl;
    }
    ofsvar << "];\n";
  }

  // Write number of rejections
  ofsvar << m_optionsObj->m_prefix << "rejected = " << (double) m_rawChainInfo.numRejections/(double) (workingChain.subSequenceSize()-1)
         << ";\n"
         << std::endl;

  if (false) { // Don't see need for code below. Let it there though, compiling, in case it is needed in the future.
    // Write names of components
    ofsvar << m_optionsObj->m_prefix << "componentNames = {";
    m_vectorSpace.printComponentsNames(ofsvar,false);
    ofsvar << "};\n";

    // Write number of out of target support
    ofsvar << m_optionsObj->m_prefix << "outTargetSupport = " << (double) m_rawChainInfo.numOutOfTargetSupport/(double) (workingChain.subSequenceSize()-1)
           << ";\n"
           << std::endl;

    // Write chain run time
    ofsvar << m_optionsObj->m_prefix << "runTime = " << m_rawChainInfo.runTime
           << ";\n"
           << std::endl;
  }

  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_ov.m_totallyMute == false)) {
    *m_env.subDisplayFile() << "\n-----------------------------------------------------"
                            << "\n Finished writing more information about the Markov chain " << workingChain.name()
                            << "\n-----------------------------------------------------"
                            << "\n"
                            << std::endl;
  }

  return iRC;
}
//--------------------------------------------------
// Operator declared outside class definition-------
//--------------------------------------------------

template<class P_V,class P_M>
std::ostream& operator<<(std::ostream& os, const uqMetropolisHastingsSGClass<P_V,P_M>& obj)
{
  obj.print(os);

  return os;
}

}  // End namespace QUESO

#endif // __UQ_MH_SG1_H__
