/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of the QUESO Library (Quantification of Uncertainty
 * for Estimation, Simulation and Optimization).
 *
 * QUESO is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QUESO is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QUESO. If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * $Id$
 *
 * Brief description of this file: 
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __UQ_MAC_SG1_H__
#define __UQ_MAC_SG1_H__

#include <uqMarkovChainSGOptions.h>
#include <uqTKGroup.h>
#include <uqVectorRV.h>
#include <uqVectorSpace.h>
#include <uqMarkovChainPositionData.h>
#include <uqScalarFunctionSynchronizer.h>
#include <uqSequenceOfVectors.h>
#include <uqArrayOfSequences.h>
#include <sys/time.h>
#include <fstream>

/*! A templated class that represents a Markov Chain generator. 'SG' stands for 'Sequence Generator'.
 */
template <class P_V,class P_M>
class uqMarkovChainSGClass
{
public:

  /*! Constructor: */
  uqMarkovChainSGClass(/*! Prefix                 */ const char*                         prefix,
                       /*! The source rv          */ const uqBaseVectorRVClass<P_V,P_M>& sourceRv,
                       /*! Initial chain position */ const P_V&                          initialPosition,
                       /*! Proposal cov. matrix   */ const P_M*                          inputProposalCovMatrix);
  uqMarkovChainSGClass(/*! Options                */ const uqMLSamplingLevelOptionsClass& inputOptions,
                       /*! The source rv          */ const uqBaseVectorRVClass<P_V,P_M>&  sourceRv,
                       /*! Initial chain position */ const P_V&                           initialPosition,
                       /*! Proposal cov. matrix   */ const P_M*                           inputProposalCovMatrix);
  /*! Destructor: */
 ~uqMarkovChainSGClass();

  /*! Operation to generate the chain */
  void   generateSequence           (uqBaseVectorSequenceClass<P_V,P_M>& workingChain,
                                     uqScalarSequenceClass<double>*      workingTargetValues);
  void   checkTheParallelEnvironment();

  void   print                      (std::ostream& os) const;


private:
  void   commonConstructor        ();
//void   proc0GenerateSequence    (uqBaseVectorSequenceClass<P_V,P_M>& workingChain); /*! */
  void   resetChainAndRelatedInfo ();

  void   generateWhiteNoiseChain  (      unsigned int                                   chainSize,
                                   uqBaseVectorSequenceClass<P_V,P_M>&                  workingChain);
  void   generateUniformChain     (      unsigned int                                   chainSize,
                                   uqBaseVectorSequenceClass<P_V,P_M>&                  workingChain);
  void   generateFullChain        (const P_V&                                           valuesOf1stPosition,
                                         unsigned int                                   chainSize,
                                   uqBaseVectorSequenceClass<P_V,P_M>&                  workingChain,
                                   uqScalarSequenceClass<double>*                       workingTargetValues);
  void   readFullChain            (const std::string&                                   inputFileName,
                                         unsigned int                                   chainSize,
                                   uqBaseVectorSequenceClass<P_V,P_M>&                  workingChain);
  void   updateAdaptedCovMatrix   (const uqBaseVectorSequenceClass<P_V,P_M>&            subChain,
                                   unsigned int                                         idOfFirstPositionInSubChain,
                                   double&                                              lastChainSize,
                                   P_V&                                                 lastMean,
                                   P_M&                                                 lastAdaptedCovMatrix);

#ifdef UQ_USES_TK_CLASS
#else
  int    computeInitialCholFactors();
  void   updateTK                 ();
  double logProposal              (const uqMarkovChainPositionDataClass<P_V>&               x,
                                   const uqMarkovChainPositionDataClass<P_V>&               y,
                                   unsigned int                                             idOfProposalCovMatrix);
  double logProposal              (const std::vector<uqMarkovChainPositionDataClass<P_V>*>& inputPositions);
#endif
  double alpha                    (const uqMarkovChainPositionDataClass<P_V>&               x,
                                   const uqMarkovChainPositionDataClass<P_V>&               y,
                                   unsigned int                                             xStageId,
                                   unsigned int                                             yStageId,
                                   double*                                                  alphaQuotientPtr = NULL);
  double alpha                    (const std::vector<uqMarkovChainPositionDataClass<P_V>*>& inputPositions,
                                   const std::vector<unsigned int                        >& inputTKStageIds);
  bool   acceptAlpha              (double                                                   alpha);

  int    writeInfo                (const uqBaseVectorSequenceClass<P_V,P_M>&                workingChain,
                                   std::ofstream&                                           ofsvar) const;
                                 //const P_M*                                               mahalanobisMatrix = NULL,
                                 //bool                                                     applyMahalanobisInvert = true) const;

  const uqBaseEnvironmentClass&                     m_env;
  const uqVectorSpaceClass <P_V,P_M>&               m_vectorSpace;
  const uqBaseJointPdfClass<P_V,P_M>&               m_targetPdf;
        P_V                                         m_initialPosition;
  const P_M*                                        m_initialProposalCovMatrix;
        bool                                        m_nullInputProposalCovMatrix;
  const uqScalarFunctionSynchronizerClass<P_V,P_M>* m_targetPdfSynchronizer;

        uqBaseTKGroupClass<P_V,P_M>*       m_tk;
#ifdef UQ_USES_TK_CLASS
#else
        bool                               m_tkIsSymmetric;
        std::vector<P_M*>                  m_lowerCholProposalCovMatrices;
        std::vector<P_M*>                  m_proposalCovMatrices;
#ifdef UQ_MAC_SG_REQUIRES_INVERTED_COV_MATRICES
        std::vector<P_M*>                  m_upperCholProposalPrecMatrices;
        std::vector<P_M*>                  m_proposalPrecMatrices;
#endif
#endif
        std::vector<unsigned int>          m_idsOfUniquePositions;
        std::vector<double>                m_logTargets;
        std::vector<double>                m_alphaQuotients;
        double                             m_rawChainRunTime;
        unsigned int                       m_numRejections;
        unsigned int                       m_numOutOfTargetSupport;
        double                             m_lastChainSize;
        P_V*                               m_lastMean;
        P_M*                               m_lastAdaptedCovMatrix;

        uqMarkovChainSGOptionsClass        m_options;
};

template<class P_V,class P_M>
std::ostream& operator<<(std::ostream& os, const uqMarkovChainSGClass<P_V,P_M>& obj);

#include <uqMarkovChainSG2.h>

template<class P_V,class P_M>
uqMarkovChainSGClass<P_V,P_M>::uqMarkovChainSGClass(
  const char*                         prefix,
  const uqBaseVectorRVClass<P_V,P_M>& sourceRv,
  const P_V&                          initialPosition,
  const P_M*                          inputProposalCovMatrix)
  :
  m_env                          (sourceRv.env()),
  m_vectorSpace                  (sourceRv.imageSet().vectorSpace()),
  m_targetPdf                    (sourceRv.pdf()),
  m_initialPosition              (initialPosition),
  m_initialProposalCovMatrix     (inputProposalCovMatrix),
  m_nullInputProposalCovMatrix   (inputProposalCovMatrix == NULL),
  m_targetPdfSynchronizer        (new uqScalarFunctionSynchronizerClass<P_V,P_M>(m_targetPdf,m_initialPosition)),
  m_tk                           (NULL),
#ifdef UQ_USES_TK_CLASS
#else
  m_tkIsSymmetric                (true),
  m_lowerCholProposalCovMatrices (1),//NULL),
  m_proposalCovMatrices          (1),//NULL),
#ifdef UQ_MAC_SG_REQUIRES_INVERTED_COV_MATRICES
  m_upperCholProposalPrecMatrices(1),//NULL),
  m_proposalPrecMatrices         (1),//NULL),
#endif
#endif
  m_idsOfUniquePositions         (0),//0.),
  m_logTargets                   (0),//0.),
  m_alphaQuotients               (0),//0.),
  m_rawChainRunTime              (0.),
  m_numRejections                (0),
  m_numOutOfTargetSupport        (0),
  m_lastChainSize                (0),
  m_lastMean                     (NULL),
  m_lastAdaptedCovMatrix         (NULL),
  m_options                      (m_env,prefix)
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Entering uqMarkovChainSGClass<P_V,P_M>::constructor(1)"
                            << std::endl;
  }

  m_options.scanOptionsValues();

  commonConstructor();

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Leaving uqMarkovChainSGClass<P_V,P_M>::constructor(1)"
                           << std::endl;
  }
}

template<class P_V,class P_M>
uqMarkovChainSGClass<P_V,P_M>::uqMarkovChainSGClass(
  const uqMLSamplingLevelOptionsClass& inputOptions,
  const uqBaseVectorRVClass<P_V,P_M>&  sourceRv,
  const P_V&                           initialPosition,
  const P_M*                           inputProposalCovMatrix)
  :
  m_env                          (sourceRv.env()),
  m_vectorSpace                  (sourceRv.imageSet().vectorSpace()),
  m_targetPdf                    (sourceRv.pdf()),
  m_initialPosition              (initialPosition),
  m_initialProposalCovMatrix     (inputProposalCovMatrix),
  m_nullInputProposalCovMatrix   (inputProposalCovMatrix == NULL),
  m_targetPdfSynchronizer        (new uqScalarFunctionSynchronizerClass<P_V,P_M>(m_targetPdf,m_initialPosition)),
  m_tk                           (NULL),
#ifdef UQ_USES_TK_CLASS
#else
  m_tkIsSymmetric                (true),
  m_lowerCholProposalCovMatrices (1),//NULL),
  m_proposalCovMatrices          (1),//NULL),
#ifdef UQ_MAC_SG_REQUIRES_INVERTED_COV_MATRICES
  m_upperCholProposalPrecMatrices(1),//NULL),
  m_proposalPrecMatrices         (1),//NULL),
#endif
#endif
  m_idsOfUniquePositions         (0),//0.),
  m_logTargets                   (0),//0.),
  m_alphaQuotients               (0),//0.),
  m_rawChainRunTime              (0.),
  m_numRejections                (0),
  m_numOutOfTargetSupport        (0),
  m_lastChainSize                (0),
  m_lastMean                     (NULL),
  m_lastAdaptedCovMatrix         (NULL),
  m_options                      (m_env,inputOptions)
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Entering uqMarkovChainSGClass<P_V,P_M>::constructor(2)"
                            << std::endl;
  }

  commonConstructor();

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Leaving uqMarkovChainSGClass<P_V,P_M>::constructor(2)"
                           << std::endl;
  }
}

template<class P_V,class P_M>
void
uqMarkovChainSGClass<P_V,P_M>::commonConstructor()
{
  /////////////////////////////////////////////////////////////////
  // Instantiate the appropriate TK
  /////////////////////////////////////////////////////////////////
#ifdef UQ_USES_TK_CLASS
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "In uqMarkovChainSGClass<P_V,P_M>::constructor()"
                           << ": running with UQ_USES_TK_CLASS flag defined"
                           << std::endl;
  }
#else
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "In uqMarkovChainSGClass<P_V,P_M>::constructor()"
                           << ": running with UQ_USES_TK_CLASS flag undefined"
                           << std::endl;
  }
#endif
  if (m_options.m_tkUseLocalHessian) {
    m_tk = new uqHessianCovMatricesTKGroupClass<P_V,P_M>(m_options.m_prefix.c_str(),
                                                         m_vectorSpace,
                                                         m_options.m_drScalesForCovMatrices,
                                                         *m_targetPdfSynchronizer);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqMarkovChainSGClass<P_V,P_M>::constructor()"
                             << ": just instantiated a 'HessianCovMatrices' TK class"
                             << std::endl;
    }
  }
  else {
    UQ_FATAL_TEST_MACRO((m_initialProposalCovMatrix == NULL),
                        m_env.fullRank(),
                        "uqMarkovChainSGClass<P_V,P_M>::constructor()",
                        "proposal cov matrix should have been passed by user, since, according to the input algorithm options, local Hessians will not be used in the proposal");

    m_tk = new uqScaledCovMatrixTKGroupClass<P_V,P_M>(m_options.m_prefix.c_str(),
                                                      m_vectorSpace,
                                                      m_options.m_drScalesForCovMatrices,
                                                      *m_initialProposalCovMatrix);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqMarkovChainSGClass<P_V,P_M>::constructor()"
                              << ": just instantiated a 'ScaledCovMatrix' TK class"
                              << std::endl;
    }
  }

  return;
}

template<class P_V,class P_M>
uqMarkovChainSGClass<P_V,P_M>::~uqMarkovChainSGClass()
{
  //if (m_env.subDisplayFile()) {
  //  *m_env.subDisplayFile() << "Entering uqMarkovChainSGClass<P_V,P_M>::destructor()"
  //                          << std::endl;
  //}

  resetChainAndRelatedInfo();
  //if (m_lowerCholProposalCovMatrices[0]) delete m_lowerCholProposalCovMatrices[0]; // Loop inside 'resetChainAndRelatedInfo()' deletes just from position '1' on

  if (m_tk                        ) delete m_tk;
  if (m_targetPdfSynchronizer     ) delete m_targetPdfSynchronizer;
  if (m_nullInputProposalCovMatrix) delete m_initialProposalCovMatrix;

  //if (m_env.subDisplayFile()) {
  //  *m_env.subDisplayFile() << "Leaving uqMarkovChainSGClass<P_V,P_M>::destructor()"
  //                          << std::endl;
  //}
}

template<class P_V,class P_M>
void
uqMarkovChainSGClass<P_V,P_M>::resetChainAndRelatedInfo()
{
  if (m_lastAdaptedCovMatrix) delete m_lastAdaptedCovMatrix;
  if (m_lastMean)             delete m_lastMean;
  m_lastChainSize         = 0;
  m_numOutOfTargetSupport = 0;
  m_rawChainRunTime       = 0.;
  m_numRejections         = 0;
  m_alphaQuotients.clear();
  m_logTargets.clear();
  m_idsOfUniquePositions.clear();

#ifdef UQ_USES_TK_CLASS
#else
#ifdef UQ_MAC_SG_REQUIRES_INVERTED_COV_MATRICES
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
#endif

  return;
}

#ifdef UQ_USES_TK_CLASS
#else
template<class P_V,class P_M>
int
uqMarkovChainSGClass<P_V,P_M>::computeInitialCholFactors()
{
//const P_M& proposalCovMatrix (*m_initialProposalCovMatrix);
//const P_M& proposalPrecMatrix(...);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqMarkovChainSGClass<P_V,P_M>::computeInitialCholFactors()..."
              << std::endl;
  }

  int iRC = UQ_OK_RC;

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "In uqMarkovChainSGClass<P_V,P_M>::computeInitialCholFactors()"
                                   << ": using supplied initialProposalCovMatrix, whose contents are:"
                                   << std::endl;
    *m_env.subDisplayFile() << *m_initialProposalCovMatrix; // FIX ME: output might need to be in parallel
    *m_env.subDisplayFile() << std::endl;
  }

#ifdef UQ_USES_TK_CLASS
#else
  m_lowerCholProposalCovMatrices[0] = new P_M(*m_initialProposalCovMatrix); 
  iRC = m_lowerCholProposalCovMatrices[0]->chol();
  UQ_FATAL_RC_MACRO(iRC,
                    m_env.fullRank(),
                    "uqMarkovChainSGClass<P_V,P_M>::computeInitialCholFactors()",
                    "initialProposalCovMatrix is not positive definite");
  m_lowerCholProposalCovMatrices[0]->zeroUpper(false);

  m_proposalCovMatrices[0] = new P_M(*m_initialProposalCovMatrix);

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "In uqMarkovChainSGClass<P_V,P_M>::computeInitialCholFactors()"
                                   << ", m_lowerCholProposalCovMatrices[0] contents are:"
                                   << std::endl;
    *m_env.subDisplayFile() << *(m_lowerCholProposalCovMatrices[0]); // FIX ME: output might need to be in parallel
    *m_env.subDisplayFile() << std::endl;
  }

#ifdef UQ_MAC_SG_REQUIRES_INVERTED_COV_MATRICES
  const P_M* internalProposalPrecMatrix = proposalPrecMatrix;
  if (proposalPrecMatrix == NULL) {
    UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                      m_env.fullRank(),
                      "uqMarkovChainSGClass<P_V,P_M>::computeInitialCholFactors()",
                      "not yet implemented for the case 'proposalPrecMatrix == NULL'");
  }

  m_upperCholProposalPrecMatrices[0] = new P_M(proposalPrecMatrix); 
  iRC = m_upperCholProposalPrecMatrices[0]->chol();
  UQ_FATAL_RC_MACRO(iRC,
                    m_env.fullRank(),
                    "uqMarkovChainSGClass<P_V,P_M>::computeInitialCholFactors()",
                    "proposalPrecMatrix is not positive definite");
  m_upperCholProposalPrecMatrices[0]->zeroLower(false);

  m_proposalPrecMatrices[0] = new P_M(proposalPrecMatrix);

  //if (m_env.subDisplayFile()) {
  //  *m_env.subDisplayFile() << "In uqMarkovChainSGClass<P_V,P_M>::computeInitialCholFactors()"
  //                                 << ", m_upperCholProposalPrecMatrices[0] contents are:"
  //                                 << std::endl;
  //  *m_env.subDisplayFile() << *(m_upperCholProposalPrecMatrices[0]); // FIX ME: output might need to be in parallel
  //  *m_env.subDisplayFile() << std::endl;
  //}

#endif
#endif

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqMarkovChainSGClass<P_V,P_M>::computeInitialCholFactors()"
                           << std::endl;
  }

  return iRC;
}

template<class P_V,class P_M>
void
uqMarkovChainSGClass<P_V,P_M>::updateTK()
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqMarkovChainSGClass<P_V,P_M>::updateTK()"
                           << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "In uqMarkovChainSGClass<P_V,P_M>::updateTK()"
                            << ": m_options.m_drMaxNumExtraStages = "           << m_options.m_drMaxNumExtraStages
                            << ", m_options.m_drScalesForCovMatrices.size() = " << m_options.m_drScalesForCovMatrices.size()
#ifdef UQ_USES_TK_CLASS
#else
                            << ", m_lowerCholProposalCovMatrices.size() = " << m_lowerCholProposalCovMatrices.size()
#endif
                            << std::endl;
  }

#ifdef UQ_USES_TK_CLASS
#else
  for (unsigned int i = 1; i < (m_options.m_drMaxNumExtraStages+1); ++i) {
    double scale = m_options.m_drScalesForCovMatrices[i];
    if (m_lowerCholProposalCovMatrices[i]) delete m_lowerCholProposalCovMatrices[i];
    m_lowerCholProposalCovMatrices [i]   = new P_M(*(m_lowerCholProposalCovMatrices[i-1]));
  *(m_lowerCholProposalCovMatrices [i]) /= scale;
    if (m_proposalCovMatrices[i]) delete m_proposalCovMatrices[i];
    m_proposalCovMatrices[i]             = new P_M(*(m_proposalCovMatrices[i-1]));
  *(m_proposalCovMatrices[i])           /= (scale*scale);
#ifdef UQ_MAC_SG_REQUIRES_INVERTED_COV_MATRICES
    m_upperCholProposalPrecMatrices[i]   = new P_M(*(m_upperCholProposalPrecMatrices[i-1]));
  *(m_upperCholProposalPrecMatrices[i]) *= scale;
    m_proposalPrecMatrices[i]            = new P_M(*(m_proposalPrecMatrices[i-1]));
  *(m_proposalPrecMatrices[i])          *= (scale*scale);
#endif
  }
#endif

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqMarkovChainSGClass<P_V,P_M>::updateTK()"
                           << std::endl;
  }

  return;
}

template<class P_V,class P_M>
double
uqMarkovChainSGClass<P_V,P_M>::logProposal(
  const uqMarkovChainPositionDataClass<P_V>& x,
  const uqMarkovChainPositionDataClass<P_V>& y,
  unsigned int                               idOfProposalCovMatrix)
{
#ifdef UQ_USES_TK_CLASS
  return 0.;
#else
  P_V diffVec(y.vecValues() - x.vecValues());
#ifdef UQ_MAC_SG_REQUIRES_INVERTED_COV_MATRICES
  double value = -0.5 * scalarProduct(diffVec, *(m_proposalPrecMatrices[idOfProposalCovMatrix]) * diffVec);
#else
  double value = -0.5 * scalarProduct(diffVec, m_proposalCovMatrices[idOfProposalCovMatrix]->invertMultiply(diffVec));
#endif
  return value;
#endif
}

template<class P_V,class P_M>
double
uqMarkovChainSGClass<P_V,P_M>::logProposal(const std::vector<uqMarkovChainPositionDataClass<P_V>*>& inputPositions)
{
  unsigned int inputSize = inputPositions.size();
  UQ_FATAL_TEST_MACRO((inputSize < 2),
                      m_env.fullRank(),
                      "uqMarkovChainSGClass<P_V,P_M>::logProposal()",
                      "inputPositions has size < 2");

  return this->logProposal(*(inputPositions[0            ]),
                           *(inputPositions[inputSize - 1]),
                           inputSize-2);
}
#endif // UQ_USES_TK_CLASS

template<class P_V,class P_M>
double
uqMarkovChainSGClass<P_V,P_M>::alpha(
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
        (isnan(x.logTarget())      )) {
      std::cerr << "WARNING In uqMarkovChainSGClass<P_V,P_M>::alpha(x,y)"
                << ", fullRank "        << m_env.fullRank()
                << ", subEnvironment "  << m_env.subId()
                << ", subRank "         << m_env.subRank()
                << ", inter0Rank "      << m_env.inter0Rank()
                << ": x.logTarget() = " << x.logTarget()
                << ", x.values() = "    << x.vecValues()
                << ", y.values() = "    << y.vecValues()
                << std::endl;
    }
    else if ((y.logTarget() == -INFINITY) ||
             (y.logTarget() ==  INFINITY) ||
             (isnan(y.logTarget())      )) {
      std::cerr << "WARNING In uqMarkovChainSGClass<P_V,P_M>::alpha(x,y)"
                << ", fullRank "        << m_env.fullRank()
                << ", subEnvironment "  << m_env.subId()
                << ", subRank "         << m_env.subRank()
                << ", inter0Rank "      << m_env.inter0Rank()
                << ": y.logTarget() = " << y.logTarget()
                << ", x.values() = "    << x.vecValues()
                << ", y.values() = "    << y.vecValues()
                << std::endl;
    }
    else {
      double yLogTargetToUse = y.logTarget();
#ifdef UQ_MAC_SG_REQUIRES_TARGET_DISTRIBUTION_ONLY
#else
      if (m_likelihoodObjComputesMisfits &&
          m_observableSpace.shouldVariancesBeUpdated()) {
        // Divide the misfitVector of 'y' by the misfitVarianceVector of 'x'
        yLogTargetToUse = -0.5 * ( y.m2lPrior() + (y.misfitVector()/x.misfitVarianceVector()).sumOfComponents() );
      }
#endif
#ifdef UQ_USES_TK_CLASS
      if (m_tk->symmetric()) {
#else
      if (m_tkIsSymmetric) {
#endif
        alphaQuotient = std::exp(yLogTargetToUse - x.logTarget());
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 10)) {
          *m_env.subDisplayFile() << "In uqMarkovChainSGClass<P_V,P_M>::alpha(x,y)"
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
#ifdef UQ_USES_TK_CLASS // AQUI
        double qyx = -.5 * m_tk->rv(yStageId).pdf().minus2LnValue(x.vecValues(),NULL,NULL,NULL,NULL);
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 10)) {
          const uqGaussianJointPdfClass<P_V,P_M>* pdfYX = dynamic_cast< const uqGaussianJointPdfClass<P_V,P_M>* >(&(m_tk->rv(yStageId).pdf()));
          *m_env.subDisplayFile() << "In uqMarkovChainSGClass<P_V,P_M>::alpha(x,y)"
                                 << ", rvYX.domainExpVector = " << pdfYX->domainExpVector()
                                 << ", rvYX.domainVarVector = " << pdfYX->domainVarVector()
                                 << ", rvYX.covMatrix = "       << pdfYX->covMatrix()
                                 << std::endl;
        }
        double qxy = -.5 * m_tk->rv(xStageId).pdf().minus2LnValue(y.vecValues(),NULL,NULL,NULL,NULL);
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 10)) {
          const uqGaussianJointPdfClass<P_V,P_M>* pdfXY = dynamic_cast< const uqGaussianJointPdfClass<P_V,P_M>* >(&(m_tk->rv(xStageId).pdf()));
          *m_env.subDisplayFile() << "In uqMarkovChainSGClass<P_V,P_M>::alpha(x,y)"
                                 << ", rvXY.domainExpVector = " << pdfXY->domainExpVector()
                                 << ", rvXY.domainVarVector = " << pdfXY->domainVarVector()
                                 << ", rvXY.covMatrix = "       << pdfXY->covMatrix()
                                 << std::endl;
        }
#else
        double qyx = logProposal(y,x,0);
        double qxy = logProposal(x,y,0);
#endif
        alphaQuotient = std::exp(yLogTargetToUse +
                                 qyx -
                                 x.logTarget() -
                                 qxy);
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 10)) {
          *m_env.subDisplayFile() << "In uqMarkovChainSGClass<P_V,P_M>::alpha(x,y)"
                                 << ": unsymmetric proposal case"
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
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 10)) {
      *m_env.subDisplayFile() << "In uqMarkovChainSGClass<P_V,P_M>::alpha(x,y)"
                             << ": x.outOfTargetSupport = " << x.outOfTargetSupport()
                             << ", y.outOfTargetSupport = " << y.outOfTargetSupport()
                             << std::endl;
    }
  }
  if (alphaQuotientPtr != NULL) *alphaQuotientPtr = alphaQuotient;

  return std::min(1.,alphaQuotient);
}

template<class P_V,class P_M>
double
uqMarkovChainSGClass<P_V,P_M>::alpha(
  const std::vector<uqMarkovChainPositionDataClass<P_V>*>& inputPositionsData,
  const std::vector<unsigned int                        >& inputTKStageIds)
{
  unsigned int inputSize = inputPositionsData.size();
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 10)) {
    *m_env.subDisplayFile() << "Entering uqMarkovChainSGClass<P_V,P_M>::alpha(vec)"
                           << ", inputSize = " << inputSize
                           << std::endl;
  }
  UQ_FATAL_TEST_MACRO((inputSize < 2),
                      m_env.fullRank(),
                      "uqMarkovChainSGClass<P_V,P_M>::alpha(vec)",
                      "inputPositionsData has size < 2");

  // If necessary, return 0. right away
  if (inputPositionsData[0          ]->outOfTargetSupport()) return 0.;
  if (inputPositionsData[inputSize-1]->outOfTargetSupport()) return 0.;

  if ((inputPositionsData[0]->logTarget() == -INFINITY) ||
      (inputPositionsData[0]->logTarget() ==  INFINITY) ||
      (isnan(inputPositionsData[0]->logTarget())      )) {
    std::cerr << "WARNING In uqMarkovChainSGClass<P_V,P_M>::alpha(vec)"
              << ", fullRank "       << m_env.fullRank()
              << ", subEnvironment " << m_env.subId()
              << ", subRank "        << m_env.subRank()
              << ", inter0Rank "     << m_env.inter0Rank()
              << ": inputSize = "    << inputSize
              << ", inputPositionsData[0]->logTarget() = " << inputPositionsData[0]->logTarget()
              << ", [0]->values() = "                      << inputPositionsData[0]->vecValues()
              << ", [inputSize - 1]->values() = "          << inputPositionsData[inputSize-1]->vecValues()
              << std::endl;
    return 0.;
  }
  else if ((inputPositionsData[inputSize - 1]->logTarget() == -INFINITY) ||
           (inputPositionsData[inputSize - 1]->logTarget() ==  INFINITY) ||
           (isnan(inputPositionsData[inputSize - 1]->logTarget())      )) {
    std::cerr << "WARNING In uqMarkovChainSGClass<P_V,P_M>::alpha(vec)"
              << ", fullRank "       << m_env.fullRank()
              << ", subEnvironment " << m_env.subId()
              << ", subRank "        << m_env.subRank()
              << ", inter0Rank "     << m_env.inter0Rank()
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
#ifdef UQ_USES_TK_CLASS // AQUI
  const P_V& _lastTKPosition         = m_tk->preComputingPosition(        tkStageIds[inputSize-1]);
  const P_V& _lastBackwardTKPosition = m_tk->preComputingPosition(backwardTKStageIds[inputSize-1]);

  double numContrib = -.5 * m_tk->rv(backwardTKStageIdsLess1).pdf().minus2LnValue(_lastBackwardTKPosition,NULL,NULL,NULL,NULL);
  double denContrib = -.5 * m_tk->rv(        tkStageIdsLess1).pdf().minus2LnValue(_lastTKPosition        ,NULL,NULL,NULL,NULL);
#else
  double numContrib = logProposal(backwardPositionsData);
  double denContrib = logProposal(        positionsData);
#endif
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 10)) {
    *m_env.subDisplayFile() << "In uqMarkovChainSGClass<P_V,P_M>::alpha(vec)"
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

#ifdef UQ_USES_TK_CLASS // AQUI
    const P_V& lastTKPosition         = m_tk->preComputingPosition(        tkStageIds[inputSize-2-i]);
    const P_V& lastBackwardTKPosition = m_tk->preComputingPosition(backwardTKStageIds[inputSize-2-i]);

            tkStageIds.pop_back();
    backwardTKStageIds.pop_back();

            tkStageIdsLess1.pop_back();
    backwardTKStageIdsLess1.pop_back();

    numContrib = -.5 * m_tk->rv(backwardTKStageIdsLess1).pdf().minus2LnValue(lastBackwardTKPosition,NULL,NULL,NULL,NULL);
    denContrib = -.5 * m_tk->rv(        tkStageIdsLess1).pdf().minus2LnValue(lastTKPosition        ,NULL,NULL,NULL,NULL);
#else
    numContrib = logProposal(backwardPositionsData);
    denContrib = logProposal(        positionsData);
#endif
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 10)) {
      *m_env.subDisplayFile() << "In uqMarkovChainSGClass<P_V,P_M>::alpha(vec)"
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
#ifdef UQ_MAC_SG_REQUIRES_TARGET_DISTRIBUTION_ONLY
#else
  if (m_likelihoodObjComputesMisfits &&
      m_observableSpace.shouldVariancesBeUpdated()) {
    // Divide the misfitVector of 'back[0]' by the misfitVarianceVector of 'pos[0]'
    numeratorLogTargetToUse = -0.5 * ( backwardPositionsData[0]->m2lPrior() +
      (backwardPositionsData[0]->misfitVector()/positionsData[0]->misfitVarianceVector()).sumOfComponents() );
  }
#endif
  numContrib = numeratorLogTargetToUse;
  denContrib = positionsData[0]->logTarget();
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 10)) {
    *m_env.subDisplayFile() << "In uqMarkovChainSGClass<P_V,P_M>::alpha(vec)"
                           << ", inputSize = "  << inputSize
                           << ", after loop"
                           << ": numContrib = " << numContrib
                           << ", denContrib = " << denContrib
                           << std::endl;
  }
  logNumerator   += numContrib;
  logDenominator += denContrib;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 10)) {
    *m_env.subDisplayFile() << "Leaving uqMarkovChainSGClass<P_V,P_M>::alpha(vec)"
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

template<class P_V,class P_M>
bool
uqMarkovChainSGClass<P_V,P_M>::acceptAlpha(double alpha)
{
  bool result = false;

  if      (alpha <= 0.                          ) result = false;
  else if (alpha >= 1.                          ) result = true;
  else if (alpha >= gsl_rng_uniform(m_env.rng())) result = true;
  else                                            result = false;

  return result;
}

template<class P_V,class P_M>
int
uqMarkovChainSGClass<P_V,P_M>::writeInfo(
  const uqBaseVectorSequenceClass<P_V,P_M>& workingChain,
  std::ofstream&                            ofsvar) const
//const P_M*                   mahalanobisMatrix,
//bool                         applyMahalanobisInvert) const
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\n"
                            << "\n-----------------------------------------------------"
                            << "\n Writing extra information about the Markov chain " << workingChain.name() << " to output file ..."
                            << "\n-----------------------------------------------------"
                            << "\n"
                            << std::endl;
  }

  int iRC = UQ_OK_RC;

  if (m_options.m_rawChainGenerateExtra) {
    // Write m_logTargets
    ofsvar << m_options.m_prefix << "logTargets_sub" << m_env.subIdString() << " = zeros(" << m_logTargets.size()
           << ","                                                           << 1
           << ");"
           << std::endl;
    ofsvar << m_options.m_prefix << "logTargets_sub" << m_env.subIdString() << " = [";
    for (unsigned int i = 0; i < m_logTargets.size(); ++i) {
      ofsvar << m_logTargets[i]
             << std::endl;
    }
    ofsvar << "];\n";

    // Write m_alphaQuotients
    ofsvar << m_options.m_prefix << "alphaQuotients_sub" << m_env.subIdString() << " = zeros(" << m_alphaQuotients.size()
           << ","                                                               << 1
           << ");"
           << std::endl;
    ofsvar << m_options.m_prefix << "alphaQuotients_sub" << m_env.subIdString() << " = [";
    for (unsigned int i = 0; i < m_alphaQuotients.size(); ++i) {
      ofsvar << m_alphaQuotients[i]
             << std::endl;
    }
    ofsvar << "];\n";
  }

  // Write names of components
  ofsvar << m_options.m_prefix << "componentNames = {";
  m_vectorSpace.printComponentsNames(ofsvar,false);
  ofsvar << "};\n";

#if 0
  // Write mahalanobis distances
  if (mahalanobisMatrix != NULL) {
    P_V diffVec(m_vectorSpace.zeroVector());
    ofsvar << m_options.m_prefix << "d = [";
    if (applyMahalanobisInvert) {
      P_V tmpVec(m_vectorSpace.zeroVector());
      P_V vec0(m_vectorSpace.zeroVector());
      workingChain.getPositionValues(0,vec0);
      for (unsigned int i = 0; i < workingChain.subSequenceSize(); ++i) {
        workingChain.getPositionValues(i,tmpVec);
        diffVec = tmpVec - vec0;
        //diffVec = *(workingChain[i]) - *(workingChain[0]);
        ofsvar << scalarProduct(diffVec, mahalanobisMatrix->invertMultiply(diffVec))
               << std::endl;
      }
    }
    else {
      P_V tmpVec(m_vectorSpace.zeroVector());
      P_V vec0(m_vectorSpace.zeroVector());
      workingChain.getPositionValues(0,vec0);
      for (unsigned int i = 0; i < workingChain.subSequenceSize(); ++i) {
        workingChain.getPositionValues(i,tmpVec);
        diffVec = tmpVec - vec0;
        //diffVec = *(workingChain[i]) - *(workingChain[0]);
        ofsvar << scalarProduct(diffVec, *mahalanobisMatrix * diffVec)
               << std::endl;
      }
    }
    ofsvar << "];\n";
  }
#endif

#if 0
  // Write prior mean values
  ofsvar << m_options.m_prefix << "priorMeanValues = ["
         << m_vectorSpace.priorMuValues()
         << "];\n";

  // Write prior sigma values
  ofsvar << m_options.m_prefix << "priorSigmaValues = ["
         << m_vectorSpace.priorSigmaValues()
         << "];\n";

#if 0
  ofsvar << m_options.m_prefix << "results.prior = [queso_priorMeanValues',queso_priorSigmaValues'];\n";
#endif

  // Write vector space lower bounds
  ofsvar << m_options.m_prefix << "minValues = ["
         << m_vectorSpace.minValues()
         << "];\n";

  // Write vector space upper bounds
  ofsvar << m_options.m_prefix << "maxValues = ["
         << m_vectorSpace.maxValues()
         << "];\n";
#endif

#if 0
  ofsvar << m_options.m_prefix << "results.limits = [queso_low',queso_upp'];\n";

  // Write out data for mcmcpred.m of MATLAB MCMC toolbox
  ofsvar << m_options.m_prefix << "results.parind = ["; // FIXME
  for (unsigned int i = 0; i < m_vectorSpace.dim(); ++i) {
    ofsvar << i+1
           << std::endl;
  }
  ofsvar << "];\n";

  ofsvar << m_options.m_prefix << "results.local = [\n"; // FIXME
  for (unsigned int i = 0; i < m_vectorSpace.dim(); ++i) {
    ofsvar << " 0";
    //<< std::endl;
  }
  ofsvar << "];\n";

  if (m_options.m_rawChainUse2) {
  }
  else {
    bool savedVectorPrintState = workingChain[workingChain.subSequenceSize()-1]->getPrintHorizontally();
    workingChain[workingChain.subSequenceSize()-1]->setPrintHorizontally(false);
    ofsvar << m_options.m_prefix << "results.theta = ["
           << *(workingChain[workingChain.subSequenceSize()-1])
           << "];\n";
    workingChain[workingChain.subSequenceSize()-1]->setPrintHorizontally(savedVectorPrintState);
  }
  
  ofsvar << m_options.m_prefix << "results.nbatch = 1.;\n"; // FIXME

  if (mahalanobisMatrix != NULL) {
    // Write covMatrix
    ofsvar << m_options.m_prefix << "mahalanobisMatrix = ["
           << *mahalanobisMatrix
           << "];\n";
  }
#endif

  // Write number of rejections
  ofsvar << m_options.m_prefix << "rejected = " << (double) m_numRejections/(double) (workingChain.subSequenceSize()-1)
         << ";\n"
         << std::endl;

  // Write number of out of target support
  ofsvar << m_options.m_prefix << "outTargetSupport = " << (double) m_numOutOfTargetSupport/(double) (workingChain.subSequenceSize()-1)
         << ";\n"
         << std::endl;

  // Write chain run time
  ofsvar << m_options.m_prefix << "runTime = " << m_rawChainRunTime
         << ";\n"
         << std::endl;

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "\n-----------------------------------------------------"
                            << "\n Finished writing extra information about the Markov chain " << workingChain.name()
                            << "\n-----------------------------------------------------"
                            << "\n"
                            << std::endl;
  }

  return iRC;
}

template<class P_V,class P_M>
std::ostream& operator<<(std::ostream& os, const uqMarkovChainSGClass<P_V,P_M>& obj)
{
  obj.print(os);

  return os;
}
#endif // __UQ_MAC_SG1_H__
