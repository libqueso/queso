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

#include <queso/asserts.h>
#include <queso/MetropolisHastingsSG.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

#include <queso/HessianCovMatricesTKGroup.h>
#include <queso/ScaledCovMatrixTKGroup.h>

#include <queso/GaussianJointPdf.h>

namespace QUESO {

// Default constructor -----------------------------
MHRawChainInfoStruct::MHRawChainInfoStruct()
{
  reset();
}
// Copy constructor----------------------------------
MHRawChainInfoStruct::MHRawChainInfoStruct(const MHRawChainInfoStruct& rhs)
{
  this->copy(rhs);
}
// Destructor ---------------------------------------
MHRawChainInfoStruct::~MHRawChainInfoStruct()
{
}
// Set methods---------------------------------------
MHRawChainInfoStruct&
MHRawChainInfoStruct::operator=(const MHRawChainInfoStruct& rhs)
{
  this->copy(rhs);
  return *this;
}
//---------------------------------------------------
MHRawChainInfoStruct&
MHRawChainInfoStruct::operator+=(const MHRawChainInfoStruct& rhs)
{
  runTime          += rhs.runTime;
  candidateRunTime += rhs.candidateRunTime;
  targetRunTime    += rhs.targetRunTime;
  mhAlphaRunTime   += rhs.mhAlphaRunTime;
  drAlphaRunTime   += rhs.drAlphaRunTime;
  drRunTime        += rhs.drRunTime;
  amRunTime        += rhs.amRunTime;

  numTargetCalls            += rhs.numTargetCalls;
  numDRs                    += rhs.numDRs;
  numOutOfTargetSupport     += rhs.numOutOfTargetSupport;
  numOutOfTargetSupportInDR += rhs.numOutOfTargetSupportInDR;
  numRejections             += rhs.numRejections;

  return *this;
}
// Misc methods--------------------------------------------------
void
MHRawChainInfoStruct::reset()
{
  runTime          = 0.;
  candidateRunTime = 0.;
  targetRunTime    = 0.;
  mhAlphaRunTime   = 0.;
  drAlphaRunTime   = 0.;
  drRunTime        = 0.;
  amRunTime        = 0.;

  numTargetCalls            = 0;
  numDRs                    = 0;
  numOutOfTargetSupport     = 0;
  numOutOfTargetSupportInDR = 0;
  numRejections             = 0;
}
//---------------------------------------------------
void
MHRawChainInfoStruct::copy(const MHRawChainInfoStruct& rhs)
{
  runTime          = rhs.runTime;
  candidateRunTime = rhs.candidateRunTime;
  targetRunTime    = rhs.targetRunTime;
  mhAlphaRunTime   = rhs.mhAlphaRunTime;
  drAlphaRunTime   = rhs.drAlphaRunTime;
  drRunTime        = rhs.drRunTime;
  amRunTime        = rhs.amRunTime;

  numTargetCalls            = rhs.numTargetCalls;
  numDRs                    = rhs.numDRs;
  numOutOfTargetSupport     = rhs.numOutOfTargetSupport;
  numOutOfTargetSupportInDR = rhs.numOutOfTargetSupportInDR;
  numRejections             = rhs.numRejections;

  return;
}
//---------------------------------------------------
void
MHRawChainInfoStruct::mpiSum(const MpiComm& comm, MHRawChainInfoStruct& sumInfo) const
{
  comm.Allreduce((void *) &runTime, (void *) &sumInfo.runTime, (int) 7, RawValue_MPI_DOUBLE, RawValue_MPI_SUM,
                 "MHRawChainInfoStruct::mpiSum()",
                 "failed MPI.Allreduce() for sum of doubles");

  comm.Allreduce((void *) &numTargetCalls, (void *) &sumInfo.numTargetCalls, (int) 5, RawValue_MPI_UNSIGNED, RawValue_MPI_SUM,
                 "MHRawChainInfoStruct::mpiSum()",
                 "failed MPI.Allreduce() for sum of unsigned ints");

  return;
}

// Default constructor -----------------------------
template<class P_V,class P_M>
MetropolisHastingsSG<P_V,P_M>::MetropolisHastingsSG(
  /*! Prefix                     */ const char*                         prefix,
  /*! Options (if no input file) */ const MhOptionsValues*       alternativeOptionsValues, // dakota
  /*! The source RV              */ const BaseVectorRV<P_V,P_M>& sourceRv,
  /*! Initial chain position     */ const P_V&                          initialPosition,
  /*! Proposal cov. matrix       */ const P_M*                          inputProposalCovMatrix)
  :
  m_env                       (sourceRv.env()),
  m_vectorSpace               (sourceRv.imageSet().vectorSpace()),
  m_targetPdf                 (sourceRv.pdf()),
  m_initialPosition           (initialPosition),
  m_initialProposalCovMatrix  (m_vectorSpace.zeroVector()),
  m_nullInputProposalCovMatrix(inputProposalCovMatrix == NULL),
  m_numDisabledParameters     (0), // gpmsa2
  m_parameterEnabledStatus    (m_vectorSpace.dimLocal(),true), // gpmsa2
  m_targetPdfSynchronizer     (new ScalarFunctionSynchronizer<P_V,P_M>(m_targetPdf,m_initialPosition)),
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
  m_optionsObj                (NULL),
  m_computeInitialPriorAndLikelihoodValues(true),
  m_initialLogPriorValue      (0.),
  m_initialLogLikelihoodValue (0.)
{
  if (inputProposalCovMatrix != NULL) {
    m_initialProposalCovMatrix = *inputProposalCovMatrix;
  }
  if (alternativeOptionsValues) m_alternativeOptionsValues = *alternativeOptionsValues;
  if (m_env.optionsInputFileName() == "") {
    m_optionsObj = new MetropolisHastingsSGOptions(m_env,prefix,m_alternativeOptionsValues);
  }
  else {
    m_optionsObj = new MetropolisHastingsSGOptions(m_env,prefix);
    m_optionsObj->scanOptionsValues();
  }

  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_ov.m_totallyMute == false)) {
    *m_env.subDisplayFile() << "Entering MetropolisHastingsSG<P_V,P_M>::constructor(1)"
                            << ": prefix = " << prefix
                            << ", alternativeOptionsValues = " << alternativeOptionsValues
                            << ", m_env.optionsInputFileName() = " << m_env.optionsInputFileName()
                            << ", m_initialProposalCovMatrix = " << m_initialProposalCovMatrix
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO(sourceRv.imageSet().vectorSpace().dimLocal() != initialPosition.sizeLocal(),
                      m_env.worldRank(),
                      "MetropolisHastingsSG<P_V,P_M>::constructor(1)",
                      "'sourceRv' and 'initialPosition' should have equal dimensions");

  if (inputProposalCovMatrix) {
    UQ_FATAL_TEST_MACRO(sourceRv.imageSet().vectorSpace().dimLocal() != inputProposalCovMatrix->numRowsLocal(),
                        m_env.worldRank(),
                        "MetropolisHastingsSG<P_V,P_M>::constructor(1)",
                        "'sourceRv' and 'inputProposalCovMatrix' should have equal dimensions");
    UQ_FATAL_TEST_MACRO(inputProposalCovMatrix->numCols() != inputProposalCovMatrix->numRowsGlobal(),
                        m_env.worldRank(),
                        "MetropolisHastingsSG<P_V,P_M>::constructor(1)",
                        "'inputProposalCovMatrix' should be a square matrix");
  }

  commonConstructor();

  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_ov.m_totallyMute == false)) {
    *m_env.subDisplayFile() << "Leaving MetropolisHastingsSG<P_V,P_M>::constructor(1)"
                            << std::endl;
  }
}
template<class P_V,class P_M>
MetropolisHastingsSG<P_V,P_M>::MetropolisHastingsSG(
  /*! Prefix                     */ const char*                         prefix,
  /*! Options (if no input file) */ const MhOptionsValues*       alternativeOptionsValues, // dakota
  /*! The source RV              */ const BaseVectorRV<P_V,P_M>& sourceRv,
  /*! Initial chain position     */ const P_V&                          initialPosition,
  /*! Initial chain prior        */ double                              initialLogPrior,
  /*! Initial chain likelihood   */ double                              initialLogLikelihood,
  /*! Proposal cov. matrix       */ const P_M*                          inputProposalCovMatrix)
  :
  m_env                       (sourceRv.env()),
  m_vectorSpace               (sourceRv.imageSet().vectorSpace()),
  m_targetPdf                 (sourceRv.pdf()),
  m_initialPosition           (initialPosition),
  m_initialProposalCovMatrix  (m_vectorSpace.zeroVector()),
  m_nullInputProposalCovMatrix(inputProposalCovMatrix == NULL),
  m_numDisabledParameters     (0), // gpmsa2
  m_parameterEnabledStatus    (m_vectorSpace.dimLocal(),true), // gpmsa2
  m_targetPdfSynchronizer     (new ScalarFunctionSynchronizer<P_V,P_M>(m_targetPdf,m_initialPosition)),
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
  m_optionsObj                (NULL),
  m_computeInitialPriorAndLikelihoodValues(false),
  m_initialLogPriorValue      (initialLogPrior),
  m_initialLogLikelihoodValue (initialLogLikelihood)
{
  if (inputProposalCovMatrix != NULL) {
    m_initialProposalCovMatrix = *inputProposalCovMatrix;
  }
  if (alternativeOptionsValues) m_alternativeOptionsValues = *alternativeOptionsValues;
  if (m_env.optionsInputFileName() == "") {
    m_optionsObj = new MetropolisHastingsSGOptions(m_env,prefix,m_alternativeOptionsValues);
  }
  else {
    m_optionsObj = new MetropolisHastingsSGOptions(m_env,prefix);
    m_optionsObj->scanOptionsValues();
  }

  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_ov.m_totallyMute == false)) {
    *m_env.subDisplayFile() << "Entering MetropolisHastingsSG<P_V,P_M>::constructor(2)"
                            << ": prefix = " << prefix
                            << ", alternativeOptionsValues = " << alternativeOptionsValues
                            << ", m_env.optionsInputFileName() = " << m_env.optionsInputFileName()
                            << ", m_initialProposalCovMatrix = " << m_initialProposalCovMatrix
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO(sourceRv.imageSet().vectorSpace().dimLocal() != initialPosition.sizeLocal(),
                      m_env.worldRank(),
                      "MetropolisHastingsSG<P_V,P_M>::constructor(2)",
                      "'sourceRv' and 'initialPosition' should have equal dimensions");

  if (inputProposalCovMatrix) {
    UQ_FATAL_TEST_MACRO(sourceRv.imageSet().vectorSpace().dimLocal() != inputProposalCovMatrix->numRowsLocal(),
                        m_env.worldRank(),
                        "MetropolisHastingsSG<P_V,P_M>::constructor(2)",
                        "'sourceRv' and 'inputProposalCovMatrix' should have equal dimensions");
    UQ_FATAL_TEST_MACRO(inputProposalCovMatrix->numCols() != inputProposalCovMatrix->numRowsGlobal(),
                        m_env.worldRank(),
                        "MetropolisHastingsSG<P_V,P_M>::constructor(2)",
                        "'inputProposalCovMatrix' should be a square matrix");
  }

  commonConstructor();

  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_ov.m_totallyMute == false)) {
    *m_env.subDisplayFile() << "Leaving MetropolisHastingsSG<P_V,P_M>::constructor(2)"
                            << std::endl;
  }
}
// Constructor -------------------------------------
template<class P_V,class P_M>
MetropolisHastingsSG<P_V,P_M>::MetropolisHastingsSG(
  const MLSamplingLevelOptions& mlOptions,
  const BaseVectorRV<P_V,P_M>&  sourceRv,
  const P_V&                           initialPosition, // KEY
  const P_M*                           inputProposalCovMatrix)
  :
  m_env                       (sourceRv.env()),
  m_vectorSpace               (sourceRv.imageSet().vectorSpace()),
  m_targetPdf                 (sourceRv.pdf()),
  m_initialPosition           (initialPosition),
  m_initialProposalCovMatrix  (m_vectorSpace.zeroVector()),
  m_nullInputProposalCovMatrix(inputProposalCovMatrix == NULL),
  m_numDisabledParameters     (0), // gpmsa2
  m_parameterEnabledStatus    (m_vectorSpace.dimLocal(),true), // gpmsa2
  m_targetPdfSynchronizer     (new ScalarFunctionSynchronizer<P_V,P_M>(m_targetPdf,m_initialPosition)),
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
  m_optionsObj                (new MetropolisHastingsSGOptions(mlOptions)),
  m_computeInitialPriorAndLikelihoodValues(true),
  m_initialLogPriorValue      (0.),
  m_initialLogLikelihoodValue (0.)
{
  if (inputProposalCovMatrix != NULL) {
    m_initialProposalCovMatrix = *inputProposalCovMatrix;
    if ((m_env.subDisplayFile()                   ) &&
        (m_optionsObj->m_ov.m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::constructor(3)"
                              << ": just set m_initialProposalCovMatrix = " << m_initialProposalCovMatrix
                              << std::endl;
    }
  }
  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_ov.m_totallyMute == false)) {
    *m_env.subDisplayFile() << "Entering MetropolisHastingsSG<P_V,P_M>::constructor(3)"
                            << std::endl;
  }

  commonConstructor();

  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_ov.m_totallyMute == false)) {
    *m_env.subDisplayFile() << "Leaving MetropolisHastingsSG<P_V,P_M>::constructor(3)"
                            << std::endl;
  }
}
// Constructor -------------------------------------
template<class P_V,class P_M>
MetropolisHastingsSG<P_V,P_M>::MetropolisHastingsSG(
  const MLSamplingLevelOptions& mlOptions,
  const BaseVectorRV<P_V,P_M>&  sourceRv,
  const P_V&                           initialPosition, // KEY
  double                               initialLogPrior,
  double                               initialLogLikelihood,
  const P_M*                           inputProposalCovMatrix)
  :
  m_env                       (sourceRv.env()),
  m_vectorSpace               (sourceRv.imageSet().vectorSpace()),
  m_targetPdf                 (sourceRv.pdf()),
  m_initialPosition           (initialPosition),
  m_initialProposalCovMatrix  (m_vectorSpace.zeroVector()),
  m_nullInputProposalCovMatrix(inputProposalCovMatrix == NULL),
  m_numDisabledParameters     (0), // gpmsa2
  m_parameterEnabledStatus    (m_vectorSpace.dimLocal(),true), // gpmsa2
  m_targetPdfSynchronizer     (new ScalarFunctionSynchronizer<P_V,P_M>(m_targetPdf,m_initialPosition)),
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
  m_optionsObj                (new MetropolisHastingsSGOptions(mlOptions)),
  m_computeInitialPriorAndLikelihoodValues(false),
  m_initialLogPriorValue      (initialLogPrior),
  m_initialLogLikelihoodValue (initialLogLikelihood)
{
  if (inputProposalCovMatrix != NULL) {
    m_initialProposalCovMatrix = *inputProposalCovMatrix;
    if ((m_env.subDisplayFile()                   ) &&
        (m_optionsObj->m_ov.m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::constructor(4)"
                              << ": just set m_initialProposalCovMatrix = " << m_initialProposalCovMatrix
                              << std::endl;
    }
  }
  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_ov.m_totallyMute == false)) {
    *m_env.subDisplayFile() << "Entering MetropolisHastingsSG<P_V,P_M>::constructor(4)"
                            << std::endl;
  }

  commonConstructor();

  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_ov.m_totallyMute == false)) {
    *m_env.subDisplayFile() << "Leaving MetropolisHastingsSG<P_V,P_M>::constructor(4)"
                            << std::endl;
  }
}
// Destructor ---------------------------------------
template<class P_V,class P_M>
MetropolisHastingsSG<P_V,P_M>::~MetropolisHastingsSG()
{
  //if (m_env.subDisplayFile()) {
  //  *m_env.subDisplayFile() << "Entering MetropolisHastingsSG<P_V,P_M>::destructor()"
  //                          << std::endl;
  //}

  if (m_lastAdaptedCovMatrix) delete m_lastAdaptedCovMatrix;
  if (m_lastMean)             delete m_lastMean;
  m_lastChainSize             = 0;
  m_rawChainInfo.reset();
  m_alphaQuotients.clear();
  m_logTargets.clear();
  m_numDisabledParameters  = 0; // gpmsa2
  m_parameterEnabledStatus.clear(); // gpmsa2
  m_positionIdForDebugging = 0;
  m_stageIdForDebugging    = 0;
  m_idsOfUniquePositions.clear();

  if (m_tk                   ) delete m_tk;
  if (m_targetPdfSynchronizer) delete m_targetPdfSynchronizer;
  if (m_optionsObj           ) delete m_optionsObj;

  //if (m_env.subDisplayFile()) {
  //  *m_env.subDisplayFile() << "Leaving MetropolisHastingsSG<P_V,P_M>::destructor()"
  //                          << std::endl;
  //}
}

// Private methods----------------------------------
template<class P_V,class P_M>
void
MetropolisHastingsSG<P_V,P_M>::commonConstructor()
{
  /////////////////////////////////////////////////////////////////
  // Instantiate the appropriate TK (transition kernel)
  /////////////////////////////////////////////////////////////////
  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_ov.m_totallyMute == false)) {
    *m_env.subDisplayFile() << "Entering MetropolisHastingsSG<P_V,P_M>::commonConstructor()"
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
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::commonConstructor()"
                              << ": just read initial position contents = " << m_initialPosition
                              << std::endl;
    }
  }

  if (m_optionsObj->m_ov.m_parameterDisabledSet.size() > 0) { // gpmsa2
    for (std::set<unsigned int>::iterator setIt = m_optionsObj->m_ov.m_parameterDisabledSet.begin(); setIt != m_optionsObj->m_ov.m_parameterDisabledSet.end(); ++setIt) {
      unsigned int paramId = *setIt;
      if (paramId < m_vectorSpace.dimLocal()) {
        m_numDisabledParameters++;
        m_parameterEnabledStatus[paramId] = false;
      }
    }
  }

  std::vector<double> drScalesAll(m_optionsObj->m_ov.m_drScalesForExtraStages.size()+1,1.);
  for (unsigned int i = 1; i < (m_optionsObj->m_ov.m_drScalesForExtraStages.size()+1); ++i) {
    drScalesAll[i] = m_optionsObj->m_ov.m_drScalesForExtraStages[i-1];
  }
  if (m_optionsObj->m_ov.m_tkUseLocalHessian) { // sep2011
    m_tk = new HessianCovMatricesTKGroup<P_V,P_M>(m_optionsObj->m_prefix.c_str(),
                                                         m_vectorSpace,
                                                         drScalesAll,
                                                         *m_targetPdfSynchronizer);
    if ((m_env.subDisplayFile()                   ) &&
        (m_optionsObj->m_ov.m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::commonConstructor()"
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
        *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::commonConstructor()"
                                << ": just read initial proposal cov matrix contents = " << m_initialProposalCovMatrix
                                << std::endl;
      }
    }
    else {
      UQ_FATAL_TEST_MACRO(m_nullInputProposalCovMatrix,
                          m_env.worldRank(),
                          "MetropolisHastingsSG<P_V,P_M>::commonConstructor()",
                          "proposal cov matrix should have been passed by user, since, according to the input algorithm options, local Hessians will not be used in the proposal");
    }

    m_tk = new ScaledCovMatrixTKGroup<P_V,P_M>(m_optionsObj->m_prefix.c_str(),
                                                      m_vectorSpace,
                                                      drScalesAll,
                                                      m_initialProposalCovMatrix);
    if ((m_env.subDisplayFile()                   ) &&
        (m_optionsObj->m_ov.m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::commonConstructor()"
                              << ": just instantiated a 'ScaledCovMatrix' TK class"
                              << std::endl;
    }
  }

  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_ov.m_totallyMute == false)) {
    *m_env.subDisplayFile() << "Leaving MetropolisHastingsSG<P_V,P_M>::commonConstructor()"
                            << std::endl;
  }
  return;
}
//--------------------------------------------------
template<class P_V,class P_M>
double
MetropolisHastingsSG<P_V,P_M>::alpha(
  const MarkovChainPositionData<P_V>& x,
  const MarkovChainPositionData<P_V>& y,
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
      std::cerr << "WARNING In MetropolisHastingsSG<P_V,P_M>::alpha(x,y)"
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
      std::cerr << "WARNING In MetropolisHastingsSG<P_V,P_M>::alpha(x,y)"
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
          *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::alpha(x,y)"
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
          const GaussianJointPdf<P_V,P_M>* pdfYX = dynamic_cast< const GaussianJointPdf<P_V,P_M>* >(&(m_tk->rv(yStageId).pdf()));
          *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::alpha(x,y)"
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
          const GaussianJointPdf<P_V,P_M>* pdfXY = dynamic_cast< const GaussianJointPdf<P_V,P_M>* >(&(m_tk->rv(xStageId).pdf()));
          *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::alpha(x,y)"
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
          *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::alpha(x,y)"
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
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::alpha(x,y)"
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
MetropolisHastingsSG<P_V,P_M>::alpha(
  const std::vector<MarkovChainPositionData<P_V>*>& inputPositionsData,
  const std::vector<unsigned int                        >& inputTKStageIds)
{
  unsigned int inputSize = inputPositionsData.size();
  if ((m_env.subDisplayFile()                   ) &&
      (m_env.displayVerbosity() >= 10           ) &&
      (m_optionsObj->m_ov.m_totallyMute == false)) {
    *m_env.subDisplayFile() << "Entering MetropolisHastingsSG<P_V,P_M>::alpha(vec)"
                           << ", inputSize = " << inputSize
                           << std::endl;
  }
  UQ_FATAL_TEST_MACRO((inputSize < 2),
                      m_env.worldRank(),
                      "MetropolisHastingsSG<P_V,P_M>::alpha(vec)",
                      "inputPositionsData has size < 2");

  // If necessary, return 0. right away
  if (inputPositionsData[0          ]->outOfTargetSupport()) return 0.;
  if (inputPositionsData[inputSize-1]->outOfTargetSupport()) return 0.;

  if ((inputPositionsData[0]->logTarget() == -INFINITY           ) ||
      (inputPositionsData[0]->logTarget() ==  INFINITY           ) ||
      ( (boost::math::isnan)(inputPositionsData[0]->logTarget()) )) {
    std::cerr << "WARNING In MetropolisHastingsSG<P_V,P_M>::alpha(vec)"
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
    std::cerr << "WARNING In MetropolisHastingsSG<P_V,P_M>::alpha(vec)"
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
  std::vector<MarkovChainPositionData<P_V>*>         positionsData  (inputSize,NULL);
  std::vector<MarkovChainPositionData<P_V>*> backwardPositionsData  (inputSize,NULL);

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
    *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::alpha(vec)"
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
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::alpha(vec)"
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
    *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::alpha(vec)"
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
    *m_env.subDisplayFile() << "Leaving MetropolisHastingsSG<P_V,P_M>::alpha(vec)"
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
MetropolisHastingsSG<P_V,P_M>::acceptAlpha(double alpha)
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
MetropolisHastingsSG<P_V,P_M>::writeInfo(
  const BaseVectorSequence<P_V,P_M>& workingChain,
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

// Statistical methods -----------------------------
/* This operation currently implements the DRAM algorithm (Heikki Haario, Marko
 * Laine, Antonietta Mira and Eero Saksman, "DRAM: Efficient Adaptive MCMC",
 * Statistics and Computing (2006), 16:339-354). It also provides support for
 * Stochastic Newton algorithm through the TK (transition kernel) class. Stochastic
 * Newton is not totally implemented yet though, since it is being researched by
 * James Martin and Omar Ghattas at ICES at the University of Texas at Austin.*/
template <class P_V,class P_M>
void
MetropolisHastingsSG<P_V,P_M>::generateSequence(
  BaseVectorSequence<P_V,P_M>& workingChain,
  ScalarSequence<double>*      workingLogLikelihoodValues, // KEY: add LogPriorValues
  ScalarSequence<double>*      workingLogTargetValues)
{
  if ((m_env.subDisplayFile()                   ) &&
      (m_env.displayVerbosity() >= 5            ) &&
      (m_optionsObj->m_ov.m_totallyMute == false)) {
    *m_env.subDisplayFile() << "Entering MetropolisHastingsSG<P_V,P_M>::generateSequence()..."
                            << std::endl;
  }

  if (m_vectorSpace.dimLocal() != workingChain.vectorSizeLocal()) {
    std::cerr << "'m_vectorSpace' and 'workingChain' are related to vector"
              << "spaces of different dimensions"
              << std::endl;
    queso_error();
  }

  // Set a flag to write out log likelihood or not
  bool writeLogLikelihood;
  if ((workingLogLikelihoodValues != NULL) &&
      (m_optionsObj->m_ov.m_outputLogLikelihood)) {
    writeLogLikelihood = true;
  }
  else {
    writeLogLikelihood = false;
  }

  // Set a flag to write out log target or not
  bool writeLogTarget;
  if ((workingLogTargetValues != NULL) &&
      (m_optionsObj->m_ov.m_outputLogTarget)) {
    writeLogTarget = true;
  }
  else {
    writeLogTarget = false;
  }

  MiscCheckTheParallelEnvironment<P_V,P_V>(m_initialPosition,
                                             m_initialPosition);

  P_V valuesOf1stPosition(m_initialPosition);
  int iRC = UQ_OK_RC;

  workingChain.setName(m_optionsObj->m_prefix + "rawChain");

  //****************************************************
  // Generate chain
  //****************************************************
  if (m_optionsObj->m_ov.m_rawChainDataInputFileName == UQ_MH_SG_FILENAME_FOR_NO_FILE) {
    generateFullChain(valuesOf1stPosition,
                      m_optionsObj->m_ov.m_rawChainSize,
                      workingChain,
                      workingLogLikelihoodValues,
                      workingLogTargetValues);
  }
  else {
    readFullChain(m_optionsObj->m_ov.m_rawChainDataInputFileName,
                  m_optionsObj->m_ov.m_rawChainDataInputFileType,
                  m_optionsObj->m_ov.m_rawChainSize,
                  workingChain);
  }

  //****************************************************
  // Open generic output file
  //****************************************************
  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_ov.m_totallyMute == false)) {
    *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateSequence()"
                            << ", prefix = "                                         << m_optionsObj->m_prefix
                            << ", chain name = "                                     << workingChain.name()
                            << ": about to try to open generic output file '"        << m_optionsObj->m_ov.m_dataOutputFileName
                            << "."                                                   << UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT // Yes, always ".m"
                            << "', subId = "                                         << m_env.subId()
                            << ", subenv is allowed to write (1/true or 0/false) = " << (m_optionsObj->m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != m_optionsObj->m_ov.m_dataOutputAllowedSet.end())
                            << "..."
                            << std::endl;
  }

  FilePtrSetStruct genericFilePtrSet;
  m_env.openOutputFile(m_optionsObj->m_ov.m_dataOutputFileName,
                       UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT, // Yes, always ".m"
                       m_optionsObj->m_ov.m_dataOutputAllowedSet,
                       false,
                       genericFilePtrSet);

  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_ov.m_totallyMute == false)) {
    *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateSequence()"
                            << ", prefix = "                                   << m_optionsObj->m_prefix
                            << ", raw chain name = "                           << workingChain.name()
                            << ": returned from opening generic output file '" << m_optionsObj->m_ov.m_dataOutputFileName
                            << "."                                             << UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT // Yes, always ".m"
                            << "', subId = "                                   << m_env.subId()
                            << std::endl;
  }

  //****************************************************************************************
  // Eventually:
  // --> write raw chain
  // --> compute statistics on it
  //****************************************************************************************
  if ((m_optionsObj->m_ov.m_rawChainDataOutputFileName != UQ_MH_SG_FILENAME_FOR_NO_FILE) &&
      (m_optionsObj->m_ov.m_totallyMute == false                                       )) {

    // Take "sub" care of raw chain
    if ((m_env.subDisplayFile()                   ) &&
        (m_optionsObj->m_ov.m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateSequence()"
                              << ", prefix = "                                         << m_optionsObj->m_prefix
                              << ", raw chain name = "                                 << workingChain.name()
                              << ": about to try to write raw sub chain output file '" << m_optionsObj->m_ov.m_rawChainDataOutputFileName
                              << "."                                                   << m_optionsObj->m_ov.m_rawChainDataOutputFileType
                              << "', subId = "                                         << m_env.subId()
                              << ", subenv is allowed to write  1/true or 0/false) = " << (m_optionsObj->m_ov.m_rawChainDataOutputAllowedSet.find(m_env.subId()) != m_optionsObj->m_ov.m_rawChainDataOutputAllowedSet.end())
                              << "..."
                              << std::endl;
    }

    if ((m_numPositionsNotSubWritten                     >  0  ) &&
        (m_optionsObj->m_ov.m_rawChainDataOutputFileName != ".")) {
      workingChain.subWriteContents(m_optionsObj->m_ov.m_rawChainSize - m_numPositionsNotSubWritten,
                                    m_numPositionsNotSubWritten,
                                    m_optionsObj->m_ov.m_rawChainDataOutputFileName,
                                    m_optionsObj->m_ov.m_rawChainDataOutputFileType,
                                    m_optionsObj->m_ov.m_rawChainDataOutputAllowedSet);
      if ((m_env.subDisplayFile()                   ) &&
          (m_optionsObj->m_ov.m_totallyMute == false)) {
        *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateSequence()"
                                << ": just wrote (per period request) remaining " << m_numPositionsNotSubWritten << " chain positions "
                                << ", " << m_optionsObj->m_ov.m_rawChainSize - m_numPositionsNotSubWritten << " <= pos <= " << m_optionsObj->m_ov.m_rawChainSize - 1
                                << std::endl;
      }

      if (writeLogLikelihood) {
        workingLogLikelihoodValues->subWriteContents(m_optionsObj->m_ov.m_rawChainSize - m_numPositionsNotSubWritten,
                                                     m_numPositionsNotSubWritten,
                                                     m_optionsObj->m_ov.m_rawChainDataOutputFileName + "_loglikelihood",
                                                     m_optionsObj->m_ov.m_rawChainDataOutputFileType,
                                                     m_optionsObj->m_ov.m_rawChainDataOutputAllowedSet);
      }

      if (writeLogTarget) {
        workingLogTargetValues->subWriteContents(m_optionsObj->m_ov.m_rawChainSize - m_numPositionsNotSubWritten,
                                                 m_numPositionsNotSubWritten,
                                                 m_optionsObj->m_ov.m_rawChainDataOutputFileName + "_logtarget",
                                                 m_optionsObj->m_ov.m_rawChainDataOutputFileType,
                                                 m_optionsObj->m_ov.m_rawChainDataOutputAllowedSet);
      }

      m_numPositionsNotSubWritten = 0;
    }

    // Compute raw sub MLE
    if (workingLogLikelihoodValues) {
      SequenceOfVectors<P_V,P_M> rawSubMLEpositions(m_vectorSpace,0,m_optionsObj->m_prefix+"rawSubMLEseq");
      double rawSubMLEvalue = workingChain.subPositionsOfMaximum(*workingLogLikelihoodValues,
                                                                 rawSubMLEpositions);
      UQ_FATAL_TEST_MACRO(rawSubMLEpositions.subSequenceSize() == 0,
                          m_env.worldRank(),
                          "MetropolisHastingsSG<P_V,P_M>::generateSequence()",
                          "rawSubMLEpositions.subSequenceSize() = 0");

      if ((m_env.subDisplayFile()                   ) &&
          (m_optionsObj->m_ov.m_totallyMute == false)) {
        P_V tmpVec(m_vectorSpace.zeroVector());
        rawSubMLEpositions.getPositionValues(0,tmpVec);
        *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateSequence()"
                                << ": just computed MLE"
                                << ", rawSubMLEvalue = "                       << rawSubMLEvalue
                                << ", rawSubMLEpositions.subSequenceSize() = " << rawSubMLEpositions.subSequenceSize()
                                << ", rawSubMLEpositions[0] = "                << tmpVec
                                << std::endl;
      }
    }

    // Compute raw sub MAP
    if (workingLogTargetValues) {
      SequenceOfVectors<P_V,P_M> rawSubMAPpositions(m_vectorSpace,0,m_optionsObj->m_prefix+"rawSubMAPseq");
      double rawSubMAPvalue = workingChain.subPositionsOfMaximum(*workingLogTargetValues,
                                                                 rawSubMAPpositions);
      UQ_FATAL_TEST_MACRO(rawSubMAPpositions.subSequenceSize() == 0,
                          m_env.worldRank(),
                          "MetropolisHastingsSG<P_V,P_M>::generateSequence()",
                          "rawSubMAPpositions.subSequenceSize() = 0");

      if ((m_env.subDisplayFile()                   ) &&
          (m_optionsObj->m_ov.m_totallyMute == false)) {
        P_V tmpVec(m_vectorSpace.zeroVector());
        rawSubMAPpositions.getPositionValues(0,tmpVec);
        *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateSequence()"
                                << ": just computed MAP"
                                << ", rawSubMAPvalue = "                       << rawSubMAPvalue
                                << ", rawSubMAPpositions.subSequenceSize() = " << rawSubMAPpositions.subSequenceSize()
                                << ", rawSubMAPpositions[0] = "                << tmpVec
                                << std::endl;
      }
    }

    if ((m_env.subDisplayFile()                   ) &&
        (m_optionsObj->m_ov.m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateSequence()"
                              << ", prefix = "                                         << m_optionsObj->m_prefix
                              << ", raw chain name = "                                 << workingChain.name()
                              << ": returned from writing raw sub chain output file '" << m_optionsObj->m_ov.m_rawChainDataOutputFileName
                              << "."                                                   << m_optionsObj->m_ov.m_rawChainDataOutputFileType
                              << "', subId = "                                         << m_env.subId()
                              << std::endl;
    }

    // Take "unified" care of raw chain
    if ((m_env.subDisplayFile()                   ) &&
        (m_optionsObj->m_ov.m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateSequence()"
                              << ", prefix = "                                             << m_optionsObj->m_prefix
                              << ", raw chain name = "                                     << workingChain.name()
                              << ": about to try to write raw unified chain output file '" << m_optionsObj->m_ov.m_rawChainDataOutputFileName
                              << "."                                                       << m_optionsObj->m_ov.m_rawChainDataOutputFileType
                              << "', subId = "                                             << m_env.subId()
                              << "..."
                              << std::endl;
    }

    workingChain.unifiedWriteContents(m_optionsObj->m_ov.m_rawChainDataOutputFileName,
                                      m_optionsObj->m_ov.m_rawChainDataOutputFileType);
    if ((m_env.subDisplayFile()                   ) &&
        (m_optionsObj->m_ov.m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateSequence()"
                              << ", prefix = "                                             << m_optionsObj->m_prefix
                              << ", raw chain name = "                                     << workingChain.name()
                              << ": returned from writing raw unified chain output file '" << m_optionsObj->m_ov.m_rawChainDataOutputFileName
                              << "."                                                       << m_optionsObj->m_ov.m_rawChainDataOutputFileType
                              << "', subId = "                                             << m_env.subId()
                              << std::endl;
    }

    if (writeLogLikelihood) {
      workingLogLikelihoodValues->unifiedWriteContents(m_optionsObj->m_ov.m_rawChainDataOutputFileName + "_loglikelihood",
                                                       m_optionsObj->m_ov.m_rawChainDataOutputFileType);
    }

    if (writeLogTarget) {
      workingLogTargetValues->unifiedWriteContents(m_optionsObj->m_ov.m_rawChainDataOutputFileName + "_logtarget",
                                                   m_optionsObj->m_ov.m_rawChainDataOutputFileType);
    }

    // Compute raw unified MLE
    if (workingLogLikelihoodValues) {
      SequenceOfVectors<P_V,P_M> rawUnifiedMLEpositions(m_vectorSpace,0,m_optionsObj->m_prefix+"rawUnifiedMLEseq");
      double rawUnifiedMLEvalue = workingChain.unifiedPositionsOfMaximum(*workingLogLikelihoodValues,
                                                                         rawUnifiedMLEpositions);
      UQ_FATAL_TEST_MACRO(rawUnifiedMLEpositions.subSequenceSize() == 0,
                          m_env.worldRank(),
                          "MetropolisHastingsSG<P_V,P_M>::generateSequence()",
                          "rawUnifiedMLEpositions.subSequenceSize() = 0");

      if ((m_env.subDisplayFile()                   ) &&
          (m_optionsObj->m_ov.m_totallyMute == false)) {
        P_V tmpVec(m_vectorSpace.zeroVector());
        rawUnifiedMLEpositions.getPositionValues(0,tmpVec);
        *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateSequence()"
                                << ": just computed MLE"
                                << ", rawUnifiedMLEvalue = "                       << rawUnifiedMLEvalue
                                << ", rawUnifiedMLEpositions.subSequenceSize() = " << rawUnifiedMLEpositions.subSequenceSize()
                                << ", rawUnifiedMLEpositions[0] = "                << tmpVec
                                << std::endl;
      }
    }

    // Compute raw unified MAP
    if (workingLogTargetValues) {
      SequenceOfVectors<P_V,P_M> rawUnifiedMAPpositions(m_vectorSpace,0,m_optionsObj->m_prefix+"rawUnifiedMAPseq");
      double rawUnifiedMAPvalue = workingChain.unifiedPositionsOfMaximum(*workingLogTargetValues,
                                                                         rawUnifiedMAPpositions);

      UQ_FATAL_TEST_MACRO(rawUnifiedMAPpositions.subSequenceSize() == 0,
                          m_env.worldRank(),
                          "MetropolisHastingsSG<P_V,P_M>::generateSequence()",
                          "rawUnifiedMAPpositions.subSequenceSize() = 0");

      if ((m_env.subDisplayFile()                   ) &&
          (m_optionsObj->m_ov.m_totallyMute == false)) {
        P_V tmpVec(m_vectorSpace.zeroVector());
        rawUnifiedMAPpositions.getPositionValues(0,tmpVec);
        *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateSequence()"
                                << ": just computed MAP"
                                << ", rawUnifiedMAPvalue = "                       << rawUnifiedMAPvalue
                                << ", rawUnifiedMAPpositions.subSequenceSize() = " << rawUnifiedMAPpositions.subSequenceSize()
                                << ", rawUnifiedMAPpositions[0] = "                << tmpVec
                                << std::endl;
      }
    }
  }

  // Take care of other aspects of raw chain
  if ((genericFilePtrSet.ofsVar                 ) &&
      (m_optionsObj->m_ov.m_totallyMute == false)) {
    // Write likelihoodValues and alphaValues, if they were requested by user
    iRC = writeInfo(workingChain,
                    *genericFilePtrSet.ofsVar);
    UQ_FATAL_RC_MACRO(iRC,
                      m_env.worldRank(),
                      "MetropolisHastingsSG<P_V,P_M>::generateSequence()",
                      "improper writeInfo() return");
  }

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  if (m_optionsObj->m_ov.m_rawChainComputeStats) {
    workingChain.computeStatistics(*m_optionsObj->m_rawChainStatisticalOptionsObj,
                                   genericFilePtrSet.ofsVar);
  }
#endif

  //****************************************************************************************
  // Eventually:
  // --> filter the raw chain
  // --> write it
  // --> compute statistics on it
  //****************************************************************************************
  if (m_optionsObj->m_ov.m_filteredChainGenerate) {
    // Compute filter parameters
    unsigned int filterInitialPos = (unsigned int) (m_optionsObj->m_ov.m_filteredChainDiscardedPortion * (double) workingChain.subSequenceSize());
    unsigned int filterSpacing    = m_optionsObj->m_ov.m_filteredChainLag;
    if (filterSpacing == 0) {
      workingChain.computeFilterParams(genericFilePtrSet.ofsVar,
                                       filterInitialPos,
                                       filterSpacing);
    }

    // Filter positions from the converged portion of the chain
    workingChain.filter(filterInitialPos,
                        filterSpacing);
    workingChain.setName(m_optionsObj->m_prefix + "filtChain");

    if (workingLogLikelihoodValues) workingLogLikelihoodValues->filter(filterInitialPos,
                                                                       filterSpacing);

    if (workingLogTargetValues) workingLogTargetValues->filter(filterInitialPos,
                                                               filterSpacing);

    // Write filtered chain
    if ((m_env.subDisplayFile()                   ) &&
        (m_optionsObj->m_ov.m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateSequence()"
                              << ", prefix = "                                                      << m_optionsObj->m_prefix
                              << ": checking necessity of opening output files for filtered chain " << workingChain.name()
                              << "..."
                              << std::endl;
    }

    // Take "sub" care of filtered chain
    if ((m_optionsObj->m_ov.m_filteredChainDataOutputFileName != UQ_MH_SG_FILENAME_FOR_NO_FILE) &&
        (m_optionsObj->m_ov.m_totallyMute == false                                            )) {
      workingChain.subWriteContents(0,
                                    workingChain.subSequenceSize(),
                                    m_optionsObj->m_ov.m_filteredChainDataOutputFileName,
                                    m_optionsObj->m_ov.m_filteredChainDataOutputFileType,
                                    m_optionsObj->m_ov.m_filteredChainDataOutputAllowedSet);
      if ((m_env.subDisplayFile()                   ) &&
          (m_optionsObj->m_ov.m_totallyMute == false)) {
        *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateSequence()"
                                << ", prefix = "                << m_optionsObj->m_prefix
                                << ": closed sub output file '" << m_optionsObj->m_ov.m_filteredChainDataOutputFileName
                                << "' for filtered chain "      << workingChain.name()
                                << std::endl;
      }

      if (writeLogLikelihood) {
        workingLogLikelihoodValues->subWriteContents(0,
                                                     workingChain.subSequenceSize(),
                                                     m_optionsObj->m_ov.m_filteredChainDataOutputFileName + "_loglikelihood",
                                                     m_optionsObj->m_ov.m_filteredChainDataOutputFileType,
                                                     m_optionsObj->m_ov.m_filteredChainDataOutputAllowedSet);
      }

      if (writeLogTarget) {
        workingLogTargetValues->subWriteContents(0,
                                                 workingChain.subSequenceSize(),
                                                 m_optionsObj->m_ov.m_filteredChainDataOutputFileName + "_logtarget",
                                                 m_optionsObj->m_ov.m_filteredChainDataOutputFileType,
                                                 m_optionsObj->m_ov.m_filteredChainDataOutputAllowedSet);
      }
    }

    // Compute sub filtered MLE and sub filtered MAP

    // Take "unified" care of filtered chain
    if ((m_optionsObj->m_ov.m_filteredChainDataOutputFileName != UQ_MH_SG_FILENAME_FOR_NO_FILE) &&
        (m_optionsObj->m_ov.m_totallyMute == false                                            )) {
      workingChain.unifiedWriteContents(m_optionsObj->m_ov.m_filteredChainDataOutputFileName,
                                        m_optionsObj->m_ov.m_filteredChainDataOutputFileType);
      if ((m_env.subDisplayFile()                   ) &&
          (m_optionsObj->m_ov.m_totallyMute == false)) {
        *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateSequence()"
                                << ", prefix = "                    << m_optionsObj->m_prefix
                                << ": closed unified output file '" << m_optionsObj->m_ov.m_filteredChainDataOutputFileName
                                << "' for filtered chain "          << workingChain.name()
                                << std::endl;
      }

      if (writeLogLikelihood) {
        workingLogLikelihoodValues->unifiedWriteContents(m_optionsObj->m_ov.m_filteredChainDataOutputFileName + "_loglikelihood",
                                                         m_optionsObj->m_ov.m_filteredChainDataOutputFileType);
      }

      if (writeLogTarget) {
        workingLogTargetValues->unifiedWriteContents(m_optionsObj->m_ov.m_filteredChainDataOutputFileName + "_logtarget",
                                                     m_optionsObj->m_ov.m_filteredChainDataOutputFileType);
      }
    }

    // Compute unified filtered MLE and unified filtered MAP

    // Compute statistics
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
    if (m_optionsObj->m_ov.m_filteredChainComputeStats) {
      workingChain.computeStatistics(*m_optionsObj->m_filteredChainStatisticalOptionsObj,
                                     genericFilePtrSet.ofsVar);
    }
#endif
  }

  //****************************************************
  // Close generic output file
  //****************************************************
  if (genericFilePtrSet.ofsVar) {
    //genericFilePtrSet.ofsVar->close();
    delete genericFilePtrSet.ofsVar;
    if ((m_env.subDisplayFile()                   ) &&
        (m_optionsObj->m_ov.m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateSequence()"
                              << ", prefix = "                    << m_optionsObj->m_prefix
                              << ": closed generic output file '" << m_optionsObj->m_ov.m_dataOutputFileName
                              << "' (chain name is "              << workingChain.name()
                              << ")"
                              << std::endl;
    }
  }

  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_ov.m_totallyMute == false)) {
    *m_env.subDisplayFile() << std::endl;
  }

  //m_env.syncPrintDebugMsg("Leaving MetropolisHastingsSG<P_V,P_M>::generateSequence()",2,3000000,m_env.fullComm()); // Dangerous to barrier on fullComm ... // KAUST

  if ((m_env.subDisplayFile()                   ) &&
      (m_env.displayVerbosity() >= 5            ) &&
      (m_optionsObj->m_ov.m_totallyMute == false)) {
    *m_env.subDisplayFile() << "Leaving MetropolisHastingsSG<P_V,P_M>::generateSequence()"
                            << std::endl;
  }

  return;
}

// -------------------------------------------------
template<class P_V,class P_M>
void
MetropolisHastingsSG<P_V,P_M>::getRawChainInfo(MHRawChainInfoStruct& info) const
{
  info = m_rawChainInfo;
  return;
}
//--------------------------------------------------
template <class P_V,class P_M>
void
MetropolisHastingsSG<P_V,P_M>::readFullChain(
  const std::string&                  inputFileName,
  const std::string&                  inputFileType,
        unsigned int                  chainSize,
  BaseVectorSequence<P_V,P_M>& workingChain)
{
  workingChain.unifiedReadContents(inputFileName,inputFileType,chainSize);
  return;
}
// Private methods ---------------------------------
template <class P_V,class P_M>
void
MetropolisHastingsSG<P_V,P_M>::generateFullChain(
  const P_V&                          valuesOf1stPosition,
        unsigned int                  chainSize,
  BaseVectorSequence<P_V,P_M>& workingChain,
  ScalarSequence<double>*      workingLogLikelihoodValues,
  ScalarSequence<double>*      workingLogTargetValues)
{
  //m_env.syncPrintDebugMsg("Entering MetropolisHastingsSG<P_V,P_M>::generateFullChain()",3,3000000,m_env.fullComm()); // Dangerous to barrier on fullComm ... // KAUST

  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_ov.m_totallyMute == false)) {
    *m_env.subDisplayFile() << "Starting the generation of Markov chain " << workingChain.name()
                            << ", with "                                  << chainSize
                            << " positions..."
                            << std::endl;
  }

  int iRC = UQ_OK_RC;
  struct timeval timevalChain;
  struct timeval timevalCandidate;
  struct timeval timevalTarget;
  struct timeval timevalMhAlpha;
  struct timeval timevalDrAlpha;
  struct timeval timevalDR;
  struct timeval timevalAM;

  m_positionIdForDebugging = 0;
  m_stageIdForDebugging    = 0;

  m_rawChainInfo.reset();

  iRC = gettimeofday(&timevalChain, NULL);

  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_ov.m_totallyMute == false)) {
    *m_env.subDisplayFile() << "\nIn MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                            << ": contents of initial position are:";
    *m_env.subDisplayFile() << valuesOf1stPosition; // FIX ME: might need parallelism
    *m_env.subDisplayFile() << "\nIn MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                            << ": targetPdf.domaintSet() info is:"
                            << m_targetPdf.domainSet();
    *m_env.subDisplayFile() << std::endl;
  }

  // Set a flag to write out log likelihood or not
  bool writeLogLikelihood;
  if ((workingLogLikelihoodValues != NULL) &&
      (m_optionsObj->m_ov.m_outputLogLikelihood)) {
    writeLogLikelihood = true;
  }
  else {
    writeLogLikelihood = false;
  }

  // Set a flag to write out log target or not
  bool writeLogTarget;
  if ((workingLogTargetValues != NULL) &&
      (m_optionsObj->m_ov.m_outputLogTarget)) {
    writeLogTarget = true;
  }
  else {
    writeLogTarget = false;
  }

  bool outOfTargetSupport = !m_targetPdf.domainSet().contains(valuesOf1stPosition);
  if ((m_env.subDisplayFile()) &&
      (outOfTargetSupport    )) {
    *m_env.subDisplayFile() << "ERROR: In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                            << ": contents of initial position are:\n";
    *m_env.subDisplayFile() << valuesOf1stPosition; // FIX ME: might need parallelism
    *m_env.subDisplayFile() << "\nERROR: In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                            << ": targetPdf.domaintSet() info is:\n"
                            << m_targetPdf.domainSet();
    *m_env.subDisplayFile() << std::endl;
  }
  UQ_FATAL_TEST_MACRO(outOfTargetSupport,
                      m_env.worldRank(),
                      "MetropolisHastingsSG<P_V,P_M>::generateFullChain()",
                      "initial position should not be out of target pdf support");
  
  double logPrior      = 0.;
  double logLikelihood = 0.;
  double logTarget     = 0.;
  if (m_computeInitialPriorAndLikelihoodValues) {
    if (m_optionsObj->m_ov.m_rawChainMeasureRunTimes) iRC = gettimeofday(&timevalTarget, NULL);
#ifdef QUESO_EXPECTS_LN_LIKELIHOOD_INSTEAD_OF_MINUS_2_LN
    logTarget =        m_targetPdfSynchronizer->callFunction(&valuesOf1stPosition,NULL,NULL,NULL,NULL,&logPrior,&logLikelihood); // Might demand parallel environment // KEY
#else
    logTarget = -0.5 * m_targetPdfSynchronizer->callFunction(&valuesOf1stPosition,NULL,NULL,NULL,NULL,&logPrior,&logLikelihood); // Might demand parallel environment
#endif
    if (m_optionsObj->m_ov.m_rawChainMeasureRunTimes) {
      m_rawChainInfo.targetRunTime += MiscGetEllapsedSeconds(&timevalTarget);
    }
    m_rawChainInfo.numTargetCalls++;
    if ((m_env.subDisplayFile()                   ) &&
        (m_env.displayVerbosity() >= 3            ) &&
        (m_optionsObj->m_ov.m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
			      << ": just returned from likelihood() for initial chain position"
			      << ", m_rawChainInfo.numTargetCalls = " << m_rawChainInfo.numTargetCalls
			      << ", logPrior = "      << logPrior
			      << ", logLikelihood = " << logLikelihood
			      << ", logTarget = "     << logTarget
			      << std::endl;
    }
  }
  else {
    logPrior       = m_initialLogPriorValue;
    logLikelihood  = m_initialLogLikelihoodValue;
    logTarget      = logPrior + logLikelihood;
    if ((m_env.subDisplayFile()                   ) &&
        (m_env.displayVerbosity() >= 3            ) &&
        (m_optionsObj->m_ov.m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
			      << ": used input prior and likelihood for initial chain position"
			      << ", m_rawChainInfo.numTargetCalls = " << m_rawChainInfo.numTargetCalls
			      << ", logPrior = "      << logPrior
			      << ", logLikelihood = " << logLikelihood
			      << ", logTarget = "     << logTarget
			      << std::endl;
    }
  }

  //*m_env.subDisplayFile() << "AQUI 001" << std::endl;
  MarkovChainPositionData<P_V> currentPositionData(m_env,
                                                          valuesOf1stPosition,
                                                          outOfTargetSupport,
                                                          logLikelihood,
                                                          logTarget);

  P_V gaussianVector(m_vectorSpace.zeroVector());
  P_V tmpVecValues(m_vectorSpace.zeroVector());
  MarkovChainPositionData<P_V> currentCandidateData(m_env);

  //****************************************************
  // Set chain position with positionId = 0
  //****************************************************
  workingChain.resizeSequence(chainSize);
  m_numPositionsNotSubWritten = 0;
  if (workingLogLikelihoodValues) workingLogLikelihoodValues->resizeSequence(chainSize);
  if (workingLogTargetValues    ) workingLogTargetValues->resizeSequence    (chainSize);
  if (true/*m_uniqueChainGenerate*/) m_idsOfUniquePositions.resize(chainSize,0);
  if (m_optionsObj->m_ov.m_rawChainGenerateExtra) {
    m_logTargets.resize    (chainSize,0.);
    m_alphaQuotients.resize(chainSize,0.);
  }

  unsigned int uniquePos = 0;
  workingChain.setPositionValues(0,currentPositionData.vecValues());
  m_numPositionsNotSubWritten++;
  if ((m_optionsObj->m_ov.m_rawChainDataOutputPeriod           >  0  ) &&
      (((0+1) % m_optionsObj->m_ov.m_rawChainDataOutputPeriod) == 0  ) &&
      (m_optionsObj->m_ov.m_rawChainDataOutputFileName         != ".")) {
    workingChain.subWriteContents(0 + 1 - m_optionsObj->m_ov.m_rawChainDataOutputPeriod,
                                  m_optionsObj->m_ov.m_rawChainDataOutputPeriod,
                                  m_optionsObj->m_ov.m_rawChainDataOutputFileName,
                                  m_optionsObj->m_ov.m_rawChainDataOutputFileType,
                                  m_optionsObj->m_ov.m_rawChainDataOutputAllowedSet);
    if ((m_env.subDisplayFile()                   ) &&
        (m_optionsObj->m_ov.m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                              << ": just wrote (per period request) " << m_numPositionsNotSubWritten << " chain positions "
                              << ", " << 0 + 1 - m_optionsObj->m_ov.m_rawChainDataOutputPeriod << " <= pos <= " << 0
                              << std::endl;
    }

    if (writeLogLikelihood) {
      workingLogLikelihoodValues->subWriteContents(0 + 1 - m_optionsObj->m_ov.m_rawChainDataOutputPeriod,
                                                   m_optionsObj->m_ov.m_rawChainDataOutputPeriod,
                                                   m_optionsObj->m_ov.m_rawChainDataOutputFileName + "_loglikelihood",
                                                   m_optionsObj->m_ov.m_rawChainDataOutputFileType,
                                                   m_optionsObj->m_ov.m_rawChainDataOutputAllowedSet);
    }

    if (writeLogTarget) {
      workingLogTargetValues->subWriteContents(0 + 1 - m_optionsObj->m_ov.m_rawChainDataOutputPeriod,
                                               m_optionsObj->m_ov.m_rawChainDataOutputPeriod,
                                               m_optionsObj->m_ov.m_rawChainDataOutputFileName + "_logtarget",
                                               m_optionsObj->m_ov.m_rawChainDataOutputFileType,
                                               m_optionsObj->m_ov.m_rawChainDataOutputAllowedSet);
    }

    m_numPositionsNotSubWritten = 0;
  }

  if (workingLogLikelihoodValues) (*workingLogLikelihoodValues)[0] = currentPositionData.logLikelihood();
  if (workingLogTargetValues    ) (*workingLogTargetValues    )[0] = currentPositionData.logTarget();
  if (true/*m_uniqueChainGenerate*/) m_idsOfUniquePositions[uniquePos++] = 0;
  if (m_optionsObj->m_ov.m_rawChainGenerateExtra) {
    m_logTargets    [0] = currentPositionData.logTarget();
    m_alphaQuotients[0] = 1.;
  }
  //*m_env.subDisplayFile() << "AQUI 002" << std::endl;

  if ((m_env.subDisplayFile()                   ) &&
      (m_env.displayVerbosity() >= 10           ) &&
      (m_optionsObj->m_ov.m_totallyMute == false)) {
    *m_env.subDisplayFile() << "\n"
                            << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                            << "\n"
                            << std::endl;
  }

  //m_env.syncPrintDebugMsg("In MetropolisHastingsSG<P_V,P_M>::generateFullChain(), right before main loop",3,3000000,m_env.fullComm()); // Dangerous to barrier on fullComm ... // KAUST

  //****************************************************
  // Begin chain loop from positionId = 1
  //****************************************************
  if ((m_env.numSubEnvironments() < (unsigned int) m_env.fullComm().NumProc()) &&
      (m_initialPosition.numOfProcsForStorage() == 1                         ) &&
      (m_env.subRank()                          != 0                         )) {
    // subRank != 0 --> Enter the barrier and wait for processor 0 to decide to call the targetPdf
    double aux = 0.;
    aux = m_targetPdfSynchronizer->callFunction(NULL,
                                                NULL,
                                                NULL,
                                                NULL,
                                                NULL,
                                                NULL,
                                                NULL);
    if (aux) {}; // just to remove compiler warning
    for (unsigned int positionId = 1; positionId < workingChain.subSequenceSize(); ++positionId) {
      // Multiply by position values by 'positionId' in order to avoid a constant sequence,
      // which would cause zero variance and eventually OVERFLOW flags raised
      workingChain.setPositionValues(positionId,((double) positionId) * currentPositionData.vecValues());
      m_rawChainInfo.numRejections++;
    }
  }
  else for (unsigned int positionId = 1; positionId < workingChain.subSequenceSize(); ++positionId) {
    //****************************************************
    // Point 1/6 of logic for new position
    // Loop: initialize variables and print some information
    //****************************************************
    m_positionIdForDebugging = positionId;
    if ((m_env.subDisplayFile()                   ) &&
        (m_env.displayVerbosity() >= 3            ) &&
        (m_optionsObj->m_ov.m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                              << ": beginning chain position of id = "  << positionId
                              << ", m_optionsObj->m_ov.m_drMaxNumExtraStages = " << m_optionsObj->m_ov.m_drMaxNumExtraStages
                              << std::endl;
    }
    unsigned int stageId  = 0;
    m_stageIdForDebugging = stageId;

    m_tk->clearPreComputingPositions();

    if ((m_env.subDisplayFile()                   ) &&
        (m_env.displayVerbosity() >= 5            ) &&
        (m_optionsObj->m_ov.m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                              << ": about to set TK pre computing position of local id " << 0
                              << ", values = " << currentPositionData.vecValues()
                              << std::endl;
    }
    bool validPreComputingPosition = m_tk->setPreComputingPosition(currentPositionData.vecValues(),0);
    if ((m_env.subDisplayFile()                   ) &&
        (m_env.displayVerbosity() >= 5            ) &&
        (m_optionsObj->m_ov.m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                              << ": returned from setting TK pre computing position of local id " << 0
                              << ", values = " << currentPositionData.vecValues()
                              << ", valid = "  << validPreComputingPosition
                              << std::endl;
    }
    UQ_FATAL_TEST_MACRO(validPreComputingPosition == false,
                        m_env.worldRank(),
                        "MetropolisHastingsSG<P_V,P_M>::generateFullChain()",
                        "initial position should not be an invalid pre computing position");

    //****************************************************
    // Point 2/6 of logic for new position
    // Loop: generate new position
    //****************************************************
    // sep2011
    bool keepGeneratingCandidates = true;
    while (keepGeneratingCandidates) {
      if (m_optionsObj->m_ov.m_rawChainMeasureRunTimes) iRC = gettimeofday(&timevalCandidate, NULL);
      m_tk->rv(0).realizer().realization(tmpVecValues);
      if (m_numDisabledParameters > 0) { // gpmsa2
        for (unsigned int paramId = 0; paramId < m_vectorSpace.dimLocal(); ++paramId) {
          if (m_parameterEnabledStatus[paramId] == false) {
            tmpVecValues[paramId] = m_initialPosition[paramId];
          }
        }
      }
      if (m_optionsObj->m_ov.m_rawChainMeasureRunTimes) m_rawChainInfo.candidateRunTime += MiscGetEllapsedSeconds(&timevalCandidate);

      outOfTargetSupport = !m_targetPdf.domainSet().contains(tmpVecValues);

      bool displayDetail = (m_env.displayVerbosity() >= 10/*99*/) || m_optionsObj->m_ov.m_displayCandidates;
      if ((m_env.subDisplayFile()                   ) &&
          (displayDetail                            ) &&
          (m_optionsObj->m_ov.m_totallyMute == false)) {
        *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                                << ": for chain position of id = " << positionId
                                << ", candidate = "                << tmpVecValues // FIX ME: might need parallelism
                                << ", outOfTargetSupport = "       << outOfTargetSupport
                                << std::endl;
      }

      if (m_optionsObj->m_ov.m_putOutOfBoundsInChain) keepGeneratingCandidates = false;
      else                                            keepGeneratingCandidates = outOfTargetSupport;
    }

    if ((m_env.subDisplayFile()                   ) &&
        (m_env.displayVerbosity() >= 5            ) &&
        (m_optionsObj->m_ov.m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                              << ": about to set TK pre computing position of local id " << stageId+1
                              << ", values = " << tmpVecValues
                              << std::endl;
    }
    validPreComputingPosition = m_tk->setPreComputingPosition(tmpVecValues,stageId+1);
    if ((m_env.subDisplayFile()                   ) &&
        (m_env.displayVerbosity() >= 5            ) &&
        (m_optionsObj->m_ov.m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                              << ": returned from setting TK pre computing position of local id " << stageId+1
                              << ", values = " << tmpVecValues
                              << ", valid = "  << validPreComputingPosition
                              << std::endl;
    }

    if (outOfTargetSupport) {
      m_rawChainInfo.numOutOfTargetSupport++;
      logPrior      = -INFINITY;
      logLikelihood = -INFINITY;
      logTarget     = -INFINITY;
    }
    else {
      if (m_optionsObj->m_ov.m_rawChainMeasureRunTimes) iRC = gettimeofday(&timevalTarget, NULL);
#ifdef QUESO_EXPECTS_LN_LIKELIHOOD_INSTEAD_OF_MINUS_2_LN
      logTarget =        m_targetPdfSynchronizer->callFunction(&tmpVecValues,NULL,NULL,NULL,NULL,&logPrior,&logLikelihood); // Might demand parallel environment
#else
      logTarget = -0.5 * m_targetPdfSynchronizer->callFunction(&tmpVecValues,NULL,NULL,NULL,NULL,&logPrior,&logLikelihood); // Might demand parallel environment
#endif
      if (m_optionsObj->m_ov.m_rawChainMeasureRunTimes) m_rawChainInfo.targetRunTime += MiscGetEllapsedSeconds(&timevalTarget);
      m_rawChainInfo.numTargetCalls++;
      if ((m_env.subDisplayFile()                   ) &&
          (m_env.displayVerbosity() >= 3            ) &&
          (m_optionsObj->m_ov.m_totallyMute == false)) {
        *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                                << ": just returned from likelihood() for chain position of id " << positionId
                                << ", m_rawChainInfo.numTargetCalls = " << m_rawChainInfo.numTargetCalls
                                << ", logPrior = "      << logPrior
                                << ", logLikelihood = " << logLikelihood
                                << ", logTarget = "     << logTarget
                                << std::endl;
      }
    }
    currentCandidateData.set(tmpVecValues,
                             outOfTargetSupport,
                             logLikelihood,
                             logTarget);

    if ((m_env.subDisplayFile()                   ) &&
        (m_env.displayVerbosity() >= 10           ) &&
        (m_optionsObj->m_ov.m_totallyMute == false)) {
      *m_env.subDisplayFile() << "\n"
                              << "\n-----------------------------------------------------------\n"
                              << "\n"
                              << std::endl;
    }
    bool accept = false;
    double alphaFirstCandidate = 0.;
    if (outOfTargetSupport) {
      if (m_optionsObj->m_ov.m_rawChainGenerateExtra) {
        m_alphaQuotients[positionId] = 0.;
      }
    }
    else {
      if (m_optionsObj->m_ov.m_rawChainMeasureRunTimes) iRC = gettimeofday(&timevalMhAlpha, NULL);
      if (m_optionsObj->m_ov.m_rawChainGenerateExtra) {
        alphaFirstCandidate = this->alpha(currentPositionData,currentCandidateData,0,1,&m_alphaQuotients[positionId]);
      }
      else {
        alphaFirstCandidate = this->alpha(currentPositionData,currentCandidateData,0,1,NULL);
      }
      if (m_optionsObj->m_ov.m_rawChainMeasureRunTimes) m_rawChainInfo.mhAlphaRunTime += MiscGetEllapsedSeconds(&timevalMhAlpha);
      if ((m_env.subDisplayFile()                   ) &&
          (m_env.displayVerbosity() >= 10           ) &&
          (m_optionsObj->m_ov.m_totallyMute == false)) {
        *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                                << ": for chain position of id = " << positionId
                                << std::endl;
      }
      accept = acceptAlpha(alphaFirstCandidate);
    }

    bool displayDetail = (m_env.displayVerbosity() >= 10/*99*/) || m_optionsObj->m_ov.m_displayCandidates;
    if ((m_env.subDisplayFile()                   ) &&
        (displayDetail                            ) &&
        (m_optionsObj->m_ov.m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                              << ": for chain position of id = " << positionId
                              << ", outOfTargetSupport = "       << outOfTargetSupport
                              << ", alpha = "                    << alphaFirstCandidate
                              << ", accept = "                   << accept
                              << ", currentCandidateData.vecValues() = ";
      *m_env.subDisplayFile() << currentCandidateData.vecValues(); // FIX ME: might need parallelism
      *m_env.subDisplayFile() << "\n"
                              << "\n curLogTarget  = "           << currentPositionData.logTarget()
                              << "\n"
                              << "\n canLogTarget  = "           << currentCandidateData.logTarget()
                              << "\n"
                              << std::endl;
    }
    if ((m_env.subDisplayFile()                   ) &&
        (m_env.displayVerbosity() >= 10           ) &&
        (m_optionsObj->m_ov.m_totallyMute == false)) {
      *m_env.subDisplayFile() << "\n"
                              << "\n-----------------------------------------------------------\n"
                              << "\n"
                              << std::endl;
    }

    //****************************************************
    // Point 3/6 of logic for new position
    // Loop: delayed rejection
    //****************************************************
    // sep2011
    std::vector<MarkovChainPositionData<P_V>*> drPositionsData(stageId+2,NULL);
    std::vector<unsigned int> tkStageIds (stageId+2,0);
    if ((accept                                   == false) &&
        (outOfTargetSupport                       == false) && // IMPORTANT
        (m_optionsObj->m_ov.m_drMaxNumExtraStages >  0    )) {
      if ((m_optionsObj->m_ov.m_drDuringAmNonAdaptiveInt  == false     ) &&
          (m_optionsObj->m_ov.m_tkUseLocalHessian         == false     ) &&
          (m_optionsObj->m_ov.m_amInitialNonAdaptInterval >  0         ) &&
          (m_optionsObj->m_ov.m_amAdaptInterval           >  0         ) &&
          (positionId <= m_optionsObj->m_ov.m_amInitialNonAdaptInterval)) {
        // Avoid DR now
      }
      else {
        if (m_optionsObj->m_ov.m_rawChainMeasureRunTimes) iRC = gettimeofday(&timevalDR, NULL);

        drPositionsData[0] = new MarkovChainPositionData<P_V>(currentPositionData );
        drPositionsData[1] = new MarkovChainPositionData<P_V>(currentCandidateData);

        tkStageIds[0] = 0;
        tkStageIds[1] = 1;

        while ((validPreComputingPosition == true                 ) &&
               (accept                    == false                ) &&
               (stageId < m_optionsObj->m_ov.m_drMaxNumExtraStages)) {
          if ((m_env.subDisplayFile()                   ) &&
              (m_env.displayVerbosity() >= 10           ) &&
              (m_optionsObj->m_ov.m_totallyMute == false)) {
            *m_env.subDisplayFile() << "\n"
                                    << "\n+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-"
                                    << "\n"
                                    << std::endl;
          }
          m_rawChainInfo.numDRs++;
          stageId++;
          m_stageIdForDebugging = stageId;
          if ((m_env.subDisplayFile()                   ) &&
              (m_env.displayVerbosity() >= 10           ) &&
              (m_optionsObj->m_ov.m_totallyMute == false)) {
            *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                                    << ": for chain position of id = " << positionId
                                    << ", beginning stageId = "        << stageId
                                    << std::endl;
          }

          keepGeneratingCandidates = true;
          while (keepGeneratingCandidates) {
            if (m_optionsObj->m_ov.m_rawChainMeasureRunTimes) iRC = gettimeofday(&timevalCandidate, NULL);
            m_tk->rv(tkStageIds).realizer().realization(tmpVecValues);
            if (m_numDisabledParameters > 0) { // gpmsa2
              for (unsigned int paramId = 0; paramId < m_vectorSpace.dimLocal(); ++paramId) {
                if (m_parameterEnabledStatus[paramId] == false) {
                  tmpVecValues[paramId] = m_initialPosition[paramId];
                }
              }
            }
            if (m_optionsObj->m_ov.m_rawChainMeasureRunTimes) m_rawChainInfo.candidateRunTime += MiscGetEllapsedSeconds(&timevalCandidate);

            outOfTargetSupport = !m_targetPdf.domainSet().contains(tmpVecValues);

            if (m_optionsObj->m_ov.m_putOutOfBoundsInChain) keepGeneratingCandidates = false;
            else                                            keepGeneratingCandidates = outOfTargetSupport;
          }

          if ((m_env.subDisplayFile()                   ) &&
              (m_env.displayVerbosity() >= 5            ) &&
              (m_optionsObj->m_ov.m_totallyMute == false)) {
            *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                                    << ": about to set TK pre computing position of local id " << stageId+1
                                    << ", values = " << tmpVecValues
                                    << std::endl;
          }
          validPreComputingPosition = m_tk->setPreComputingPosition(tmpVecValues,stageId+1);
          if ((m_env.subDisplayFile()                   ) &&
              (m_env.displayVerbosity() >= 5            ) &&
              (m_optionsObj->m_ov.m_totallyMute == false)) {
            *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                                    << ": returned from setting TK pre computing position of local id " << stageId+1
                                    << ", values = " << tmpVecValues
                                    << ", valid = "  << validPreComputingPosition
                                    << std::endl;
          }

          if (outOfTargetSupport) {
            m_rawChainInfo.numOutOfTargetSupportInDR++; // new 2010/May/12
            logPrior      = -INFINITY;
            logLikelihood = -INFINITY;
            logTarget     = -INFINITY;
          }
          else {
            if (m_optionsObj->m_ov.m_rawChainMeasureRunTimes) iRC = gettimeofday(&timevalTarget, NULL);
#ifdef QUESO_EXPECTS_LN_LIKELIHOOD_INSTEAD_OF_MINUS_2_LN
            logTarget =        m_targetPdfSynchronizer->callFunction(&tmpVecValues,NULL,NULL,NULL,NULL,&logPrior,&logLikelihood); // Might demand parallel environment
#else
            logTarget = -0.5 * m_targetPdfSynchronizer->callFunction(&tmpVecValues,NULL,NULL,NULL,NULL,&logPrior,&logLikelihood); // Might demand parallel environment
#endif
            if (m_optionsObj->m_ov.m_rawChainMeasureRunTimes) m_rawChainInfo.targetRunTime += MiscGetEllapsedSeconds(&timevalTarget);
            m_rawChainInfo.numTargetCalls++;
            if ((m_env.subDisplayFile()                   ) &&
                (m_env.displayVerbosity() >= 3            ) &&
                (m_optionsObj->m_ov.m_totallyMute == false)) {
              *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                                      << ": just returned from likelihood() for chain position of id " << positionId
                                      << ", m_rawChainInfo.numTargetCalls = " << m_rawChainInfo.numTargetCalls
                                      << ", stageId = "       << stageId
                                      << ", logPrior = "      << logPrior
                                      << ", logLikelihood = " << logLikelihood
                                      << ", logTarget = "     << logTarget
                                      << std::endl;
            }
          }
          currentCandidateData.set(tmpVecValues,
                                   outOfTargetSupport,
                                   logLikelihood,
                                   logTarget);

          drPositionsData.push_back(new MarkovChainPositionData<P_V>(currentCandidateData));
          tkStageIds.push_back     (stageId+1);

          double alphaDR = 0.;
          if (outOfTargetSupport == false) {
            if (m_optionsObj->m_ov.m_rawChainMeasureRunTimes) iRC = gettimeofday(&timevalDrAlpha, NULL);
            alphaDR = this->alpha(drPositionsData,tkStageIds);
            if (m_optionsObj->m_ov.m_rawChainMeasureRunTimes) m_rawChainInfo.drAlphaRunTime += MiscGetEllapsedSeconds(&timevalDrAlpha);
            accept = acceptAlpha(alphaDR);
          }

          displayDetail = (m_env.displayVerbosity() >= 10/*99*/) || m_optionsObj->m_ov.m_displayCandidates;
          if ((m_env.subDisplayFile()                   ) &&
              (displayDetail                            ) &&
              (m_optionsObj->m_ov.m_totallyMute == false)) {
            *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                                    << ": for chain position of id = " << positionId
                                    << " and stageId = "               << stageId
                                    << ", outOfTargetSupport = "       << outOfTargetSupport
                                    << ", alpha = "                    << alphaDR
                                    << ", accept = "                   << accept
                                    << ", currentCandidateData.vecValues() = ";
            *m_env.subDisplayFile() << currentCandidateData.vecValues(); // FIX ME: might need parallelism
            *m_env.subDisplayFile() << std::endl;
          }
        } // while

        if (m_optionsObj->m_ov.m_rawChainMeasureRunTimes) m_rawChainInfo.drRunTime += MiscGetEllapsedSeconds(&timevalDR);
      } // if-else "Avoid DR now"
    } // end of 'delayed rejection' logic

    for (unsigned int i = 0; i < drPositionsData.size(); ++i) {
      if (drPositionsData[i]) delete drPositionsData[i];
    }

    //****************************************************
    // Point 4/6 of logic for new position
    // Loop: update chain
    //****************************************************
    if (accept) {
      workingChain.setPositionValues(positionId,currentCandidateData.vecValues());
      if (true/*m_uniqueChainGenerate*/) m_idsOfUniquePositions[uniquePos++] = positionId;
      currentPositionData = currentCandidateData;
    }
    else {
      workingChain.setPositionValues(positionId,currentPositionData.vecValues());
      m_rawChainInfo.numRejections++;
    }
    m_numPositionsNotSubWritten++;
    if ((m_optionsObj->m_ov.m_rawChainDataOutputPeriod                    >  0  ) &&
        (((positionId+1) % m_optionsObj->m_ov.m_rawChainDataOutputPeriod) == 0  ) &&
        (m_optionsObj->m_ov.m_rawChainDataOutputFileName                  != ".")) {
      if ((m_env.subDisplayFile()                   ) &&
          (m_env.displayVerbosity()         >= 10   ) &&
          (m_optionsObj->m_ov.m_totallyMute == false)) {
        *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                                << ", for chain position of id = " << positionId
                                << ": about to write (per period request) " << m_numPositionsNotSubWritten << " chain positions "
                                << ", " << positionId + 1 - m_optionsObj->m_ov.m_rawChainDataOutputPeriod << " <= pos <= " << positionId
                                << std::endl;
      }
      workingChain.subWriteContents(positionId + 1 - m_optionsObj->m_ov.m_rawChainDataOutputPeriod,
                                    m_optionsObj->m_ov.m_rawChainDataOutputPeriod,
                                    m_optionsObj->m_ov.m_rawChainDataOutputFileName,
                                    m_optionsObj->m_ov.m_rawChainDataOutputFileType,
                                    m_optionsObj->m_ov.m_rawChainDataOutputAllowedSet);
      if ((m_env.subDisplayFile()                   ) &&
          (m_optionsObj->m_ov.m_totallyMute == false)) {
        *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                                << ", for chain position of id = " << positionId
                                << ": just wrote (per period request) " << m_numPositionsNotSubWritten << " chain positions "
                                << ", " << positionId + 1 - m_optionsObj->m_ov.m_rawChainDataOutputPeriod << " <= pos <= " << positionId
                                << std::endl;
      }

      if (writeLogLikelihood) {
        workingLogLikelihoodValues->subWriteContents(0 + 1 - m_optionsObj->m_ov.m_rawChainDataOutputPeriod,
                                                     m_optionsObj->m_ov.m_rawChainDataOutputPeriod,
                                                     m_optionsObj->m_ov.m_rawChainDataOutputFileName + "_loglikelihood",
                                                     m_optionsObj->m_ov.m_rawChainDataOutputFileType,
                                                     m_optionsObj->m_ov.m_rawChainDataOutputAllowedSet);
      }

      if (writeLogTarget) {
        workingLogTargetValues->subWriteContents(0 + 1 - m_optionsObj->m_ov.m_rawChainDataOutputPeriod,
                                                 m_optionsObj->m_ov.m_rawChainDataOutputPeriod,
                                                 m_optionsObj->m_ov.m_rawChainDataOutputFileName + "_logtarget",
                                                 m_optionsObj->m_ov.m_rawChainDataOutputFileType,
                                                 m_optionsObj->m_ov.m_rawChainDataOutputAllowedSet);
      }

      m_numPositionsNotSubWritten = 0;
    }


    if (workingLogLikelihoodValues) (*workingLogLikelihoodValues)[positionId] = currentPositionData.logLikelihood();
    if (workingLogTargetValues    ) (*workingLogTargetValues    )[positionId] = currentPositionData.logTarget();

    if (m_optionsObj->m_ov.m_rawChainGenerateExtra) {
      m_logTargets[positionId] = currentPositionData.logTarget();
    }

    if( m_optionsObj->m_ov.m_enableBrooksGelmanConvMonitor > 0 ) {
      if( positionId%m_optionsObj->m_ov.m_enableBrooksGelmanConvMonitor == 0 &&
    positionId > m_optionsObj->m_ov.m_BrooksGelmanLag+1 ) { //+1 to help ensure there are at least 2 samples to use

  double conv_est = workingChain.estimateConvBrooksGelman( m_optionsObj->m_ov.m_BrooksGelmanLag,
                 positionId - m_optionsObj->m_ov.m_BrooksGelmanLag );

  if ( m_env.subDisplayFile() ) {
      *m_env.subDisplayFile() << "positionId = " << positionId
            << ", conv_est = " << conv_est << std::endl;
      (*m_env.subDisplayFile()).flush();
  }
      }
    }

    //****************************************************
    // Point 5/6 of logic for new position
    // Loop: adaptive Metropolis (adaptation of covariance matrix)
    //****************************************************
    // sep2011
    if ((m_optionsObj->m_ov.m_tkUseLocalHessian ==    false) && // IMPORTANT
        (m_optionsObj->m_ov.m_amInitialNonAdaptInterval > 0) &&
        (m_optionsObj->m_ov.m_amAdaptInterval           > 0)) {
      if (m_optionsObj->m_ov.m_rawChainMeasureRunTimes) iRC = gettimeofday(&timevalAM, NULL);

      // Now might be the moment to adapt
      unsigned int idOfFirstPositionInSubChain = 0;
      SequenceOfVectors<P_V,P_M> partialChain(m_vectorSpace,0,m_optionsObj->m_prefix+"partialChain");

      // Check if now is indeed the moment to adapt
      bool printAdaptedMatrix = false;
      if (positionId < m_optionsObj->m_ov.m_amInitialNonAdaptInterval) {
        // Do nothing
      }
      else if (positionId == m_optionsObj->m_ov.m_amInitialNonAdaptInterval) {
        idOfFirstPositionInSubChain = 0;
        partialChain.resizeSequence(m_optionsObj->m_ov.m_amInitialNonAdaptInterval+1);
        m_lastMean             = m_vectorSpace.newVector();
        m_lastAdaptedCovMatrix = m_vectorSpace.newMatrix();
        printAdaptedMatrix = true;
      }
      else {
        unsigned int interval = positionId - m_optionsObj->m_ov.m_amInitialNonAdaptInterval;
        if ((interval % m_optionsObj->m_ov.m_amAdaptInterval) == 0) {
          idOfFirstPositionInSubChain = positionId - m_optionsObj->m_ov.m_amAdaptInterval;
          partialChain.resizeSequence(m_optionsObj->m_ov.m_amAdaptInterval);

          if (m_optionsObj->m_ov.m_amAdaptedMatricesDataOutputPeriod > 0) {
            if ((interval % m_optionsObj->m_ov.m_amAdaptedMatricesDataOutputPeriod) == 0) {
              printAdaptedMatrix = true;
            }
          }
        }
      }

      // If now is indeed the moment to adapt, then do it!
      if (partialChain.subSequenceSize() > 0) {
        P_V transporterVec(m_vectorSpace.zeroVector());
        for (unsigned int i = 0; i < partialChain.subSequenceSize(); ++i) {
          workingChain.getPositionValues(idOfFirstPositionInSubChain+i,transporterVec);
          partialChain.setPositionValues(i,transporterVec);
        }
        updateAdaptedCovMatrix(partialChain,
                               idOfFirstPositionInSubChain,
                               m_lastChainSize,
                               *m_lastMean,
                               *m_lastAdaptedCovMatrix);

        if ((printAdaptedMatrix                                       == true) &&
            (m_optionsObj->m_ov.m_amAdaptedMatricesDataOutputFileName != "." )) { // palms
          char varNamePrefix[64];
          sprintf(varNamePrefix,"mat_am%d",positionId);

          char tmpChar[64];
          sprintf(tmpChar,"_am%d",positionId);

          std::set<unsigned int> tmpSet;
          tmpSet.insert(m_env.subId());

          m_lastAdaptedCovMatrix->subWriteContents(varNamePrefix,
                                                   (m_optionsObj->m_ov.m_amAdaptedMatricesDataOutputFileName+tmpChar),
                                                   m_optionsObj->m_ov.m_amAdaptedMatricesDataOutputFileType,
                                                   tmpSet);
          if ((m_env.subDisplayFile()                   ) &&
              (m_optionsObj->m_ov.m_totallyMute == false)) {
            *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                                    << ": just wrote last adapted proposal cov matrix contents = " << *m_lastAdaptedCovMatrix
                                    << std::endl;
          }
        } // if (printAdaptedMatrix && ...)

        bool tmpCholIsPositiveDefinite = false;
        P_M tmpChol(*m_lastAdaptedCovMatrix);
        P_M attemptedMatrix(tmpChol);
        if ((m_env.subDisplayFile()        ) &&
            (m_env.displayVerbosity() >= 10)) {
    //(m_optionsObj->m_ov.m_totallyMute == false)) {
          *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                                  << ", positionId = "  << positionId
                                  << ": 'am' calling first tmpChol.chol()"
                                  << std::endl;
        }
        iRC = tmpChol.chol();
        if (iRC) {
          std::cerr << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain(): first chol failed\n";
        }
        if ((m_env.subDisplayFile()        ) &&
            (m_env.displayVerbosity() >= 10)) {
    //(m_optionsObj->m_ov.m_totallyMute == false)) {
          *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                                  << ", positionId = "  << positionId
                                  << ": 'am' got first tmpChol.chol() with iRC = " << iRC
                                  << std::endl;
          if (iRC == 0) {
            double diagMult = 1.;
            for (unsigned int j = 0; j < tmpChol.numRowsLocal(); ++j) {
              diagMult *= tmpChol(j,j);
            }
            *m_env.subDisplayFile() << "diagMult = " << diagMult
                                    << std::endl;
          }
        }
#if 0 // tentative logic
        if (iRC == 0) {
          double diagMult = 1.;
          for (unsigned int j = 0; j < tmpChol.numRowsLocal(); ++j) {
            diagMult *= tmpChol(j,j);
          }
          if (diagMult < 1.e-40) {
            iRC = UQ_MATRIX_IS_NOT_POS_DEFINITE_RC;
          }
        }
#endif

        if (iRC) {
          UQ_FATAL_TEST_MACRO(iRC != UQ_MATRIX_IS_NOT_POS_DEFINITE_RC,
                              m_env.worldRank(),
                              "MetropolisHastingsSG<P_V,P_M>::generateFullChain()",
                              "invalid iRC returned from first chol()");
          // Matrix is not positive definite
          P_M* tmpDiag = m_vectorSpace.newDiagMatrix(m_optionsObj->m_ov.m_amEpsilon);
          if (m_numDisabledParameters > 0) { // gpmsa2
            for (unsigned int paramId = 0; paramId < m_vectorSpace.dimLocal(); ++paramId) {
              if (m_parameterEnabledStatus[paramId] == false) {
                (*tmpDiag)(paramId,paramId) = 0.;
              }
            }
          }
          tmpChol = *m_lastAdaptedCovMatrix + *tmpDiag;
          attemptedMatrix = tmpChol;
          delete tmpDiag;
          if ((m_env.subDisplayFile()        ) &&
              (m_env.displayVerbosity() >= 10)) {
      //(m_optionsObj->m_ov.m_totallyMute == false)) {
            *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                                    << ", positionId = "  << positionId
                                    << ": 'am' calling second tmpChol.chol()"
                                    << std::endl;
          }
          iRC = tmpChol.chol();
          if (iRC) {
            std::cerr << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain(): second chol failed\n";
          }
          if ((m_env.subDisplayFile()        ) &&
              (m_env.displayVerbosity() >= 10)) {
      //(m_optionsObj->m_ov.m_totallyMute == false)) {
            *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                                    << ", positionId = " << positionId
                                    << ": 'am' got second tmpChol.chol() with iRC = " << iRC
                                    << std::endl;
            if (iRC == 0) {
              double diagMult = 1.;
              for (unsigned int j = 0; j < tmpChol.numRowsLocal(); ++j) {
                diagMult *= tmpChol(j,j);
              }
              *m_env.subDisplayFile() << "diagMult = " << diagMult
                                      << std::endl;
            }
            else {
              *m_env.subDisplayFile() << "attemptedMatrix = " << attemptedMatrix // FIX ME: might demand parallelism
                                      << std::endl;
            }
          }
          if (iRC) {
            UQ_FATAL_TEST_MACRO(iRC != UQ_MATRIX_IS_NOT_POS_DEFINITE_RC,
                                m_env.worldRank(),
                                "MetropolisHastingsSG<P_V,P_M>::generateFullChain()",
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
          ScaledCovMatrixTKGroup<P_V,P_M>* tempTK = dynamic_cast<ScaledCovMatrixTKGroup<P_V,P_M>* >(m_tk);
          P_M tmpMatrix(m_optionsObj->m_ov.m_amEta*attemptedMatrix);
          if (m_numDisabledParameters > 0) { // gpmsa2
            for (unsigned int paramId = 0; paramId < m_vectorSpace.dimLocal(); ++paramId) {
              if (m_parameterEnabledStatus[paramId] == false) {
                for (unsigned int i = 0; i < m_vectorSpace.dimLocal(); ++i) {
                  tmpMatrix(i,paramId) = 0.;
                }
                for (unsigned int j = 0; j < m_vectorSpace.dimLocal(); ++j) {
                  tmpMatrix(paramId,j) = 0.;
                }
                tmpMatrix(paramId,paramId) = 1.;
              }
            }
          }
          tempTK->updateLawCovMatrix(tmpMatrix);

#ifdef UQ_DRAM_MCG_REQUIRES_INVERTED_COV_MATRICES
          UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                            m_env.worldRank(),
                            "MetropolisHastingsSG<P_V,P_M>::generateFullChain()",
                            "need to code the update of m_upperCholProposalPrecMatrices");
#endif
        }

        //for (unsigned int i = 0; i < partialChain.subSequenceSize(); ++i) {
        //  if (partialChain[i]) delete partialChain[i];
        //}
      } // if (partialChain.subSequenceSize() > 0)

      if (m_optionsObj->m_ov.m_rawChainMeasureRunTimes) m_rawChainInfo.amRunTime += MiscGetEllapsedSeconds(&timevalAM);
    } // End of 'adaptive Metropolis' logic

    //****************************************************
    // Point 6/6 of logic for new position
    // Loop: print some information before going to the next chain position
    //****************************************************
    if ((m_env.subDisplayFile()                   ) &&
        (m_env.displayVerbosity() >= 3            ) &&
        (m_optionsObj->m_ov.m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                              << ": finishing chain position of id = " << positionId
                              << ", accept = "                         << accept
                              << ", curLogTarget  = "                  << currentPositionData.logTarget()
                              << ", canLogTarget  = "                  << currentCandidateData.logTarget()
                              << std::endl;
    }

    if ((m_optionsObj->m_ov.m_rawChainDisplayPeriod                     > 0) &&
        (((positionId+1) % m_optionsObj->m_ov.m_rawChainDisplayPeriod) == 0)) {
      if ((m_env.subDisplayFile()                   ) &&
          (m_optionsObj->m_ov.m_totallyMute == false)) {
        *m_env.subDisplayFile() << "Finished generating " << positionId+1
                                << " positions"  // root
                                << ", curret rejection percentage = " << (100. * ((double) m_rawChainInfo.numRejections)/((double) (positionId+1)))
                                << " %"
                                << std::endl;
      }
    }

    if ((m_env.subDisplayFile()                   ) &&
        (m_env.displayVerbosity() >= 10           ) &&
        (m_optionsObj->m_ov.m_totallyMute == false)) {
      *m_env.subDisplayFile() << "\n"
                              << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                              << "\n"
                              << std::endl;
    }
  } // end chain loop [for (unsigned int positionId = 1; positionId < workingChain.subSequenceSize(); ++positionId) {]

  if ((m_env.numSubEnvironments() < (unsigned int) m_env.fullComm().NumProc()) &&
      (m_initialPosition.numOfProcsForStorage() == 1                         ) &&
      (m_env.subRank()                          == 0                         )) {
    // subRank == 0 --> Tell all other processors to exit barrier now that the chain has been fully generated
    double aux = 0.;
    aux = m_targetPdfSynchronizer->callFunction(NULL,
                                                NULL,
                                                NULL,
                                                NULL,
                                                NULL,
                                                NULL,
                                                NULL);
    if (aux) {}; // just to remove compiler warning
  }

  //****************************************************
  // Print basic information about the chain
  //****************************************************
  m_rawChainInfo.runTime += MiscGetEllapsedSeconds(&timevalChain);
  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_ov.m_totallyMute == false)) {
    *m_env.subDisplayFile() << "Finished the generation of Markov chain " << workingChain.name()
                            << ", with sub "                              << workingChain.subSequenceSize()
                            << " positions";
    *m_env.subDisplayFile() << "\nSome information about this chain:"
                            << "\n  Chain run time       = " << m_rawChainInfo.runTime
                            << " seconds";
    if (m_optionsObj->m_ov.m_rawChainMeasureRunTimes) {
      *m_env.subDisplayFile() << "\n\n Breaking of the chain run time:\n";
      *m_env.subDisplayFile() << "\n  Candidate run time   = " << m_rawChainInfo.candidateRunTime
                              << " seconds ("                  << 100.*m_rawChainInfo.candidateRunTime/m_rawChainInfo.runTime
                              << "%)";
      *m_env.subDisplayFile() << "\n  Num target calls  = "    << m_rawChainInfo.numTargetCalls;
      *m_env.subDisplayFile() << "\n  Target d. run time   = " << m_rawChainInfo.targetRunTime
                              << " seconds ("                  << 100.*m_rawChainInfo.targetRunTime/m_rawChainInfo.runTime
                              << "%)";
      *m_env.subDisplayFile() << "\n  Avg target run time   = " << m_rawChainInfo.targetRunTime/((double) m_rawChainInfo.numTargetCalls)
                              << " seconds";
      *m_env.subDisplayFile() << "\n  Mh alpha run time    = " << m_rawChainInfo.mhAlphaRunTime
                              << " seconds ("                  << 100.*m_rawChainInfo.mhAlphaRunTime/m_rawChainInfo.runTime
                              << "%)";
      *m_env.subDisplayFile() << "\n  Dr alpha run time    = " << m_rawChainInfo.drAlphaRunTime
                              << " seconds ("                  << 100.*m_rawChainInfo.drAlphaRunTime/m_rawChainInfo.runTime
                              << "%)";
      *m_env.subDisplayFile() << "\n----------------------   --------------";
      double sumRunTime = m_rawChainInfo.candidateRunTime + m_rawChainInfo.targetRunTime + m_rawChainInfo.mhAlphaRunTime + m_rawChainInfo.drAlphaRunTime;
      *m_env.subDisplayFile() << "\n  Sum                  = " << sumRunTime
                              << " seconds ("                  << 100.*sumRunTime/m_rawChainInfo.runTime
                              << "%)";
      *m_env.subDisplayFile() << "\n\n Other run times:";
      *m_env.subDisplayFile() << "\n  DR run time          = " << m_rawChainInfo.drRunTime
                              << " seconds ("                  << 100.*m_rawChainInfo.drRunTime/m_rawChainInfo.runTime
                              << "%)";
      *m_env.subDisplayFile() << "\n  AM run time          = " << m_rawChainInfo.amRunTime
                              << " seconds ("                  << 100.*m_rawChainInfo.amRunTime/m_rawChainInfo.runTime
                              << "%)";
    }
    *m_env.subDisplayFile() << "\n  Number of DRs = "  << m_rawChainInfo.numDRs << "(num_DRs/chain_size = " << (double) m_rawChainInfo.numDRs/(double) workingChain.subSequenceSize()
                            << ")";
    *m_env.subDisplayFile() << "\n  Out of target support in DR = " << m_rawChainInfo.numOutOfTargetSupportInDR;
    *m_env.subDisplayFile() << "\n  Rejection percentage = "        << 100. * (double) m_rawChainInfo.numRejections/(double) workingChain.subSequenceSize()
                            << " %";
    *m_env.subDisplayFile() << "\n  Out of target support percentage = " << 100. * (double) m_rawChainInfo.numOutOfTargetSupport/(double) workingChain.subSequenceSize()
                            << " %";
    *m_env.subDisplayFile() << std::endl;
  }

  //****************************************************
  // Release memory before leaving routine
  //****************************************************
  //m_env.syncPrintDebugMsg("Leaving MetropolisHastingsSG<P_V,P_M>::generateFullChain()",3,3000000,m_env.fullComm()); // Dangerous to barrier on fullComm ... // KAUST

  return;
}
//--------------------------------------------------
template <class P_V,class P_M>
void
MetropolisHastingsSG<P_V,P_M>::updateAdaptedCovMatrix(
  const BaseVectorSequence<P_V,P_M>& partialChain,
  unsigned int                              idOfFirstPositionInSubChain,
  double&                                   lastChainSize,
  P_V&                                      lastMean,
  P_M&                                      lastAdaptedCovMatrix)
{
  double doubleSubChainSize = (double) partialChain.subSequenceSize();
  if (lastChainSize == 0) {
    UQ_FATAL_TEST_MACRO(partialChain.subSequenceSize() < 2,
                        m_env.worldRank(),
                        "MetropolisHastingsSG<P_V,P_M>::updateAdaptedCovMatrix()",
                        "'partialChain.subSequenceSize()' should be >= 2");

#if 1 // prudenci-2012-07-06
    lastMean = partialChain.subMeanPlain();
#else
    partialChain.subMeanExtra(0,partialChain.subSequenceSize(),lastMean);
#endif

    P_V tmpVec(m_vectorSpace.zeroVector());
    lastAdaptedCovMatrix = -doubleSubChainSize * matrixProduct(lastMean,lastMean);
    for (unsigned int i = 0; i < partialChain.subSequenceSize(); ++i) {
      partialChain.getPositionValues(i,tmpVec);
      lastAdaptedCovMatrix += matrixProduct(tmpVec,tmpVec);
    }
    lastAdaptedCovMatrix /= (doubleSubChainSize - 1.); // That is why partialChain size must be >= 2
  }
  else {
    UQ_FATAL_TEST_MACRO(partialChain.subSequenceSize() < 1,
                        m_env.worldRank(),
                        "MetropolisHastingsSG<P_V,P_M>::updateAdaptedCovMatrix()",
                        "'partialChain.subSequenceSize()' should be >= 1");

    UQ_FATAL_TEST_MACRO(idOfFirstPositionInSubChain < 1,
                        m_env.worldRank(),
                        "MetropolisHastingsSG<P_V,P_M>::updateAdaptedCovMatrix()",
                        "'idOfFirstPositionInSubChain' should be >= 1");

    P_V tmpVec (m_vectorSpace.zeroVector());
    P_V diffVec(m_vectorSpace.zeroVector());
    for (unsigned int i = 0; i < partialChain.subSequenceSize(); ++i) {
      double doubleCurrentId  = (double) (idOfFirstPositionInSubChain+i);
      partialChain.getPositionValues(i,tmpVec);
      diffVec = tmpVec - lastMean;

      double ratio1         = (1. - 1./doubleCurrentId); // That is why idOfFirstPositionInSubChain must be >= 1
      double ratio2         = (1./(1.+doubleCurrentId));
      lastAdaptedCovMatrix  = ratio1 * lastAdaptedCovMatrix + ratio2 * matrixProduct(diffVec,diffVec);
      lastMean             += ratio2 * diffVec;
    }
  }
  lastChainSize += doubleSubChainSize;

  if (m_numDisabledParameters > 0) { // gpmsa2
    for (unsigned int paramId = 0; paramId < m_vectorSpace.dimLocal(); ++paramId) {
      if (m_parameterEnabledStatus[paramId] == false) {
        for (unsigned int i = 0; i < m_vectorSpace.dimLocal(); ++i) {
          lastAdaptedCovMatrix(i,paramId) = 0.;
        }
        for (unsigned int j = 0; j < m_vectorSpace.dimLocal(); ++j) {
          lastAdaptedCovMatrix(paramId,j) = 0.;
        }
        lastAdaptedCovMatrix(paramId,paramId) = 1.;
      }
    }
  }

  return;
}

}  // End namespace QUESO

template class QUESO::MetropolisHastingsSG<QUESO::GslVector, QUESO::GslMatrix>;
