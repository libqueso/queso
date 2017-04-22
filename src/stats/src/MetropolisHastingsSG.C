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

#include <queso/asserts.h>
#include <queso/MetropolisHastingsSG.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

#include <queso/HessianCovMatricesTKGroup.h>
#include <queso/ScaledCovMatrixTKGroup.h>
#include <queso/TransformedScaledCovMatrixTKGroup.h>

#include <queso/InvLogitGaussianJointPdf.h>

#include <queso/GaussianJointPdf.h>

#include <queso/TKFactoryInitializer.h>
#include <queso/TransitionKernelFactory.h>
#include <queso/AlgorithmFactoryInitializer.h>
#include <queso/AlgorithmFactory.h>

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
MHRawChainInfoStruct::mpiSum(const MpiComm& comm, MHRawChainInfoStruct& sumInfo)
{
  comm.Allreduce<double>(&runTime, &sumInfo.runTime, (int) 7, RawValue_MPI_SUM,
                 "MHRawChainInfoStruct::mpiSum()",
                 "failed MPI.Allreduce() for sum of doubles");

  comm.Allreduce<unsigned int>(&numTargetCalls, &sumInfo.numTargetCalls, (int) 5, RawValue_MPI_SUM,
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
  m_tk                        (),
  m_algorithm                 (),
  m_positionIdForDebugging    (0),
  m_stageIdForDebugging       (0),
  m_idsOfUniquePositions      (0),//0.),
  m_logTargets                (0),//0.),
  m_alphaQuotients            (0),//0.),
  m_lastChainSize             (0),
  m_lastMean                  (),
  m_lastAdaptedCovMatrix      (),
  m_numPositionsNotSubWritten (0),
  m_optionsObj                (),
  m_computeInitialPriorAndLikelihoodValues(true),
  m_initialLogPriorValue      (0.),
  m_initialLogLikelihoodValue (0.),
  m_userDidNotProvideOptions(false),
  m_latestDirtyCovMatrixIteration(0)
{
  if (inputProposalCovMatrix != NULL) {
    m_initialProposalCovMatrix = *inputProposalCovMatrix;
  }

  // If user provided options, copy their object
  if (alternativeOptionsValues != NULL) {
    m_optionsObj.reset(new MhOptionsValues(*alternativeOptionsValues));
  }
  else {
    // Otherwise, we create one with the default
    m_optionsObj.reset(new MhOptionsValues(&m_env, prefix));
  }

  if (m_optionsObj->m_help != "") {
    if (m_env.subDisplayFile() && !m_optionsObj->m_totallyMute) {
      *m_env.subDisplayFile() << (*m_optionsObj) << std::endl;
    }
  }

  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_totallyMute == false)) {
    *m_env.subDisplayFile() << "Entering MetropolisHastingsSG<P_V,P_M>::constructor(1)"
                            << ": prefix = " << prefix
                            << ", alternativeOptionsValues = " << alternativeOptionsValues
                            << ", m_env.optionsInputFileName() = " << m_env.optionsInputFileName()
                            << ", m_initialProposalCovMatrix = " << m_initialProposalCovMatrix
                            << std::endl;
  }

  queso_require_equal_to_msg(sourceRv.imageSet().vectorSpace().dimLocal(), initialPosition.sizeLocal(), "'sourceRv' and 'initialPosition' should have equal dimensions");

  if (inputProposalCovMatrix) {
    queso_require_equal_to_msg(sourceRv.imageSet().vectorSpace().dimLocal(), inputProposalCovMatrix->numRowsLocal(), "'sourceRv' and 'inputProposalCovMatrix' should have equal dimensions");
    queso_require_equal_to_msg(inputProposalCovMatrix->numCols(), inputProposalCovMatrix->numRowsGlobal(), "'inputProposalCovMatrix' should be a square matrix");
  }

  commonConstructor();

  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_totallyMute == false)) {
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
  m_tk                        (),
  m_algorithm                 (),
  m_positionIdForDebugging    (0),
  m_stageIdForDebugging       (0),
  m_idsOfUniquePositions      (0),//0.),
  m_logTargets                (0),//0.),
  m_alphaQuotients            (0),//0.),
  m_lastChainSize             (0),
  m_lastMean                  (),
  m_lastAdaptedCovMatrix      (),
  m_numPositionsNotSubWritten (0),
  m_optionsObj                (),
  m_computeInitialPriorAndLikelihoodValues(false),
  m_initialLogPriorValue      (initialLogPrior),
  m_initialLogLikelihoodValue (initialLogLikelihood),
  m_userDidNotProvideOptions(false),
  m_latestDirtyCovMatrixIteration(0)
{
  if (inputProposalCovMatrix != NULL) {
    m_initialProposalCovMatrix = *inputProposalCovMatrix;
  }

  // If user provided options, copy their object
  if (alternativeOptionsValues != NULL) {
    m_optionsObj.reset(new MhOptionsValues(*alternativeOptionsValues));
  }
  else {
    // Otherwise we create one
    m_optionsObj.reset(new MhOptionsValues(&m_env, prefix));
  }

  if (m_optionsObj->m_help != "") {
    if (m_env.subDisplayFile() && !m_optionsObj->m_totallyMute) {
      *m_env.subDisplayFile() << (*m_optionsObj) << std::endl;
    }
  }

  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_totallyMute == false)) {
    *m_env.subDisplayFile() << "Entering MetropolisHastingsSG<P_V,P_M>::constructor(2)"
                            << ": prefix = " << prefix
                            << ", alternativeOptionsValues = " << alternativeOptionsValues
                            << ", m_env.optionsInputFileName() = " << m_env.optionsInputFileName()
                            << ", m_initialProposalCovMatrix = " << m_initialProposalCovMatrix
                            << std::endl;
  }

  queso_require_equal_to_msg(sourceRv.imageSet().vectorSpace().dimLocal(), initialPosition.sizeLocal(), "'sourceRv' and 'initialPosition' should have equal dimensions");

  if (inputProposalCovMatrix) {
    queso_require_equal_to_msg(sourceRv.imageSet().vectorSpace().dimLocal(), inputProposalCovMatrix->numRowsLocal(), "'sourceRv' and 'inputProposalCovMatrix' should have equal dimensions");
    queso_require_equal_to_msg(inputProposalCovMatrix->numCols(), inputProposalCovMatrix->numRowsGlobal(), "'inputProposalCovMatrix' should be a square matrix");
  }

  commonConstructor();

  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_totallyMute == false)) {
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
  m_tk                        (),
  m_algorithm                 (),
  m_positionIdForDebugging    (0),
  m_stageIdForDebugging       (0),
  m_idsOfUniquePositions      (0),//0.),
  m_logTargets                (0),//0.),
  m_alphaQuotients            (0),//0.),
  m_lastChainSize             (0),
  m_lastMean                  (),
  m_lastAdaptedCovMatrix      (),
  m_computeInitialPriorAndLikelihoodValues(true),
  m_initialLogPriorValue      (0.),
  m_initialLogLikelihoodValue (0.),
  m_userDidNotProvideOptions(true),
  m_latestDirtyCovMatrixIteration(0)
{
  // We do this dance because one of the MetropolisHastingsSGOptions
  // constructors takes one of the old-style MLSamplingLevelOptions options
  // objects.
  //
  // We do a copy and then pull out the raw values from m_ov.  We also need it
  // as a member (m_oldOptions) because otherwise m_ov will die when the
  // MetropolisHastingsSGOptions instance dies.
  m_oldOptions.reset(new MetropolisHastingsSGOptions(mlOptions));
  m_optionsObj.reset(new MhOptionsValues(m_oldOptions->m_ov));

  if (inputProposalCovMatrix != NULL) {
    m_initialProposalCovMatrix = *inputProposalCovMatrix;
    if ((m_env.subDisplayFile()                   ) &&
        (m_optionsObj->m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::constructor(3)"
                              << ": just set m_initialProposalCovMatrix = " << m_initialProposalCovMatrix
                              << std::endl;
    }
  }
  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_totallyMute == false)) {
    *m_env.subDisplayFile() << "Entering MetropolisHastingsSG<P_V,P_M>::constructor(3)"
                            << std::endl;
  }

  commonConstructor();

  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_totallyMute == false)) {
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
  m_tk                        (),
  m_algorithm                 (),
  m_positionIdForDebugging    (0),
  m_stageIdForDebugging       (0),
  m_idsOfUniquePositions      (0),//0.),
  m_logTargets                (0),//0.),
  m_alphaQuotients            (0),//0.),
  m_lastChainSize             (0),
  m_lastMean                  (),
  m_lastAdaptedCovMatrix      (),
  m_computeInitialPriorAndLikelihoodValues(false),
  m_initialLogPriorValue      (initialLogPrior),
  m_initialLogLikelihoodValue (initialLogLikelihood),
  m_userDidNotProvideOptions(true),
  m_latestDirtyCovMatrixIteration(0)
{
  // We do this dance because one of the MetropolisHastingsSGOptions
  // constructors takes one of the old-style MLSamplingLevelOptions options
  // objects.
  //
  // We do a copy and then pull out the raw values from m_ov.  We also need it
  // as a member (m_oldOptions) because otherwise m_ov will die when the
  // MetropolisHastingsSGOptions instance dies.
  m_oldOptions.reset(new MetropolisHastingsSGOptions(mlOptions));
  m_optionsObj.reset(new MhOptionsValues(m_oldOptions->m_ov));

  if (inputProposalCovMatrix != NULL) {
    m_initialProposalCovMatrix = *inputProposalCovMatrix;
    if ((m_env.subDisplayFile()                   ) &&
        (m_optionsObj->m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::constructor(4)"
                              << ": just set m_initialProposalCovMatrix = " << m_initialProposalCovMatrix
                              << std::endl;
    }
  }
  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_totallyMute == false)) {
    *m_env.subDisplayFile() << "Entering MetropolisHastingsSG<P_V,P_M>::constructor(4)"
                            << std::endl;
  }

  commonConstructor();

  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_totallyMute == false)) {
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

  m_lastChainSize             = 0;
  m_rawChainInfo.reset();
  m_alphaQuotients.clear();
  m_logTargets.clear();
  m_numDisabledParameters  = 0; // gpmsa2
  m_parameterEnabledStatus.clear(); // gpmsa2
  m_positionIdForDebugging = 0;
  m_stageIdForDebugging    = 0;
  m_idsOfUniquePositions.clear();

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
      (m_optionsObj->m_totallyMute == false)) {
    *m_env.subDisplayFile() << "Entering MetropolisHastingsSG<P_V,P_M>::commonConstructor()"
                            << std::endl;
  }

  if (m_optionsObj->m_initialPositionDataInputFileName != ".") { // palms
    std::set<unsigned int> tmpSet;
    tmpSet.insert(m_env.subId());
    m_initialPosition.subReadContents((m_optionsObj->m_initialPositionDataInputFileName+"_sub"+m_env.subIdString()),
                                      m_optionsObj->m_initialPositionDataInputFileType,
                                      tmpSet);
    if ((m_env.subDisplayFile()                   ) &&
        (m_optionsObj->m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::commonConstructor()"
                              << ": just read initial position contents = " << m_initialPosition
                              << std::endl;
    }
  }

  if (m_optionsObj->m_parameterDisabledSet.size() > 0) { // gpmsa2
    for (std::set<unsigned int>::iterator setIt = m_optionsObj->m_parameterDisabledSet.begin(); setIt != m_optionsObj->m_parameterDisabledSet.end(); ++setIt) {
      unsigned int paramId = *setIt;
      if (paramId < m_vectorSpace.dimLocal()) {
        m_numDisabledParameters++;
        m_parameterEnabledStatus[paramId] = false;
      }
    }
  }

  std::vector<double> drScalesAll(m_optionsObj->m_drScalesForExtraStages.size()+1,1.);
  for (unsigned int i = 1; i < (m_optionsObj->m_drScalesForExtraStages.size()+1); ++i) {
    drScalesAll[i] = m_optionsObj->m_drScalesForExtraStages[i-1];
  }

  // Deprecate m_doLogitTransform
  if (m_optionsObj->m_doLogitTransform != UQ_MH_SG_DO_LOGIT_TRANSFORM) {
    std::string msg;
    msg = "The doLogitTransform option is deprecated.  ";
    msg += "Use both ip_mh_algorithm and ip_mh_tk instead.";
    queso_warning(msg.c_str());
  }

  if (m_optionsObj->m_initialProposalCovMatrixDataInputFileName != ".") { // palms
    std::set<unsigned int> tmpSet;
    tmpSet.insert(m_env.subId());
    m_initialProposalCovMatrix.subReadContents((m_optionsObj->m_initialProposalCovMatrixDataInputFileName+"_sub"+m_env.subIdString()),
                                               m_optionsObj->m_initialProposalCovMatrixDataInputFileType,
                                               tmpSet);
    if ((m_env.subDisplayFile()                   ) &&
        (m_optionsObj->m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::commonConstructor()"
                              << ": just read initial proposal cov matrix contents = " << m_initialProposalCovMatrix
                              << std::endl;
    }
  }
  else {
    queso_require_msg(!(m_nullInputProposalCovMatrix), "proposal cov matrix should have been passed by user, since, according to the input algorithm options, local Hessians will not be used in the proposal");
  }

  // Only transform prop cov matrix if we're doing a logit random walk.
  // Also note we're transforming *after* we potentially read it from the input
  // file.
  if ((m_optionsObj->m_algorithm == "logit_random_walk") &&
      (m_optionsObj->m_tk        == "logit_random_walk") &&
      m_optionsObj->m_doLogitTransform) {
    // Variable transform initial proposal cov matrix
    transformInitialCovMatrixToGaussianSpace(
        dynamic_cast<const BoxSubset<P_V, P_M> & >(m_targetPdf.domainSet()));
  }

  // This instantiates all the transition kernels with their associated
  // factories
  TKFactoryInitializer tk_factory_initializer;

  TransitionKernelFactory::set_vectorspace(m_vectorSpace);
  TransitionKernelFactory::set_options(*m_optionsObj);
  TransitionKernelFactory::set_pdf_synchronizer(*m_targetPdfSynchronizer);
  TransitionKernelFactory::set_initial_cov_matrix(m_initialProposalCovMatrix);
  TransitionKernelFactory::set_dr_scales(drScalesAll);
  TransitionKernelFactory::set_target_pdf(m_targetPdf);
  m_tk = TransitionKernelFactory::build(m_optionsObj->m_tk);

  // This instantiates all the algorithms with their associated factories
  AlgorithmFactoryInitializer algorithm_factory_initializer;

  AlgorithmFactory::set_environment(m_env);
  AlgorithmFactory::set_tk(*m_tk);
  m_algorithm = AlgorithmFactory::build(m_optionsObj->m_algorithm);
}
//--------------------------------------------------
template<class P_V,class P_M>
double
MetropolisHastingsSG<P_V,P_M>::alpha(
  const std::vector<MarkovChainPositionData<P_V>*>& inputPositionsData,
  const std::vector<unsigned int                 >& inputTKStageIds)
{
  // inputPositionsData is all the DR position data, except for the first
  // two elements.  The first element is the current state, and the second
  // element is the would-be candidate before DR.
  //
  // inputTKStageIds is a vector containing 0, 1, 2, ..., n

  unsigned int inputSize = inputPositionsData.size();
  if ((m_env.subDisplayFile()                   ) &&
      (m_env.displayVerbosity() >= 10           ) &&
      (m_optionsObj->m_totallyMute == false)) {
    *m_env.subDisplayFile() << "Entering MetropolisHastingsSG<P_V,P_M>::alpha(vec)"
                           << ", inputSize = " << inputSize
                           << std::endl;
  }
  queso_require_greater_equal_msg(inputSize, 2, "inputPositionsData has size < 2");
  queso_require_equal_to_msg(inputSize, inputPositionsData.size(), "inputPositionsData and inputTKStageIds have different lengths");

  // If necessary, return 0. right away
  if (inputPositionsData[0          ]->outOfTargetSupport()) return 0.;
  if (inputPositionsData[inputSize-1]->outOfTargetSupport()) return 0.;

  if ((inputPositionsData[0]->logTarget() == -INFINITY           ) ||
      (inputPositionsData[0]->logTarget() ==  INFINITY           ) ||
      ( queso_isnan(inputPositionsData[0]->logTarget()) )) {
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
           ( queso_isnan(inputPositionsData[inputSize - 1]->logTarget()) )) {
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
  if (inputSize == 2) {
    const P_V & tk_pos_x = m_tk->preComputingPosition(inputTKStageIds[inputSize-1]);
    const P_V & tk_pos_y = m_tk->preComputingPosition(inputTKStageIds[0]);
    return this->m_algorithm->acceptance_ratio(
        *(inputPositionsData[0]),
        *(inputPositionsData[inputSize - 1]),
        tk_pos_x,
        tk_pos_y);
  }

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

  double numContrib = m_tk->rv(backwardTKStageIdsLess1).pdf().lnValue(_lastBackwardTKPosition);

  double denContrib = m_tk->rv(tkStageIdsLess1).pdf().lnValue(_lastTKPosition);
  if ((m_env.subDisplayFile()                   ) &&
      (m_env.displayVerbosity() >= 10           ) &&
      (m_optionsObj->m_totallyMute == false)) {
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

    numContrib = m_tk->rv(backwardTKStageIdsLess1).pdf().lnValue(lastBackwardTKPosition);

    denContrib = m_tk->rv(tkStageIdsLess1).pdf().lnValue(lastTKPosition);
    if ((m_env.subDisplayFile()                   ) &&
        (m_env.displayVerbosity() >= 10           ) &&
        (m_optionsObj->m_totallyMute == false)) {
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
      (m_optionsObj->m_totallyMute == false)) {
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
      (m_optionsObj->m_totallyMute == false)) {
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
      (m_optionsObj->m_totallyMute == false)) {
    *m_env.subDisplayFile() << "\n"
                            << "\n-----------------------------------------------------"
                            << "\n Writing more information about the Markov chain " << workingChain.name() << " to output file ..."
                            << "\n-----------------------------------------------------"
                            << "\n"
                            << std::endl;
  }

  int iRC = UQ_OK_RC;

  if (m_optionsObj->m_rawChainGenerateExtra) {
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
      (m_optionsObj->m_totallyMute == false)) {
    *m_env.subDisplayFile() << "\n-----------------------------------------------------"
                            << "\n Finished writing more information about the Markov chain " << workingChain.name()
                            << "\n-----------------------------------------------------"
                            << "\n"
                            << std::endl;
  }

  return iRC;
}

template <class P_V, class P_M>
const BaseTKGroup<P_V, P_M> &
MetropolisHastingsSG<P_V, P_M>::transitionKernel() const
{
  return *m_tk;
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
      (m_optionsObj->m_totallyMute == false)) {
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
      (m_optionsObj->m_outputLogLikelihood)) {
    writeLogLikelihood = true;
  }
  else {
    writeLogLikelihood = false;
  }

  // Set a flag to write out log target or not
  bool writeLogTarget;
  if ((workingLogTargetValues != NULL) &&
      (m_optionsObj->m_outputLogTarget)) {
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
  if (m_optionsObj->m_rawChainDataInputFileName == UQ_MH_SG_FILENAME_FOR_NO_FILE) {
    generateFullChain(valuesOf1stPosition,
                      m_optionsObj->m_rawChainSize,
                      workingChain,
                      workingLogLikelihoodValues,
                      workingLogTargetValues);
  }
  else {
    readFullChain(m_optionsObj->m_rawChainDataInputFileName,
                  m_optionsObj->m_rawChainDataInputFileType,
                  m_optionsObj->m_rawChainSize,
                  workingChain);
  }

  //****************************************************
  // Open generic output file
  //****************************************************
  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_totallyMute == false)) {
    *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateSequence()"
                            << ", prefix = "                                         << m_optionsObj->m_prefix
                            << ", chain name = "                                     << workingChain.name()
                            << ": about to try to open generic output file '"        << m_optionsObj->m_dataOutputFileName
                            << "."                                                   << UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT // Yes, always ".m"
                            << "', subId = "                                         << m_env.subId()
                            << ", subenv is allowed to write (1/true or 0/false) = " << (m_optionsObj->m_dataOutputAllowedSet.find(m_env.subId()) != m_optionsObj->m_dataOutputAllowedSet.end())
                            << "..."
                            << std::endl;
  }

  FilePtrSetStruct genericFilePtrSet;
  m_env.openOutputFile(m_optionsObj->m_dataOutputFileName,
                       UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT, // Yes, always ".m"
                       m_optionsObj->m_dataOutputAllowedSet,
                       false,
                       genericFilePtrSet);

  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_totallyMute == false)) {
    *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateSequence()"
                            << ", prefix = "                                   << m_optionsObj->m_prefix
                            << ", raw chain name = "                           << workingChain.name()
                            << ": returned from opening generic output file '" << m_optionsObj->m_dataOutputFileName
                            << "."                                             << UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT // Yes, always ".m"
                            << "', subId = "                                   << m_env.subId()
                            << std::endl;
  }

  //****************************************************************************************
  // Eventually:
  // --> write raw chain
  // --> compute statistics on it
  //****************************************************************************************
  if ((m_optionsObj->m_rawChainDataOutputFileName != UQ_MH_SG_FILENAME_FOR_NO_FILE) &&
      (m_optionsObj->m_totallyMute == false                                       )) {

    // Take "sub" care of raw chain
    if ((m_env.subDisplayFile()                   ) &&
        (m_optionsObj->m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateSequence()"
                              << ", prefix = "                                         << m_optionsObj->m_prefix
                              << ", raw chain name = "                                 << workingChain.name()
                              << ": about to try to write raw sub chain output file '" << m_optionsObj->m_rawChainDataOutputFileName
                              << "."                                                   << m_optionsObj->m_rawChainDataOutputFileType
                              << "', subId = "                                         << m_env.subId()
                              << ", subenv is allowed to write  1/true or 0/false) = " << (m_optionsObj->m_rawChainDataOutputAllowedSet.find(m_env.subId()) != m_optionsObj->m_rawChainDataOutputAllowedSet.end())
                              << "..."
                              << std::endl;
    }

    if ((m_numPositionsNotSubWritten                     >  0  ) &&
        (m_optionsObj->m_rawChainDataOutputFileName != ".")) {
      workingChain.subWriteContents(m_optionsObj->m_rawChainSize - m_numPositionsNotSubWritten,
                                    m_numPositionsNotSubWritten,
                                    m_optionsObj->m_rawChainDataOutputFileName,
                                    m_optionsObj->m_rawChainDataOutputFileType,
                                    m_optionsObj->m_rawChainDataOutputAllowedSet);
      if ((m_env.subDisplayFile()                   ) &&
          (m_optionsObj->m_totallyMute == false)) {
        *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateSequence()"
                                << ": just wrote (per period request) remaining " << m_numPositionsNotSubWritten << " chain positions "
                                << ", " << m_optionsObj->m_rawChainSize - m_numPositionsNotSubWritten << " <= pos <= " << m_optionsObj->m_rawChainSize - 1
                                << std::endl;
      }

      if (writeLogLikelihood) {
        workingLogLikelihoodValues->subWriteContents(m_optionsObj->m_rawChainSize - m_numPositionsNotSubWritten,
                                                     m_numPositionsNotSubWritten,
                                                     m_optionsObj->m_rawChainDataOutputFileName + "_loglikelihood",
                                                     m_optionsObj->m_rawChainDataOutputFileType,
                                                     m_optionsObj->m_rawChainDataOutputAllowedSet);
      }

      if (writeLogTarget) {
        workingLogTargetValues->subWriteContents(m_optionsObj->m_rawChainSize - m_numPositionsNotSubWritten,
                                                 m_numPositionsNotSubWritten,
                                                 m_optionsObj->m_rawChainDataOutputFileName + "_logtarget",
                                                 m_optionsObj->m_rawChainDataOutputFileType,
                                                 m_optionsObj->m_rawChainDataOutputAllowedSet);
      }

      m_numPositionsNotSubWritten = 0;
    }

    // Compute raw sub MLE
    if (workingLogLikelihoodValues) {
      SequenceOfVectors<P_V,P_M> rawSubMLEpositions(m_vectorSpace,0,m_optionsObj->m_prefix+"rawSubMLEseq");
      double rawSubMLEvalue = workingChain.subPositionsOfMaximum(*workingLogLikelihoodValues,
                                                                 rawSubMLEpositions);
      queso_require_not_equal_to_msg(rawSubMLEpositions.subSequenceSize(), 0, "rawSubMLEpositions.subSequenceSize() = 0");

      if ((m_env.subDisplayFile()                   ) &&
          (m_optionsObj->m_totallyMute == false)) {
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
      queso_require_not_equal_to_msg(rawSubMAPpositions.subSequenceSize(), 0, "rawSubMAPpositions.subSequenceSize() = 0");

      if ((m_env.subDisplayFile()                   ) &&
          (m_optionsObj->m_totallyMute == false)) {
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
        (m_optionsObj->m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateSequence()"
                              << ", prefix = "                                         << m_optionsObj->m_prefix
                              << ", raw chain name = "                                 << workingChain.name()
                              << ": returned from writing raw sub chain output file '" << m_optionsObj->m_rawChainDataOutputFileName
                              << "."                                                   << m_optionsObj->m_rawChainDataOutputFileType
                              << "', subId = "                                         << m_env.subId()
                              << std::endl;
    }

    // Take "unified" care of raw chain
    if ((m_env.subDisplayFile()                   ) &&
        (m_optionsObj->m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateSequence()"
                              << ", prefix = "                                             << m_optionsObj->m_prefix
                              << ", raw chain name = "                                     << workingChain.name()
                              << ": about to try to write raw unified chain output file '" << m_optionsObj->m_rawChainDataOutputFileName
                              << "."                                                       << m_optionsObj->m_rawChainDataOutputFileType
                              << "', subId = "                                             << m_env.subId()
                              << "..."
                              << std::endl;
    }

    workingChain.unifiedWriteContents(m_optionsObj->m_rawChainDataOutputFileName,
                                      m_optionsObj->m_rawChainDataOutputFileType);
    if ((m_env.subDisplayFile()                   ) &&
        (m_optionsObj->m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateSequence()"
                              << ", prefix = "                                             << m_optionsObj->m_prefix
                              << ", raw chain name = "                                     << workingChain.name()
                              << ": returned from writing raw unified chain output file '" << m_optionsObj->m_rawChainDataOutputFileName
                              << "."                                                       << m_optionsObj->m_rawChainDataOutputFileType
                              << "', subId = "                                             << m_env.subId()
                              << std::endl;
    }

    if (writeLogLikelihood) {
      workingLogLikelihoodValues->unifiedWriteContents(m_optionsObj->m_rawChainDataOutputFileName + "_loglikelihood",
                                                       m_optionsObj->m_rawChainDataOutputFileType);
    }

    if (writeLogTarget) {
      workingLogTargetValues->unifiedWriteContents(m_optionsObj->m_rawChainDataOutputFileName + "_logtarget",
                                                   m_optionsObj->m_rawChainDataOutputFileType);
    }

    // Compute raw unified MLE only in inter0Comm
    if (workingLogLikelihoodValues && (m_env.subRank() == 0)) {
      SequenceOfVectors<P_V,P_M> rawUnifiedMLEpositions(m_vectorSpace,0,m_optionsObj->m_prefix+"rawUnifiedMLEseq");

      double rawUnifiedMLEvalue = workingChain.unifiedPositionsOfMaximum(*workingLogLikelihoodValues,
                                                                         rawUnifiedMLEpositions);

      if ((m_env.subDisplayFile()                   ) &&
          (m_optionsObj->m_totallyMute == false)) {
        P_V tmpVec(m_vectorSpace.zeroVector());

        // Make sure the positions vector (which only contains stuff on process
        // zero) actually contains positions
        if (rawUnifiedMLEpositions.subSequenceSize() > 0) {
          rawUnifiedMLEpositions.getPositionValues(0,tmpVec);
          *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateSequence()"
                                  << ": just computed MLE"
                                  << ", rawUnifiedMLEvalue = "                       << rawUnifiedMLEvalue
                                  << ", rawUnifiedMLEpositions.subSequenceSize() = " << rawUnifiedMLEpositions.subSequenceSize()
                                  << ", rawUnifiedMLEpositions[0] = "                << tmpVec
                                  << std::endl;
        }
      }
    }

    // Compute raw unified MAP only in inter0Comm
    if (workingLogTargetValues && (m_env.subRank() == 0)) {
      SequenceOfVectors<P_V,P_M> rawUnifiedMAPpositions(m_vectorSpace,0,m_optionsObj->m_prefix+"rawUnifiedMAPseq");
      double rawUnifiedMAPvalue = workingChain.unifiedPositionsOfMaximum(*workingLogTargetValues,
                                                                         rawUnifiedMAPpositions);

      if ((m_env.subDisplayFile()                   ) &&
          (m_optionsObj->m_totallyMute == false)) {
        P_V tmpVec(m_vectorSpace.zeroVector());

        // Make sure the positions vector (which only contains stuff on process
        // zero) actually contains positions
        if (rawUnifiedMAPpositions.subSequenceSize() > 0) {
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
  }

  // Take care of other aspects of raw chain
  if ((genericFilePtrSet.ofsVar                 ) &&
      (m_optionsObj->m_totallyMute == false)) {
    // Write likelihoodValues and alphaValues, if they were requested by user
    iRC = writeInfo(workingChain,
                    *genericFilePtrSet.ofsVar);
    queso_require_msg(!(iRC), "improper writeInfo() return");
  }

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  if (m_optionsObj->m_rawChainComputeStats) {
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
  if (m_optionsObj->m_filteredChainGenerate) {
    // Compute filter parameters
    unsigned int filterInitialPos = (unsigned int) (m_optionsObj->m_filteredChainDiscardedPortion * (double) workingChain.subSequenceSize());
    unsigned int filterSpacing    = m_optionsObj->m_filteredChainLag;
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
        (m_optionsObj->m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateSequence()"
                              << ", prefix = "                                                      << m_optionsObj->m_prefix
                              << ": checking necessity of opening output files for filtered chain " << workingChain.name()
                              << "..."
                              << std::endl;
    }

    // Take "sub" care of filtered chain
    if ((m_optionsObj->m_filteredChainDataOutputFileName != UQ_MH_SG_FILENAME_FOR_NO_FILE) &&
        (m_optionsObj->m_totallyMute == false                                            )) {
      workingChain.subWriteContents(0,
                                    workingChain.subSequenceSize(),
                                    m_optionsObj->m_filteredChainDataOutputFileName,
                                    m_optionsObj->m_filteredChainDataOutputFileType,
                                    m_optionsObj->m_filteredChainDataOutputAllowedSet);
      if ((m_env.subDisplayFile()                   ) &&
          (m_optionsObj->m_totallyMute == false)) {
        *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateSequence()"
                                << ", prefix = "                << m_optionsObj->m_prefix
                                << ": closed sub output file '" << m_optionsObj->m_filteredChainDataOutputFileName
                                << "' for filtered chain "      << workingChain.name()
                                << std::endl;
      }

      if (writeLogLikelihood) {
        workingLogLikelihoodValues->subWriteContents(0,
                                                     workingChain.subSequenceSize(),
                                                     m_optionsObj->m_filteredChainDataOutputFileName + "_loglikelihood",
                                                     m_optionsObj->m_filteredChainDataOutputFileType,
                                                     m_optionsObj->m_filteredChainDataOutputAllowedSet);
      }

      if (writeLogTarget) {
        workingLogTargetValues->subWriteContents(0,
                                                 workingChain.subSequenceSize(),
                                                 m_optionsObj->m_filteredChainDataOutputFileName + "_logtarget",
                                                 m_optionsObj->m_filteredChainDataOutputFileType,
                                                 m_optionsObj->m_filteredChainDataOutputAllowedSet);
      }
    }

    // Compute sub filtered MLE and sub filtered MAP

    // Take "unified" care of filtered chain
    if ((m_optionsObj->m_filteredChainDataOutputFileName != UQ_MH_SG_FILENAME_FOR_NO_FILE) &&
        (m_optionsObj->m_totallyMute == false                                            )) {
      workingChain.unifiedWriteContents(m_optionsObj->m_filteredChainDataOutputFileName,
                                        m_optionsObj->m_filteredChainDataOutputFileType);
      if ((m_env.subDisplayFile()                   ) &&
          (m_optionsObj->m_totallyMute == false)) {
        *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateSequence()"
                                << ", prefix = "                    << m_optionsObj->m_prefix
                                << ": closed unified output file '" << m_optionsObj->m_filteredChainDataOutputFileName
                                << "' for filtered chain "          << workingChain.name()
                                << std::endl;
      }

      if (writeLogLikelihood) {
        workingLogLikelihoodValues->unifiedWriteContents(m_optionsObj->m_filteredChainDataOutputFileName + "_loglikelihood",
                                                         m_optionsObj->m_filteredChainDataOutputFileType);
      }

      if (writeLogTarget) {
        workingLogTargetValues->unifiedWriteContents(m_optionsObj->m_filteredChainDataOutputFileName + "_logtarget",
                                                     m_optionsObj->m_filteredChainDataOutputFileType);
      }
    }

    // Compute unified filtered MLE and unified filtered MAP

    // Compute statistics
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
    if (m_optionsObj->m_filteredChainComputeStats) {
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
        (m_optionsObj->m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateSequence()"
                              << ", prefix = "                    << m_optionsObj->m_prefix
                              << ": closed generic output file '" << m_optionsObj->m_dataOutputFileName
                              << "' (chain name is "              << workingChain.name()
                              << ")"
                              << std::endl;
    }
  }

  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_totallyMute == false)) {
    *m_env.subDisplayFile() << std::endl;
  }

  //m_env.syncPrintDebugMsg("Leaving MetropolisHastingsSG<P_V,P_M>::generateSequence()",2,3000000,m_env.fullComm()); // Dangerous to barrier on fullComm ... // KAUST

  if ((m_env.subDisplayFile()                   ) &&
      (m_env.displayVerbosity() >= 5            ) &&
      (m_optionsObj->m_totallyMute == false)) {
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
      (m_optionsObj->m_totallyMute == false)) {
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

  m_positionIdForDebugging = 0;
  m_stageIdForDebugging    = 0;

  m_rawChainInfo.reset();

  iRC = gettimeofday(&timevalChain, NULL);
  queso_require_equal_to_msg(iRC, 0, "gettimeofday called failed");

  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_totallyMute == false)) {
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
      (m_optionsObj->m_outputLogLikelihood)) {
    writeLogLikelihood = true;
  }
  else {
    writeLogLikelihood = false;
  }

  // Set a flag to write out log target or not
  bool writeLogTarget;
  if ((workingLogTargetValues != NULL) &&
      (m_optionsObj->m_outputLogTarget)) {
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
  queso_require_msg(!(outOfTargetSupport), "initial position should not be out of target pdf support");

  double logPrior      = 0.;
  double logLikelihood = 0.;
  double logTarget     = 0.;
  if (m_computeInitialPriorAndLikelihoodValues) {
    if (m_optionsObj->m_rawChainMeasureRunTimes) {
      iRC = gettimeofday(&timevalTarget, NULL);
      queso_require_equal_to_msg(iRC, 0, "gettimeofday called failed");
    }
    logTarget = m_targetPdfSynchronizer->callFunction(&valuesOf1stPosition,&logPrior,&logLikelihood); // Might demand parallel environment // KEY
    if (m_optionsObj->m_rawChainMeasureRunTimes) {
      m_rawChainInfo.targetRunTime += MiscGetEllapsedSeconds(&timevalTarget);
    }
    m_rawChainInfo.numTargetCalls++;
    if ((m_env.subDisplayFile()                   ) &&
        (m_env.displayVerbosity() >= 3            ) &&
        (m_optionsObj->m_totallyMute == false)) {
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
        (m_optionsObj->m_totallyMute == false)) {
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
  if (m_optionsObj->m_rawChainGenerateExtra) {
    m_logTargets.resize    (chainSize,0.);
    m_alphaQuotients.resize(chainSize,0.);
  }

  unsigned int uniquePos = 0;
  workingChain.setPositionValues(0,currentPositionData.vecValues());
  m_numPositionsNotSubWritten++;
  if ((m_optionsObj->m_rawChainDataOutputPeriod           >  0  ) &&
      (((0+1) % m_optionsObj->m_rawChainDataOutputPeriod) == 0  ) &&
      (m_optionsObj->m_rawChainDataOutputFileName         != ".")) {
    workingChain.subWriteContents(0 + 1 - m_optionsObj->m_rawChainDataOutputPeriod,
                                  m_optionsObj->m_rawChainDataOutputPeriod,
                                  m_optionsObj->m_rawChainDataOutputFileName,
                                  m_optionsObj->m_rawChainDataOutputFileType,
                                  m_optionsObj->m_rawChainDataOutputAllowedSet);
    if ((m_env.subDisplayFile()                   ) &&
        (m_optionsObj->m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                              << ": just wrote (per period request) " << m_numPositionsNotSubWritten << " chain positions "
                              << ", " << 0 + 1 - m_optionsObj->m_rawChainDataOutputPeriod << " <= pos <= " << 0
                              << std::endl;
    }

    if (writeLogLikelihood) {
      workingLogLikelihoodValues->subWriteContents(0 + 1 - m_optionsObj->m_rawChainDataOutputPeriod,
                                                   m_optionsObj->m_rawChainDataOutputPeriod,
                                                   m_optionsObj->m_rawChainDataOutputFileName + "_loglikelihood",
                                                   m_optionsObj->m_rawChainDataOutputFileType,
                                                   m_optionsObj->m_rawChainDataOutputAllowedSet);
    }

    if (writeLogTarget) {
      workingLogTargetValues->subWriteContents(0 + 1 - m_optionsObj->m_rawChainDataOutputPeriod,
                                               m_optionsObj->m_rawChainDataOutputPeriod,
                                               m_optionsObj->m_rawChainDataOutputFileName + "_logtarget",
                                               m_optionsObj->m_rawChainDataOutputFileType,
                                               m_optionsObj->m_rawChainDataOutputAllowedSet);
    }

    m_numPositionsNotSubWritten = 0;
  }

  if (workingLogLikelihoodValues) (*workingLogLikelihoodValues)[0] = currentPositionData.logLikelihood();
  if (workingLogTargetValues    ) (*workingLogTargetValues    )[0] = currentPositionData.logTarget();
  if (true/*m_uniqueChainGenerate*/) m_idsOfUniquePositions[uniquePos++] = 0;
  if (m_optionsObj->m_rawChainGenerateExtra) {
    m_logTargets    [0] = currentPositionData.logTarget();
    m_alphaQuotients[0] = 1.;
  }
  //*m_env.subDisplayFile() << "AQUI 002" << std::endl;

  if ((m_env.subDisplayFile()                   ) &&
      (m_env.displayVerbosity() >= 10           ) &&
      (m_optionsObj->m_totallyMute == false)) {
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
        (m_optionsObj->m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                              << ": beginning chain position of id = "  << positionId
                              << ", m_optionsObj->m_drMaxNumExtraStages = " << m_optionsObj->m_drMaxNumExtraStages
                              << std::endl;
    }
    unsigned int stageId  = 0;
    m_stageIdForDebugging = stageId;

    m_tk->clearPreComputingPositions();

    if ((m_env.subDisplayFile()                   ) &&
        (m_env.displayVerbosity() >= 5            ) &&
        (m_optionsObj->m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                              << ": about to set TK pre computing position of local id " << 0
                              << ", values = " << currentPositionData.vecValues()
                              << std::endl;
    }
    bool validPreComputingPosition = m_tk->setPreComputingPosition(currentPositionData.vecValues(),0);
    if ((m_env.subDisplayFile()                   ) &&
        (m_env.displayVerbosity() >= 5            ) &&
        (m_optionsObj->m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                              << ": returned from setting TK pre computing position of local id " << 0
                              << ", values = " << currentPositionData.vecValues()
                              << ", valid = "  << validPreComputingPosition
                              << std::endl;
    }
    queso_require_msg(validPreComputingPosition, "initial position should not be an invalid pre computing position");

    //****************************************************
    // Point 2/6 of logic for new position
    // Loop: generate new position
    //****************************************************
    bool keepGeneratingCandidates = true;
    while (keepGeneratingCandidates) {
      if (m_optionsObj->m_rawChainMeasureRunTimes) {
        iRC = gettimeofday(&timevalCandidate, NULL);
        queso_require_equal_to_msg(iRC, 0, "gettimeofday called failed");
      }

      m_tk->rv(currentPositionData.vecValues()).realizer().realization(tmpVecValues);

      if (m_numDisabledParameters > 0) { // gpmsa2
        for (unsigned int paramId = 0; paramId < m_vectorSpace.dimLocal(); ++paramId) {
          if (m_parameterEnabledStatus[paramId] == false) {
            tmpVecValues[paramId] = m_initialPosition[paramId];
          }
        }
      }
      if (m_optionsObj->m_rawChainMeasureRunTimes) m_rawChainInfo.candidateRunTime += MiscGetEllapsedSeconds(&timevalCandidate);

      outOfTargetSupport = !m_targetPdf.domainSet().contains(tmpVecValues);

      bool displayDetail = (m_env.displayVerbosity() >= 10/*99*/) || m_optionsObj->m_displayCandidates;
      if ((m_env.subDisplayFile()                   ) &&
          (displayDetail                            ) &&
          (m_optionsObj->m_totallyMute == false)) {
        *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                                << ": for chain position of id = " << positionId
                                << ", candidate = "                << tmpVecValues // FIX ME: might need parallelism
                                << ", outOfTargetSupport = "       << outOfTargetSupport
                                << std::endl;
      }

      if (m_optionsObj->m_putOutOfBoundsInChain) keepGeneratingCandidates = false;
      else                                            keepGeneratingCandidates = outOfTargetSupport;
    }

    if ((m_env.subDisplayFile()                   ) &&
        (m_env.displayVerbosity() >= 5            ) &&
        (m_optionsObj->m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                              << ": about to set TK pre computing position of local id " << stageId+1
                              << ", values = " << tmpVecValues
                              << std::endl;
    }
    validPreComputingPosition = m_tk->setPreComputingPosition(tmpVecValues,stageId+1);
    if ((m_env.subDisplayFile()                   ) &&
        (m_env.displayVerbosity() >= 5            ) &&
        (m_optionsObj->m_totallyMute == false)) {
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
      if (m_optionsObj->m_rawChainMeasureRunTimes) {
        iRC = gettimeofday(&timevalTarget, NULL);
        queso_require_equal_to_msg(iRC, 0, "gettimeofday called failed");
      }
      logTarget = m_targetPdfSynchronizer->callFunction(&tmpVecValues,&logPrior,&logLikelihood); // Might demand parallel environment
      if (m_optionsObj->m_rawChainMeasureRunTimes) m_rawChainInfo.targetRunTime += MiscGetEllapsedSeconds(&timevalTarget);
      m_rawChainInfo.numTargetCalls++;
      if ((m_env.subDisplayFile()                   ) &&
          (m_env.displayVerbosity() >= 3            ) &&
          (m_optionsObj->m_totallyMute == false)) {
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
        (m_optionsObj->m_totallyMute == false)) {
      *m_env.subDisplayFile() << "\n"
                              << "\n-----------------------------------------------------------\n"
                              << "\n"
                              << std::endl;
    }
    bool accept = false;
    double alphaFirstCandidate = 0.;
    if (outOfTargetSupport) {
      if (m_optionsObj->m_rawChainGenerateExtra) {
        m_alphaQuotients[positionId] = 0.;
      }
    }
    else {
      if (m_optionsObj->m_rawChainMeasureRunTimes) {
        iRC = gettimeofday(&timevalMhAlpha, NULL);
        queso_require_equal_to_msg(iRC, 0, "gettimeofday called failed");
      }
      if (m_optionsObj->m_rawChainGenerateExtra) {
        alphaFirstCandidate = m_algorithm->acceptance_ratio(
            currentPositionData,
            currentCandidateData,
            currentCandidateData.vecValues(),
            currentPositionData.vecValues());
      }
      else {
        alphaFirstCandidate = m_algorithm->acceptance_ratio(
            currentPositionData,
            currentCandidateData,
            currentCandidateData.vecValues(),
            currentPositionData.vecValues());
      }
      if (m_optionsObj->m_rawChainMeasureRunTimes) m_rawChainInfo.mhAlphaRunTime += MiscGetEllapsedSeconds(&timevalMhAlpha);
      if ((m_env.subDisplayFile()                   ) &&
          (m_env.displayVerbosity() >= 10           ) &&
          (m_optionsObj->m_totallyMute == false)) {
        *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                                << ": for chain position of id = " << positionId
                                << std::endl;
      }
      accept = acceptAlpha(alphaFirstCandidate);
    }

    bool displayDetail = (m_env.displayVerbosity() >= 10/*99*/) || m_optionsObj->m_displayCandidates;
    if ((m_env.subDisplayFile()                   ) &&
        (displayDetail                            ) &&
        (m_optionsObj->m_totallyMute == false)) {
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
        (m_optionsObj->m_totallyMute == false)) {
      *m_env.subDisplayFile() << "\n"
                              << "\n-----------------------------------------------------------\n"
                              << "\n"
                              << std::endl;
    }

    //****************************************************
    // Point 3/6 of logic for new position
    // Loop: delayed rejection
    //****************************************************
    if ((accept == false) &&
        (outOfTargetSupport == false) &&
        (m_optionsObj->m_drMaxNumExtraStages > 0)) {
      accept = delayedRejection(positionId,
                                currentPositionData,
                                currentCandidateData);
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
    if ((m_optionsObj->m_rawChainDataOutputPeriod                    >  0  ) &&
        (((positionId+1) % m_optionsObj->m_rawChainDataOutputPeriod) == 0  ) &&
        (m_optionsObj->m_rawChainDataOutputFileName                  != ".")) {
      if ((m_env.subDisplayFile()                   ) &&
          (m_env.displayVerbosity()         >= 10   ) &&
          (m_optionsObj->m_totallyMute == false)) {
        *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                                << ", for chain position of id = " << positionId
                                << ": about to write (per period request) " << m_numPositionsNotSubWritten << " chain positions "
                                << ", " << positionId + 1 - m_optionsObj->m_rawChainDataOutputPeriod << " <= pos <= " << positionId
                                << std::endl;
      }
      workingChain.subWriteContents(positionId + 1 - m_optionsObj->m_rawChainDataOutputPeriod,
                                    m_optionsObj->m_rawChainDataOutputPeriod,
                                    m_optionsObj->m_rawChainDataOutputFileName,
                                    m_optionsObj->m_rawChainDataOutputFileType,
                                    m_optionsObj->m_rawChainDataOutputAllowedSet);
      if ((m_env.subDisplayFile()                   ) &&
          (m_optionsObj->m_totallyMute == false)) {
        *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                                << ", for chain position of id = " << positionId
                                << ": just wrote (per period request) " << m_numPositionsNotSubWritten << " chain positions "
                                << ", " << positionId + 1 - m_optionsObj->m_rawChainDataOutputPeriod << " <= pos <= " << positionId
                                << std::endl;
      }

      if (writeLogLikelihood) {
        workingLogLikelihoodValues->subWriteContents(0 + 1 - m_optionsObj->m_rawChainDataOutputPeriod,
                                                     m_optionsObj->m_rawChainDataOutputPeriod,
                                                     m_optionsObj->m_rawChainDataOutputFileName + "_loglikelihood",
                                                     m_optionsObj->m_rawChainDataOutputFileType,
                                                     m_optionsObj->m_rawChainDataOutputAllowedSet);
      }

      if (writeLogTarget) {
        workingLogTargetValues->subWriteContents(0 + 1 - m_optionsObj->m_rawChainDataOutputPeriod,
                                                 m_optionsObj->m_rawChainDataOutputPeriod,
                                                 m_optionsObj->m_rawChainDataOutputFileName + "_logtarget",
                                                 m_optionsObj->m_rawChainDataOutputFileType,
                                                 m_optionsObj->m_rawChainDataOutputAllowedSet);
      }

      m_numPositionsNotSubWritten = 0;
    }


    if (workingLogLikelihoodValues) (*workingLogLikelihoodValues)[positionId] = currentPositionData.logLikelihood();
    if (workingLogTargetValues    ) (*workingLogTargetValues    )[positionId] = currentPositionData.logTarget();

    if (m_optionsObj->m_rawChainGenerateExtra) {
      m_logTargets[positionId] = currentPositionData.logTarget();
    }

    if (m_optionsObj->m_enableBrooksGelmanConvMonitor > 0) {
      if (positionId % m_optionsObj->m_enableBrooksGelmanConvMonitor == 0 &&
          positionId > m_optionsObj->m_BrooksGelmanLag + 1) {  // +1 to help ensure there are at least 2 samples to use

        double conv_est = workingChain.estimateConvBrooksGelman(
            m_optionsObj->m_BrooksGelmanLag,
            positionId - m_optionsObj->m_BrooksGelmanLag);

        if (m_env.subDisplayFile()) {
          *m_env.subDisplayFile() << "positionId = " << positionId
                                  << ", conv_est = " << conv_est
                                  << std::endl;
          (*m_env.subDisplayFile()).flush();
        }
      }
    }

    // Possibly user-overridden to implement strange things, but we allow it.
    if (positionId % m_optionsObj->m_updateInterval == 0) {
      m_tk->updateTK();
    }

    // If the user dirtied the cov matrix, keep track of the latest iteration
    // number it happened at.
    if (m_tk->covMatrixIsDirty()) {
      m_latestDirtyCovMatrixIteration = positionId;

      // Clean the covariance matrix so that the last dirty iteration tracker
      // doesn't get wiped in the next iteration.
      m_tk->cleanCovMatrix();
    }

    //****************************************************
    // Point 5/6 of logic for new position
    // Adaptive Metropolis calculation
    //****************************************************
    this->adapt(positionId, workingChain);

    //****************************************************
    // Point 6/6 of logic for new position
    // Loop: print some information before going to the next chain position
    //****************************************************
    if ((m_env.subDisplayFile()                   ) &&
        (m_env.displayVerbosity() >= 3            ) &&
        (m_optionsObj->m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                              << ": finishing chain position of id = " << positionId
                              << ", accept = "                         << accept
                              << ", curLogTarget  = "                  << currentPositionData.logTarget()
                              << ", canLogTarget  = "                  << currentCandidateData.logTarget()
                              << std::endl;
    }

    if ((m_optionsObj->m_rawChainDisplayPeriod                     > 0) &&
        (((positionId+1) % m_optionsObj->m_rawChainDisplayPeriod) == 0)) {
      if ((m_env.subDisplayFile()                   ) &&
          (m_optionsObj->m_totallyMute == false)) {
        *m_env.subDisplayFile() << "Finished generating " << positionId+1
                                << " positions"  // root
                                << ", current rejection percentage = " << (100. * ((double) m_rawChainInfo.numRejections)/((double) (positionId+1)))
                                << " %"
                                << std::endl;
      }
    }

    if ((m_env.subDisplayFile()                   ) &&
        (m_env.displayVerbosity() >= 10           ) &&
        (m_optionsObj->m_totallyMute == false)) {
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
                                                NULL);
    if (aux) {}; // just to remove compiler warning
  }

  //****************************************************
  // Print basic information about the chain
  //****************************************************
  m_rawChainInfo.runTime += MiscGetEllapsedSeconds(&timevalChain);
  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_totallyMute == false)) {
    *m_env.subDisplayFile() << "Finished the generation of Markov chain " << workingChain.name()
                            << ", with sub "                              << workingChain.subSequenceSize()
                            << " positions";
    *m_env.subDisplayFile() << "\nSome information about this chain:"
                            << "\n  Chain run time       = " << m_rawChainInfo.runTime
                            << " seconds";
    if (m_optionsObj->m_rawChainMeasureRunTimes) {
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

template <class P_V, class P_M>
void
MetropolisHastingsSG<P_V, P_M>::adapt(unsigned int positionId,
    BaseVectorSequence<P_V, P_M> & workingChain)
{
  int iRC = UQ_OK_RC;
  struct timeval timevalAM;

  // Bail early if we don't satisfy conditions needed to do adaptation
  if ((m_optionsObj->m_tkUseLocalHessian         == true) || // IMPORTANT
      (m_optionsObj->m_amInitialNonAdaptInterval == 0) ||
      (m_optionsObj->m_amAdaptInterval           == 0)) {
    return;
  }

  // Get timing info if we're measuring run times
  if (m_optionsObj->m_rawChainMeasureRunTimes) {
    iRC = gettimeofday(&timevalAM, NULL);
    queso_require_equal_to_msg(iRC, 0, "gettimeofday called failed");
  }

  unsigned int idOfFirstPositionInSubChain = 0;
  SequenceOfVectors<P_V,P_M> partialChain(m_vectorSpace,0,m_optionsObj->m_prefix+"partialChain");

  // Check if now is indeed the moment to adapt
  bool printAdaptedMatrix = false;
  if (positionId < m_optionsObj->m_amInitialNonAdaptInterval) {
    // Do nothing
  }
  else if (positionId == m_optionsObj->m_amInitialNonAdaptInterval) {
    idOfFirstPositionInSubChain = 0;
    partialChain.resizeSequence(m_optionsObj->m_amInitialNonAdaptInterval+1);
    m_lastMean.reset(m_vectorSpace.newVector());
    m_lastAdaptedCovMatrix.reset(m_vectorSpace.newMatrix());
    printAdaptedMatrix = true;
  }
  else {
    unsigned int interval = positionId - m_optionsObj->m_amInitialNonAdaptInterval;
    if ((interval % m_optionsObj->m_amAdaptInterval) == 0) {

      // If the user has set the proposal cov matrix to 'dirty', already
      // recorded the positionId at which that happened in m_latestDirtyCovMatrixIteration
      //
      // If the user didn't dirty it, we're good.
      if (m_latestDirtyCovMatrixIteration > 0) {
        m_lastMean->cwSet(0.0);
        m_lastAdaptedCovMatrix->cwSet(0.0);

        // We'll adapt over the states from when the user dirtied the matrix
        // until the current one
        unsigned int iter_diff = positionId - m_latestDirtyCovMatrixIteration;
        idOfFirstPositionInSubChain = iter_diff;
        partialChain.resizeSequence(iter_diff);

        // Finally set the latest dirty iteration back to zero.  If the user
        // sets the dirty flag again, then this will change.
        m_latestDirtyCovMatrixIteration = 0;
      }
      else {
        idOfFirstPositionInSubChain = positionId - m_optionsObj->m_amAdaptInterval;
        partialChain.resizeSequence(m_optionsObj->m_amAdaptInterval);
      }

      if (m_optionsObj->m_amAdaptedMatricesDataOutputPeriod > 0) {
        if ((interval % m_optionsObj->m_amAdaptedMatricesDataOutputPeriod) == 0) {
          printAdaptedMatrix = true;
        }
      }
    }
  }

  // Bail out if we don't have the samples to adapt
  if (partialChain.subSequenceSize() == 0) {
    // Save timings and bail
    if (m_optionsObj->m_rawChainMeasureRunTimes) {
      m_rawChainInfo.amRunTime += MiscGetEllapsedSeconds(&timevalAM);
    }

    return;
  }

  // If now is indeed the moment to adapt, then do it!
  P_V transporterVec(m_vectorSpace.zeroVector());
  for (unsigned int i = 0; i < partialChain.subSequenceSize(); ++i) {
    workingChain.getPositionValues(idOfFirstPositionInSubChain+i,transporterVec);

    // Transform to the space without boundaries.  This is the space
    // where the proposal distribution is Gaussian
    if (this->m_optionsObj->m_tk == "logit_random_walk") {
      // Only do this when we don't use the Hessian (this may change in
      // future, but transformToGaussianSpace() is only implemented in
      // TransformedScaledCovMatrixTKGroup
      P_V transformedTransporterVec(m_vectorSpace.zeroVector());
      dynamic_cast<TransformedScaledCovMatrixTKGroup<P_V, P_M>* >(
          m_tk.get())->transformToGaussianSpace(transporterVec,
            transformedTransporterVec);
      partialChain.setPositionValues(i, transformedTransporterVec);
    }
    else {
      partialChain.setPositionValues(i, transporterVec);
    }
  }

  updateAdaptedCovMatrix(partialChain,
                         idOfFirstPositionInSubChain,
                         m_lastChainSize,
                         *m_lastMean,
                         *m_lastAdaptedCovMatrix);

  // Print adapted matrix info
  if ((printAdaptedMatrix == true) &&
      (m_optionsObj->m_amAdaptedMatricesDataOutputFileName != "." )) {
    char varNamePrefix[64];
    sprintf(varNamePrefix,"mat_am%d",positionId);

    char tmpChar[64];
    sprintf(tmpChar,"_am%d",positionId);

    std::set<unsigned int> tmpSet;
    tmpSet.insert(m_env.subId());

    m_lastAdaptedCovMatrix->subWriteContents(varNamePrefix,
                                             (m_optionsObj->m_amAdaptedMatricesDataOutputFileName+tmpChar),
                                             m_optionsObj->m_amAdaptedMatricesDataOutputFileType,
                                             tmpSet);
    if ((m_env.subDisplayFile()                   ) &&
        (m_optionsObj->m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                              << ": just wrote last adapted proposal cov matrix contents = " << *m_lastAdaptedCovMatrix
                              << std::endl;
    }
  }

  // Check if adapted matrix is positive definite
  bool tmpCholIsPositiveDefinite = false;
  P_M tmpChol(*m_lastAdaptedCovMatrix);
  P_M attemptedMatrix(tmpChol);
  if ((m_env.subDisplayFile()        ) &&
      (m_env.displayVerbosity() >= 10)) {
//(m_optionsObj->m_totallyMute == false)) {
    *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                            << ", positionId = "  << positionId
                            << ": 'am' calling first tmpChol.chol()"
                            << std::endl;
  }
  iRC = tmpChol.chol();
  if (iRC) {
    std::string err1 = "In MetropolisHastingsSG<P_V,P_M>::adapt(): first ";
    err1 += "Cholesky factorisation of proposal covariance matrix ";
    err1 += "failed.  QUESO will now attempt to regularise the ";
    err1 += "matrix before proceeding.  This is not a fatal error.";
    std::cerr << err1 << std::endl;
  }

  // Print some infor about the Cholesky factorisation
  if ((m_env.subDisplayFile()        ) &&
      (m_env.displayVerbosity() >= 10)) {
//(m_optionsObj->m_totallyMute == false)) {
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

  // If the Cholesky factorisation failed, add a regularisation to the
  // diagonal components (of size m_amEpsilon) and try the factorisation again
  if (iRC) {
    queso_require_equal_to_msg(iRC, UQ_MATRIX_IS_NOT_POS_DEFINITE_RC, "invalid iRC returned from first chol()");
    // Matrix is not positive definite
    P_M* tmpDiag = m_vectorSpace.newDiagMatrix(m_optionsObj->m_amEpsilon);
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
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                              << ", positionId = "  << positionId
                              << ": 'am' calling second tmpChol.chol()"
                              << std::endl;
    }

    // Trying the second (regularised Cholesky factorisation)
    iRC = tmpChol.chol();
    if (iRC) {
      std::string err2 = "In MetropolisHastingsSG<P_V,P_M>::adapt(): second ";
      err2 += "Cholesky factorisation of (regularised) proposal ";
      err2 += "covariance matrix failed.  QUESO is falling back to ";
      err2 += "the last factorisable proposal covariance matrix.  ";
      err2 += "This is not a fatal error.";
      std::cerr << err2 << std::endl;
    }

    // Print some diagnostic info
    if ((m_env.subDisplayFile()        ) &&
        (m_env.displayVerbosity() >= 10)) {
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

    // If the second (regularised) Cholesky factorisation failed, do nothing
    if (iRC) {
      queso_require_equal_to_msg(iRC, UQ_MATRIX_IS_NOT_POS_DEFINITE_RC, "invalid iRC returned from second chol()");
      // Do nothing
    }
    else {
      tmpCholIsPositiveDefinite = true;
    }
  }
  else {
    tmpCholIsPositiveDefinite = true;
  }

  // If the adapted matrix is pos. def., scale by \eta (s_d in Haario paper)
  if (tmpCholIsPositiveDefinite) {
    P_M tmpMatrix(m_optionsObj->m_amEta*attemptedMatrix);
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

    // Transform the proposal covariance matrix if we have Logit transforms
    // turned on.
    //
    // This logic should really be done inside the kernel itself
    if (this->m_optionsObj->m_tk == "logit_random_walk") {
      (dynamic_cast<TransformedScaledCovMatrixTKGroup<P_V,P_M>* >(m_tk.get()))
        ->updateLawCovMatrix(tmpMatrix);
    }
    else {
      // Assume that if we're not doing a logit transform, we're doing a random
      // something sensible
      (dynamic_cast<ScaledCovMatrixTKGroup<P_V,P_M>* >(m_tk.get()))
        ->updateLawCovMatrix(tmpMatrix);
    }

#ifdef UQ_DRAM_MCG_REQUIRES_INVERTED_COV_MATRICES
    queso_require_msg(!(UQ_INCOMPLETE_IMPLEMENTATION_RC), "need to code the update of m_upperCholProposalPrecMatrices");
#endif
  }

  // FIXME: What is this for?  Check the destructor frees the memory and nuke
  //        the commented code below
  //for (unsigned int i = 0; i < partialChain.subSequenceSize(); ++i) {
  //  if (partialChain[i]) delete partialChain[i];
  //}

  if (m_optionsObj->m_rawChainMeasureRunTimes) {
    m_rawChainInfo.amRunTime += MiscGetEllapsedSeconds(&timevalAM);
  }
}

template <class P_V, class P_M>
bool
MetropolisHastingsSG<P_V, P_M>::delayedRejection(unsigned int positionId,
    const MarkovChainPositionData<P_V> & currentPositionData,
    MarkovChainPositionData<P_V> & currentCandidateData)
{
  if ((m_optionsObj->m_drDuringAmNonAdaptiveInt  == false     ) &&
      (m_optionsObj->m_tkUseLocalHessian         == false     ) &&
      (m_optionsObj->m_amInitialNonAdaptInterval >  0         ) &&
      (m_optionsObj->m_amAdaptInterval           >  0         ) &&
      (positionId <= m_optionsObj->m_amInitialNonAdaptInterval)) {
    return false;
  }

  unsigned int stageId = 0;

  bool validPreComputingPosition;

  m_tk->clearPreComputingPositions();

  validPreComputingPosition = m_tk->setPreComputingPosition(
      currentPositionData.vecValues(), 0);

  validPreComputingPosition = m_tk->setPreComputingPosition(
      currentCandidateData.vecValues(), stageId + 1);

  std::vector<MarkovChainPositionData<P_V>*> drPositionsData(stageId+2,NULL);
  std::vector<unsigned int> tkStageIds (stageId+2,0);

  int iRC = UQ_OK_RC;
  struct timeval timevalDR;
  struct timeval timevalDrAlpha;
  struct timeval timevalCandidate;
  struct timeval timevalTarget;

  if (m_optionsObj->m_rawChainMeasureRunTimes) {
    iRC = gettimeofday(&timevalDR, NULL);
    queso_require_equal_to_msg(iRC, 0, "gettimeofday call failed");
  }

  drPositionsData[0] = new MarkovChainPositionData<P_V>(currentPositionData );
  drPositionsData[1] = new MarkovChainPositionData<P_V>(currentCandidateData);

  tkStageIds[0] = 0;
  tkStageIds[1] = 1;

  bool accept = false;
  while ((validPreComputingPosition == true                 ) &&
         (accept                    == false                ) &&
         (stageId < m_optionsObj->m_drMaxNumExtraStages)) {
    if ((m_env.subDisplayFile()                   ) &&
        (m_env.displayVerbosity() >= 10           ) &&
        (m_optionsObj->m_totallyMute == false)) {
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
        (m_optionsObj->m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                              << ": for chain position of id = " << positionId
                              << ", beginning stageId = "        << stageId
                              << std::endl;
    }

    P_V tmpVecValues(currentCandidateData.vecValues());
    bool keepGeneratingCandidates = true;
    bool outOfTargetSupport = false;
    while (keepGeneratingCandidates) {
      if (m_optionsObj->m_rawChainMeasureRunTimes) {
        iRC = gettimeofday(&timevalCandidate, NULL);
        queso_require_equal_to_msg(iRC, 0, "gettimeofday call failed");
      }
      m_tk->rv(tkStageIds).realizer().realization(tmpVecValues);
      if (m_numDisabledParameters > 0) { // gpmsa2
        for (unsigned int paramId = 0; paramId < m_vectorSpace.dimLocal(); ++paramId) {
          if (m_parameterEnabledStatus[paramId] == false) {
            tmpVecValues[paramId] = m_initialPosition[paramId];
          }
        }
      }
      if (m_optionsObj->m_rawChainMeasureRunTimes) m_rawChainInfo.candidateRunTime += MiscGetEllapsedSeconds(&timevalCandidate);

      outOfTargetSupport = !m_targetPdf.domainSet().contains(tmpVecValues);

      if (m_optionsObj->m_putOutOfBoundsInChain) keepGeneratingCandidates = false;
      else                                            keepGeneratingCandidates = outOfTargetSupport;
    }

    if ((m_env.subDisplayFile()                   ) &&
        (m_env.displayVerbosity() >= 5            ) &&
        (m_optionsObj->m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                              << ": about to set TK pre computing position of local id " << stageId+1
                              << ", values = " << tmpVecValues
                              << std::endl;
    }
    validPreComputingPosition = m_tk->setPreComputingPosition(tmpVecValues,stageId+1);
    if ((m_env.subDisplayFile()                   ) &&
        (m_env.displayVerbosity() >= 5            ) &&
        (m_optionsObj->m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In MetropolisHastingsSG<P_V,P_M>::generateFullChain()"
                              << ": returned from setting TK pre computing position of local id " << stageId+1
                              << ", values = " << tmpVecValues
                              << ", valid = "  << validPreComputingPosition
                              << std::endl;
    }

    double logPrior;
    double logLikelihood;
    double logTarget;
    if (outOfTargetSupport) {
      m_rawChainInfo.numOutOfTargetSupportInDR++; // new 2010/May/12
      logPrior      = -INFINITY;
      logLikelihood = -INFINITY;
      logTarget     = -INFINITY;
    }
    else {
      if (m_optionsObj->m_rawChainMeasureRunTimes) {
        iRC = gettimeofday(&timevalTarget, NULL);
        queso_require_equal_to_msg(iRC, 0, "gettimeofday call failed");
      }
      logTarget = m_targetPdfSynchronizer->callFunction(&tmpVecValues,&logPrior,&logLikelihood); // Might demand parallel environment
      if (m_optionsObj->m_rawChainMeasureRunTimes) m_rawChainInfo.targetRunTime += MiscGetEllapsedSeconds(&timevalTarget);
      m_rawChainInfo.numTargetCalls++;
      if ((m_env.subDisplayFile()                   ) &&
          (m_env.displayVerbosity() >= 3            ) &&
          (m_optionsObj->m_totallyMute == false)) {
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

    // Ok, so we almost don't need setPreComputingPosition.  All the DR
    // position information we needed was generated in this while loop.
    drPositionsData.push_back(new MarkovChainPositionData<P_V>(currentCandidateData));
    tkStageIds.push_back     (stageId+1);

    double alphaDR = 0.;
    if (outOfTargetSupport == false) {
      if (m_optionsObj->m_rawChainMeasureRunTimes) {
        iRC = gettimeofday(&timevalDrAlpha, NULL);
        queso_require_equal_to_msg(iRC, 0, "gettimeofday call failed");
      }
      alphaDR = this->alpha(drPositionsData,tkStageIds);
      if (m_optionsObj->m_rawChainMeasureRunTimes) m_rawChainInfo.drAlphaRunTime += MiscGetEllapsedSeconds(&timevalDrAlpha);
      accept = acceptAlpha(alphaDR);
    }

    bool displayDetail = (m_env.displayVerbosity() >= 10/*99*/) || m_optionsObj->m_displayCandidates;
    if ((m_env.subDisplayFile()                   ) &&
        (displayDetail                            ) &&
        (m_optionsObj->m_totallyMute == false)) {
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

  if (m_optionsObj->m_rawChainMeasureRunTimes) m_rawChainInfo.drRunTime += MiscGetEllapsedSeconds(&timevalDR);

  for (unsigned int i = 0; i < drPositionsData.size(); ++i) {
    if (drPositionsData[i]) delete drPositionsData[i];
  }

  return accept;
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
    queso_require_greater_equal_msg(partialChain.subSequenceSize(), 2, "'partialChain.subSequenceSize()' should be >= 2");

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
    queso_require_greater_equal_msg(partialChain.subSequenceSize(), 1, "'partialChain.subSequenceSize()' should be >= 1");
    queso_require_greater_equal_msg(idOfFirstPositionInSubChain, 1, "'idOfFirstPositionInSubChain' should be >= 1");

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

template <class P_V, class P_M>
void
MetropolisHastingsSG<P_V, P_M>::transformInitialCovMatrixToGaussianSpace(
    const BoxSubset<P_V, P_M> & boxSubset)
{
  P_V min_domain_bounds(boxSubset.minValues());
  P_V max_domain_bounds(boxSubset.maxValues());

  for (unsigned int i = 0; i < min_domain_bounds.sizeLocal(); i++) {
    double min_val = min_domain_bounds[i];
    double max_val = max_domain_bounds[i];

    if (queso_isfinite(min_val) && queso_isfinite(max_val)) {
      if (m_initialProposalCovMatrix(i, i) >= max_val - min_val) {
        // User is trying to specify a uniform proposal distribution, which
        // is unsupported.  Throw an error for now.
        std::cerr << "Proposal variance element "
                  << i
                  << " is "
                  << m_initialProposalCovMatrix(i, i)
                  << " but domain is of size "
                  << max_val - min_val
                  << std::endl;
        std::cerr << "QUESO does not support uniform-like proposal "
                  << "distributions.  Try making the proposal variance smaller"
                  << std::endl;
      }

      // The jacobian at the midpoint of the domain
      double transformJacobian = 4.0 / (max_val - min_val);

      // Just do the multiplication by hand for now.  There's no method in
      // Gsl(Vector|Matrix) to do this for me.
      for (unsigned int j = 0; j < min_domain_bounds.sizeLocal(); j++) {
        // Multiply column j by element j
        m_initialProposalCovMatrix(j, i) *= transformJacobian;
      }
      for (unsigned int j = 0; j < min_domain_bounds.sizeLocal(); j++) {
        // Multiply row j by element j
        m_initialProposalCovMatrix(i, j) *= transformJacobian;
      }
    }
  }
}

template class MetropolisHastingsSG<GslVector, GslMatrix>;

}  // End namespace QUESO
