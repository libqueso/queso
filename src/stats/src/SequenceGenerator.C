//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2015 The PECOS Development Team
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

#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/VectorRV.h>
#include <queso/TKGroup.h>
#include <queso/ScalarFunctionSynchronizer.h>
#include <queso/MetropolisHastingsSGOptions.h>
#include <queso/MarkovChainPositionData.h>

#include <queso/SequenceGenerator.h>

namespace QUESO {

template<class V, class M>
SequenceGenerator<V,M>::SequenceGenerator(const char * prefix,
                                          const MhOptionsValues * alternativeOptionsValues,
                                          const BaseVectorRV<V, M> & sourceRv,
                                          const V & initialPosition,
                                          const M * inputProposalCovMatrix)
: m_env(sourceRv.env()),
  m_vectorSpace(sourceRv.imageSet().vectorSpace()),
  m_targetPdf(sourceRv.pdf()),
  m_initialPosition(initialPosition),
  m_initialProposalCovMatrix(m_vectorSpace.zeroVector()),
  m_nullInputProposalCovMatrix(inputProposalCovMatrix == NULL),
  m_targetPdfSynchronizer(new ScalarFunctionSynchronizer<V,M>(m_targetPdf,m_initialPosition)),
  m_tk(NULL),
  m_idsOfUniquePositions(0),
  m_logTargets(0),
  m_alphaQuotients(0),
  m_lastChainSize(0),
  m_lastMean(NULL),
  m_lastAdaptedCovMatrix(NULL),
  m_optionsObj(alternativeOptionsValues),
  m_computeInitialPriorAndLikelihoodValues(true),
  m_initialLogPriorValue(0.),
  m_initialLogLikelihoodValue(0.),
  m_userDidNotProvideOptions(false)
{
  if (inputProposalCovMatrix != NULL) {
    m_initialProposalCovMatrix = *inputProposalCovMatrix;
  }

  queso_require_equal_to_msg(sourceRv.imageSet().vectorSpace().dimLocal(), initialPosition.sizeLocal(), "'sourceRv' and 'initialPosition' should have equal dimensions");

  if (inputProposalCovMatrix) {
    queso_require_equal_to_msg(sourceRv.imageSet().vectorSpace().dimLocal(), inputProposalCovMatrix->numRowsLocal(), "'sourceRv' and 'inputProposalCovMatrix' should have equal dimensions");
    queso_require_equal_to_msg(inputProposalCovMatrix->numCols(), inputProposalCovMatrix->numRowsGlobal(), "'inputProposalCovMatrix' should be a square matrix");
  }
}

template<class V, class M>
SequenceGenerator<V, M>::~SequenceGenerator()
{
  if (m_lastAdaptedCovMatrix) {
    delete m_lastAdaptedCovMatrix;
  }
  if (m_lastMean) {
    delete m_lastMean;
  }
  m_lastChainSize = 0;
  m_alphaQuotients.clear();
  m_logTargets.clear();
  m_idsOfUniquePositions.clear();

  if (m_tk) {
    delete m_tk;
  }
  if (m_targetPdfSynchronizer) {
    delete m_targetPdfSynchronizer;
  }

  // Only delete if the user didn't provide the options
  // I.e., if the user created their options object, then they are resonsible
  // for freeing it.
  if (m_optionsObj && m_userDidNotProvideOptions) {
    delete m_optionsObj;
  }
}

template<class V, class M>
double
SequenceGenerator<V,M>::alpha(const MarkovChainPositionData<V> & x,
                              const MarkovChainPositionData<V> & y,
                              unsigned int xStageId,
                              unsigned int yStageId,
                              double * alphaQuotientPtr)
{
  double alphaQuotient = 0.;
  if ((x.outOfTargetSupport() == false) &&
      (y.outOfTargetSupport() == false)) {
    if ((x.logTarget() == -INFINITY) ||
        (x.logTarget() ==  INFINITY) ||
        ( (boost::math::isnan)(x.logTarget())      )) {
      std::cerr << "WARNING In SequenceGenerator<V,M>::alpha(x,y)"
                << ", worldRank "       << m_env.worldRank()
                << ", fullRank "        << m_env.fullRank()
                << ", subEnvironment "  << m_env.subId()
                << ", subRank "         << m_env.subRank()
                << ", inter0Rank "      << m_env.inter0Rank()
                << ": x.logTarget() = " << x.logTarget()
                << ", x.values() = "    << x.vecValues()
                << ", y.values() = "    << y.vecValues()
                << std::endl;
    }
    else if ((y.logTarget() == -INFINITY           ) ||
             (y.logTarget() ==  INFINITY           ) ||
             ( (boost::math::isnan)(y.logTarget()) )) {
      std::cerr << "WARNING In SequenceGenerator<V,M>::alpha(x,y)"
                << ", worldRank "       << m_env.worldRank()
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

      if (m_tk->symmetric()) {
        alphaQuotient = std::exp(yLogTargetToUse - x.logTarget());

        if ((m_env.subDisplayFile()                   ) &&
            (m_env.displayVerbosity() >= 3            ) &&
            (m_optionsObj->m_totallyMute == false)) {
          *m_env.subDisplayFile() << "In SequenceGenerator<V,M>::alpha(x,y)"
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
        double qyx = m_tk->rv(yStageId).pdf().lnValue(x.vecValues(),NULL,NULL,NULL,NULL);
        double qxy = m_tk->rv(xStageId).pdf().lnValue(y.vecValues(),NULL,NULL,NULL,NULL);
        alphaQuotient = std::exp(yLogTargetToUse +
                                 qyx -
                                 x.logTarget() -
                                 qxy);
        if ((m_env.subDisplayFile()                   ) &&
            (m_env.displayVerbosity() >= 3            ) &&
            (m_optionsObj->m_totallyMute == false)) {
          *m_env.subDisplayFile() << "In SequenceGenerator<V,M>::alpha(x,y)"
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
        (m_optionsObj->m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In SequenceGenerator<V,M>::alpha(x,y)"
                             << ": x.outOfTargetSupport = " << x.outOfTargetSupport()
                             << ", y.outOfTargetSupport = " << y.outOfTargetSupport()
                             << std::endl;
    }
  }
  if (alphaQuotientPtr != NULL) *alphaQuotientPtr = alphaQuotient;

  return std::min(1.,alphaQuotient);
}

template<class V, class M>
bool
SequenceGenerator<V, M>::acceptAlpha(double alpha)
{
  bool result = false;

  if      (alpha <= 0.                                ) result = false;
  else if (alpha >= 1.                                ) result = true;
  else if (alpha >= m_env.rngObject()->uniformSample()) result = true;
  else                                                  result = false;

  return result;
}

template<class V, class M>
int
SequenceGenerator<V,M>::writeInfo(const BaseVectorSequence<V,M> & workingChain,
                                  std::ofstream & ofsvar) const
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

  if (false) { // Don't see need for code below. Let it there though, compiling, in case it is needed in the future.
    // Write names of components
    ofsvar << m_optionsObj->m_prefix << "componentNames = {";
    m_vectorSpace.printComponentsNames(ofsvar,false);
    ofsvar << "};\n";
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

template <class V, class M>
const BaseTKGroup<V, M> &
SequenceGenerator<V, M>::transitionKernel() const
{
  return *m_tk;
}

template <class V, class M>
void
SequenceGenerator<V, M>::generateSequence(
  BaseVectorSequence<V, M> & workingChain,
  ScalarSequence<double> * workingLogLikelihoodValues,
  ScalarSequence<double> * workingLogTargetValues)
{
  if ((m_env.subDisplayFile()                   ) &&
      (m_env.displayVerbosity() >= 5            ) &&
      (m_optionsObj->m_totallyMute == false)) {
    *m_env.subDisplayFile() << "Entering SequenceGenerator<V,M>::generateSequence()..."
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

  MiscCheckTheParallelEnvironment<V,V>(m_initialPosition,
                                             m_initialPosition);

  V valuesOf1stPosition(m_initialPosition);
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
    *m_env.subDisplayFile() << "In SequenceGenerator<V,M>::generateSequence()"
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
    *m_env.subDisplayFile() << "In SequenceGenerator<V,M>::generateSequence()"
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
      *m_env.subDisplayFile() << "In SequenceGenerator<V,M>::generateSequence()"
                              << ", prefix = "                                         << m_optionsObj->m_prefix
                              << ", raw chain name = "                                 << workingChain.name()
                              << ": about to try to write raw sub chain output file '" << m_optionsObj->m_rawChainDataOutputFileName
                              << "."                                                   << m_optionsObj->m_rawChainDataOutputFileType
                              << "', subId = "                                         << m_env.subId()
                              << ", subenv is allowed to write  1/true or 0/false) = " << (m_optionsObj->m_rawChainDataOutputAllowedSet.find(m_env.subId()) != m_optionsObj->m_rawChainDataOutputAllowedSet.end())
                              << "..."
                              << std::endl;
    }

    // Compute raw sub MLE
    if (workingLogLikelihoodValues) {
      SequenceOfVectors<V,M> rawSubMLEpositions(m_vectorSpace,0,m_optionsObj->m_prefix+"rawSubMLEseq");
      double rawSubMLEvalue = workingChain.subPositionsOfMaximum(*workingLogLikelihoodValues,
                                                                 rawSubMLEpositions);
      queso_require_not_equal_to_msg(rawSubMLEpositions.subSequenceSize(), 0, "rawSubMLEpositions.subSequenceSize() = 0");

      if ((m_env.subDisplayFile()                   ) &&
          (m_optionsObj->m_totallyMute == false)) {
        V tmpVec(m_vectorSpace.zeroVector());
        rawSubMLEpositions.getPositionValues(0,tmpVec);
        *m_env.subDisplayFile() << "In SequenceGenerator<V,M>::generateSequence()"
                                << ": just computed MLE"
                                << ", rawSubMLEvalue = "                       << rawSubMLEvalue
                                << ", rawSubMLEpositions.subSequenceSize() = " << rawSubMLEpositions.subSequenceSize()
                                << ", rawSubMLEpositions[0] = "                << tmpVec
                                << std::endl;
      }
    }

    // Compute raw sub MAP
    if (workingLogTargetValues) {
      SequenceOfVectors<V,M> rawSubMAPpositions(m_vectorSpace,0,m_optionsObj->m_prefix+"rawSubMAPseq");
      double rawSubMAPvalue = workingChain.subPositionsOfMaximum(*workingLogTargetValues,
                                                                 rawSubMAPpositions);
      queso_require_not_equal_to_msg(rawSubMAPpositions.subSequenceSize(), 0, "rawSubMAPpositions.subSequenceSize() = 0");

      if ((m_env.subDisplayFile()                   ) &&
          (m_optionsObj->m_totallyMute == false)) {
        V tmpVec(m_vectorSpace.zeroVector());
        rawSubMAPpositions.getPositionValues(0,tmpVec);
        *m_env.subDisplayFile() << "In SequenceGenerator<V,M>::generateSequence()"
                                << ": just computed MAP"
                                << ", rawSubMAPvalue = "                       << rawSubMAPvalue
                                << ", rawSubMAPpositions.subSequenceSize() = " << rawSubMAPpositions.subSequenceSize()
                                << ", rawSubMAPpositions[0] = "                << tmpVec
                                << std::endl;
      }
    }

    if ((m_env.subDisplayFile()                   ) &&
        (m_optionsObj->m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In SequenceGenerator<V,M>::generateSequence()"
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
      *m_env.subDisplayFile() << "In SequenceGenerator<V,M>::generateSequence()"
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
      *m_env.subDisplayFile() << "In SequenceGenerator<V,M>::generateSequence()"
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

    // Compute raw unified MLE
    if (workingLogLikelihoodValues) {
      SequenceOfVectors<V,M> rawUnifiedMLEpositions(m_vectorSpace,0,m_optionsObj->m_prefix+"rawUnifiedMLEseq");

      double rawUnifiedMLEvalue = workingChain.unifiedPositionsOfMaximum(*workingLogLikelihoodValues,
                                                                         rawUnifiedMLEpositions);

      if ((m_env.subDisplayFile()                   ) &&
          (m_optionsObj->m_totallyMute == false)) {
        V tmpVec(m_vectorSpace.zeroVector());

        // Make sure the positions vector (which only contains stuff on process
        // zero) actually contains positions
        if (rawUnifiedMLEpositions.subSequenceSize() > 0) {
          rawUnifiedMLEpositions.getPositionValues(0,tmpVec);
          *m_env.subDisplayFile() << "In SequenceGenerator<V,M>::generateSequence()"
                                  << ": just computed MLE"
                                  << ", rawUnifiedMLEvalue = "                       << rawUnifiedMLEvalue
                                  << ", rawUnifiedMLEpositions.subSequenceSize() = " << rawUnifiedMLEpositions.subSequenceSize()
                                  << ", rawUnifiedMLEpositions[0] = "                << tmpVec
                                  << std::endl;
        }
      }
    }

    // Compute raw unified MAP
    if (workingLogTargetValues) {
      SequenceOfVectors<V,M> rawUnifiedMAPpositions(m_vectorSpace,0,m_optionsObj->m_prefix+"rawUnifiedMAPseq");
      double rawUnifiedMAPvalue = workingChain.unifiedPositionsOfMaximum(*workingLogTargetValues,
                                                                         rawUnifiedMAPpositions);

      if ((m_env.subDisplayFile()                   ) &&
          (m_optionsObj->m_totallyMute == false)) {
        V tmpVec(m_vectorSpace.zeroVector());

        // Make sure the positions vector (which only contains stuff on process
        // zero) actually contains positions
        if (rawUnifiedMAPpositions.subSequenceSize() > 0) {
          rawUnifiedMAPpositions.getPositionValues(0,tmpVec);
          *m_env.subDisplayFile() << "In SequenceGenerator<V,M>::generateSequence()"
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
      *m_env.subDisplayFile() << "In SequenceGenerator<V,M>::generateSequence()"
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
        *m_env.subDisplayFile() << "In SequenceGenerator<V,M>::generateSequence()"
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
        *m_env.subDisplayFile() << "In SequenceGenerator<V,M>::generateSequence()"
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
      *m_env.subDisplayFile() << "In SequenceGenerator<V,M>::generateSequence()"
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

  //m_env.syncPrintDebugMsg("Leaving SequenceGenerator<V,M>::generateSequence()",2,3000000,m_env.fullComm()); // Dangerous to barrier on fullComm ... // KAUST

  if ((m_env.subDisplayFile()                   ) &&
      (m_env.displayVerbosity() >= 5            ) &&
      (m_optionsObj->m_totallyMute == false)) {
    *m_env.subDisplayFile() << "Leaving SequenceGenerator<V,M>::generateSequence()"
                            << std::endl;
  }

  return;
}


template <class V, class M>
void
SequenceGenerator<V, M>::readFullChain(const std::string & inputFileName,
                                       const std::string & inputFileType,
                                       unsigned int chainSize,
                                       BaseVectorSequence<V, M> & workingChain)
{
  workingChain.unifiedReadContents(inputFileName, inputFileType, chainSize);
}

template <class V, class M>
void
SequenceGenerator<V, M>::generateFullChain(
    const V & valuesOf1stPosition,
    unsigned int chainSize,
    BaseVectorSequence<V, M> & workingChain,
    ScalarSequence<double> * workingLogLikelihoodValues,
    ScalarSequence<double> * workingLogTargetValues)
{
  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_totallyMute == false)) {
    *m_env.subDisplayFile() << "Starting the generation of Markov chain " << workingChain.name()
                            << ", with "                                  << chainSize
                            << " positions..."
                            << std::endl;
  }

  int iRC = UQ_OK_RC;

  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_totallyMute == false)) {
    *m_env.subDisplayFile() << "\nIn SequenceGenerator<V,M>::generateFullChain()"
                            << ": contents of initial position are:";
    *m_env.subDisplayFile() << valuesOf1stPosition; // FIX ME: might need parallelism
    *m_env.subDisplayFile() << "\nIn SequenceGenerator<V,M>::generateFullChain()"
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
    *m_env.subDisplayFile() << "ERROR: In SequenceGenerator<V,M>::generateFullChain()"
                            << ": contents of initial position are:\n";
    *m_env.subDisplayFile() << valuesOf1stPosition; // FIX ME: might need parallelism
    *m_env.subDisplayFile() << "\nERROR: In SequenceGenerator<V,M>::generateFullChain()"
                            << ": targetPdf.domaintSet() info is:\n"
                            << m_targetPdf.domainSet();
    *m_env.subDisplayFile() << std::endl;
  }
  queso_require_msg(!(outOfTargetSupport), "initial position should not be out of target pdf support");

  double logPrior      = 0.;
  double logLikelihood = 0.;
  double logTarget     = 0.;
  if (m_computeInitialPriorAndLikelihoodValues) {
    logTarget = m_targetPdfSynchronizer->callFunction(&valuesOf1stPosition,NULL,NULL,NULL,NULL,&logPrior,&logLikelihood); // Might demand parallel environment // KEY
    if ((m_env.subDisplayFile()                   ) &&
        (m_env.displayVerbosity() >= 3            ) &&
        (m_optionsObj->m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In SequenceGenerator<V,M>::generateFullChain()"
            << ": just returned from likelihood() for initial chain position"
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
      *m_env.subDisplayFile() << "In SequenceGenerator<V,M>::generateFullChain()"
            << ": used input prior and likelihood for initial chain position"
            << ", logPrior = "      << logPrior
            << ", logLikelihood = " << logLikelihood
            << ", logTarget = "     << logTarget
            << std::endl;
    }
  }

  MarkovChainPositionData<V> currentPositionData(m_env,
                                                 valuesOf1stPosition,
                                                 outOfTargetSupport,
                                                 logLikelihood,
                                                 logTarget);

  V gaussianVector(m_vectorSpace.zeroVector());
  V tmpVecValues(m_vectorSpace.zeroVector());
  MarkovChainPositionData<V> currentCandidateData(m_env);

  //****************************************************
  // Set chain position with positionId = 0
  //****************************************************
  workingChain.resizeSequence(chainSize);
  if (workingLogLikelihoodValues) workingLogLikelihoodValues->resizeSequence(chainSize);
  if (workingLogTargetValues    ) workingLogTargetValues->resizeSequence    (chainSize);
  if (true/*m_uniqueChainGenerate*/) m_idsOfUniquePositions.resize(chainSize,0);
  if (m_optionsObj->m_rawChainGenerateExtra) {
    m_logTargets.resize    (chainSize,0.);
    m_alphaQuotients.resize(chainSize,0.);
  }

  unsigned int uniquePos = 0;
  workingChain.setPositionValues(0,currentPositionData.vecValues());
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
      *m_env.subDisplayFile() << "In SequenceGenerator<V,M>::generateFullChain()"
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
  }

  if (workingLogLikelihoodValues) (*workingLogLikelihoodValues)[0] = currentPositionData.logLikelihood();
  if (workingLogTargetValues    ) (*workingLogTargetValues    )[0] = currentPositionData.logTarget();
  if (true/*m_uniqueChainGenerate*/) m_idsOfUniquePositions[uniquePos++] = 0;
  if (m_optionsObj->m_rawChainGenerateExtra) {
    m_logTargets    [0] = currentPositionData.logTarget();
    m_alphaQuotients[0] = 1.;
  }

  if ((m_env.subDisplayFile()                   ) &&
      (m_env.displayVerbosity() >= 10           ) &&
      (m_optionsObj->m_totallyMute == false)) {
    *m_env.subDisplayFile() << "\n"
                            << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                            << "\n"
                            << std::endl;
  }

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
    }
  }
  else {
    for (unsigned int positionId = 1; positionId < workingChain.subSequenceSize(); ++positionId) {
      m_tk->clearPreComputingPositions();

      //****************************************************
      // Point 1/6 of logic for new position
      // Loop: generate new position
      //****************************************************

      // positionId is the index of the proposed state, so we subtract 1 to get
      // the index of the current state
      propose(positionId - 1, workingChain, tmpVecValues);

      queso_require_msg(m_targetPdf.domainSet().contains(tmpVecValues),
          "generated proposal lies outside of the support of the posterior");

      outOfTargetSupport = false;

      logTarget = m_targetPdfSynchronizer->callFunction(&tmpVecValues,NULL,NULL,NULL,NULL,&logPrior,&logLikelihood); // Might demand parallel environment

      if ((m_env.subDisplayFile()                   ) &&
          (m_env.displayVerbosity() >= 3            ) &&
          (m_optionsObj->m_totallyMute == false)) {
        *m_env.subDisplayFile() << "In SequenceGenerator<V,M>::generateFullChain()"
                                << ": just returned from likelihood() for chain position of id " << positionId
                                << ", logPrior = "      << logPrior
                                << ", logLikelihood = " << logLikelihood
                                << ", logTarget = "     << logTarget
                                << std::endl;
      }

      currentCandidateData.set(tmpVecValues,
                               outOfTargetSupport,
                               logLikelihood,
                               logTarget);

      bool accept = false;
      double alphaFirstCandidate = 0.;
      if (m_optionsObj->m_rawChainGenerateExtra) {
        alphaFirstCandidate = this->alpha(currentPositionData,currentCandidateData,0,1,&m_alphaQuotients[positionId]);
      }
      else {
        alphaFirstCandidate = this->alpha(currentPositionData,currentCandidateData,0,1,NULL);
      }
      accept = acceptAlpha(alphaFirstCandidate);

      bool displayDetail = (m_env.displayVerbosity() >= 10/*99*/) || m_optionsObj->m_displayCandidates;
      if ((m_env.subDisplayFile()                   ) &&
          (displayDetail                            ) &&
          (m_optionsObj->m_totallyMute == false)) {
        *m_env.subDisplayFile() << "In SequenceGenerator<V,M>::generateFullChain()"
                                << ": for chain position of id = " << positionId
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

      //****************************************************
      // Point 2/6 of logic for new position
      // Loop: update chain
      //****************************************************
      if (accept) {
        workingChain.setPositionValues(positionId,currentCandidateData.vecValues());
        if (true/*m_uniqueChainGenerate*/) m_idsOfUniquePositions[uniquePos++] = positionId;
        currentPositionData = currentCandidateData;
      }
      else {
        workingChain.setPositionValues(positionId,currentPositionData.vecValues());
      }
      if ((m_optionsObj->m_rawChainDataOutputPeriod                    >  0  ) &&
          (((positionId+1) % m_optionsObj->m_rawChainDataOutputPeriod) == 0  ) &&
          (m_optionsObj->m_rawChainDataOutputFileName                  != ".")) {
        if ((m_env.subDisplayFile()                   ) &&
            (m_env.displayVerbosity()         >= 10   ) &&
            (m_optionsObj->m_totallyMute == false)) {
          *m_env.subDisplayFile() << "In SequenceGenerator<V,M>::generateFullChain()"
                                  << ", for chain position of id = " << positionId
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
          *m_env.subDisplayFile() << "In SequenceGenerator<V,M>::generateFullChain()"
                                  << ", for chain position of id = " << positionId
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

      //****************************************************
      // Point 3/6 of logic for new position
      // Loop: print some information before going to the next chain position
      //****************************************************
      if ((m_env.subDisplayFile()                   ) &&
          (m_env.displayVerbosity() >= 3            ) &&
          (m_optionsObj->m_totallyMute == false)) {
        *m_env.subDisplayFile() << "In SequenceGenerator<V,M>::generateFullChain()"
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
  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_totallyMute == false)) {
    *m_env.subDisplayFile() << "Finished the generation of Markov chain " << workingChain.name()
                            << ", with sub "                              << workingChain.subSequenceSize()
                            << " positions";
    *m_env.subDisplayFile() << std::endl;
  }

  //****************************************************
  // Release memory before leaving routine
  //****************************************************

  return;
}

}  // End namespace QUESO

template class QUESO::SequenceGenerator<QUESO::GslVector, QUESO::GslMatrix>;
