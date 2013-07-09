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

#ifndef __UQ_MH_SG2_H__
#define __UQ_MH_SG2_H__

// Statistical methods -----------------------------
/* This operation currently implements the DRAM algorithm (Heikki Haario, Marko
 * Laine, Antonietta Mira and Eero Saksman, "DRAM: Efficient Adaptive MCMC", 
 * Statistics and Computing (2006), 16:339-354). It also provides support for 
 * Stochastic Newton algorithm through the TK (transition kernel) class. Stochastic
 * Newton is not totally implemented yet though, since it is being researched by
 * James Martin and Omar Ghattas at ICES at the University of Texas at Austin.*/    
template <class P_V,class P_M>
void
uqMetropolisHastingsSGClass<P_V,P_M>::generateSequence(
  uqBaseVectorSequenceClass<P_V,P_M>& workingChain,
  uqScalarSequenceClass<double>*      workingLogLikelihoodValues, // KEY: add LogPriorValues
  uqScalarSequenceClass<double>*      workingLogTargetValues)
{
  if ((m_env.subDisplayFile()                   ) &&
      (m_env.displayVerbosity() >= 5            ) &&
      (m_optionsObj->m_ov.m_totallyMute == false)) {
    *m_env.subDisplayFile() << "Entering uqMetropolisHastingsSGClass<P_V,P_M>::generateSequence()..."
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO(m_vectorSpace.dimLocal() != workingChain.vectorSizeLocal(),
                      m_env.worldRank(),
                      "uqMetropolisHastingsSGClass<P_V,P_M>::generateSequence()",
                      "'m_vectorSpace' and 'workingChain' are related to vector spaces of different dimensions");

  //m_env.syncPrintDebugMsg("Entering uqMetropolisHastingsSGClass<P_V,P_M>::generateSequence()",2,3000000,m_env.fullComm());  // Dangerous to barrier on fullComm ... // KAUST
  uqMiscCheckTheParallelEnvironment<P_V,P_V>(m_initialPosition,
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
    *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateSequence()"
                            << ", prefix = "                                         << m_optionsObj->m_prefix
                            << ", chain name = "                                     << workingChain.name()
                            << ": about to try to open generic output file '"        << m_optionsObj->m_ov.m_dataOutputFileName
                            << "."                                                   << UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT // Yes, always ".m"
                            << "', subId = "                                         << m_env.subId()
                            << ", subenv is allowed to write (1/true or 0/false) = " << (m_optionsObj->m_ov.m_dataOutputAllowedSet.find(m_env.subId()) != m_optionsObj->m_ov.m_dataOutputAllowedSet.end())
                            << "..."
                            << std::endl;
  }

  uqFilePtrSetStruct genericFilePtrSet;
  m_env.openOutputFile(m_optionsObj->m_ov.m_dataOutputFileName,
                       UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT, // Yes, always ".m"
                       m_optionsObj->m_ov.m_dataOutputAllowedSet,
                       false,
                       genericFilePtrSet);

  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_ov.m_totallyMute == false)) {
    *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateSequence()"
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
      *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateSequence()"
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
        *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateSequence()"
                                << ": just wrote (per period request) remaining " << m_numPositionsNotSubWritten << " chain positions "
                                << ", " << m_optionsObj->m_ov.m_rawChainSize - m_numPositionsNotSubWritten << " <= pos <= " << m_optionsObj->m_ov.m_rawChainSize - 1
                                << std::endl;
      }

      if (workingLogLikelihoodValues) {
        workingLogLikelihoodValues->subWriteContents(m_optionsObj->m_ov.m_rawChainSize - m_numPositionsNotSubWritten,
                                                     m_numPositionsNotSubWritten,
                                                     m_optionsObj->m_ov.m_rawChainDataOutputFileName + "_likelihood",
                                                     m_optionsObj->m_ov.m_rawChainDataOutputFileType,
                                                     m_optionsObj->m_ov.m_rawChainDataOutputAllowedSet);
      }

      if (workingLogTargetValues) {
        workingLogTargetValues->subWriteContents(m_optionsObj->m_ov.m_rawChainSize - m_numPositionsNotSubWritten,
                                                 m_numPositionsNotSubWritten,
                                                 m_optionsObj->m_ov.m_rawChainDataOutputFileName + "_target",
                                                 m_optionsObj->m_ov.m_rawChainDataOutputFileType,
                                                 m_optionsObj->m_ov.m_rawChainDataOutputAllowedSet);
      }

      m_numPositionsNotSubWritten = 0;
    }

    // Compute raw sub MLE
    if (workingLogLikelihoodValues) {
      uqSequenceOfVectorsClass<P_V,P_M> rawSubMLEpositions(m_vectorSpace,0,m_optionsObj->m_prefix+"rawSubMLEseq");
      double rawSubMLEvalue = workingChain.subPositionsOfMaximum(*workingLogLikelihoodValues,
                                                                 rawSubMLEpositions);
      UQ_FATAL_TEST_MACRO(rawSubMLEpositions.subSequenceSize() == 0,
                          m_env.worldRank(),
                          "uqMetropolisHastingsSGClass<P_V,P_M>::generateSequence()",
                          "rawSubMLEpositions.subSequenceSize() = 0");

      if ((m_env.subDisplayFile()                   ) &&
          (m_optionsObj->m_ov.m_totallyMute == false)) {
        P_V tmpVec(m_vectorSpace.zeroVector());
        rawSubMLEpositions.getPositionValues(0,tmpVec);
        *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateSequence()"
                                << ": just computed MLE"
                                << ", rawSubMLEvalue = "                       << rawSubMLEvalue
                                << ", rawSubMLEpositions.subSequenceSize() = " << rawSubMLEpositions.subSequenceSize()
                                << ", rawSubMLEpositions[0] = "                << tmpVec
                                << std::endl;
      }
    }

    // Compute raw sub MAP
    if (workingLogTargetValues) {
      uqSequenceOfVectorsClass<P_V,P_M> rawSubMAPpositions(m_vectorSpace,0,m_optionsObj->m_prefix+"rawSubMAPseq");
      double rawSubMAPvalue = workingChain.subPositionsOfMaximum(*workingLogTargetValues,
                                                                 rawSubMAPpositions);
      UQ_FATAL_TEST_MACRO(rawSubMAPpositions.subSequenceSize() == 0,
                          m_env.worldRank(),
                          "uqMetropolisHastingsSGClass<P_V,P_M>::generateSequence()",
                          "rawSubMAPpositions.subSequenceSize() = 0");

      if ((m_env.subDisplayFile()                   ) &&
          (m_optionsObj->m_ov.m_totallyMute == false)) {
        P_V tmpVec(m_vectorSpace.zeroVector());
        rawSubMAPpositions.getPositionValues(0,tmpVec);
        *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateSequence()"
                                << ": just computed MAP"
                                << ", rawSubMAPvalue = "                       << rawSubMAPvalue
                                << ", rawSubMAPpositions.subSequenceSize() = " << rawSubMAPpositions.subSequenceSize()
                                << ", rawSubMAPpositions[0] = "                << tmpVec
                                << std::endl;
      }
    }

    if ((m_env.subDisplayFile()                   ) &&
        (m_optionsObj->m_ov.m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateSequence()"
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
      *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateSequence()"
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
      *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateSequence()"
                              << ", prefix = "                                             << m_optionsObj->m_prefix
                              << ", raw chain name = "                                     << workingChain.name()
                              << ": returned from writing raw unified chain output file '" << m_optionsObj->m_ov.m_rawChainDataOutputFileName
                              << "."                                                       << m_optionsObj->m_ov.m_rawChainDataOutputFileType
                              << "', subId = "                                             << m_env.subId()
                              << std::endl;
    }

    if (workingLogLikelihoodValues) {
      workingLogLikelihoodValues->unifiedWriteContents(m_optionsObj->m_ov.m_rawChainDataOutputFileName + "_likelihood",
                                                       m_optionsObj->m_ov.m_rawChainDataOutputFileType);
    }

    if (workingLogTargetValues) {
      workingLogTargetValues->unifiedWriteContents(m_optionsObj->m_ov.m_rawChainDataOutputFileName + "_target",
                                                   m_optionsObj->m_ov.m_rawChainDataOutputFileType);
    }

    // Compute raw unified MLE
    if (workingLogLikelihoodValues) {
      uqSequenceOfVectorsClass<P_V,P_M> rawUnifiedMLEpositions(m_vectorSpace,0,m_optionsObj->m_prefix+"rawUnifiedMLEseq");
      double rawUnifiedMLEvalue = workingChain.unifiedPositionsOfMaximum(*workingLogLikelihoodValues,
                                                                         rawUnifiedMLEpositions);
      UQ_FATAL_TEST_MACRO(rawUnifiedMLEpositions.subSequenceSize() == 0,
                          m_env.worldRank(),
                          "uqMetropolisHastingsSGClass<P_V,P_M>::generateSequence()",
                          "rawUnifiedMLEpositions.subSequenceSize() = 0");

      if ((m_env.subDisplayFile()                   ) &&
          (m_optionsObj->m_ov.m_totallyMute == false)) {
        P_V tmpVec(m_vectorSpace.zeroVector());
        rawUnifiedMLEpositions.getPositionValues(0,tmpVec);
        *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateSequence()"
                                << ": just computed MLE"
                                << ", rawUnifiedMLEvalue = "                       << rawUnifiedMLEvalue
                                << ", rawUnifiedMLEpositions.subSequenceSize() = " << rawUnifiedMLEpositions.subSequenceSize()
                                << ", rawUnifiedMLEpositions[0] = "                << tmpVec
                                << std::endl;
      }
    }

    // Compute raw unified MAP
    if (workingLogTargetValues) {
      uqSequenceOfVectorsClass<P_V,P_M> rawUnifiedMAPpositions(m_vectorSpace,0,m_optionsObj->m_prefix+"rawUnifiedMAPseq");
      double rawUnifiedMAPvalue = workingChain.unifiedPositionsOfMaximum(*workingLogTargetValues,
                                                                         rawUnifiedMAPpositions);

      UQ_FATAL_TEST_MACRO(rawUnifiedMAPpositions.subSequenceSize() == 0,
                          m_env.worldRank(),
                          "uqMetropolisHastingsSGClass<P_V,P_M>::generateSequence()",
                          "rawUnifiedMAPpositions.subSequenceSize() = 0");

      if ((m_env.subDisplayFile()                   ) &&
          (m_optionsObj->m_ov.m_totallyMute == false)) {
        P_V tmpVec(m_vectorSpace.zeroVector());
        rawUnifiedMAPpositions.getPositionValues(0,tmpVec);
        *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateSequence()"
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
                      "uqMetropolisHastingsSGClass<P_V,P_M>::generateSequence()",
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
      *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateSequence()"
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
        *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateSequence()"
                                << ", prefix = "                << m_optionsObj->m_prefix
                                << ": closed sub output file '" << m_optionsObj->m_ov.m_filteredChainDataOutputFileName
                                << "' for filtered chain "      << workingChain.name()
                                << std::endl;
      }

      if (workingLogLikelihoodValues) {
        workingLogLikelihoodValues->subWriteContents(0,
                                                     workingChain.subSequenceSize(),
                                                     m_optionsObj->m_ov.m_filteredChainDataOutputFileName + "_likelihood",
                                                     m_optionsObj->m_ov.m_filteredChainDataOutputFileType,
                                                     m_optionsObj->m_ov.m_filteredChainDataOutputAllowedSet);
      }

      if (workingLogTargetValues) {
        workingLogTargetValues->subWriteContents(0,
                                                 workingChain.subSequenceSize(),
                                                 m_optionsObj->m_ov.m_filteredChainDataOutputFileName + "_target",
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
        *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateSequence()"
                                << ", prefix = "                    << m_optionsObj->m_prefix
                                << ": closed unified output file '" << m_optionsObj->m_ov.m_filteredChainDataOutputFileName
                                << "' for filtered chain "          << workingChain.name()
                                << std::endl;
      }

      if (workingLogLikelihoodValues) {
        workingLogLikelihoodValues->unifiedWriteContents(m_optionsObj->m_ov.m_filteredChainDataOutputFileName + "_likelihood",
                                                         m_optionsObj->m_ov.m_filteredChainDataOutputFileType);
      }

      if (workingLogTargetValues) {
        workingLogTargetValues->unifiedWriteContents(m_optionsObj->m_ov.m_filteredChainDataOutputFileName + "_target",
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
      *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateSequence()"
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

  //m_env.syncPrintDebugMsg("Leaving uqMetropolisHastingsSGClass<P_V,P_M>::generateSequence()",2,3000000,m_env.fullComm()); // Dangerous to barrier on fullComm ... // KAUST

  if ((m_env.subDisplayFile()                   ) &&
      (m_env.displayVerbosity() >= 5            ) &&
      (m_optionsObj->m_ov.m_totallyMute == false)) {
    *m_env.subDisplayFile() << "Leaving uqMetropolisHastingsSGClass<P_V,P_M>::generateSequence()"
                            << std::endl;
  }

  return;
}

// -------------------------------------------------
template<class P_V,class P_M>
void
uqMetropolisHastingsSGClass<P_V,P_M>::getRawChainInfo(uqMHRawChainInfoStruct& info) const
{
  info = m_rawChainInfo;
  return;
}
//--------------------------------------------------
template <class P_V,class P_M>
void
uqMetropolisHastingsSGClass<P_V,P_M>::readFullChain(
  const std::string&                  inputFileName,
  const std::string&                  inputFileType,
        unsigned int                  chainSize,
  uqBaseVectorSequenceClass<P_V,P_M>& workingChain)
{
  workingChain.unifiedReadContents(inputFileName,inputFileType,chainSize);
  return;
}
// Private methods ---------------------------------
template <class P_V,class P_M>
void
uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain(
  const P_V&                          valuesOf1stPosition,
        unsigned int                  chainSize,
  uqBaseVectorSequenceClass<P_V,P_M>& workingChain,
  uqScalarSequenceClass<double>*      workingLogLikelihoodValues,
  uqScalarSequenceClass<double>*      workingLogTargetValues)
{
  //m_env.syncPrintDebugMsg("Entering uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain()",3,3000000,m_env.fullComm()); // Dangerous to barrier on fullComm ... // KAUST

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
    *m_env.subDisplayFile() << "\nIn uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain()"
                            << ": contents of initial position are:";
    *m_env.subDisplayFile() << valuesOf1stPosition; // FIX ME: might need parallelism
    *m_env.subDisplayFile() << "\nIn uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain()"
                            << ": targetPdf.domaintSet() info is:"
                            << m_targetPdf.domainSet();
    *m_env.subDisplayFile() << std::endl;
  }

  bool outOfTargetSupport = !m_targetPdf.domainSet().contains(valuesOf1stPosition);
  if ((m_env.subDisplayFile()) &&
      (outOfTargetSupport    )) {
    *m_env.subDisplayFile() << "ERROR: In uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain()"
                            << ": contents of initial position are:\n";
    *m_env.subDisplayFile() << valuesOf1stPosition; // FIX ME: might need parallelism
    *m_env.subDisplayFile() << "\nERROR: In uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain()"
                            << ": targetPdf.domaintSet() info is:\n"
                            << m_targetPdf.domainSet();
    *m_env.subDisplayFile() << std::endl;
  }
  UQ_FATAL_TEST_MACRO(outOfTargetSupport,
                      m_env.worldRank(),
                      "uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain()",
                      "initial position should not be out of target pdf support");
  if (m_optionsObj->m_ov.m_rawChainMeasureRunTimes) iRC = gettimeofday(&timevalTarget, NULL);
  double logPrior      = 0.;
  double logLikelihood = 0.;
#ifdef QUESO_EXPECTS_LN_LIKELIHOOD_INSTEAD_OF_MINUS_2_LN
  double logTarget =        m_targetPdfSynchronizer->callFunction(&valuesOf1stPosition,NULL,NULL,NULL,NULL,&logPrior,&logLikelihood); // Might demand parallel environment // KEY
#else
  double logTarget = -0.5 * m_targetPdfSynchronizer->callFunction(&valuesOf1stPosition,NULL,NULL,NULL,NULL,&logPrior,&logLikelihood); // Might demand parallel environment
#endif
  if (m_optionsObj->m_ov.m_rawChainMeasureRunTimes) m_rawChainInfo.targetRunTime += uqMiscGetEllapsedSeconds(&timevalTarget);
  m_rawChainInfo.numTargetCalls++;
  if ((m_env.subDisplayFile()                   ) &&
      (m_env.displayVerbosity() >= 3            ) &&
      (m_optionsObj->m_ov.m_totallyMute == false)) {
    *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain()"
                            << ": just returned from likelihood() for initial chain position"
                            << ", m_rawChainInfo.numTargetCalls = " << m_rawChainInfo.numTargetCalls
                            << ", logPrior = "      << logPrior
                            << ", logLikelihood = " << logLikelihood
                            << ", logTarget = "     << logTarget
                            << std::endl;
  }

  //*m_env.subDisplayFile() << "AQUI 001" << std::endl;
  uqMarkovChainPositionDataClass<P_V> currentPositionData(m_env,
                                                          valuesOf1stPosition,
                                                          outOfTargetSupport,
                                                          logLikelihood,
                                                          logTarget);

  P_V gaussianVector(m_vectorSpace.zeroVector());
  P_V tmpVecValues(m_vectorSpace.zeroVector());
  uqMarkovChainPositionDataClass<P_V> currentCandidateData(m_env);

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
      *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain()"
                              << ": just wrote (per period request) " << m_numPositionsNotSubWritten << " chain positions "
                              << ", " << 0 + 1 - m_optionsObj->m_ov.m_rawChainDataOutputPeriod << " <= pos <= " << 0
                              << std::endl;
    }

    if (workingLogLikelihoodValues) {
      workingLogLikelihoodValues->subWriteContents(0 + 1 - m_optionsObj->m_ov.m_rawChainDataOutputPeriod,
                                                   m_optionsObj->m_ov.m_rawChainDataOutputPeriod, 
                                                   m_optionsObj->m_ov.m_rawChainDataOutputFileName + "_likelihood",
                                                   m_optionsObj->m_ov.m_rawChainDataOutputFileType,
                                                   m_optionsObj->m_ov.m_rawChainDataOutputAllowedSet);
    }

    if (workingLogTargetValues) {
      workingLogTargetValues->subWriteContents(0 + 1 - m_optionsObj->m_ov.m_rawChainDataOutputPeriod,
                                               m_optionsObj->m_ov.m_rawChainDataOutputPeriod, 
                                               m_optionsObj->m_ov.m_rawChainDataOutputFileName + "_target",
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

  //m_env.syncPrintDebugMsg("In uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain(), right before main loop",3,3000000,m_env.fullComm()); // Dangerous to barrier on fullComm ... // KAUST

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
      *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain()"
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
      *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain()"
                              << ": about to set TK pre computing position of local id " << 0
                              << ", values = " << currentPositionData.vecValues()
                              << std::endl;
    }
    bool validPreComputingPosition = m_tk->setPreComputingPosition(currentPositionData.vecValues(),0);
    if ((m_env.subDisplayFile()                   ) &&
        (m_env.displayVerbosity() >= 5            ) &&
        (m_optionsObj->m_ov.m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain()"
                              << ": returned from setting TK pre computing position of local id " << 0
                              << ", values = " << currentPositionData.vecValues()
                              << ", valid = "  << validPreComputingPosition
                              << std::endl;
    }
    UQ_FATAL_TEST_MACRO(validPreComputingPosition == false,
                        m_env.worldRank(),
                        "uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain()",
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
      if (m_optionsObj->m_ov.m_rawChainMeasureRunTimes) m_rawChainInfo.candidateRunTime += uqMiscGetEllapsedSeconds(&timevalCandidate);

      outOfTargetSupport = !m_targetPdf.domainSet().contains(tmpVecValues);

      bool displayDetail = (m_env.displayVerbosity() >= 10/*99*/) || m_optionsObj->m_ov.m_displayCandidates;
      if ((m_env.subDisplayFile()                   ) &&
          (displayDetail                            ) &&
          (m_optionsObj->m_ov.m_totallyMute == false)) {
        *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain()"
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
      *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain()"
                              << ": about to set TK pre computing position of local id " << stageId+1
                              << ", values = " << tmpVecValues
                              << std::endl;
    }
    validPreComputingPosition = m_tk->setPreComputingPosition(tmpVecValues,stageId+1);
    if ((m_env.subDisplayFile()                   ) &&
        (m_env.displayVerbosity() >= 5            ) &&
        (m_optionsObj->m_ov.m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain()"
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
      if (m_optionsObj->m_ov.m_rawChainMeasureRunTimes) m_rawChainInfo.targetRunTime += uqMiscGetEllapsedSeconds(&timevalTarget);
      m_rawChainInfo.numTargetCalls++;
      if ((m_env.subDisplayFile()                   ) &&
          (m_env.displayVerbosity() >= 3            ) &&
          (m_optionsObj->m_ov.m_totallyMute == false)) {
        *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain()"
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
      if (m_optionsObj->m_ov.m_rawChainMeasureRunTimes) m_rawChainInfo.mhAlphaRunTime += uqMiscGetEllapsedSeconds(&timevalMhAlpha);
      if ((m_env.subDisplayFile()                   ) &&
          (m_env.displayVerbosity() >= 10           ) &&
          (m_optionsObj->m_ov.m_totallyMute == false)) {
        *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain()"
                                << ": for chain position of id = " << positionId
                                << std::endl;
      }
      accept = acceptAlpha(alphaFirstCandidate);
    }

    bool displayDetail = (m_env.displayVerbosity() >= 10/*99*/) || m_optionsObj->m_ov.m_displayCandidates;
    if ((m_env.subDisplayFile()                   ) &&
        (displayDetail                            ) &&
        (m_optionsObj->m_ov.m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain()"
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
    std::vector<uqMarkovChainPositionDataClass<P_V>*> drPositionsData(stageId+2,NULL);
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

        drPositionsData[0] = new uqMarkovChainPositionDataClass<P_V>(currentPositionData );
        drPositionsData[1] = new uqMarkovChainPositionDataClass<P_V>(currentCandidateData);

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
            *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain()"
                                    << ": for chain position of id = " << positionId
                                    << ", beginning stageId = "        << stageId
                                    << std::endl;
          }

          keepGeneratingCandidates = true;
          while (keepGeneratingCandidates) {
            if (m_optionsObj->m_ov.m_rawChainMeasureRunTimes) iRC = gettimeofday(&timevalCandidate, NULL);
            m_tk->rv(tkStageIds).realizer().realization(tmpVecValues);
            if (m_optionsObj->m_ov.m_rawChainMeasureRunTimes) m_rawChainInfo.candidateRunTime += uqMiscGetEllapsedSeconds(&timevalCandidate);

            outOfTargetSupport = !m_targetPdf.domainSet().contains(tmpVecValues);

            if (m_optionsObj->m_ov.m_putOutOfBoundsInChain) keepGeneratingCandidates = false;
            else                                            keepGeneratingCandidates = outOfTargetSupport;
          }
  
          if ((m_env.subDisplayFile()                   ) &&
              (m_env.displayVerbosity() >= 5            ) &&
              (m_optionsObj->m_ov.m_totallyMute == false)) {
            *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain()"
                                    << ": about to set TK pre computing position of local id " << stageId+1
                                    << ", values = " << tmpVecValues
                                    << std::endl;
          }
          validPreComputingPosition = m_tk->setPreComputingPosition(tmpVecValues,stageId+1);
          if ((m_env.subDisplayFile()                   ) &&
              (m_env.displayVerbosity() >= 5            ) &&
              (m_optionsObj->m_ov.m_totallyMute == false)) {
            *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain()"
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
            if (m_optionsObj->m_ov.m_rawChainMeasureRunTimes) m_rawChainInfo.targetRunTime += uqMiscGetEllapsedSeconds(&timevalTarget);
            m_rawChainInfo.numTargetCalls++;
            if ((m_env.subDisplayFile()                   ) &&
                (m_env.displayVerbosity() >= 3            ) &&
                (m_optionsObj->m_ov.m_totallyMute == false)) {
              *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain()"
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

          drPositionsData.push_back(new uqMarkovChainPositionDataClass<P_V>(currentCandidateData));
          tkStageIds.push_back     (stageId+1);

          double alphaDR = 0.;
          if (outOfTargetSupport == false) {
            if (m_optionsObj->m_ov.m_rawChainMeasureRunTimes) iRC = gettimeofday(&timevalDrAlpha, NULL);
            alphaDR = this->alpha(drPositionsData,tkStageIds);
            if (m_optionsObj->m_ov.m_rawChainMeasureRunTimes) m_rawChainInfo.drAlphaRunTime += uqMiscGetEllapsedSeconds(&timevalDrAlpha);
            accept = acceptAlpha(alphaDR);
          }

          displayDetail = (m_env.displayVerbosity() >= 10/*99*/) || m_optionsObj->m_ov.m_displayCandidates;
          if ((m_env.subDisplayFile()                   ) &&
              (displayDetail                            ) &&
              (m_optionsObj->m_ov.m_totallyMute == false)) {
            *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain()"
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

        if (m_optionsObj->m_ov.m_rawChainMeasureRunTimes) m_rawChainInfo.drRunTime += uqMiscGetEllapsedSeconds(&timevalDR);
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
        *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain()"
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
        *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain()"
                                << ", for chain position of id = " << positionId
                                << ": just wrote (per period request) " << m_numPositionsNotSubWritten << " chain positions "
                                << ", " << positionId + 1 - m_optionsObj->m_ov.m_rawChainDataOutputPeriod << " <= pos <= " << positionId
                                << std::endl;
      }

      if (workingLogLikelihoodValues) {
        workingLogLikelihoodValues->subWriteContents(0 + 1 - m_optionsObj->m_ov.m_rawChainDataOutputPeriod,
                                                     m_optionsObj->m_ov.m_rawChainDataOutputPeriod, 
                                                     m_optionsObj->m_ov.m_rawChainDataOutputFileName + "_likelihood",
                                                     m_optionsObj->m_ov.m_rawChainDataOutputFileType,
                                                     m_optionsObj->m_ov.m_rawChainDataOutputAllowedSet);
      }

      if (workingLogTargetValues) {
        workingLogTargetValues->subWriteContents(0 + 1 - m_optionsObj->m_ov.m_rawChainDataOutputPeriod,
                                                 m_optionsObj->m_ov.m_rawChainDataOutputPeriod, 
                                                 m_optionsObj->m_ov.m_rawChainDataOutputFileName + "_target",
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
      uqSequenceOfVectorsClass<P_V,P_M> partialChain(m_vectorSpace,0,m_optionsObj->m_prefix+"partialChain");

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
            *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain()"
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
          *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain()"
                                  << ", positionId = "  << positionId
                                  << ": 'am' calling first tmpChol.chol()"
                                  << std::endl;
        }
        iRC = tmpChol.chol();
        if (iRC) {
          std::cerr << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain(): first chol failed\n";
        }
        if ((m_env.subDisplayFile()        ) &&
            (m_env.displayVerbosity() >= 10)) {
	  //(m_optionsObj->m_ov.m_totallyMute == false)) {
          *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain()"
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
                              "uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain()",
                              "invalid iRC returned from first chol()");
          // Matrix is not positive definite
          P_M* tmpDiag = m_vectorSpace.newDiagMatrix(m_optionsObj->m_ov.m_amEpsilon);
          tmpChol = *m_lastAdaptedCovMatrix + *tmpDiag;
          attemptedMatrix = tmpChol;
          delete tmpDiag;
          if ((m_env.subDisplayFile()        ) &&
              (m_env.displayVerbosity() >= 10)) {
	    //(m_optionsObj->m_ov.m_totallyMute == false)) {
            *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain()"
                                    << ", positionId = "  << positionId
                                    << ": 'am' calling second tmpChol.chol()"
                                    << std::endl;
          }
          iRC = tmpChol.chol();
          if (iRC) {
            std::cerr << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain(): second chol failed\n";
          }
          if ((m_env.subDisplayFile()        ) &&
              (m_env.displayVerbosity() >= 10)) {
	    //(m_optionsObj->m_ov.m_totallyMute == false)) {
            *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain()"
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
                                "uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain()",
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
          uqScaledCovMatrixTKGroupClass<P_V,P_M>* tempTK = dynamic_cast<uqScaledCovMatrixTKGroupClass<P_V,P_M>* >(m_tk);
          tempTK->updateLawCovMatrix(m_optionsObj->m_ov.m_amEta*attemptedMatrix);

#ifdef UQ_DRAM_MCG_REQUIRES_INVERTED_COV_MATRICES
          UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                            m_env.worldRank(),
                            "uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain()",
                            "need to code the update of m_upperCholProposalPrecMatrices");
#endif
        }

        //for (unsigned int i = 0; i < partialChain.subSequenceSize(); ++i) {
        //  if (partialChain[i]) delete partialChain[i];
        //}
      } // if (partialChain.subSequenceSize() > 0)

      if (m_optionsObj->m_ov.m_rawChainMeasureRunTimes) m_rawChainInfo.amRunTime += uqMiscGetEllapsedSeconds(&timevalAM);
    } // End of 'adaptive Metropolis' logic

    //****************************************************
    // Point 6/6 of logic for new position
    // Loop: print some information before going to the next chain position
    //****************************************************
    if ((m_env.subDisplayFile()                   ) &&
        (m_env.displayVerbosity() >= 3            ) &&
        (m_optionsObj->m_ov.m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain()"
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
                                << " positions"
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
  m_rawChainInfo.runTime += uqMiscGetEllapsedSeconds(&timevalChain);
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
  //m_env.syncPrintDebugMsg("Leaving uqMetropolisHastingsSGClass<P_V,P_M>::generateFullChain()",3,3000000,m_env.fullComm()); // Dangerous to barrier on fullComm ... // KAUST

  return;
}
//--------------------------------------------------
template <class P_V,class P_M>
void
uqMetropolisHastingsSGClass<P_V,P_M>::updateAdaptedCovMatrix(
  const uqBaseVectorSequenceClass<P_V,P_M>& partialChain,
  unsigned int                              idOfFirstPositionInSubChain,
  double&                                   lastChainSize,
  P_V&                                      lastMean,
  P_M&                                      lastAdaptedCovMatrix)
{
  double doubleSubChainSize = (double) partialChain.subSequenceSize();
  if (lastChainSize == 0) {
    UQ_FATAL_TEST_MACRO(partialChain.subSequenceSize() < 2,
                        m_env.worldRank(),
                        "uqMetropolisHastingsSGClass<P_V,P_M>::updateAdaptedCovMatrix()",
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
                        "uqMetropolisHastingsSGClass<P_V,P_M>::updateAdaptedCovMatrix()",
                        "'partialChain.subSequenceSize()' should be >= 1");

    UQ_FATAL_TEST_MACRO(idOfFirstPositionInSubChain < 1,
                        m_env.worldRank(),
                        "uqMetropolisHastingsSGClass<P_V,P_M>::updateAdaptedCovMatrix()",
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

  return;
}
#endif // __UQ_MH_SG2_H__
