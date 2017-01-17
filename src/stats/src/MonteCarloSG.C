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

#include <queso/MonteCarloSG.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Default constructor -----------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
MonteCarloSG<P_V,P_M,Q_V,Q_M>::MonteCarloSG(
  /*! Prefix                     */ const char*                                       prefix,
  /*! Options (if no input file) */ const McOptionsValues*                     alternativeOptionsValues, // dakota
  /*! The parameter RV           */ const BaseVectorRV      <P_V,P_M>&         paramRv,
  /*! The QoI function           */ const BaseVectorFunction<P_V,P_M,Q_V,Q_M>& qoiFunction)
  :
  m_env                     (paramRv.env()),
  m_paramRv                 (paramRv),
  m_qoiFunction             (qoiFunction),
  m_paramSpace              (m_paramRv.imageSet().vectorSpace()),
  m_qoiSpace                (m_qoiFunction.imageSet().vectorSpace()),
  m_qoiFunctionSynchronizer (new VectorFunctionSynchronizer<P_V,P_M,Q_V,Q_M>(m_qoiFunction,m_paramRv.imageSet().vectorSpace().zeroVector(),m_qoiFunction.imageSet().vectorSpace().zeroVector())),
  m_numPsNotSubWritten      (0),
  m_numQsNotSubWritten      (0),
  m_optionsObj              (alternativeOptionsValues),
  m_userDidNotProvideOptions(false)
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Entering MonteCarloSG<P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": prefix = " << prefix
                            << ", alternativeOptionsValues = " << alternativeOptionsValues
                            << ", m_env.optionsInputFileName() = " << m_env.optionsInputFileName()
                            << std::endl;
  }

  // If NULL, we create one
  if (m_optionsObj == NULL) {
    McOptionsValues * tempOptions = new McOptionsValues(&m_env, prefix);

    // We did this dance because scanOptionsValues is not a const method, but
    // m_optionsObj is a pointer to const
    m_optionsObj = tempOptions;

    // We do this so we don't delete the user's object in the dtor
    m_userDidNotProvideOptions = true;
  }

  if (m_optionsObj->m_help != "") {
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << (*m_optionsObj) << std::endl;
    }
  }

  queso_require_equal_to_msg(paramRv.imageSet().vectorSpace().dimLocal(), qoiFunction.domainSet().vectorSpace().dimLocal(), "'paramRv' and 'qoiFunction' are related to vector spaces of different dimensions");

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Leaving MonteCarloSG<P_V,P_M,Q_V,Q_M>::constructor()"
                            << std::endl;
  }
}
// Destructor ---------------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
MonteCarloSG<P_V,P_M,Q_V,Q_M>::~MonteCarloSG()
{
  if (m_optionsObj && m_userDidNotProvideOptions) {
    delete m_optionsObj;
  }

  if (m_qoiFunctionSynchronizer) delete m_qoiFunctionSynchronizer;
}
// Statistical methods ------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
void
MonteCarloSG<P_V,P_M,Q_V,Q_M>::generateSequence(
  BaseVectorSequence<P_V,P_M>& workingPSeq,
  BaseVectorSequence<Q_V,Q_M>& workingQSeq)
{
  queso_require_equal_to_msg(m_qoiFunction.domainSet().vectorSpace().dimLocal(), workingPSeq.vectorSizeLocal(), "'m_qoiFunction.domainSet' and 'workingPSeq' are related to vector spaces of different dimensions");

  queso_require_equal_to_msg(m_qoiFunction.imageSet().vectorSpace().dimLocal(), workingQSeq.vectorSizeLocal(), "'m_qoiFunction.imageSet' and 'workingQSeq' are related to vector spaces of different dimensions");

  MiscCheckTheParallelEnvironment<P_V,Q_V>(m_paramRv.imageSet().vectorSpace().zeroVector(),
                                             m_qoiFunction.imageSet().vectorSpace().zeroVector());
  internGenerateSequence(m_paramRv,workingPSeq,workingQSeq);

  return;
}
// I/O methods---------------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
void
MonteCarloSG<P_V,P_M,Q_V,Q_M>::print(std::ostream& os) const
{
  return;
}
// Private methods----------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
void
MonteCarloSG<P_V,P_M,Q_V,Q_M>::internGenerateSequence(
  const BaseVectorRV      <P_V,P_M>& paramRv,
        BaseVectorSequence<P_V,P_M>& workingPSeq,
        BaseVectorSequence<Q_V,Q_M>& workingQSeq)
{
  workingPSeq.setName(m_optionsObj->m_prefix+"ParamSeq");
  workingQSeq.setName(m_optionsObj->m_prefix+"QoiSeq");

  //****************************************************
  // Generate sequence of QoI values
  //****************************************************
  unsigned int subActualSizeBeforeGeneration = std::min(m_optionsObj->m_qseqSize,paramRv.realizer().subPeriod());
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "In MonteCarloSG<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                            << ": m_optionsObj->m_qseqSize = "                             << m_optionsObj->m_qseqSize
                            << ", paramRv.realizer().subPeriod() = "                            << paramRv.realizer().subPeriod()
                            << ", about to call actualGenerateSequence() with subActualSize = " << subActualSizeBeforeGeneration
                            << std::endl;
  }
  if (m_optionsObj->m_qseqDataInputFileName == UQ_MOC_SG_FILENAME_FOR_NO_FILE) {
    actualGenerateSequence(paramRv,
                           workingPSeq,
                           workingQSeq,
                           subActualSizeBeforeGeneration);
  }
  else {
    actualReadSequence(paramRv,
                       m_optionsObj->m_qseqDataInputFileName,
                       m_optionsObj->m_qseqDataInputFileType,
                       workingPSeq,
                       workingQSeq,
                       subActualSizeBeforeGeneration);
  }
  unsigned int subActualSizeAfterGeneration = workingPSeq.subSequenceSize();
  queso_require_equal_to_msg(subActualSizeAfterGeneration, workingQSeq.subSequenceSize(), "P and Q sequences should have the same size!");

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "In MonteCarloSG<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                            << ": returned from call to actualGenerateSequence() with subActualSize = " << subActualSizeAfterGeneration
                            << std::endl;
  }

  //****************************************************
  // Open generic output file
  //****************************************************
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "In MonteCarloSG<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                            << ", prefix = "                                                        << m_optionsObj->m_prefix
                            << ": checking necessity of opening generic output file (qseq name is " << workingQSeq.name()
                            << ") ..."
                            << std::endl;
  }
  FilePtrSetStruct genericFilePtrSet;
  m_env.openOutputFile(m_optionsObj->m_dataOutputFileName,
                       UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT, // Yes, always ".m"
                       m_optionsObj->m_dataOutputAllowedSet,
                       false,
                       genericFilePtrSet);

  //****************************************************
  // Eventually:
  // --> write parameter sequence
  // --> compute statistics on it
  //****************************************************
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "In MonteCarloSG<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                            << ", prefix = "                                            << m_optionsObj->m_prefix
                            << ": checking necessity of opening output files for pseq " << workingPSeq.name()
                            << "..."
                            << std::endl;
  }

  // Take "sub" care of pseq
  if ((m_numPsNotSubWritten                        >  0                             ) &&
      (m_optionsObj->m_pseqDataOutputFileName != UQ_MOC_SG_FILENAME_FOR_NO_FILE)) {
    workingPSeq.subWriteContents(subActualSizeBeforeGeneration - m_numPsNotSubWritten,
                                 m_numPsNotSubWritten,
                                 m_optionsObj->m_pseqDataOutputFileName,
                                 m_optionsObj->m_pseqDataOutputFileType,
                                 m_optionsObj->m_pseqDataOutputAllowedSet);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In MonteCarloG<P_V,P_M>::internGenerateSequence()"
                              << ": just wrote remaining pseq positions (per period request)"
                              << std::endl;
    }
    m_numPsNotSubWritten = 0;
    //if (m_env.subDisplayFile()) {
    //  *m_env.subDisplayFile() << "In MonteCarloSG<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
    //                          << ", prefix = "                 << m_optionsObj->m_prefix
    //                          << ": closed data output file '" << m_optionsObj->m_pseqDataOutputFileName
    //                          << "' for pseq "                 << workingPSeq.name()
    //                          << std::endl;
    //}
  }

  // Take "unified" care of pseq
  if (m_optionsObj->m_pseqDataOutputFileName != UQ_MOC_SG_FILENAME_FOR_NO_FILE) {
    workingPSeq.unifiedWriteContents(m_optionsObj->m_pseqDataOutputFileName,m_optionsObj->m_pseqDataOutputFileType);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In MonteCarloSG<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                              << ", prefix = "                         << m_optionsObj->m_prefix
                              << ": closed unified data output file '" << m_optionsObj->m_pseqDataOutputFileName
                              << "' for pseq "                         << workingPSeq.name()
                              << std::endl;
    }
  }

  // Take case of other aspects of pseq
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  if (m_optionsObj->m_pseqComputeStats) {
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In MonteCarloSG<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                              << ": about to call 'workingPSeq.computeStatistics()'"
                              << std::endl;
    }
    workingPSeq.computeStatistics(*m_optionsObj->m_pseqStatisticalOptionsObj,
                                  genericFilePtrSet.ofsVar);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In MonteCarloSG<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                              << ": returned from call to 'workingPSeq.computeStatistics()'"
                              << std::endl;
    }
  }
#endif
  //****************************************************
  // Eventually:
  // --> write QoI sequence
  // --> compute statistics on it
  //****************************************************
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "In MonteCarloSG<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                            << ", prefix = "                                            << m_optionsObj->m_prefix
                            << ": checking necessity of opening output files for qseq " << workingQSeq.name()
                            << "..."
                            << std::endl;
  }

  // Take "sub" care of qseq
  if ((m_numQsNotSubWritten                        >  0                             ) &&
      (m_optionsObj->m_qseqDataOutputFileName != UQ_MOC_SG_FILENAME_FOR_NO_FILE)) {
    workingQSeq.subWriteContents(subActualSizeBeforeGeneration - m_numQsNotSubWritten,
                                 m_numQsNotSubWritten,
                                 m_optionsObj->m_qseqDataOutputFileName,
                                 m_optionsObj->m_qseqDataOutputFileType,
                                 m_optionsObj->m_qseqDataOutputAllowedSet);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In MonteCarloG<P_V,P_M>::internGenerateSequence()"
                              << ": just wrote remaining qseq positions (per period request)"
                              << std::endl;
    }
    m_numQsNotSubWritten = 0;
    //if (m_env.subDisplayFile()) {
    //  *m_env.subDisplayFile() << "In MonteCarloSG<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
    //                          << ", prefix = "                 << m_optionsObj->m_prefix
    //                          << ": closed data output file '" << m_optionsObj->m_qseqDataOutputFileName
    //                          << "' for qseq "                 << workingQSeq.name()
    //                          << std::endl;
    //}
  }

  // Take "unified" care of qseq
  if (m_optionsObj->m_qseqDataOutputFileName != UQ_MOC_SG_FILENAME_FOR_NO_FILE) {
    workingQSeq.unifiedWriteContents(m_optionsObj->m_qseqDataOutputFileName,m_optionsObj->m_qseqDataOutputFileType);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In MonteCarloSG<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                              << ", prefix = "                         << m_optionsObj->m_prefix
                              << ": closed unified data output file '" << m_optionsObj->m_qseqDataOutputFileName
                              << "' for qseq "                         << workingQSeq.name()
                              << std::endl;
    }
  }

  // Take case of other aspects of qseq
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  if (m_optionsObj->m_qseqComputeStats) {
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In MonteCarloSG<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                              << ": about to call 'workingQSeq.computeStatistics()'"
                              << std::endl;
    }
    workingQSeq.computeStatistics(*m_optionsObj->m_qseqStatisticalOptionsObj,
                                  genericFilePtrSet.ofsVar);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In MonteCarloSG<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                              << ": returned from call to 'workingQSeq.computeStatistics()'"
                              << std::endl;
    }
  }
#endif
  //****************************************************
  // Close generic output file
  //****************************************************
  if (genericFilePtrSet.ofsVar) {
    //std::cout << "TODAY 000" << std::endl;
    delete genericFilePtrSet.ofsVar;
    //genericFilePtrSet.ofsVar->close();
    //std::cout << "TODAY 001" << std::endl;
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In MonteCarloSG<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                              << ", prefix = "                         << m_optionsObj->m_prefix
                              << ": closed generic data output file '" << m_optionsObj->m_dataOutputFileName
                              << "' for QoI sequence "                 << workingQSeq.name()
                              << std::endl;
    }
  }
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << std::endl;
  }

  return;
}
// --------------------------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
void
MonteCarloSG<P_V,P_M,Q_V,Q_M>::actualGenerateSequence(
  const BaseVectorRV      <P_V,P_M>& paramRv,
        BaseVectorSequence<P_V,P_M>& workingPSeq,
        BaseVectorSequence<Q_V,Q_M>& workingQSeq,
        unsigned int                        requestedSeqSize)
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Starting the generation of qoi sequence " << workingQSeq.name()
                            << ", with "                                  << requestedSeqSize
                            << " samples..."
                            << std::endl;
  }

  int iRC = UQ_OK_RC;
  struct timeval timevalSeq;
  struct timeval timevalQoIFunction;

  double seqRunTime         = 0;
  double qoiFunctionRunTime = 0;

  iRC = gettimeofday(&timevalSeq, NULL);
  if (iRC) {}; // just to remover compiler warning

  workingPSeq.resizeSequence(requestedSeqSize);
  m_numPsNotSubWritten = 0;
  workingQSeq.resizeSequence(requestedSeqSize);
  m_numQsNotSubWritten = 0;

  P_V tmpP(m_paramSpace.zeroVector());
  Q_V tmpQ(m_qoiSpace.zeroVector());

  unsigned int actualSeqSize = 0;
  for (unsigned int i = 0; i < requestedSeqSize; ++i) {
    paramRv.realizer().realization(tmpP);

    if (m_optionsObj->m_qseqMeasureRunTimes) iRC = gettimeofday(&timevalQoIFunction, NULL);
    m_qoiFunctionSynchronizer->callFunction(&tmpP,NULL,&tmpQ,NULL,NULL,NULL); // Might demand parallel environment
    if (m_optionsObj->m_qseqMeasureRunTimes) qoiFunctionRunTime += MiscGetEllapsedSeconds(&timevalQoIFunction);

    bool allQsAreFinite = true;
    for (unsigned int j = 0; j < tmpQ.sizeLocal(); ++j) {
      if ((tmpQ[j] == INFINITY) || (tmpQ[j] == -INFINITY)) {
  std::cerr << "WARNING In MonteCarloSG<P_V,P_M,Q_V,Q_M>::actualGenerateSequence()"
                  << ", worldRank "      << m_env.worldRank()
                  << ", fullRank "       << m_env.fullRank()
                  << ", subEnvironment " << m_env.subId()
                  << ", subRank "        << m_env.subRank()
                  << ", inter0Rank "     << m_env.inter0Rank()
                  << ": i = "            << i
                  << ", tmpQ[" << j << "] = " << tmpQ[j]
                  << ", tmpP = "         << tmpP
                  << ", tmpQ = "         << tmpQ
                  << std::endl;
        allQsAreFinite = false;

        if (i > 0) {
          workingPSeq.getPositionValues(i-1,tmpP); // FIXME: temporary code
          workingQSeq.getPositionValues(i-1,tmpQ); // FIXME: temporary code
        }

        break;
      }
    }
    if (allQsAreFinite) {}; // just to remover compiler warning

    //if (allQsAreFinite) { // FIXME: this will cause different processors to have sequences of different sizes
      workingPSeq.setPositionValues(i,tmpP);
      m_numPsNotSubWritten++;
      if ((m_optionsObj->m_pseqDataOutputPeriod           >  0  ) &&
          (((i+1) % m_optionsObj->m_pseqDataOutputPeriod) == 0  ) &&
          (m_optionsObj->m_pseqDataOutputFileName         != ".")) {
        workingPSeq.subWriteContents(i + 1 - m_optionsObj->m_pseqDataOutputPeriod,
                                     m_optionsObj->m_pseqDataOutputPeriod,
                                     m_optionsObj->m_pseqDataOutputFileName,
                                     m_optionsObj->m_pseqDataOutputFileType,
                                     m_optionsObj->m_pseqDataOutputAllowedSet);
        if (m_env.subDisplayFile()) {
          *m_env.subDisplayFile() << "In MonteCarloG<P_V,P_M>::actualGenerateSequence()"
                                  << ": just wrote pseq positions (per period request)"
                                  << std::endl;
        }
        m_numPsNotSubWritten = 0;
      }

      workingQSeq.setPositionValues(i,tmpQ);
      m_numQsNotSubWritten++;
      if ((m_optionsObj->m_qseqDataOutputPeriod           >  0  ) &&
          (((i+1) % m_optionsObj->m_qseqDataOutputPeriod) == 0  ) &&
          (m_optionsObj->m_qseqDataOutputFileName         != ".")) {
        workingQSeq.subWriteContents(i + 1 - m_optionsObj->m_qseqDataOutputPeriod,
                                     m_optionsObj->m_qseqDataOutputPeriod,
                                     m_optionsObj->m_qseqDataOutputFileName,
                                     m_optionsObj->m_qseqDataOutputFileType,
                                     m_optionsObj->m_qseqDataOutputAllowedSet);
        if (m_env.subDisplayFile()) {
          *m_env.subDisplayFile() << "In MonteCarloG<P_V,P_M>::actualGenerateSequence()"
                                  << ": just wrote qseq positions (per period request)"
                                  << std::endl;
        }
        m_numQsNotSubWritten = 0;
      }

      actualSeqSize++;

    //}

    if ((m_optionsObj->m_qseqDisplayPeriod            > 0) &&
        (((i+1) % m_optionsObj->m_qseqDisplayPeriod) == 0)) {
      if (m_env.subDisplayFile()) {
        *m_env.subDisplayFile() << "Finished generating " << i+1
                                << " qoi samples"
                                << std::endl;
      }
    }
  }

  //if (actualSeqSize != requestedSeqSize) {
  //  workingPSeq.resizeSequence(actualSeqSize);
  //  workingQSeq.resizeSequence(actualSeqSize);
  //}

  seqRunTime = MiscGetEllapsedSeconds(&timevalSeq);

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Finished the generation of qoi sequence " << workingQSeq.name()
                            << ", with sub "                              << workingQSeq.subSequenceSize()
                            << " samples"
                            << "\nSome information about this sequence:"
                            << "\n  Sequence run time = " << seqRunTime
                            << " seconds"
                            << "\n\n Breaking of the seq run time:\n"
                            << "\n  QoI function run time   = " << qoiFunctionRunTime
                            << " seconds ("                     << 100.*qoiFunctionRunTime/seqRunTime
                            << "%)"
                            << std::endl;
  }

  return;
}
// --------------------------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
void
MonteCarloSG<P_V,P_M,Q_V,Q_M>::actualReadSequence(
  const BaseVectorRV      <P_V,P_M>& paramRv,
  const std::string&                        dataInputFileName,
  const std::string&                        dataInputFileType,
        BaseVectorSequence<P_V,P_M>& workingPSeq,
        BaseVectorSequence<Q_V,Q_M>& workingQSeq,
        unsigned int                        requestedSeqSize)
{
  workingPSeq.resizeSequence(requestedSeqSize);
  P_V tmpP(m_paramSpace.zeroVector());
  for (unsigned int i = 0; i < requestedSeqSize; ++i) {
    paramRv.realizer().realization(tmpP);
    workingPSeq.setPositionValues(i,tmpP);
  }

  workingQSeq.unifiedReadContents(dataInputFileName,dataInputFileType,requestedSeqSize);

  return;
}

}  // End namespace QUESO

template class QUESO::MonteCarloSG<QUESO::GslVector, QUESO::GslMatrix, QUESO::GslVector, QUESO::GslMatrix>;
