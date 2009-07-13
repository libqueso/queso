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

#ifndef __UQ_MOC_SG_H__
#define __UQ_MOC_SG_H__

#include <uqVectorRV.h>
#include <uqVectorFunction.h>
#include <uqVectorFunctionSynchronizer.h>
#include <uqMonteCarloSGOptions.h>

/*! A templated class that implements a Monte Carlo Distribution Calculator
 */
template <class P_V,class P_M,class Q_V,class Q_M>
class uqMonteCarloSGClass
{
public:

  /*! Constructor: */
  uqMonteCarloSGClass(/*! Prefix           */ const char*                                       prefix,         
                      /*! The parameter rv */ const uqBaseVectorRVClass      <P_V,P_M>&         paramRv,        
                      /*! The qoi function */ const uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M>& qoiFunctionObj, 
                      /*! The qoi rv       */       uqBaseVectorRVClass      <Q_V,Q_M>&         qoiRv);         
  /*! Destructor: */
 ~uqMonteCarloSGClass();

  /*! Operation to generate the chain */
  void generateSequence           (uqBaseVectorSequenceClass<P_V,P_M>& workingPSeq,
                                   uqBaseVectorSequenceClass<Q_V,Q_M>& workingQSeq);
  void checkTheParallelEnvironment();

  void print                      (std::ostream& os) const;

private:
  void internGenerateSequence(const uqBaseVectorRVClass      <P_V,P_M>& paramRv,
                                    uqBaseVectorSequenceClass<P_V,P_M>& workingPSeq,
                                    uqBaseVectorSequenceClass<Q_V,Q_M>& workingQSeq);

  void actualGenerateSequence(const uqBaseVectorRVClass      <P_V,P_M>& paramRv,
                                    uqBaseVectorSequenceClass<P_V,P_M>& workingPSeq,
                                    uqBaseVectorSequenceClass<Q_V,Q_M>& workingQSeq,
                                    unsigned int                        seqSize);
  void actualReadSequence    (const uqBaseVectorRVClass      <P_V,P_M>& paramRv,
                              const std::string&                        dataInputFileName,
                                    uqBaseVectorSequenceClass<P_V,P_M>& workingPSeq,
                                    uqBaseVectorSequenceClass<Q_V,Q_M>& workingQSeq,
                                    unsigned int                        seqSize);

  const uqBaseEnvironmentClass&                             m_env;
  const uqBaseVectorRVClass              <P_V,P_M>&         m_paramRv;
  const uqBaseVectorFunctionClass        <P_V,P_M,Q_V,Q_M>& m_qoiFunctionObj;
  const uqBaseVectorRVClass              <Q_V,Q_M>&         m_qoiRv;
  const uqVectorSpaceClass               <P_V,P_M>&         m_paramSpace;
  const uqVectorSpaceClass               <Q_V,Q_M>&         m_qoiSpace;
  const uqVectorFunctionSynchronizerClass<P_V,P_M,Q_V,Q_M>* m_qoiFunctionSynchronizer;

        uqMonteCarloSGOptionsClass                          m_options;
};

template<class P_V,class P_M,class Q_V,class Q_M>
std::ostream& operator<<(std::ostream& os, const uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>& obj);

template <class P_V,class P_M,class Q_V,class Q_M>
uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::uqMonteCarloSGClass(
  const char*                                       prefix,
  const uqBaseVectorRVClass      <P_V,P_M>&         paramRv,
  const uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M>& qoiFunctionObj,
        uqBaseVectorRVClass      <Q_V,Q_M>&         qoiRv)
  :
  m_env                             (paramRv.env()),
  m_paramRv                         (paramRv),
  m_qoiFunctionObj                  (qoiFunctionObj),
  m_qoiRv                           (qoiRv),
  m_paramSpace                      (m_paramRv.imageSet().vectorSpace()),
  m_qoiSpace                        (m_qoiRv.imageSet().vectorSpace()),
  m_qoiFunctionSynchronizer         (new uqVectorFunctionSynchronizerClass<P_V,P_M,Q_V,Q_M>(m_qoiFunctionObj,m_paramRv.imageSet().vectorSpace().zeroVector(),m_qoiRv.imageSet().vectorSpace().zeroVector())),
  m_options                         (m_env,prefix)
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Entering uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::constructor()"
                            << std::endl;
  }

  m_options.scanOptionsValues();

  if (m_options.m_pseqComputeStats) m_options.m_pseqStatisticalOptions = new uqSequenceStatisticalOptionsClass(m_env,m_options.m_prefix + "pseq_");
  if (m_options.m_qseqComputeStats) m_options.m_qseqStatisticalOptions = new uqSequenceStatisticalOptionsClass(m_env,m_options.m_prefix + "qseq_");

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Leaving uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::constructor()"
                            << std::endl;
  }
}

template <class P_V,class P_M,class Q_V,class Q_M>
uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::~uqMonteCarloSGClass()
{
  if (m_options.m_qseqStatisticalOptions ) delete m_options.m_qseqStatisticalOptions;
  if (m_options.m_pseqStatisticalOptions ) delete m_options.m_pseqStatisticalOptions;
  if (m_qoiFunctionSynchronizer) delete m_qoiFunctionSynchronizer;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::generateSequence(
  uqBaseVectorSequenceClass<P_V,P_M>& workingPSeq,
  uqBaseVectorSequenceClass<Q_V,Q_M>& workingQSeq)
{
  checkTheParallelEnvironment();
  internGenerateSequence(m_paramRv,workingPSeq,workingQSeq);

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence(
  const uqBaseVectorRVClass      <P_V,P_M>& paramRv,
        uqBaseVectorSequenceClass<P_V,P_M>& workingPSeq,
        uqBaseVectorSequenceClass<Q_V,Q_M>& workingQSeq)
{
  workingPSeq.setName(m_options.m_prefix+"ParamSeq");
  workingQSeq.setName(m_options.m_prefix+"QoiSeq");

  //****************************************************
  // Generate sequence of qoi values
  //****************************************************
  unsigned int subActualSizeBeforeGeneration = std::min(m_options.m_qseqSize,paramRv.realizer().subPeriod());
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                                  << ": m_options.m_qseqSize = "                                                << m_options.m_qseqSize
                                  << ", paramRv.realizer().subPeriod() = "                            << paramRv.realizer().subPeriod()
                                  << ", about to call actualGenerateSequence() with subActualSize = " << subActualSizeBeforeGeneration
                                  << std::endl;
  }
  if (m_options.m_qseqDataInputFileName == UQ_MOC_SG_FILENAME_FOR_NO_FILE) {
    actualGenerateSequence(paramRv,
                           workingPSeq,
                           workingQSeq,
                           subActualSizeBeforeGeneration);
  }
  else {
    actualReadSequence(paramRv,
                       m_options.m_qseqDataInputFileName,
                       workingPSeq,
                       workingQSeq,
                       subActualSizeBeforeGeneration);
  }
  unsigned int subActualSizeAfterGeneration = workingPSeq.subSequenceSize();
  UQ_FATAL_TEST_MACRO(subActualSizeAfterGeneration != workingQSeq.subSequenceSize(),
                      m_env.fullRank(),
                      "uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()",
                      "P and Q sequences should have the same size!");

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                                  << ": returned from call to actualGenerateSequence() with subActualSize = " << subActualSizeAfterGeneration
                                  << std::endl;
  }

  //****************************************************
  // Open generic output file      
  //****************************************************
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                                  << ", prefix = "                                                        << m_options.m_prefix
                                  << ": checking necessity of opening generic output file (qseq name is " << workingQSeq.name()
                                  << ") ..."
                                  << std::endl;
  }
  std::ofstream* genericOfsVar = NULL;
  m_env.openOutputFile(m_options.m_dataOutputFileName,
                       UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT,
                       m_options.m_dataOutputAllowedSet,
                       false,
                       genericOfsVar);

  //****************************************************
  // Eventually:
  // --> write parameter sequence
  // --> compute statistics on it
  //****************************************************
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                           << ", prefix = "                                            << m_options.m_prefix
                           << ": checking necessity of opening output files for pseq " << workingPSeq.name()
                           << "..."
                           << std::endl;
  }

  // Take "sub" care of pseq
  std::ofstream* pseqOfsVar = NULL;
  m_env.openOutputFile(m_options.m_pseqDataOutputFileName,
                       UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT,
                       m_options.m_pseqDataOutputAllowedSet,
                       false, // A 'true' causes problems when the user chooses (via options
                              // in the input file) to use just one file for all outputs.
                       pseqOfsVar);

  if (pseqOfsVar) {
    workingPSeq.subWriteContents(*pseqOfsVar);
  }

  if (pseqOfsVar) {
    pseqOfsVar->close();
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                             << ", prefix = "                 << m_options.m_prefix
                             << ": closed data output file '" << m_options.m_pseqDataOutputFileName
                             << "' for pseq "                 << workingPSeq.name()
                             << std::endl;
    }
  }

  // Take "unified" care of pseq
#if 0
  std::ofstream* unifiedPSeqOfsVar = NULL;
  m_env.openUnifiedOutputFile(m_options.m_pseqDataOutputFileName,
                              UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT,
                              false, // A 'true' causes problems when the user chooses (via options
                                     // in the input file) to use just one file for all outputs.
                              unifiedPSeqOfsVar);

  if (unifiedPSeqOfsVar) {
    workingPSeq.unifiedWriteContents(*unifiedPSeqOfsVar);
  }

  if (unifiedPSeqOfsVar) {
    unifiedPSeqOfsVar->close();
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                             << ", prefix = "                         << m_options.m_prefix
                             << ": closed unified data output file '" << m_options.m_pseqDataOutputFileName
                             << "' for pseq "                         << workingPSeq.name()
                             << std::endl;
    }
  }
#else
  if (m_options.m_pseqDataOutputFileName != UQ_MOC_SG_FILENAME_FOR_NO_FILE) {
    workingPSeq.unifiedWriteContents(m_options.m_pseqDataOutputFileName);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                             << ", prefix = "                         << m_options.m_prefix
                             << ": closed unified data output file '" << m_options.m_pseqDataOutputFileName
                             << "' for pseq "                         << workingPSeq.name()
                             << std::endl;
    }
  }
#endif

  // Take case of other aspects of pseq
  if (m_options.m_pseqComputeStats) {
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                             << ": about to call 'workingPSeq.computeStatistics()'"
                             << std::endl;
    }
    workingPSeq.computeStatistics(*m_options.m_pseqStatisticalOptions,
                                  genericOfsVar);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                             << ": returned from call to 'workingPSeq.computeStatistics()'"
                             << std::endl;
    }
  }

  //****************************************************
  // Eventually:
  // --> write qoi sequence
  // --> compute statistics on it
  //****************************************************
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                           << ", prefix = "                                            << m_options.m_prefix
                           << ": checking necessity of opening output files for qseq " << workingQSeq.name()
                           << "..."
                           << std::endl;
  }

  // Take "sub" care of qseq
  std::ofstream* qseqOfsVar = NULL;
  m_env.openOutputFile(m_options.m_qseqDataOutputFileName,
                       UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT,
                       m_options.m_qseqDataOutputAllowedSet,
                       false, // A 'true' causes problems when the user chooses (via options
                              // in the input file) to use just one file for all outputs.
                       qseqOfsVar);

  if (qseqOfsVar) {
    workingQSeq.subWriteContents(*qseqOfsVar);
  }

  if (qseqOfsVar) {
    qseqOfsVar->close();
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                             << ", prefix = "                 << m_options.m_prefix
                             << ": closed data output file '" << m_options.m_qseqDataOutputFileName
                             << "' for qseq "                 << workingQSeq.name()
                             << std::endl;
    }
  }

  // Take "unified" care of qseq
  if (m_options.m_qseqDataOutputFileName != UQ_MOC_SG_FILENAME_FOR_NO_FILE) {
    workingQSeq.unifiedWriteContents(m_options.m_qseqDataOutputFileName);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                             << ", prefix = "                         << m_options.m_prefix
                             << ": closed unified data output file '" << m_options.m_qseqDataOutputFileName
                             << "' for qseq "                         << workingQSeq.name()
                             << std::endl;
    }
  }

  // Take case of other aspects of qseq
  if (m_options.m_qseqComputeStats) {
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                              << ": about to call 'workingQSeq.computeStatistics()'"
                              << std::endl;
    }
    workingQSeq.computeStatistics(*m_options.m_qseqStatisticalOptions,
                                  genericOfsVar);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                              << ": returned from call to 'workingQSeq.computeStatistics()'"
                              << std::endl;
    }
  }

  //****************************************************
  // Close generic output file      
  //****************************************************
  if (genericOfsVar) {
    genericOfsVar->close();
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                             << ", prefix = "                         << m_options.m_prefix
                             << ": closed generic data output file '" << m_options.m_dataOutputFileName
                             << "' for qoi sequence "                 << workingQSeq.name()
                             << std::endl;
    }
  }
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << std::endl;
  }

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::actualGenerateSequence(
  const uqBaseVectorRVClass      <P_V,P_M>& paramRv,
        uqBaseVectorSequenceClass<P_V,P_M>& workingPSeq,
        uqBaseVectorSequenceClass<Q_V,Q_M>& workingQSeq,
        unsigned int                        requestedSeqSize)
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Starting the generation of qoi sequence " << workingQSeq.name()
                           << ", with "                                   << requestedSeqSize
                           << " samples..."
                           << std::endl;
  }

  int iRC = UQ_OK_RC;
  struct timeval timevalSeq;
  struct timeval timevalQoIFunction;

  double seqRunTime         = 0;
  double qoiFunctionRunTime = 0;

  iRC = gettimeofday(&timevalSeq, NULL);

  workingPSeq.resizeSequence(requestedSeqSize);
  workingQSeq.resizeSequence(requestedSeqSize);
  P_V tmpP(m_paramSpace.zeroVector());
  Q_V tmpQ(m_qoiSpace.zeroVector());

  unsigned int actualSeqSize = 0;
  for (unsigned int i = 0; i < requestedSeqSize; ++i) {
    paramRv.realizer().realization(tmpP);

    if (m_options.m_qseqMeasureRunTimes) iRC = gettimeofday(&timevalQoIFunction, NULL);
    m_qoiFunctionSynchronizer->callFunction(&tmpP,NULL,&tmpQ,NULL,NULL,NULL); // Might demand parallel environment
    if (m_options.m_qseqMeasureRunTimes) qoiFunctionRunTime += uqMiscGetEllapsedSeconds(&timevalQoIFunction);

    bool allQsAreFinite = true;
    for (unsigned int j = 0; j < tmpQ.sizeLocal(); ++j) {
      if ((tmpQ[j] == INFINITY) || (tmpQ[j] == -INFINITY)) {
	std::cerr << "WARNING In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::actualGenerateSequence()"
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

    //if (allQsAreFinite) { // FIXME: this will cause different processors to have sequences of different sizes
      workingPSeq.setPositionValues(i,tmpP);
      workingQSeq.setPositionValues(i,tmpQ);
      actualSeqSize++;
    //}

    if ((m_options.m_qseqDisplayPeriod            > 0) && 
        (((i+1) % m_options.m_qseqDisplayPeriod) == 0)) {
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

  seqRunTime = uqMiscGetEllapsedSeconds(&timevalSeq);

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

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::actualReadSequence(
  const uqBaseVectorRVClass      <P_V,P_M>& paramRv,
  const std::string&                        dataInputFileName,
        uqBaseVectorSequenceClass<P_V,P_M>& workingPSeq,
        uqBaseVectorSequenceClass<Q_V,Q_M>& workingQSeq,
        unsigned int                        requestedSeqSize)
{
  workingPSeq.resizeSequence(requestedSeqSize);
  P_V tmpP(m_paramSpace.zeroVector());
  for (unsigned int i = 0; i < requestedSeqSize; ++i) {
    paramRv.realizer().realization(tmpP);
    workingPSeq.setPositionValues(i,tmpP);
  }

  workingQSeq.unifiedReadContents(dataInputFileName,requestedSeqSize);

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::checkTheParallelEnvironment()
{
  if (m_env.numSubEnvironments() == (unsigned int) m_env.fullComm().NumProc()) {
    UQ_FATAL_TEST_MACRO(m_env.subRank() != 0,
                        m_env.fullRank(),
                        "uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::checkTheParallelEnvironment()",
                        "there should exist only one processor per sub environment");
    UQ_FATAL_TEST_MACRO(m_paramRv.imageSet().vectorSpace().zeroVector().numberOfProcessorsRequiredForStorage() != 1,
                        m_env.fullRank(),
                        "uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::checkTheParallelEnvironment()",
                        "only 1 processor (per sub environment) should be necessary for the storage of a parameter vector");
  }
  else if (m_env.numSubEnvironments() < (unsigned int) m_env.fullComm().NumProc()) {
    UQ_FATAL_TEST_MACRO(m_env.fullComm().NumProc()%m_env.numSubEnvironments() != 0,
                        m_env.fullRank(),
                        "uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::checkTheParallelEnvironment()",
                        "total number of processors should be a multiple of the number of sub environments");
    unsigned int numProcsPerSubEnvironment = m_env.fullComm().NumProc()/m_env.numSubEnvironments();
    UQ_FATAL_TEST_MACRO(m_env.subComm().NumProc() != (int) numProcsPerSubEnvironment,
                        m_env.fullRank(),
                        "uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::checkTheParallelEnvironment()",
                        "inconsistent number of processors per sub environment");
    if ((m_paramRv.imageSet().vectorSpace().zeroVector().numberOfProcessorsRequiredForStorage() == 1) &&
        (m_qoiRv.imageSet().vectorSpace().zeroVector().numberOfProcessorsRequiredForStorage()   == 1)) {
      // Ok
    }
    else if ((m_paramRv.imageSet().vectorSpace().zeroVector().numberOfProcessorsRequiredForStorage() == numProcsPerSubEnvironment) &&
             (m_qoiRv.imageSet().vectorSpace().zeroVector().numberOfProcessorsRequiredForStorage()   == numProcsPerSubEnvironment)) {
      UQ_FATAL_TEST_MACRO(true,
                          m_env.fullRank(),
                          "uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::checkTheParallelEnvironment()",
                          "parallel vectors are not supported yet");
    }
    else {
      UQ_FATAL_TEST_MACRO(true,
                          m_env.fullRank(),
                          "uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::checkTheParallelEnvironment()",
                          "number of processors required for a vector storage should be equal to the number of processors in the sub environment");
    }
  }
  else {
    UQ_FATAL_TEST_MACRO(true,
                        m_env.fullRank(),
                        "uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::checkTheParallelEnvironment()",
                        "number of processors per sub environment is too large");
  }

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::print(std::ostream& os) const
{
  return;
}

template<class P_V,class P_M,class Q_V,class Q_M> 
std::ostream& operator<<(std::ostream& os, const uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>& obj)
{
  obj.print(os);

  return os;
}
#endif // __UQ_MOC_SG_H__
