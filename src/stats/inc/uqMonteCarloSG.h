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

#ifndef __UQ_MOC_SG_H__
#define __UQ_MOC_SG_H__

#include <uqVectorRV.h>
#include <uqVectorFunction.h>
#include <uqVectorFunctionSynchronizer.h>
#include <uqMonteCarloSGOptions.h>

namespace QUESO {

/*! 
 * \file uqMonteCarloSG.h
 * \brief A templated class that implements a Monte Carlo generator of samples.
 * 
 * \class uqMonteCarloSGClass
 * \brief A templated class that implements a Monte Carlo generator of samples.
 *
 * This class implements a Monte Carlo generator of samples. 'SG' stands for 'Sequence Generator'.
 * Options reading is handled by class 'uqMonteCarloOptionsClass'. If options request data to be 
 * written in the output file (MATLAB .m format only, for now), the user can check which MATLAB 
 * variables are defined and set by running 'grep zeros <OUTPUT FILE NAME>' after the solution
 * procedures ends. The names of the variables are self explanatory. */

template <class P_V,class P_M,class Q_V,class Q_M>
class uqMonteCarloSGClass
{
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Constructor.
 /*! Requirements: 1) the image set of the vector random variable 'paramRv' and 2) the domain set of the
 * QoI function 'qoiFunction' should belong to vector spaces of equal dimensions. If the requirements 
 * are satisfied, the constructor then reads input options that begin with the string '\<prefix\>_mc_'.
 * For instance, if 'prefix' is 'pROblem_775_fp_', then the constructor will read all options that begin
 * with 'pROblem_775_fp_mc_'. Options reading is handled by class 'uqMonteCarloOptionsClass'.*/
  uqMonteCarloSGClass(const char*                                       prefix,         
                      const uqMcOptionsValuesClass*                     alternativeOptionsValues, // dakota
                      const uqBaseVectorRVClass      <P_V,P_M>&         paramRv,        
                      const uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M>& qoiFunction);
  
  //! Destructor.
  ~uqMonteCarloSGClass();
  //@}
 
  //! @name Statistical methods
  //@{
  //! Generates the QoI (output) sequence, it calls internGenerateSequence().
  /*! This method checks for a parallel environment (and uses it if available) and calls the private
   * member function internGenerateSequence() to generate a sequence of values of the Quantity of
   * interest (QoI).*/
  void generateSequence(uqBaseVectorSequenceClass<P_V,P_M>& workingPSeq,
                        uqBaseVectorSequenceClass<Q_V,Q_M>& workingQSeq);
  //@}
  
  //! @name I/O methods
  //@{ 
  //! Prints the sequence.
  void print                      (std::ostream& os) const;
  //@}
private:
  //! Generates the QoI (output) sequence; it calls actualGenerateSequence().
  /*!  This method calls the private member actualGenerateSequence() to \b actually generate a 
   * sequence of values of the Quantity of interest (QoI). It writes both the QoI and the parameter 
   * sequences to output files. If a parallel environment is available, the user has the choice of 
   * (via the input options file) which sub-environment will write data to files. The unified 
   * sequence of data (QoI and/or Parameter sequences) of all the sub-environments (or of the single
   * environment, in case of a serial run) is also written to files.*/
  void internGenerateSequence(const uqBaseVectorRVClass      <P_V,P_M>& paramRv,
                                    uqBaseVectorSequenceClass<P_V,P_M>& workingPSeq,
                                    uqBaseVectorSequenceClass<Q_V,Q_M>& workingQSeq);

  
  //!  This method \b actually generates the QoI sequence.
  /*! Requirements: 1) the vector space containing the domain set of the QoI function 'm_qoiFunction' 
   * should have dimension equal to the size of a vector in 'workingPSeq' and 2) the vector space 
   * containing the image set of the qoi function 'm_qoiFunction' should have dimension equal to the 
   * size of a vector in 'workingQSeq'. If the requirements are satisfied, this operation sets the 
   * size and the contents of 'workingPSeq' and 'workingQSeq' using the algorithm options set in the
   * constructor. */
  void actualGenerateSequence(const uqBaseVectorRVClass      <P_V,P_M>& paramRv,
                                    uqBaseVectorSequenceClass<P_V,P_M>& workingPSeq,
                                    uqBaseVectorSequenceClass<Q_V,Q_M>& workingQSeq,
                                    unsigned int                        seqSize);
  //! Reads the sequence.
  void actualReadSequence    (const uqBaseVectorRVClass      <P_V,P_M>& paramRv,
                              const std::string&                        dataInputFileName,
                              const std::string&                        dataInputFileType,
                                    uqBaseVectorSequenceClass<P_V,P_M>& workingPSeq,
                                    uqBaseVectorSequenceClass<Q_V,Q_M>& workingQSeq,
                                    unsigned int                        seqSize);

  const uqBaseEnvironmentClass&                             m_env;
  const uqBaseVectorRVClass              <P_V,P_M>&         m_paramRv;
  const uqBaseVectorFunctionClass        <P_V,P_M,Q_V,Q_M>& m_qoiFunction;
  const uqVectorSpaceClass               <P_V,P_M>&         m_paramSpace;
  const uqVectorSpaceClass               <Q_V,Q_M>&         m_qoiSpace;
  const uqVectorFunctionSynchronizerClass<P_V,P_M,Q_V,Q_M>* m_qoiFunctionSynchronizer;
  unsigned int                                              m_numPsNotSubWritten;
  unsigned int                                              m_numQsNotSubWritten;

  uqMcOptionsValuesClass                                    m_alternativeOptionsValues;
  uqMonteCarloSGOptionsClass*                               m_optionsObj;
};

//! Prints the object \c obj, overloading an operator.
template<class P_V,class P_M,class Q_V,class Q_M>
std::ostream& operator<<(std::ostream& os, const uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>& obj);

// Default constructor -----------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::uqMonteCarloSGClass(
  /*! Prefix                     */ const char*                                       prefix,
  /*! Options (if no input file) */ const uqMcOptionsValuesClass*                     alternativeOptionsValues, // dakota
  /*! The parameter RV           */ const uqBaseVectorRVClass      <P_V,P_M>&         paramRv,
  /*! The QoI function           */ const uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M>& qoiFunction)
  :
  m_env                     (paramRv.env()),
  m_paramRv                 (paramRv),
  m_qoiFunction             (qoiFunction),
  m_paramSpace              (m_paramRv.imageSet().vectorSpace()),
  m_qoiSpace                (m_qoiFunction.imageSet().vectorSpace()),
  m_qoiFunctionSynchronizer (new uqVectorFunctionSynchronizerClass<P_V,P_M,Q_V,Q_M>(m_qoiFunction,m_paramRv.imageSet().vectorSpace().zeroVector(),m_qoiFunction.imageSet().vectorSpace().zeroVector())),
  m_numPsNotSubWritten      (0),
  m_numQsNotSubWritten      (0),
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  m_alternativeOptionsValues(NULL,NULL),
#else
  m_alternativeOptionsValues(),
#endif
  m_optionsObj              (NULL)
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Entering uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": prefix = " << prefix
                            << ", alternativeOptionsValues = " << alternativeOptionsValues
                            << ", m_env.optionsInputFileName() = " << m_env.optionsInputFileName()
                            << std::endl;
  }

  if (alternativeOptionsValues) m_alternativeOptionsValues = *alternativeOptionsValues;
  if (m_env.optionsInputFileName() == "") {
    m_optionsObj = new uqMonteCarloSGOptionsClass(m_env,prefix,m_alternativeOptionsValues);
  }
  else {
    m_optionsObj = new uqMonteCarloSGOptionsClass(m_env,prefix);
    m_optionsObj->scanOptionsValues();
  }

  UQ_FATAL_TEST_MACRO(paramRv.imageSet().vectorSpace().dimLocal() != qoiFunction.domainSet().vectorSpace().dimLocal(),
                      m_env.worldRank(),
                      "uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::constructor()",
                      "'paramRv' and 'qoiFunction' are related to vector spaces of different dimensions");

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Leaving uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::constructor()"
                            << std::endl;
  }
}
// Destructor ---------------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::~uqMonteCarloSGClass()
{
  if (m_optionsObj             ) delete m_optionsObj;
  if (m_qoiFunctionSynchronizer) delete m_qoiFunctionSynchronizer;
}
// Statistical methods ------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
void
uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::generateSequence(
  uqBaseVectorSequenceClass<P_V,P_M>& workingPSeq,
  uqBaseVectorSequenceClass<Q_V,Q_M>& workingQSeq)
{
  UQ_FATAL_TEST_MACRO(m_qoiFunction.domainSet().vectorSpace().dimLocal() != workingPSeq.vectorSizeLocal(),
                      m_env.worldRank(),
                      "uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::generateSequence()",
                      "'m_qoiFunction.domainSet' and 'workingPSeq' are related to vector spaces of different dimensions");

  UQ_FATAL_TEST_MACRO(m_qoiFunction.imageSet().vectorSpace().dimLocal() != workingQSeq.vectorSizeLocal(),
                      m_env.worldRank(),
                      "uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::generateSequence()",
                      "'m_qoiFunction.imageSet' and 'workingQSeq' are related to vector spaces of different dimensions");

  uqMiscCheckTheParallelEnvironment<P_V,Q_V>(m_paramRv.imageSet().vectorSpace().zeroVector(),
                                             m_qoiFunction.imageSet().vectorSpace().zeroVector());
  internGenerateSequence(m_paramRv,workingPSeq,workingQSeq);

  return;
}
// I/O methods---------------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
void
uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::print(std::ostream& os) const
{
  return;
}
// Private methods----------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
void
uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence(
  const uqBaseVectorRVClass      <P_V,P_M>& paramRv,
        uqBaseVectorSequenceClass<P_V,P_M>& workingPSeq,
        uqBaseVectorSequenceClass<Q_V,Q_M>& workingQSeq)
{
  workingPSeq.setName(m_optionsObj->m_prefix+"ParamSeq");
  workingQSeq.setName(m_optionsObj->m_prefix+"QoiSeq");

  //****************************************************
  // Generate sequence of QoI values
  //****************************************************
  unsigned int subActualSizeBeforeGeneration = std::min(m_optionsObj->m_ov.m_qseqSize,paramRv.realizer().subPeriod());
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
    *m_env.subDisplayFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                            << ": m_optionsObj->m_ov.m_qseqSize = "                             << m_optionsObj->m_ov.m_qseqSize
                            << ", paramRv.realizer().subPeriod() = "                            << paramRv.realizer().subPeriod()
                            << ", about to call actualGenerateSequence() with subActualSize = " << subActualSizeBeforeGeneration
                            << std::endl;
  }
  if (m_optionsObj->m_ov.m_qseqDataInputFileName == UQ_MOC_SG_FILENAME_FOR_NO_FILE) {
    actualGenerateSequence(paramRv,
                           workingPSeq,
                           workingQSeq,
                           subActualSizeBeforeGeneration);
  }
  else {
    actualReadSequence(paramRv,
                       m_optionsObj->m_ov.m_qseqDataInputFileName,
                       m_optionsObj->m_ov.m_qseqDataInputFileType,
                       workingPSeq,
                       workingQSeq,
                       subActualSizeBeforeGeneration);
  }
  unsigned int subActualSizeAfterGeneration = workingPSeq.subSequenceSize();
  UQ_FATAL_TEST_MACRO(subActualSizeAfterGeneration != workingQSeq.subSequenceSize(),
                      m_env.worldRank(),
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
                            << ", prefix = "                                                        << m_optionsObj->m_prefix
                            << ": checking necessity of opening generic output file (qseq name is " << workingQSeq.name()
                            << ") ..."
                            << std::endl;
  }
  uqFilePtrSetStruct genericFilePtrSet;
  m_env.openOutputFile(m_optionsObj->m_ov.m_dataOutputFileName,
                       UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT, // Yes, always ".m"
                       m_optionsObj->m_ov.m_dataOutputAllowedSet,
                       false,
                       genericFilePtrSet);

  //****************************************************
  // Eventually:
  // --> write parameter sequence
  // --> compute statistics on it
  //****************************************************
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                            << ", prefix = "                                            << m_optionsObj->m_prefix
                            << ": checking necessity of opening output files for pseq " << workingPSeq.name()
                            << "..."
                            << std::endl;
  }

  // Take "sub" care of pseq
  if ((m_numPsNotSubWritten                        >  0                             ) &&
      (m_optionsObj->m_ov.m_pseqDataOutputFileName != UQ_MOC_SG_FILENAME_FOR_NO_FILE)) {
    workingPSeq.subWriteContents(subActualSizeBeforeGeneration - m_numPsNotSubWritten,
                                 m_numPsNotSubWritten, 
                                 m_optionsObj->m_ov.m_pseqDataOutputFileName,
                                 m_optionsObj->m_ov.m_pseqDataOutputFileType,
                                 m_optionsObj->m_ov.m_pseqDataOutputAllowedSet);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqMonteCarloGClass<P_V,P_M>::internGenerateSequence()"
                              << ": just wrote remaining pseq positions (per period request)"
                              << std::endl;
    }
    m_numPsNotSubWritten = 0;
    //if (m_env.subDisplayFile()) {
    //  *m_env.subDisplayFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
    //                          << ", prefix = "                 << m_optionsObj->m_prefix
    //                          << ": closed data output file '" << m_optionsObj->m_ov.m_pseqDataOutputFileName
    //                          << "' for pseq "                 << workingPSeq.name()
    //                          << std::endl;
    //}
  }

  // Take "unified" care of pseq
  if (m_optionsObj->m_ov.m_pseqDataOutputFileName != UQ_MOC_SG_FILENAME_FOR_NO_FILE) {
    workingPSeq.unifiedWriteContents(m_optionsObj->m_ov.m_pseqDataOutputFileName,m_optionsObj->m_ov.m_pseqDataOutputFileType);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                              << ", prefix = "                         << m_optionsObj->m_prefix
                              << ": closed unified data output file '" << m_optionsObj->m_ov.m_pseqDataOutputFileName
                              << "' for pseq "                         << workingPSeq.name()
                              << std::endl;
    }
  }

  // Take case of other aspects of pseq
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  if (m_optionsObj->m_ov.m_pseqComputeStats) {
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                              << ": about to call 'workingPSeq.computeStatistics()'"
                              << std::endl;
    }
    workingPSeq.computeStatistics(*m_optionsObj->m_pseqStatisticalOptionsObj,
                                  genericFilePtrSet.ofsVar);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
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
    *m_env.subDisplayFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                            << ", prefix = "                                            << m_optionsObj->m_prefix
                            << ": checking necessity of opening output files for qseq " << workingQSeq.name()
                            << "..."
                            << std::endl;
  }

  // Take "sub" care of qseq
  if ((m_numQsNotSubWritten                        >  0                             ) &&
      (m_optionsObj->m_ov.m_qseqDataOutputFileName != UQ_MOC_SG_FILENAME_FOR_NO_FILE)) {
    workingQSeq.subWriteContents(subActualSizeBeforeGeneration - m_numQsNotSubWritten,
                                 m_numQsNotSubWritten, 
                                 m_optionsObj->m_ov.m_qseqDataOutputFileName,
                                 m_optionsObj->m_ov.m_qseqDataOutputFileType,
                                 m_optionsObj->m_ov.m_qseqDataOutputAllowedSet);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqMonteCarloGClass<P_V,P_M>::internGenerateSequence()"
                              << ": just wrote remaining qseq positions (per period request)"
                              << std::endl;
    }
    m_numQsNotSubWritten = 0;
    //if (m_env.subDisplayFile()) {
    //  *m_env.subDisplayFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
    //                          << ", prefix = "                 << m_optionsObj->m_prefix
    //                          << ": closed data output file '" << m_optionsObj->m_ov.m_qseqDataOutputFileName
    //                          << "' for qseq "                 << workingQSeq.name()
    //                          << std::endl;
    //}
  }

  // Take "unified" care of qseq
  if (m_optionsObj->m_ov.m_qseqDataOutputFileName != UQ_MOC_SG_FILENAME_FOR_NO_FILE) {
    workingQSeq.unifiedWriteContents(m_optionsObj->m_ov.m_qseqDataOutputFileName,m_optionsObj->m_ov.m_qseqDataOutputFileType);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                              << ", prefix = "                         << m_optionsObj->m_prefix
                              << ": closed unified data output file '" << m_optionsObj->m_ov.m_qseqDataOutputFileName
                              << "' for qseq "                         << workingQSeq.name()
                              << std::endl;
    }
  }

  // Take case of other aspects of qseq
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  if (m_optionsObj->m_ov.m_qseqComputeStats) {
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                              << ": about to call 'workingQSeq.computeStatistics()'"
                              << std::endl;
    }
    workingQSeq.computeStatistics(*m_optionsObj->m_qseqStatisticalOptionsObj,
                                  genericFilePtrSet.ofsVar);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 0)) {
      *m_env.subDisplayFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
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
      *m_env.subDisplayFile() << "In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::internGenerateSequence()"
                              << ", prefix = "                         << m_optionsObj->m_prefix
                              << ": closed generic data output file '" << m_optionsObj->m_ov.m_dataOutputFileName
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
uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::actualGenerateSequence(
  const uqBaseVectorRVClass      <P_V,P_M>& paramRv,
        uqBaseVectorSequenceClass<P_V,P_M>& workingPSeq,
        uqBaseVectorSequenceClass<Q_V,Q_M>& workingQSeq,
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

    if (m_optionsObj->m_ov.m_qseqMeasureRunTimes) iRC = gettimeofday(&timevalQoIFunction, NULL);
    m_qoiFunctionSynchronizer->callFunction(&tmpP,NULL,&tmpQ,NULL,NULL,NULL); // Might demand parallel environment
    if (m_optionsObj->m_ov.m_qseqMeasureRunTimes) qoiFunctionRunTime += uqMiscGetEllapsedSeconds(&timevalQoIFunction);

    bool allQsAreFinite = true;
    for (unsigned int j = 0; j < tmpQ.sizeLocal(); ++j) {
      if ((tmpQ[j] == INFINITY) || (tmpQ[j] == -INFINITY)) {
	std::cerr << "WARNING In uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::actualGenerateSequence()"
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
      if ((m_optionsObj->m_ov.m_pseqDataOutputPeriod           >  0  ) && 
          (((i+1) % m_optionsObj->m_ov.m_pseqDataOutputPeriod) == 0  ) &&
          (m_optionsObj->m_ov.m_pseqDataOutputFileName         != ".")) {
        workingPSeq.subWriteContents(i + 1 - m_optionsObj->m_ov.m_pseqDataOutputPeriod,
                                     m_optionsObj->m_ov.m_pseqDataOutputPeriod, 
                                     m_optionsObj->m_ov.m_pseqDataOutputFileName,
                                     m_optionsObj->m_ov.m_pseqDataOutputFileType,
                                     m_optionsObj->m_ov.m_pseqDataOutputAllowedSet);
        if (m_env.subDisplayFile()) {
          *m_env.subDisplayFile() << "In uqMonteCarloGClass<P_V,P_M>::actualGenerateSequence()"
                                  << ": just wrote pseq positions (per period request)"
                                  << std::endl;
        }
        m_numPsNotSubWritten = 0;
      }

      workingQSeq.setPositionValues(i,tmpQ);
      m_numQsNotSubWritten++;
      if ((m_optionsObj->m_ov.m_qseqDataOutputPeriod           >  0  ) && 
          (((i+1) % m_optionsObj->m_ov.m_qseqDataOutputPeriod) == 0  ) &&
          (m_optionsObj->m_ov.m_qseqDataOutputFileName         != ".")) {
        workingQSeq.subWriteContents(i + 1 - m_optionsObj->m_ov.m_qseqDataOutputPeriod,
                                     m_optionsObj->m_ov.m_qseqDataOutputPeriod, 
                                     m_optionsObj->m_ov.m_qseqDataOutputFileName,
                                     m_optionsObj->m_ov.m_qseqDataOutputFileType,
                                     m_optionsObj->m_ov.m_qseqDataOutputAllowedSet);
        if (m_env.subDisplayFile()) {
          *m_env.subDisplayFile() << "In uqMonteCarloGClass<P_V,P_M>::actualGenerateSequence()"
                                  << ": just wrote qseq positions (per period request)"
                                  << std::endl;
        }
        m_numQsNotSubWritten = 0;
      }

      actualSeqSize++;

    //}

    if ((m_optionsObj->m_ov.m_qseqDisplayPeriod            > 0) && 
        (((i+1) % m_optionsObj->m_ov.m_qseqDisplayPeriod) == 0)) {
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
// --------------------------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
void
uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>::actualReadSequence(
  const uqBaseVectorRVClass      <P_V,P_M>& paramRv,
  const std::string&                        dataInputFileName,
  const std::string&                        dataInputFileType,
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

  workingQSeq.unifiedReadContents(dataInputFileName,dataInputFileType,requestedSeqSize);

  return;
}
//---------------------------------------------------
// Operator declared outside class definition -------
//---------------------------------------------------
template<class P_V,class P_M,class Q_V,class Q_M> 
std::ostream& operator<<(std::ostream& os, const uqMonteCarloSGClass<P_V,P_M,Q_V,Q_M>& obj)
{
  obj.print(os);

  return os;
}

}  // End namespace QUESO

#endif // __UQ_MOC_SG_H__
