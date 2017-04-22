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

#include <queso/SequenceOfVectors.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Default constructor -----------------------------
template <class V, class M>
SequenceOfVectors<V,M>::SequenceOfVectors(
  const VectorSpace<V,M>& vectorSpace,
  unsigned int                   subSequenceSize,
  const std::string&             name)
  :
  BaseVectorSequence<V,M>(vectorSpace,subSequenceSize,name),
  m_seq                         (subSequenceSize,NULL)
#ifdef UQ_CODE_HAS_MONITORS
  ,
  m_subMeanMonitorPosSeq        (NULL),
  m_subMeanVecSeq               (NULL),
  m_subMeanCltStdSeq            (NULL),
  m_subMeanInter0MonitorPosSeq  (NULL),
  m_subMeanInter0Mean           (NULL),
  m_subMeanInter0Clt95          (NULL),
  m_subMeanInter0Empirical90    (NULL),
  m_subMeanInter0Min            (NULL),
  m_subMeanInter0Max            (NULL),
  m_unifiedMeanMonitorPosSeq    (NULL),
  m_unifiedMeanVecSeq           (NULL),
  m_unifiedMeanCltStdSeq        (NULL)
#endif
{
  //if (m_env.subDisplayFile()) {
  //  *m_env.subDisplayFile() << "Entering SequenceOfVectors<V,M>::constructor()"
  //                          << std::endl;
  //}

  //if (m_env.subDisplayFile()) {
  //  *m_env.subDisplayFile() << "Leaving SequenceOfVectors<V,M>::constructor()"
  //                          << std::endl;
  //}
}
// Destructor ---------------------------------------
template <class V, class M>
SequenceOfVectors<V,M>::~SequenceOfVectors()
{
#ifdef UQ_CODE_HAS_MONITORS
  if (m_subMeanMonitorPosSeq) delete m_subMeanMonitorPosSeq;
  if (m_subMeanVecSeq       ) delete m_subMeanVecSeq;
  if (m_subMeanCltStdSeq    ) delete m_subMeanCltStdSeq;

  if (m_subMeanInter0MonitorPosSeq) delete m_subMeanInter0MonitorPosSeq;
  if (m_subMeanInter0Mean         ) delete m_subMeanInter0Mean;
  if (m_subMeanInter0Clt95        ) delete m_subMeanInter0Clt95;
  if (m_subMeanInter0Empirical90  ) delete m_subMeanInter0Empirical90;
  if (m_subMeanInter0Min          ) delete m_subMeanInter0Min;
  if (m_subMeanInter0Max          ) delete m_subMeanInter0Max;

  if (m_unifiedMeanMonitorPosSeq) delete m_unifiedMeanMonitorPosSeq;
  if (m_unifiedMeanVecSeq       ) delete m_unifiedMeanVecSeq;
  if (m_unifiedMeanCltStdSeq    ) delete m_unifiedMeanCltStdSeq;
#endif

  for (unsigned int i = 0; i < (unsigned int) m_seq.size(); ++i) {
    if (m_seq[i]) delete m_seq[i];
  }
}
// Set methods --------------------------------------
template <class V, class M>
SequenceOfVectors<V,M>&
SequenceOfVectors<V,M>::operator= (const SequenceOfVectors<V,M>& rhs)
{
  this->copy(rhs);
  return *this;
}

// Sequence methods ---------------------------------
template <class V, class M>
unsigned int
SequenceOfVectors<V,M>::subSequenceSize() const
{
  return m_seq.size();
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::resizeSequence(unsigned int newSubSequenceSize)
{
  if (newSubSequenceSize != this->subSequenceSize()) {
    if (newSubSequenceSize < this->subSequenceSize()) {
      this->resetValues(newSubSequenceSize,this->subSequenceSize()-newSubSequenceSize);
    }
    m_seq.resize(newSubSequenceSize,NULL);
    std::vector<const V*>(m_seq).swap(m_seq);
    BaseVectorSequence<V,M>::deleteStoredVectors();
  }

 return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::resetValues(unsigned int initialPos, unsigned int numPos)
{
  bool bRC = ((initialPos          <  this->subSequenceSize()) &&
              (0                   <  numPos                 ) &&
              ((initialPos+numPos) <= this->subSequenceSize()));
  if ((bRC == false) && (m_env.subDisplayFile())) {
    *m_env.subDisplayFile() << "In SequenceOfVectors<V,M>::resetValues()"
                           << ", initialPos = "              << initialPos
                           << ", this->subSequenceSize() = " << this->subSequenceSize()
                           << ", numPos = "                  << numPos
                           << std::endl;
  }
  queso_require_msg(bRC, "invalid input data");

  for (unsigned int j = 0; j < numPos; ++j) {
    if (m_seq[initialPos+j] != NULL) {
      delete m_seq[initialPos+j];
      m_seq[initialPos+j] = NULL;
    }
  }

  BaseVectorSequence<V,M>::deleteStoredVectors();

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::erasePositions(unsigned int initialPos, unsigned int numPos)
{
  bool bRC = ((initialPos          <  this->subSequenceSize()) &&
              (0                   <  numPos                 ) &&
              ((initialPos+numPos) <= this->subSequenceSize()));
  queso_require_msg(bRC, "invalid input data");

  for (unsigned int j = 0; j < numPos; ++j) {
    if (m_seq[initialPos+j] != NULL) {
      delete m_seq[initialPos+j];
      m_seq[initialPos+j] = NULL;
    }
  }

  seqVectorPositionIteratorTypedef posIteratorBegin = m_seq.begin();
  if (initialPos < this->subSequenceSize()) std::advance(posIteratorBegin,initialPos);
  else                                      posIteratorBegin = m_seq.end();

  unsigned int posEnd = initialPos + numPos;
  seqVectorPositionIteratorTypedef posIteratorEnd = m_seq.begin();
  if (posEnd < this->subSequenceSize()) std::advance(posIteratorEnd,posEnd);
  else                                  posIteratorEnd = m_seq.end();

  unsigned int oldSubSequenceSize = this->subSequenceSize();
  m_seq.erase(posIteratorBegin,posIteratorEnd);
  queso_require_equal_to_msg((oldSubSequenceSize - numPos), this->subSequenceSize(), "(oldSubSequenceSize - numPos) != this->subSequenceSize()");

  BaseVectorSequence<V,M>::deleteStoredVectors();

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::getPositionValues(unsigned int posId, V& vec) const
{
  queso_require_less_msg(posId, this->subSequenceSize(), "posId > subSequenceSize()");

  queso_require_msg(m_seq[posId], "posId is NULL");

  //if (posId == 0) { // mox
  //  std::cout << "In SequenceOfVectors<V,M>::getPositionValues(): m_seq[0] = " << m_seq[0] << ", *(m_seq[0]) = " << *(m_seq[0])
  //            << std::endl;
  //}

  vec = *(m_seq[posId]); // *(const_cast<V*>(m_seq[posId])); // prudenci 2010-06-17 mox

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::setPositionValues(unsigned int posId, const V& vec)
{
  queso_require_less_msg(posId, this->subSequenceSize(), "posId > subSequenceSize()");

  queso_require_equal_to_msg(vec.sizeLocal(), m_vectorSpace.zeroVector().sizeLocal(), "invalid vec");

  if (m_seq[posId] != NULL) delete m_seq[posId];
  m_seq[posId] = new V(vec);

  //if (posId == 0) { // mox
  //  std::cout << "In SequenceOfVectors<V,M>::setPositionValues(): m_seq[0] = " << m_seq[0] << ", *(m_seq[0]) = " << *(m_seq[0])
  //            << std::endl;
  //}

  BaseVectorSequence<V,M>::deleteStoredVectors();

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::subUniformlySampledCdf(
  const V&                       numEvaluationPointsVec,
  ArrayOfOneDGrids <V,M>& cdfGrids,
  ArrayOfOneDTables<V,M>& cdfValues) const
{
  V minDomainValues(m_vectorSpace.zeroVector());
  V maxDomainValues(m_vectorSpace.zeroVector());

  ScalarSequence<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(0,                 // initialPos
                           1,                 // spacing
                           subSequenceSize(), // numPos
                           i,
                           data);

    std::vector<double> aCdf(0);
    data.subUniformlySampledCdf((unsigned int) numEvaluationPointsVec[i],
                                minDomainValues[i],
                                maxDomainValues[i],
                                aCdf);
    cdfValues.setOneDTable(i,aCdf);
  }

  cdfGrids.setUniformGrids(numEvaluationPointsVec,
                           minDomainValues,
                           maxDomainValues);

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::unifiedUniformlySampledCdf(
  const V&                       numEvaluationPointsVec,
  ArrayOfOneDGrids <V,M>& unifiedCdfGrids,
  ArrayOfOneDTables<V,M>& unifiedCdfValues) const
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 10)) {
    *m_env.subDisplayFile() << "Entering SequenceOfVectors<V,M>::unifiedUniformlySampledCdf()"
                           << std::endl;
  }

  V unifiedMinDomainValues(m_vectorSpace.zeroVector());
  V unifiedMaxDomainValues(m_vectorSpace.zeroVector());

  ScalarSequence<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(0,                 // initialPos
                           1,                 // spacing
                           subSequenceSize(), // numPos
                           i,
                           data);

    std::vector<double> aCdf(0);
    data.unifiedUniformlySampledCdf(m_vectorSpace.numOfProcsForStorage() == 1,
                                    (unsigned int) numEvaluationPointsVec[i],
                                    unifiedMinDomainValues[i],
                                    unifiedMaxDomainValues[i],
                                    aCdf);
    unifiedCdfValues.setOneDTable(i,aCdf);
  }

  unifiedCdfGrids.setUniformGrids(numEvaluationPointsVec,
                                  unifiedMinDomainValues,
                                  unifiedMaxDomainValues);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 10)) {
    *m_env.subDisplayFile() << "Leaving SequenceOfVectors<V,M>::unifiedUniformlySampledCdf()"
                           << std::endl;
  }

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::subMeanExtra(
  unsigned int initialPos,
  unsigned int numPos,
  V&           meanVec) const
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering SequenceOfVectors<V,M>::subMeanExtra()"
                           << ": initialPos = "        << initialPos
                           << ", numPos = "            << numPos
                           << ", sub sequence size = " << this->subSequenceSize()
                           << std::endl;
  }

  bool bRC = ((initialPos          <  this->subSequenceSize()) &&
              (0                   <  numPos                 ) &&
              ((initialPos+numPos) <= this->subSequenceSize()) &&
              (this->vectorSizeLocal()  == meanVec.sizeLocal()         ));
  if ((bRC == false) && (m_env.subDisplayFile())) {
    *m_env.subDisplayFile() << "In SequenceOfVectors<V,M>::subMeanExtra()"
                           << ", initialPos = "              << initialPos
                           << ", this->subSequenceSize() = " << this->subSequenceSize()
                           << ", numPos = "                  << numPos
                           << ", this->vectorSizeLocal() = " << this->vectorSizeLocal()
                           << ", meanVec.sizeLocal() = "     << meanVec.sizeLocal()
                           << std::endl;
  }
  queso_require_msg(bRC, "invalid input data");

  ScalarSequence<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    meanVec[i] = data.subMeanExtra(0,
                                   numPos);
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving SequenceOfVectors<V,M>::subMeanExtra()"
                           << ": initialPos = "        << initialPos
                           << ", numPos = "            << numPos
                           << ", sub sequence size = " << this->subSequenceSize()
                           << ", meanVec = "           << meanVec
                           << std::endl;
  }

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::unifiedMeanExtra(
  unsigned int initialPos,
  unsigned int numPos,
  V&           unifiedMeanVec) const
{
  unsigned int tmpUnif = this->unifiedSequenceSize();
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering SequenceOfVectors<V,M>::unifiedMeanExtra()"
                            << ": initialPos = "            << initialPos
                            << ", numPos = "                << numPos
                            << ", sub sequence size = "     << this->subSequenceSize()
                            << ", unified sequence size = " << tmpUnif
                            << std::endl;
  }

  bool bRC = ((initialPos              <  this->subSequenceSize()   ) &&
              (0                       <  numPos                    ) &&
              ((initialPos+numPos)     <= this->subSequenceSize()   ) &&
              (this->vectorSizeLocal() == unifiedMeanVec.sizeLocal()));
  if ((bRC == false) && (m_env.subDisplayFile())) {
    *m_env.subDisplayFile() << "In SequenceOfVectors<V,M>::unifiedMeanExtra()"
                            << ", initialPos = "                 << initialPos
                            << ", this->subSequenceSize() = "    << this->subSequenceSize()
                            << ", numPos = "                     << numPos
                            << ", this->vectorSizeLocal() = "    << this->vectorSizeLocal()
                            << ", unifiedMeanVec.sizeLocal() = " << unifiedMeanVec.sizeLocal()
                            << std::endl;
  }
  queso_require_msg(bRC, "invalid input data");

  ScalarSequence<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    unifiedMeanVec[i] = data.unifiedMeanExtra(m_vectorSpace.numOfProcsForStorage() == 1,
                                              0,
                                              numPos);
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving SequenceOfVectors<V,M>::unifiedMeanExtra()"
                           << ": initialPos = "            << initialPos
                           << ", numPos = "                << numPos
                           << ", sub sequence size = "     << this->subSequenceSize()
                           << ", unified sequence size = " << tmpUnif
                           << ", unifiedMeanVec = "        << unifiedMeanVec
                           << std::endl;
  }

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::subMedianExtra(
  unsigned int initialPos,
  unsigned int numPos,
  V&           medianVec) const
{
  if (this->subSequenceSize() == 0) return;

  bool bRC = ((initialPos          <  this->subSequenceSize()) &&
              (0                   <  numPos                 ) &&
              ((initialPos+numPos) <= this->subSequenceSize()));
  if (bRC == false) {
    std::cerr << "In SequenceOfVectors<V,M>::subMedianExtra()"
              << ": ERROR at fullRank "         << m_env.fullRank()
              << ", initialPos = "              << initialPos
              << ", numPos = "                  << numPos
              << ", this->subSequenceSize() = " << this->subSequenceSize()
              << std::endl;
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In SequenceOfVectors<V,M>::subMedianExtra()"
                              << ": ERROR at fullRank "         << m_env.fullRank()
                              << ", initialPos = "              << initialPos
                              << ", numPos = "                  << numPos
                              << ", this->subSequenceSize() = " << this->subSequenceSize()
                              << std::endl;
    }
  }
  queso_require_msg(bRC, "invalid input data");

  ScalarSequence<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    medianVec[i] = data.subMedianExtra(0,
                                       numPos);
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving SequenceOfVectors<V,M>::subMedianExtra()"
                           << ": initialPos = "        << initialPos
                           << ", numPos = "            << numPos
                           << ", sub sequence size = " << this->subSequenceSize()
                           << ", medianVec = "         << medianVec
                           << std::endl;
  }

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::unifiedMedianExtra(
  unsigned int initialPos,
  unsigned int numPos,
  V&           unifiedMedianVec) const
{
  unsigned int tmpUnif = this->unifiedSequenceSize();
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering SequenceOfVectors<V,M>::unifiedMedianExtra()"
                            << ": initialPos = "            << initialPos
                            << ", numPos = "                << numPos
                            << ", sub sequence size = "     << this->subSequenceSize()
                            << ", unified sequence size = " << tmpUnif
                            << std::endl;
  }

  bool bRC = ((initialPos              <  this->subSequenceSize()     ) &&
              (0                       <  numPos                      ) &&
              ((initialPos+numPos)     <= this->subSequenceSize()     ) &&
              (this->vectorSizeLocal() == unifiedMedianVec.sizeLocal()));
  if ((bRC == false) && (m_env.subDisplayFile())) {
    *m_env.subDisplayFile() << "In SequenceOfVectors<V,M>::unifiedMedianExtra()"
                            << ", initialPos = "                   << initialPos
                            << ", this->subSequenceSize() = "      << this->subSequenceSize()
                            << ", numPos = "                       << numPos
                            << ", this->vectorSizeLocal() = "      << this->vectorSizeLocal()
                            << ", unifiedMedianVec.sizeLocal() = " << unifiedMedianVec.sizeLocal()
                            << std::endl;
  }
  queso_require_msg(bRC, "invalid input data");

  ScalarSequence<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    unifiedMedianVec[i] = data.unifiedMedianExtra(m_vectorSpace.numOfProcsForStorage() == 1,
                                                  0,
                                                  numPos);
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving SequenceOfVectors<V,M>::unifiedMedianExtra()"
                            << ": initialPos = "            << initialPos
                            << ", numPos = "                << numPos
                            << ", sub sequence size = "     << this->subSequenceSize()
                            << ", unified sequence size = " << tmpUnif
                            << ", unifiedMedianVec = "      << unifiedMedianVec
                            << std::endl;
  }

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::subSampleVarianceExtra(
  unsigned int initialPos,
  unsigned int numPos,
  const V&     meanVec,
  V&           samVec) const
{
  bool bRC = ((initialPos              <  this->subSequenceSize()) &&
              (0                       <  numPos                 ) &&
              ((initialPos+numPos)     <= this->subSequenceSize()) &&
              (this->vectorSizeLocal() == meanVec.sizeLocal()    ) &&
              (this->vectorSizeLocal() == samVec.sizeLocal()     ));
  queso_require_msg(bRC, "invalid input data");

  ScalarSequence<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    samVec[i] = data.subSampleVarianceExtra(0,
                                            numPos,
                                            meanVec[i]);
  }

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::unifiedSampleVarianceExtra(
  unsigned int initialPos,
  unsigned int numPos,
  const V&     unifiedMeanVec,
  V&           unifiedSamVec) const
{
  bool bRC = ((initialPos              <  this->subSequenceSize()   ) &&
              (0                       <  numPos                    ) &&
              ((initialPos+numPos)     <= this->subSequenceSize()   ) &&
              (this->vectorSizeLocal() == unifiedMeanVec.sizeLocal()) &&
              (this->vectorSizeLocal() == unifiedSamVec.sizeLocal() ));
  queso_require_msg(bRC, "invalid input data");

  ScalarSequence<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    unifiedSamVec[i] = data.unifiedSampleVarianceExtra(m_vectorSpace.numOfProcsForStorage() == 1,
                                                       0,
                                                       numPos,
                                                       unifiedMeanVec[i]);
  }

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::subSampleStd(
  unsigned int initialPos,
  unsigned int numPos,
  const V&     meanVec,
  V&           stdvec) const
{
  bool bRC = ((initialPos              <  this->subSequenceSize()) &&
              (0                       <  numPos                 ) &&
              ((initialPos+numPos)     <= this->subSequenceSize()) &&
              (this->vectorSizeLocal() == meanVec.sizeLocal()    ) &&
              (this->vectorSizeLocal() == stdvec.sizeLocal()     ));
  queso_require_msg(bRC, "invalid input data");

  ScalarSequence<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    stdvec[i] = data.subSampleStd(0,
                                  numPos,
                                  meanVec[i]);
  }

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::unifiedSampleStd(
  unsigned int initialPos,
  unsigned int numPos,
  const V&     unifiedMeanVec,
  V&           unifiedStdVec) const
{
  bool bRC = ((initialPos              <  this->subSequenceSize()   ) &&
              (0                       <  numPos                    ) &&
              ((initialPos+numPos)     <= this->subSequenceSize()   ) &&
              (this->vectorSizeLocal() == unifiedMeanVec.sizeLocal()) &&
              (this->vectorSizeLocal() == unifiedStdVec.sizeLocal() ));
  queso_require_msg(bRC, "invalid input data");

  ScalarSequence<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    unifiedStdVec[i] = data.unifiedSampleStd(m_vectorSpace.numOfProcsForStorage() == 1,
                                             0,
                                             numPos,
                                             unifiedMeanVec[i]);
  }

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::subPopulationVariance(
  unsigned int initialPos,
  unsigned int numPos,
  const V&     meanVec,
  V&           popVec) const
{
  bool bRC = ((initialPos              <  this->subSequenceSize()) &&
              (0                       <  numPos                 ) &&
              ((initialPos+numPos)     <= this->subSequenceSize()) &&
              (this->vectorSizeLocal() == meanVec.sizeLocal()    ) &&
              (this->vectorSizeLocal() == popVec.sizeLocal()     ));
  queso_require_msg(bRC, "invalid input data");

  ScalarSequence<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    popVec[i] = data.subPopulationVariance(0,
                                           numPos,
                                           meanVec[i]);
  }

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::unifiedPopulationVariance(
  unsigned int initialPos,
  unsigned int numPos,
  const V&     unifiedMeanVec,
  V&           unifiedPopVec) const
{
  bool bRC = ((initialPos              <  this->subSequenceSize()   ) &&
              (0                       <  numPos                    ) &&
              ((initialPos+numPos)     <= this->subSequenceSize()   ) &&
              (this->vectorSizeLocal() == unifiedMeanVec.sizeLocal()) &&
              (this->vectorSizeLocal() == unifiedPopVec.sizeLocal() ));
  queso_require_msg(bRC, "invalid input data");

  ScalarSequence<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    unifiedPopVec[i] = data.unifiedPopulationVariance(m_vectorSpace.numOfProcsForStorage() == 1,
                                                      0,
                                                      numPos,
                                                      unifiedMeanVec[i]);
  }

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::autoCovariance(
  unsigned int initialPos,
  unsigned int numPos,
  const V&     meanVec,
  unsigned int lag,
  V&           covVec) const
{
  bool bRC = ((initialPos              <  this->subSequenceSize()) &&
              (0                       <  numPos                 ) &&
              ((initialPos+numPos)     <= this->subSequenceSize()) &&
              (this->vectorSizeLocal() == meanVec.sizeLocal()    ) &&
              (lag                     <  numPos                 ) && // lag should not be too large
              (this->vectorSizeLocal() == covVec.sizeLocal()     ));
  queso_require_msg(bRC, "invalid input data");

  ScalarSequence<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    covVec[i] = data.autoCovariance(0,
                                    numPos,
                                    meanVec[i],
                                    lag);
  }

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::autoCorrViaDef(
  unsigned int initialPos,
  unsigned int numPos,
  unsigned int lag,
  V&           corrVec) const
{
  bool bRC = ((initialPos              <  this->subSequenceSize()) &&
              (0                       <  numPos                 ) &&
              ((initialPos+numPos)     <= this->subSequenceSize()) &&
              (lag                     <  numPos                 ) && // lag should not be too large
              (this->vectorSizeLocal() == corrVec.sizeLocal()    ));
  queso_require_msg(bRC, "invalid input data");

  ScalarSequence<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    corrVec[i] = data.autoCorrViaDef(0,
                                     numPos,
                                     lag);
  }

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::autoCorrViaFft(
  unsigned int                     initialPos,
  unsigned int                     numPos,
  const std::vector<unsigned int>& lags,
  std::vector<V*>&                 corrVecs) const
{
  bool bRC = ((initialPos          <  this->subSequenceSize()) &&
              (0                   <  numPos                 ) &&
              ((initialPos+numPos) <= this->subSequenceSize()) &&
              (0                   <  lags.size()            ) &&
              (lags[lags.size()-1] <  numPos                 )); // lag should not be too large
  queso_require_msg(bRC, "invalid input data");

  for (unsigned int j = lags.size(); j < corrVecs.size(); ++j) {
    if (corrVecs[j] != NULL) {
      delete corrVecs[j];
      corrVecs[j] = NULL;
    }
  }
  corrVecs.resize(lags.size(),NULL);
  for (unsigned int j = 0;           j < corrVecs.size(); ++j) {
    if (corrVecs[j] == NULL) corrVecs[j] = new V(m_vectorSpace.zeroVector());
  }

  ScalarSequence<double> data(m_env,0,"");
  unsigned int maxLag = lags[lags.size()-1];
  std::vector<double> autoCorrs(maxLag+1,0.); // Yes, +1

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);

    //if (m_env.subDisplayFile()) {
    //  *m_env.subDisplayFile() << "In SequenceOfVectors<V,M>::autoCorrViaFft()"
    //                         << ": about to call data.autoCorrViaFft() for paramId = " << i
    //                         << ", with numPos = "      << numPos
    //                         << ", maxLag = "           << maxLag
    //                         << ", autoCorrs.size() = " << autoCorrs.size()
    //                         << std::endl;
    //}
    data.autoCorrViaFft(0,
                        numPos,
                        maxLag,
                        autoCorrs);

    for (unsigned int j = 0; j < lags.size(); ++j) {
      (*(corrVecs[j]))[i] = autoCorrs[lags[j]];
    }
  }

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::autoCorrViaFft(
  unsigned int initialPos,
  unsigned int numPos,
  unsigned int numSum,
  V&           autoCorrsSumVec) const
{
  bool bRC = ((initialPos             <  this->subSequenceSize()) &&
              (0                      <  numPos                 ) &&
              ((initialPos+numPos)    <= this->subSequenceSize()) &&
              (0                      <  numSum                 ) &&
              (numSum                 <= numPos                 ) &&
              (autoCorrsSumVec.sizeLocal() == this->vectorSizeLocal()));
  queso_require_msg(bRC, "invalid input data");

  ScalarSequence<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);

    data.autoCorrViaFft(0,
                        numPos,
                        numSum,
                        autoCorrsSumVec[i]);
  }

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::subMinMaxExtra(
  unsigned int initialPos,
  unsigned int numPos,
  V&           minVec,
  V&           maxVec) const
{
  bool bRC = ((0                       <  numPos                 ) &&
              ((initialPos+numPos)     <= this->subSequenceSize()) &&
              (this->vectorSizeLocal() == minVec.sizeLocal()     ) &&
              (this->vectorSizeLocal() == maxVec.sizeLocal()     ));
  queso_require_msg(bRC, "invalid input data");

  //unsigned int numPos = this->subSequenceSize() - initialPos;
  unsigned int numParams = this->vectorSizeLocal();
  ScalarSequence<double> data(m_env,0,"");

  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    data.subMinMaxExtra(0,numPos,minVec[i],maxVec[i]);
  }

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::unifiedMinMaxExtra(
  unsigned int initialPos,
  unsigned int numPos,
  V&           unifiedMinVec,
  V&           unifiedMaxVec) const
{
  bool bRC = ((0                       <  numPos                   ) &&
              ((initialPos+numPos)     <= this->subSequenceSize()  ) &&
              (this->vectorSizeLocal() == unifiedMinVec.sizeLocal()) &&
              (this->vectorSizeLocal() == unifiedMaxVec.sizeLocal()));
  queso_require_msg(bRC, "invalid input data");

  //unsigned int numPos = this->subSequenceSize() - initialPos;
  unsigned int numParams = this->vectorSizeLocal();
  ScalarSequence<double> data(m_env,0,"");

  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    data.unifiedMinMaxExtra(m_vectorSpace.numOfProcsForStorage() == 1,
                            0,
                            numPos,
                            unifiedMinVec[i],
                            unifiedMaxVec[i]);
  }

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::subHistogram(
  unsigned int     initialPos,
  const V&         minVec,
  const V&         maxVec,
  std::vector<V*>& centersForAllBins,
  std::vector<V*>& quanttsForAllBins) const
{
  bool bRC = ((initialPos               <  this->subSequenceSize() ) &&
              (this->vectorSizeLocal()  == minVec.sizeLocal()      ) &&
              (this->vectorSizeLocal()  == maxVec.sizeLocal()      ) &&
              (0                        <  centersForAllBins.size()) &&
              (centersForAllBins.size() == quanttsForAllBins.size()));
  queso_require_msg(bRC, "invalid input data");

  for (unsigned int j = 0; j < quanttsForAllBins.size(); ++j) {
    centersForAllBins[j] = new V(m_vectorSpace.zeroVector());
    quanttsForAllBins [j] = new V(m_vectorSpace.zeroVector());
  }

  unsigned int dataSize = this->subSequenceSize() - initialPos;
  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    ScalarSequence<double> data(m_env,dataSize,"");
    for (unsigned int j = 0; j < dataSize; ++j) {
      data[j] = (*(m_seq[initialPos+j]))[i];
    }

    std::vector<double      > centers(centersForAllBins.size(),0.);
    std::vector<unsigned int> quantts(quanttsForAllBins.size(), 0 );
    data.subHistogram(0,
                      minVec[i],
                      maxVec[i],
                      centers,
                      quantts);

    for (unsigned int j = 0; j < quantts.size(); ++j) {
      (*(centersForAllBins[j]))[i] = centers[j];
      (*(quanttsForAllBins[j]))[i] = (double) quantts[j];
    }
  }

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::unifiedHistogram(
  unsigned int     initialPos,
  const V&         unifiedMinVec,
  const V&         unifiedMaxVec,
  std::vector<V*>& unifiedCentersForAllBins,
  std::vector<V*>& unifiedQuanttsForAllBins) const
{
  bool bRC = ((initialPos                      <  this->subSequenceSize()        ) &&
              (this->vectorSizeLocal()         == unifiedMinVec.sizeLocal()      ) &&
              (this->vectorSizeLocal()         == unifiedMaxVec.sizeLocal()      ) &&
              (0                               <  unifiedCentersForAllBins.size()) &&
              (unifiedCentersForAllBins.size() == unifiedQuanttsForAllBins.size()));
  queso_require_msg(bRC, "invalid input data");

  for (unsigned int j = 0; j < unifiedQuanttsForAllBins.size(); ++j) {
    unifiedCentersForAllBins[j] = new V(m_vectorSpace.zeroVector());
    unifiedQuanttsForAllBins [j] = new V(m_vectorSpace.zeroVector());
  }

  unsigned int dataSize = this->subSequenceSize() - initialPos;
  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    ScalarSequence<double> data(m_env,dataSize,"");
    for (unsigned int j = 0; j < dataSize; ++j) {
      data[j] = (*(m_seq[initialPos+j]))[i];
    }

    std::vector<double      > unifiedCenters(unifiedCentersForAllBins.size(),0.);
    std::vector<unsigned int> unifiedQuantts(unifiedQuanttsForAllBins.size(), 0 );
    data.unifiedHistogram(m_vectorSpace.numOfProcsForStorage() == 1,
                          0,
                          unifiedMinVec[i],
                          unifiedMaxVec[i],
                          unifiedCenters,
                          unifiedQuantts);

    for (unsigned int j = 0; j < unifiedQuantts.size(); ++j) {
      (*(unifiedCentersForAllBins[j]))[i] = unifiedCenters[j];
      (*(unifiedQuanttsForAllBins[j]))[i] = (double) unifiedQuantts[j];
    }
  }

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::subInterQuantileRange(
  unsigned int initialPos,
  V&           iqrVec) const
{
  bool bRC = ((initialPos              <  this->subSequenceSize()) &&
              (this->vectorSizeLocal() == iqrVec.sizeLocal()     ));
  queso_require_msg(bRC, "invalid input data");

  unsigned int numPos = this->subSequenceSize() - initialPos;
  ScalarSequence<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    iqrVec[i] = data.subInterQuantileRange(0);
  }

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::unifiedInterQuantileRange(
  unsigned int initialPos,
  V&           unifiedIqrVec) const
{
  bool bRC = ((initialPos              <  this->subSequenceSize()  ) &&
              (this->vectorSizeLocal() == unifiedIqrVec.sizeLocal()));
  queso_require_msg(bRC, "invalid input data");

  unsigned int numPos = this->subSequenceSize() - initialPos;
  ScalarSequence<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    unifiedIqrVec[i] = data.unifiedInterQuantileRange(m_vectorSpace.numOfProcsForStorage() == 1,
                                                      0);
  }

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::subScalesForKde(
  unsigned int initialPos,
  const V&     iqrVec,
  unsigned int kdeDimension,
  V&           scaleVec) const
{
  bool bRC = ((initialPos              <  this->subSequenceSize()) &&
              (this->vectorSizeLocal() == iqrVec.sizeLocal()     ) &&
              (this->vectorSizeLocal() == scaleVec.sizeLocal()   ));
  queso_require_msg(bRC, "invalid input data");

  unsigned int numPos = this->subSequenceSize() - initialPos;
  ScalarSequence<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    scaleVec[i] = data.subScaleForKde(0,
                                      iqrVec[i],
                                      kdeDimension);
  }

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::unifiedScalesForKde(
  unsigned int initialPos,
  const V&     unifiedIqrVec,
  unsigned int kdeDimension,
  V&           unifiedScaleVec) const
{
  bool bRC = ((initialPos              <  this->subSequenceSize()    ) &&
              (this->vectorSizeLocal() == unifiedIqrVec.sizeLocal()  ) &&
              (this->vectorSizeLocal() == unifiedScaleVec.sizeLocal()));
  queso_require_msg(bRC, "invalid input data");

  unsigned int numPos = this->subSequenceSize() - initialPos;
  ScalarSequence<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    unifiedScaleVec[i] = data.unifiedScaleForKde(m_vectorSpace.numOfProcsForStorage() == 1,
                                                 0,
                                                 unifiedIqrVec[i],
                                                 kdeDimension);
  }

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::subGaussian1dKde(
  unsigned int           initialPos,
  const V&               scaleVec,
  const std::vector<V*>& evalParamVecs,
  std::vector<V*>&       densityVecs) const
{
  bool bRC = ((initialPos              <  this->subSequenceSize()) &&
              (this->vectorSizeLocal() == scaleVec.sizeLocal()   ) &&
              (0                       <  evalParamVecs.size()   ) &&
              (evalParamVecs.size()    == densityVecs.size()     ));
  queso_require_msg(bRC, "invalid input data");

  unsigned int numPos = this->subSequenceSize() - initialPos;
  ScalarSequence<double> data(m_env,0,"");

  unsigned int numEvals = evalParamVecs.size();
  for (unsigned int j = 0; j < numEvals; ++j) {
    densityVecs[j] = new V(m_vectorSpace.zeroVector());
  }
  std::vector<double> evalParams(numEvals,0.);
  std::vector<double> densities  (numEvals,0.);

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);

    for (unsigned int j = 0; j < numEvals; ++j) {
      evalParams[j] = (*evalParamVecs[j])[i];
    }

    data.subGaussian1dKde(0,
                          scaleVec[i],
                          evalParams,
                          densities);

    for (unsigned int j = 0; j < numEvals; ++j) {
      (*densityVecs[j])[i] = densities[j];
    }
  }

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::unifiedGaussian1dKde(
  unsigned int           initialPos,
  const V&               unifiedScaleVec,
  const std::vector<V*>& unifiedEvalParamVecs,
  std::vector<V*>&       unifiedDensityVecs) const
{
  bool bRC = ((initialPos                  <  this->subSequenceSize()    ) &&
              (this->vectorSizeLocal()     == unifiedScaleVec.sizeLocal()) &&
              (0                           <  unifiedEvalParamVecs.size()) &&
              (unifiedEvalParamVecs.size() == unifiedDensityVecs.size()  ));
  queso_require_msg(bRC, "invalid input data");

  unsigned int numPos = this->subSequenceSize() - initialPos;
  ScalarSequence<double> data(m_env,0,"");

  unsigned int numEvals = unifiedEvalParamVecs.size();
  for (unsigned int j = 0; j < numEvals; ++j) {
    unifiedDensityVecs[j] = new V(m_vectorSpace.zeroVector());
  }
  std::vector<double> unifiedEvalParams(numEvals,0.);
  std::vector<double> unifiedDensities (numEvals,0.);

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);

    for (unsigned int j = 0; j < numEvals; ++j) {
      unifiedEvalParams[j] = (*unifiedEvalParamVecs[j])[i];
    }

    data.unifiedGaussian1dKde(m_vectorSpace.numOfProcsForStorage() == 1,
                              0,
                              unifiedScaleVec[i],
                              unifiedEvalParams,
                              unifiedDensities);

    for (unsigned int j = 0; j < numEvals; ++j) {
      (*unifiedDensityVecs[j])[i] = unifiedDensities[j];
    }
  }

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::subWriteContents(
  unsigned int                  initialPos,
  unsigned int                  numPos,
  const std::string&            fileName,
  const std::string&            fileType,
  const std::set<unsigned int>& allowedSubEnvIds) const
{
  queso_require_greater_equal_msg(m_env.subRank(), 0, "unexpected subRank");

  FilePtrSetStruct filePtrSet;
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
    *m_env.subDisplayFile() << "In SequenceOfVectors<V,M>::subWriteContents()"
                            << ": about to try to open file '" << fileName << "." << fileType
                            << "'"
                            << ", initialPos = " << initialPos
                            << ", numPos = "     << numPos
                            << std::endl;
  }
  if (m_env.openOutputFile(fileName,
                           fileType,
                           allowedSubEnvIds,
                           false, // A 'true' causes problems when the user chooses (via options
                                  // in the input file) to use just one file for all outputs.
                           filePtrSet)) {
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
      *m_env.subDisplayFile() << "In SequenceOfVectors<V,M>::subWriteContents()"
                              << ": successfully opened file '" << fileName << "." << fileType
                              << "'"
                              << std::endl;
    }
    this->subWriteContents(initialPos,
                           numPos,
                           filePtrSet,
                           fileType);
    m_env.closeFile(filePtrSet,fileType);
  }
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
    *m_env.subDisplayFile() << "In SequenceOfVectors<V,M>::subWriteContents()"
                            << ": before Barrier()"
                            << std::endl;
  }
  m_env.subComm().Barrier();

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::subWriteContents(
  unsigned int        initialPos,
  unsigned int        numPos,
  FilePtrSetStruct& filePtrSet,
  const std::string&  fileType) const
{
  if (fileType == UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT ||
      fileType == UQ_FILE_EXTENSION_FOR_TXT_FORMAT) {
    queso_require_msg(filePtrSet.ofsVar, "filePtrSet.ofsVar should not be NULL");
    this->subWriteContents(initialPos,
                           numPos,
                           *filePtrSet.ofsVar,
                           fileType);
  }
#ifdef QUESO_HAS_HDF5
  else if (fileType == UQ_FILE_EXTENSION_FOR_HDF_FORMAT) {

    // Check the file is still legit
    queso_require_greater_equal_msg(
        filePtrSet.h5Var,
        0,
        "filePtrSet.h5Var should not be non-negative");

    // Sanity check extent
    queso_require_less_equal_msg(
        (initialPos+numPos),
        this->subSequenceSize(),
        "invalid routine input parameters");

    unsigned int numParams = m_vectorSpace.dimLocal();
    unsigned int chainSize = this->subSequenceSize();
    hsize_t dims[2] = { chainSize, numParams };

    // Create dataspace
    hid_t dataspace_id = H5Screate_simple(2, dims, dims);

    // Sanity check dataspace
    queso_require_greater_equal_msg(
        dataspace_id,
        0,
        "error creating dataspace with id: " << dataspace_id);

    // Create dataset
    hid_t dataset_id = H5Dcreate(filePtrSet.h5Var,
                                 "data",
                                 H5T_IEEE_F64LE,
                                 dataspace_id,
                                 H5P_DEFAULT,
                                 H5P_DEFAULT,
                                 H5P_DEFAULT);

    // Sanity check dataset
    queso_require_greater_equal_msg(
        dataset_id,
        0,
        "error creating dataset with id: " << dataset_id);

    // This is so egregiously awfully badly terribly horrific I want to die.
    //
    // And, of course, if any of the subsequent H5* sanity check fail we
    // throw an exception and leak a metric fuckton of memory.
    double * data = (double *)malloc(numParams * chainSize * sizeof(double));

    for (unsigned int i = 0; i < chainSize; i++) {
      V tmpVec(*(m_seq[i]));
      for (unsigned int j = 0; j < numParams; j++) {
        data[numParams*i+j] = tmpVec[j];
      }
    }

    // Write the dataset
    herr_t status = H5Dwrite(
        dataset_id,
        H5T_NATIVE_DOUBLE,  // The type in memory
        H5S_ALL,  // The dataspace in memory
        dataspace_id,  // The file dataspace
        H5P_DEFAULT,  // Xfer property list
        data);

    // Sanity check the write
    queso_require_greater_equal_msg(
        status,
        0,
        "error writing dataset to file with id: " << filePtrSet.h5Var);

    // Clean up
    free(data);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    // Should we close the file too?  It's unclear.

  }
#endif
  else {
    queso_error_msg("invalid file type");
  }

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::subWriteContents(
  unsigned int       initialPos,
  unsigned int       numPos,
  std::ofstream&     ofs,
  const std::string& fileType) const
{
  queso_require_less_equal_msg((initialPos+numPos), this->subSequenceSize(), "invalid routine input parameters");

  if (initialPos == 0) {
    if (fileType == UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT) {
      // Need sub matlab header here since this is subWriteContents
      this->writeSubMatlabHeader(ofs,
                                 this->subSequenceSize(),
                                 this->vectorSizeLocal());
    }
    else if (fileType == UQ_FILE_EXTENSION_FOR_TXT_FORMAT) {
      this->writeTxtHeader(ofs,
                           this->subSequenceSize(),
                           this->vectorSizeLocal());
    }
  }

  for (unsigned int j = initialPos; j < initialPos+numPos; ++j) {
    bool savedVectorPrintScientific = m_seq[j]->getPrintScientific();
    bool savedVectorPrintState      = m_seq[j]->getPrintHorizontally();
    m_seq[j]->setPrintScientific  (true);
    m_seq[j]->setPrintHorizontally(true);

    ofs << *(m_seq[j])
        << std::endl;

    m_seq[j]->setPrintHorizontally(savedVectorPrintState);
    m_seq[j]->setPrintScientific  (savedVectorPrintScientific);
  }

  // Write Matlab-specific ending if desired
  if ((fileType == UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT) &&
      ((initialPos + numPos) == this->subSequenceSize())) {
    ofs << "];\n";
  }
}

template <class V, class M>
void
SequenceOfVectors<V,M>::writeSubMatlabHeader(std::ofstream & ofs,
    double sequenceSize, double vectorSizeLocal) const
{
  ofs << m_name << "_sub" << m_env.subIdString() << " = zeros(" << sequenceSize
      << ","                                                    << vectorSizeLocal
      << ");"
      << std::endl;
  ofs << m_name << "_sub" << m_env.subIdString() << " = [";
}

template <class V, class M>
void
SequenceOfVectors<V,M>::writeUnifiedMatlabHeader(std::ofstream & ofs,
    double sequenceSize, double vectorSizeLocal) const
{
  ofs << m_name << "_unified" << " = zeros(" << sequenceSize
                            << ","           << vectorSizeLocal
                            << ");"
                            << std::endl;
  ofs<< m_name << "_unified" << " = [";
}

template <class V, class M>
void
SequenceOfVectors<V,M>::writeTxtHeader(std::ofstream & ofs,
    double sequenceSize, double vectorSizeLocal) const
{
  ofs << sequenceSize << " " << vectorSizeLocal
      << std::endl;
}

//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::unifiedWriteContents(
  const std::string& fileName,
  const std::string& inputFileType) const
{
  std::string fileType(inputFileType);
#ifdef QUESO_HAS_HDF5
  // Do nothing
#else
  if (fileType == UQ_FILE_EXTENSION_FOR_HDF_FORMAT) {
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "WARNING in SequenceOfVectors<V,M>::unifiedWriteContents()"
                              << ": file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                              << "' has been requested, but this QUESO library has not been built with 'hdf5'"
                              << ". Code will therefore process the file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                              << "' instead..."
                              << std::endl;
    }
    if (m_env.subRank() == 0) {
      std::cerr << "WARNING in SequenceOfVectors<V,M>::unifiedWriteContents()"
                << ": file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                << "' has been requested, but this QUESO library has not been built with 'hdf5'"
                << ". Code will therefore process the file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                << "' instead..."
                << std::endl;
    }
    fileType = UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT;
  }
#endif

  // All processors in 'fullComm' should call this routine...

  //m_env.fullComm().Barrier(); // Dangerous to barrier on fullComm ... // prudenci-2011-01-17
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 10)) {
    *m_env.subDisplayFile() << "Entering SequenceOfVectors<V,M>::unifiedWriteContents()"
                            << ": worldRank "      << m_env.worldRank()
                            << ", fullRank "       << m_env.fullRank()
                            << ", subEnvironment " << m_env.subId()
                            << ", subRank "        << m_env.subRank()
                            << ", inter0Rank "     << m_env.inter0Rank()
                          //<< ", m_env.inter0Comm().NumProc() = " << m_env.inter0Comm().NumProc()
                            << ", fileName = "     << fileName
                            << std::endl;
  }

  if (m_env.inter0Rank() >= 0) {
    for (unsigned int r = 0; r < (unsigned int) m_env.inter0Comm().NumProc(); ++r) {
      if (m_env.inter0Rank() == (int) r) {
        // My turn
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 10)) {
          *m_env.subDisplayFile() << "In SequenceOfVectors<V,M>::unifiedWriteContents()"
                                  << ": worldRank "      << m_env.worldRank()
                                  << ", fullRank "       << m_env.fullRank()
                                  << ", subEnvironment " << m_env.subId()
                                  << ", subRank "        << m_env.subRank()
                                  << ", inter0Rank "     << m_env.inter0Rank()
                                //<< ", m_env.inter0Comm().NumProc() = " << m_env.inter0Comm().NumProc()
                                  << ", fileName = "     << fileName
                                  << ", about to open file for r = " << r
                                  << std::endl;
        }

        // bool writeOver = (r == 0);
        bool writeOver = false; // A 'true' causes problems when the user chooses (via options
                                // in the input file) to use just one file for all outputs.
        FilePtrSetStruct unifiedFilePtrSet;
        if (m_env.openUnifiedOutputFile(fileName,
                                        fileType, // "m or hdf"
                                        writeOver,
                                        unifiedFilePtrSet)) {
          if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 10)) { // 2013-02-23
            *m_env.subDisplayFile() << "In SequenceOfVectors<V,M>::unifiedWriteContents()"
                                    << ": worldRank "      << m_env.worldRank()
                                    << ", fullRank "       << m_env.fullRank()
                                    << ", subEnvironment " << m_env.subId()
                                    << ", subRank "        << m_env.subRank()
                                    << ", inter0Rank "     << m_env.inter0Rank()
                                  //<< ", m_env.inter0Comm().NumProc() = " << m_env.inter0Comm().NumProc()
                                    << ", fileName = "     << fileName
                                    << ", just opened file for r = " << r
                                    << std::endl;
          }

          unsigned int chainSize = this->subSequenceSize();
          if ((fileType == UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT) ||
              (fileType == UQ_FILE_EXTENSION_FOR_TXT_FORMAT)) {
            if (r == 0) {
              if (fileType == UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT) {
                // Need unified matlab header here since this is unifiedWriteContents
                writeUnifiedMatlabHeader(*unifiedFilePtrSet.ofsVar,
                    this->subSequenceSize()*m_env.inter0Comm().NumProc(),
                    this->vectorSizeLocal());
              }
              else {  // If we get here it's definitely a txt file not matlab
                writeTxtHeader(*unifiedFilePtrSet.ofsVar,
                    this->subSequenceSize()*m_env.inter0Comm().NumProc(),
                    this->vectorSizeLocal());
              }
            }

            for (unsigned int j = 0; j < chainSize; ++j) { // 2013-02-23
        //std::cout << "In SequenceOfVectors<V,M>::unifiedWriteContents(): m_seq[" << j << "] = " << m_seq[j]
              //          << std::endl;
            //std::cout << "In SequenceOfVectors<V,M>::unifiedWriteContents(): &(m_seq[" << j << "].map()) = " << &(m_seq[j]->map())
              //          << std::endl;
              //std::cout << "In SequenceOfVectors<V,M>::unifiedWriteContents(): (m_seq[" << j << "].map().NumMyElements = " << m_seq[j]->map().NumMyElements()
              //          << std::endl;
              //V tmpVec(*(m_seq[j]));
        //std::cout << "*(m_seq[" << j << "]) = " << tmpVec
              //          << std::endl;
        //std::cout << "*(m_seq[" << j << "]) = " << *(m_seq[j])
              //          << std::endl;

              bool savedVectorPrintScientific = m_seq[j]->getPrintScientific();
              bool savedVectorPrintState      = m_seq[j]->getPrintHorizontally();
              m_seq[j]->setPrintScientific  (true);
              m_seq[j]->setPrintHorizontally(true);

              *unifiedFilePtrSet.ofsVar << *(m_seq[j])
                                        << std::endl;

              m_seq[j]->setPrintHorizontally(savedVectorPrintState);
              m_seq[j]->setPrintScientific  (savedVectorPrintScientific);
            }
          }
#ifdef QUESO_HAS_HDF5
          else if (fileType == UQ_FILE_EXTENSION_FOR_HDF_FORMAT) {
            unsigned int numParams = m_vectorSpace.dimLocal();
            if (r == 0) {
              hid_t datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
              //std::cout << "In SequenceOfVectors<V,M>::unifiedWriteContents(): h5 case, data type created" << std::endl;
              hsize_t dimsf[2];
              dimsf[0] = chainSize;
              dimsf[1] = numParams;
              hid_t dataspace = H5Screate_simple(2, dimsf, NULL); // HDF5_rank = 2
              //std::cout << "In SequenceOfVectors<V,M>::unifiedWriteContents(): h5 case, data space created" << std::endl;
              hid_t dataset = H5Dcreate2(unifiedFilePtrSet.h5Var,
                                         "data",
                                         datatype,
                                         dataspace,
                                         H5P_DEFAULT,  // Link creation property list
                                         H5P_DEFAULT,  // Dataset creation property list
                                         H5P_DEFAULT); // Dataset access property list
              //std::cout << "In SequenceOfVectors<V,M>::unifiedWriteContents(): h5 case, data set created" << std::endl;

              struct timeval timevalBegin;
              int iRC = UQ_OK_RC;
              iRC = gettimeofday(&timevalBegin,NULL);
              if (iRC) {}; // just to remover compiler warning

              double * data;
              data = (double *)malloc(numParams * chainSize * sizeof(double));

              for (unsigned int i = 0; i < chainSize; ++i) {
                V tmpVec(*(m_seq[i]));
                for (unsigned int j = 0; j < numParams; ++j) {
                  data[numParams*i+j] = tmpVec[j];
                }
              }

              herr_t status;
              status = H5Dwrite(dataset,
                                H5T_NATIVE_DOUBLE,
                                H5S_ALL,
                                H5S_ALL,
                                H5P_DEFAULT,
                                data);
              if (status) {}; // just to remover compiler warning
              //std::cout << "In SequenceOfVectors<V,M>::unifiedWriteContents(): h5 case, data written" << std::endl;

              double writeTime = MiscGetEllapsedSeconds(&timevalBegin);
              if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
                *m_env.subDisplayFile() << "In SequenceOfVectors<V,M>::unifiedWriteContents()"
                                        << ": worldRank "      << m_env.worldRank()
                                        << ", fullRank "       << m_env.fullRank()
                                        << ", subEnvironment " << m_env.subId()
                                        << ", subRank "        << m_env.subRank()
                                        << ", inter0Rank "     << m_env.inter0Rank()
                                        << ", fileName = "     << fileName
                                        << ", numParams = "    << numParams
                                        << ", chainSize = "    << chainSize
                                        << ", writeTime = "    << writeTime << " seconds"
                                        << std::endl;
              }

              H5Dclose(dataset);
              //std::cout << "In SequenceOfVectors<V,M>::unifiedWriteContents(): h5 case, data set closed" << std::endl;
              H5Sclose(dataspace);
              //std::cout << "In SequenceOfVectors<V,M>::unifiedWriteContents(): h5 case, data space closed" << std::endl;
              H5Tclose(datatype);
              //std::cout << "In SequenceOfVectors<V,M>::unifiedWriteContents(): h5 case, data type closed" << std::endl;
              free(data);
            }
            else {
              queso_error_msg("hdf file type not supported for multiple sub-environments yet");
            }
          }
#endif
          else {
            queso_error_msg("invalid file type");
          }

          if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 10)) {
            *m_env.subDisplayFile() << "In SequenceOfVectors<V,M>::unifiedWriteContents()"
                                    << ": worldRank "      << m_env.worldRank()
                                    << ", fullRank "       << m_env.fullRank()
                                    << ", subEnvironment " << m_env.subId()
                                    << ", subRank "        << m_env.subRank()
                                    << ", inter0Rank "     << m_env.inter0Rank()
                                  //<< ", m_env.inter0Comm().NumProc() = " << m_env.inter0Comm().NumProc()
                                    << ", fileName = "     << fileName
                                    << ", about to close file for r = " << r
                                    << std::endl;
          }

          m_env.closeFile(unifiedFilePtrSet,fileType);

          if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 10)) {
            *m_env.subDisplayFile() << "In SequenceOfVectors<V,M>::unifiedWriteContents()"
                                    << ": worldRank "      << m_env.worldRank()
                                    << ", fullRank "       << m_env.fullRank()
                                    << ", subEnvironment " << m_env.subId()
                                    << ", subRank "        << m_env.subRank()
                                    << ", inter0Rank "     << m_env.inter0Rank()
                                  //<< ", m_env.inter0Comm().NumProc() = " << m_env.inter0Comm().NumProc()
                                    << ", fileName = "     << fileName
                                    << ", just closed file for r = " << r
                                    << std::endl;
          }
        } // if (m_env.openUnifiedOutputFile())
      } // if (m_env.inter0Rank() == (int) r)
      m_env.inter0Comm().Barrier();
    } // for r

    if (m_env.inter0Rank() == 0) {
      if ((fileType == UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT) ||
          (fileType == UQ_FILE_EXTENSION_FOR_TXT_FORMAT)) {
        FilePtrSetStruct unifiedFilePtrSet;
        if (m_env.openUnifiedOutputFile(fileName,
                                        fileType,
                                        false, // Yes, 'writeOver = false' in order to close the array for matlab
                                        unifiedFilePtrSet)) {

          if (fileType == UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT) {
            *unifiedFilePtrSet.ofsVar << "];\n";
          }

          m_env.closeFile(unifiedFilePtrSet,fileType);
        }
      }
      else if (fileType == UQ_FILE_EXTENSION_FOR_HDF_FORMAT) {
        // Do nothing
      }
      else {
        queso_error_msg("invalid file type");
      }
    }
  } // if (m_env.inter0Rank() >= 0)

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 10)) {
    *m_env.subDisplayFile() << "Leaving SequenceOfVectors<V,M>::unifiedWriteContents()"
                            << ", fileName = " << fileName
                            << std::endl;
  }
  //m_env.fullComm().Barrier(); // Dangerous to barrier on fullComm ... // prudenci-2011-01-17

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::unifiedReadContents(
  const std::string& fileName,
  const std::string& inputFileType,
  const unsigned int subReadSize)
{
  std::string fileType(inputFileType);
#ifdef QUESO_HAS_HDF5
  // Do nothing
#else
  if (fileType == UQ_FILE_EXTENSION_FOR_HDF_FORMAT) {
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "WARNING in SequenceOfVectors<V,M>::unifiedReadContents()"
                              << ": file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                              << "' has been requested, but this QUESO library has not been built with 'hdf5'"
                              << ". Code will therefore process the file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                              << "' instead..."
                              << std::endl;
    }
    if (m_env.subRank() == 0) {
      std::cerr << "WARNING in SequenceOfVectors<V,M>::unifiedReadContents()"
                << ": file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                << "' has been requested, but this QUESO library has not been built with 'hdf5'"
                << ". Code will therefore process the file format '" << UQ_FILE_EXTENSION_FOR_HDF_FORMAT
                << "' instead..."
                << std::endl;
    }
    fileType = UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT;
  }
#endif

  //m_env.fullComm().Barrier(); // Dangerous to barrier on fullComm ...
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Entering SequenceOfVectors<V,M>::unifiedReadContents()"
                            << ": worldRank "                      << m_env.worldRank()
                            << ", fullRank "                       << m_env.fullRank()
                            << ", subEnvironment "                 << m_env.subId()
                            << ", subRank "                        << m_env.subRank()
                            << ", inter0Rank "                     << m_env.inter0Rank()
      //<< ", m_env.inter0Comm().NumProc() = " << m_env.inter0Comm().NumProc()
                            << ", fileName = "                     << fileName
                            << ", subReadSize = "                  << subReadSize
      //<< ", unifiedReadSize = "              << unifiedReadSize
                            << std::endl;
  }

  this->resizeSequence(subReadSize);

  if (m_env.inter0Rank() >= 0) {
    double unifiedReadSize = subReadSize*m_env.inter0Comm().NumProc();

    // In the logic below, the id of a line' begins with value 0 (zero)
    unsigned int idOfMyFirstLine = 1 + m_env.inter0Rank()*subReadSize;
    unsigned int idOfMyLastLine = (1 + m_env.inter0Rank())*subReadSize;
    unsigned int numParams = this->vectorSizeLocal();

    for (unsigned int r = 0; r < (unsigned int) m_env.inter0Comm().NumProc(); ++r) { // "m or hdf"
      if (m_env.inter0Rank() == (int) r) {
        // My turn
        FilePtrSetStruct unifiedFilePtrSet;
        if (m_env.openUnifiedInputFile(fileName,
                                       fileType,
                                       unifiedFilePtrSet)) {
          if ((fileType == UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT) ||
              (fileType == UQ_FILE_EXTENSION_FOR_TXT_FORMAT)) {
            if (r == 0) {
              // Read number of chain positions in the file by taking care of the first line,
              // which resembles something like 'variable_name = zeros(n_positions,m_params);'
        std::string tmpString;

              // Read 'variable name' string
              *unifiedFilePtrSet.ifsVar >> tmpString;
        //std::cout << "Just read '" << tmpString << "'" << std::endl;

              // Read '=' sign
              *unifiedFilePtrSet.ifsVar >> tmpString;
          //std::cout << "Just read '" << tmpString << "'" << std::endl;
              queso_require_equal_to_msg(tmpString, std::string("="), std::string("string should be the '=' sign"));

              // Read 'zeros(n_positions,n_params)' string
              *unifiedFilePtrSet.ifsVar >> tmpString;
        //std::cout << "Just read '" << tmpString << "'" << std::endl;
              unsigned int posInTmpString = 6;

              // Isolate 'n_positions' in a string
              //char nPositionsString[tmpString.size()-posInTmpString+1]; // avoid compiler warning
        std::string nPositionsString((size_t) (tmpString.size()-posInTmpString+1),' ');
              unsigned int posInPositionsString = 0;
              do {
                queso_require_less_msg(posInTmpString, tmpString.size(), "symbol ',' not found in first line of file");
                nPositionsString[posInPositionsString++] = tmpString[posInTmpString++];
              } while (tmpString[posInTmpString] != ',');
              nPositionsString[posInPositionsString] = '\0';

              // Isolate 'n_params' in a string
              posInTmpString++; // Avoid reading ',' char
              //char nParamsString[tmpString.size()-posInTmpString+1]; // avoid compiler warning
        std::string nParamsString((size_t) (tmpString.size()-posInTmpString+1),' ');
              unsigned int posInParamsString = 0;
              do {
                queso_require_less_msg(posInTmpString, tmpString.size(), "symbol ')' not found in first line of file");
                nParamsString[posInParamsString++] = tmpString[posInTmpString++];
              } while (tmpString[posInTmpString] != ')');
              nParamsString[posInParamsString] = '\0';

              // Convert 'n_positions' and 'n_params' strings to numbers
              unsigned int sizeOfChainInFile = (unsigned int) strtod(nPositionsString.c_str(),NULL);
              unsigned int numParamsInFile   = (unsigned int) strtod(nParamsString.c_str(),   NULL);
              if (m_env.subDisplayFile()) {
                *m_env.subDisplayFile() << "In SequenceOfVectors<V,M>::unifiedReadContents()"
                                        << ": worldRank "           << m_env.worldRank()
                                        << ", fullRank "            << m_env.fullRank()
                                        << ", sizeOfChainInFile = " << sizeOfChainInFile
                                        << ", numParamsInFile = "   << numParamsInFile
                                        << std::endl;
              }

              // Check if [size of chain in file] >= [requested unified sequence size]
              queso_require_greater_equal_msg(sizeOfChainInFile, unifiedReadSize, "size of chain in file is not big enough");

              // Check if [num params in file] == [num params in current chain]
              queso_require_equal_to_msg(numParamsInFile, numParams, "number of parameters of chain in file is different than number of parameters in this chain object");
            } // if (r == 0)

            // Code common to any core in 'inter0Comm', including core of rank 0
            unsigned int maxCharsPerLine = 64*numParams; // Up to about 60 characters to represent each parameter value

            unsigned int lineId = 0;
            while (lineId < idOfMyFirstLine) {
              unifiedFilePtrSet.ifsVar->ignore(maxCharsPerLine,'\n');
              lineId++;
            };

            if (r == 0) {
              // Take care of initial part of the first data line,
              // which resembles something like 'variable_name = [value1 value2 ...'
        std::string tmpString;

              // Read 'variable name' string
              *unifiedFilePtrSet.ifsVar >> tmpString;
          //std::cout << "Core 0 just read '" << tmpString << "'" << std::endl;

              // Read '=' sign
              *unifiedFilePtrSet.ifsVar >> tmpString;
        //std::cout << "Core 0 just read '" << tmpString << "'" << std::endl;
              queso_require_equal_to_msg(tmpString, std::string("="), std::string("in core 0, string should be the '=' sign"));

              // Take into account the ' [' portion
        std::streampos tmpPos = unifiedFilePtrSet.ifsVar->tellg();
              unifiedFilePtrSet.ifsVar->seekg(tmpPos+(std::streampos)2);
            }

            V tmpVec(m_vectorSpace.zeroVector());
            while (lineId <= idOfMyLastLine) {
              for (unsigned int i = 0; i < numParams; ++i) {
                *unifiedFilePtrSet.ifsVar >> tmpVec[i];
              }
              this->setPositionValues(lineId - idOfMyFirstLine, tmpVec);
              lineId++;
            };
          }
#ifdef QUESO_HAS_HDF5
          else if (fileType == UQ_FILE_EXTENSION_FOR_HDF_FORMAT) {
            if (r == 0) {
              hid_t dataset = H5Dopen2(unifiedFilePtrSet.h5Var,
                                       "data",
                                       H5P_DEFAULT); // Dataset access property list
              hid_t datatype  = H5Dget_type(dataset);
              H5T_class_t t_class = H5Tget_class(datatype);
              queso_require_equal_to_msg(t_class, H5T_FLOAT, "t_class is not H5T_DOUBLE");
              hid_t dataspace = H5Dget_space(dataset);
              int   rank      = H5Sget_simple_extent_ndims(dataspace);
              queso_require_equal_to_msg(rank, 2, "hdf rank is not 2");
              hsize_t dims_in[2];
              int     status_n;
              status_n  = H5Sget_simple_extent_dims(dataspace, dims_in, NULL);
              if (status_n) {}; // just to remover compiler warning
        //std::cout << "In SequenceOfVectors<V,M>::unifiedReadContents()"
              //          << ": dims_in[0] = " << dims_in[0]
              //          << ", dims_in[1] = " << dims_in[1]
              //          << std::endl;
              queso_require_equal_to_msg(dims_in[0], numParams, "dims_in[0] is not equal to 'numParams'");
              queso_require_greater_equal_msg(dims_in[1], subReadSize, "dims_in[1] is smaller that requested 'subReadSize'");

              struct timeval timevalBegin;
              int iRC = UQ_OK_RC;
              iRC = gettimeofday(&timevalBegin,NULL);
              if (iRC) {}; // just to remover compiler warning

              unsigned int chainSizeIn = (unsigned int) dims_in[1];
              //double* dataIn[numParams]; // avoid compiler warning
        std::vector<double*> dataIn((size_t) numParams,NULL);
              dataIn[0] = (double*) malloc(numParams*chainSizeIn*sizeof(double));
              for (unsigned int i = 1; i < numParams; ++i) { // Yes, from '1'
                dataIn[i] = dataIn[i-1] + chainSizeIn; // Yes, just 'chainSizeIn', not 'chainSizeIn*sizeof(double)'
              }
              //std::cout << "In SequenceOfVectors<V,M>::unifiedReadContents(): h5 case, memory allocated" << std::endl;
              herr_t status;
              status = H5Dread(dataset,
                               H5T_NATIVE_DOUBLE,
                               H5S_ALL,
                               dataspace,
                               H5P_DEFAULT,
                               dataIn[0]);
              if (status) {}; // just to remover compiler warning
              //std::cout << "In SequenceOfVectors<V,M>::unifiedReadContents(): h5 case, data read" << std::endl;
              V tmpVec(m_vectorSpace.zeroVector());
              for (unsigned int j = 0; j < subReadSize; ++j) { // Yes, 'subReadSize', not 'chainSizeIn'
                for (unsigned int i = 0; i < numParams; ++i) {
                  tmpVec[i] = dataIn[i][j];
                }
                this->setPositionValues(j, tmpVec);
              }

              double readTime = MiscGetEllapsedSeconds(&timevalBegin);
              if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
                *m_env.subDisplayFile() << "In SequenceOfVectors<V,M>::unifiedReadContents()"
                                        << ": worldRank "      << m_env.worldRank()
                                        << ", fullRank "       << m_env.fullRank()
                                        << ", subEnvironment " << m_env.subId()
                                        << ", subRank "        << m_env.subRank()
                                        << ", inter0Rank "     << m_env.inter0Rank()
                                        << ", fileName = "     << fileName
                                        << ", numParams = "    << numParams
                                        << ", chainSizeIn = "  << chainSizeIn
                                        << ", subReadSize = "  << subReadSize
                                        << ", readTime = "     << readTime << " seconds"
                                        << std::endl;
              }

              H5Sclose(dataspace);
              H5Tclose(datatype);
              H5Dclose(dataset);
              //free(dataIn[0]); // related to the changes above for compiler warning
              for (unsigned int tmpIndex = 0; tmpIndex < dataIn.size(); tmpIndex++) {
                free (dataIn[tmpIndex]);
              }
            }
            else {
              queso_error_msg("hdf file type not supported for multiple sub-environments yet");
            }
          }
#endif
          else {
            queso_error_msg("invalid file type");
          }
          m_env.closeFile(unifiedFilePtrSet,fileType);
        } // if (m_env.openUnifiedInputFile())
      } // if (m_env.inter0Rank() == (int) r)
      m_env.inter0Comm().Barrier();
    } // for r
  } // if (m_env.inter0Rank() >= 0)
  else {
    V tmpVec(m_vectorSpace.zeroVector());
    for (unsigned int i = 1; i < subReadSize; ++i) {
      this->setPositionValues(i,tmpVec);
    }
  }

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Leaving SequenceOfVectors<V,M>::unifiedReadContents()"
                            << ", fileName = " << fileName
                            << std::endl;
  }
  //m_env.fullComm().Barrier(); // Dangerous to barrier on fullComm ...

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::select(const std::vector<unsigned int>& /* idsOfUniquePositions */)
{
  queso_error_msg("Code is not complete yet");

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::filter(
  unsigned int initialPos,
  unsigned int spacing)
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Entering SequenceOfVectors<V,M>::filter()"
                           << ": initialPos = "      << initialPos
                           << ", spacing = "         << spacing
                           << ", subSequenceSize = " << this->subSequenceSize()
                           << std::endl;
  }

  unsigned int i = 0;
  unsigned int j = initialPos;
  unsigned int originalSubSequenceSize = this->subSequenceSize();
  while (j < originalSubSequenceSize) {
    if (i != j) {
      //*m_env.subDisplayFile() << i << "--" << j << " ";
      delete m_seq[i];
      m_seq[i] = new V(*(m_seq[j]));
    }
    i++;
    j += spacing;
  }

  this->resetValues(i,originalSubSequenceSize-i);
  this->resizeSequence(i);

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Leaving SequenceOfVectors<V,M>::filter()"
                           << ": initialPos = "      << initialPos
                           << ", spacing = "         << spacing
                           << ", subSequenceSize = " << this->subSequenceSize()
                           << std::endl;
  }

  return;
}
//---------------------------------------------------
template <class V, class M>
double
SequenceOfVectors<V,M>::estimateConvBrooksGelman(
  unsigned int initialPos,
  unsigned int numPos) const
{
  // This method requires *at least* two sequences. Error if there is only one.
  queso_require_greater_equal_msg(m_env.numSubEnvironments(), 2, "At least two sequences required for Brooks-Gelman convergence test.");

  // TODO: Need special case for 1-dimensional parameter space.

  // Initialize with garbage to give the user a clue something is funky.
  double convMeasure = -1.0;

  // We only do the work on the subenvironment where the sequence data
  // is stashed.
  if( m_env.inter0Rank() >= 0 )
    {
      // Sanity Checks

      // REMEMBER: \psi is a *vector* of parameters
      // Get quantities we will use several times
      V psi_j_dot = m_vectorSpace.zeroVector();
      V psi_dot_dot = m_vectorSpace.zeroVector();
      V work = m_vectorSpace.zeroVector();

      // m = number of chains > 1
      // n = number of steps for which we are computing the metric
      int m = m_env.numSubEnvironments();
      int n = numPos;

      this->subMeanExtra    ( initialPos, numPos, psi_j_dot   );
      this->unifiedMeanExtra( initialPos, numPos, psi_dot_dot );

#if 0
      std::cout << "psi_j_dot = " << psi_j_dot << std::endl;
      std::cout << "psi_dot_dot = " << psi_dot_dot << std::endl;
#endif

      /* W = \frac{1}{m*(n-1)}*\sum_{j=1}^m \sum{t=1}^n
   (\psi_{jt} - \overline{\psi_{j\cdot}})*(\psi_{jt} - \overline{\psi_{j\cdot}})^T
   This corresponds to the "within-sequence" covariance matrix. */
      M* W_local = m_vectorSpace.newDiagMatrix( m_vectorSpace.zeroVector() );
      M* W = m_vectorSpace.newDiagMatrix( m_vectorSpace.zeroVector() );
      V  psi_j_t = m_vectorSpace.zeroVector();

      // Sum within the chain
      for( unsigned int t = initialPos; t < initialPos+numPos; ++t )
  {
    psi_j_t = *(m_seq[t]);

    work = psi_j_t - psi_j_dot;

    (*W_local) += matrixProduct( work, work );
  }

      // Now do the sum over the chains
      // W will be available on all inter0 processors
      W_local->mpiSum( m_env.inter0Comm(), (*W) );

      (*W) = 1.0/(double(m)*(double(n)-1.0)) * (*W);

#if 0
      std::cout << "n, m = " << n << ", " << m << std::endl;
      std::cout << "W_local = " << *W_local << std::endl;
      std::cout << "W = " << *W << std::endl;
#endif

      // Need to delete pointers to temporary covariance matrices
      delete W_local;

      /* B/n = \frac{1}{m-1}\sum_{j=1}^m
   (\overline{\psi_{j\cdot}} - \overline{\psi_{\cdot \cdot}})*
   (\overline{\psi_{j\cdot}} - \overline{\psi_{\cdot \cdot}})^T
   This corresponds to the "between-sequence" covariance matrix. */
      M* B_over_n_local = m_vectorSpace.newDiagMatrix( m_vectorSpace.zeroVector() );
      M* B_over_n = m_vectorSpace.newDiagMatrix( m_vectorSpace.zeroVector() );

      work = psi_j_dot - psi_dot_dot;
      (*B_over_n_local) = matrixProduct( work, work );

      B_over_n_local->mpiSum( m_env.inter0Comm(), (*B_over_n) );

      // Need to delete pointers to temporary covariance matrices
      delete B_over_n_local;

      (*B_over_n) = 1.0/(double(m)-1.0) * (*B_over_n);

#if 0
      std::cout << "B_over_n = " << *B_over_n << std::endl;
#endif


      /* R_p = (n-1)/n + (m+1)/m * \lambda
   \lambda = largest eigenvalue of W^{-1}*B/n */
      M* A = m_vectorSpace.newDiagMatrix( m_vectorSpace.zeroVector() );

      W->invertMultiply( *B_over_n, *A );

#if 0
      std::cout << "A = " << *A << std::endl;
      std::cout.flush();
#endif
      // Need to delete pointers to temporary covariance matrices
      delete W;
      delete B_over_n;

      double eigenValue;
      V eigenVector = m_vectorSpace.zeroVector();

      A->largestEigen( eigenValue, eigenVector );

      // Need to delete pointers to temporary covariance matrices
      delete A;

      // Now, finally compute the final convMeasure
      convMeasure = (double(n)-1.0)/double(n) + (double(m)+1.0)/double(m)*eigenValue;

    } // End of check on inter0Rank

  //TODO: Do we need a Barrier here?
  //TODO: Error checking on MPI.
  //m_env.fullComm().Barrier(); // Dangerous to barrier on fullComm ...

  return convMeasure;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::extractScalarSeq(
  unsigned int                   initialPos,
  unsigned int                   spacing,
  unsigned int                   numPos,
  unsigned int                   paramId,
  ScalarSequence<double>& scalarSeq) const
{
  scalarSeq.resizeSequence(numPos);
  if (spacing == 1) {
    for (unsigned int j = 0; j < numPos; ++j) {
      scalarSeq[j] = (*(m_seq[initialPos+j        ]))[paramId];
    }
  }
  else {
    for (unsigned int j = 0; j < numPos; ++j) {
      scalarSeq[j] = (*(m_seq[initialPos+j*spacing]))[paramId];
    }
  }

  return;
}
// Private methods ------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::copy(const SequenceOfVectors<V,M>& src)
{
  BaseVectorSequence<V,M>::copy(src);
  for (unsigned int i = 0; i < (unsigned int) m_seq.size(); ++i) {
    if (m_seq[i]) {
      delete m_seq[i];
      m_seq[i] = NULL;
    }
  }
  m_seq.resize(src.subSequenceSize(),NULL);
  for (unsigned int i = 0; i < m_seq.size(); ++i) {
    m_seq[i] = new V(*(src.m_seq[i]));
  }

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::extractRawData(
  unsigned int         initialPos,
  unsigned int         spacing,
  unsigned int         numPos,
  unsigned int         paramId,
  std::vector<double>& rawData) const
{
  rawData.resize(numPos);
  if (spacing == 1) {
    for (unsigned int j = 0; j < numPos; ++j) {
      rawData[j] = (*(m_seq[initialPos+j        ]))[paramId];
    }
  }
  else {
    for (unsigned int j = 0; j < numPos; ++j) {
      rawData[j] = (*(m_seq[initialPos+j*spacing]))[paramId];
    }
  }

  return;
}

// --------------------------------------------------
// Methods conditionally available ------------------
// --------------------------------------------------
// --------------------------------------------------
#ifdef UQ_SEQ_VEC_USES_OPERATOR
template <class V, class M>
const V*
SequenceOfVectors<V,M>::operator[](unsigned int posId) const
{
  queso_require_less_msg(posId, this->subSequenceSize(), "posId > subSequenceSize()");

  return (const V*) (m_seq[posId]);
}
// --------------------------------------------------
template <class V, class M>
const V*&
SequenceOfVectors<V,M>::operator[](unsigned int posId)
{
  queso_require_less_msg(posId, this->subSequenceSize(), "posId > subSequenceSize()");

  return m_seq[posId];
}
#endif

// --------------------------------------------------
// --------------------------------------------------
// --------------------------------------------------

#ifdef UQ_ALSO_COMPUTE_MDFS_WITHOUT_KDE
template <class V, class M>
void
SequenceOfVectors<V,M>::subUniformlySampledMdf(
  const V&                       numEvaluationPointsVec,
  ArrayOfOneDGrids <V,M>& mdfGrids,
  ArrayOfOneDTables<V,M>& mdfValues) const
{
  V minDomainValues(m_vectorSpace.zeroVector());
  V maxDomainValues(m_vectorSpace.zeroVector());

  ScalarSequence<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(0,                 // initialPos
                           1,                 // spacing
                           subSequenceSize(), // numPos
                           i,
                           data);

    std::vector<double> aMdf(0);
    data.subUniformlySampledMdf((unsigned int) numEvaluationPointsVec[i],
                                minDomainValues[i],
                                maxDomainValues[i],
                                aMdf);
    mdfValues.setOneDTable(i,aMdf);
  }

  mdfGrids.setUniformGrids(numEvaluationPointsVec,
                           minDomainValues,
                           maxDomainValues);

  return;
}
#endif

// --------------------------------------------------
// --------------------------------------------------
// --------------------------------------------------

#ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS
template <class V, class M>
void
SequenceOfVectors<V,M>::subMeanCltStd(
  unsigned int initialPos,
  unsigned int numPos,
  const V&     meanVec,
  V&           stdVec) const
{
  bool bRC = ((initialPos              <  this->subSequenceSize()) &&
              (0                       <  numPos                 ) &&
              ((initialPos+numPos)     <= this->subSequenceSize()) &&
              (this->vectorSizeLocal() == meanVec.sizeLocal()    ) &&
              (this->vectorSizeLocal() == stdVec.sizeLocal()     ));
  queso_require_msg(bRC, "invalid input data");

  ScalarSequence<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    stdVec[i] = data.subMeanCltStd(0,
                                   numPos,
                                   meanVec[i]);
  }

  return;
}

template <class V, class M>
void
SequenceOfVectors<V,M>::unifiedMeanCltStd(
  unsigned int initialPos,
  unsigned int numPos,
  const V&     unifiedMeanVec,
  V&           unifiedSamVec) const
{
  bool bRC = ((initialPos              <  this->subSequenceSize()   ) &&
              (0                       <  numPos                    ) &&
              ((initialPos+numPos)     <= this->subSequenceSize()   ) &&
              (this->vectorSizeLocal() == unifiedMeanVec.sizeLocal()) &&
              (this->vectorSizeLocal() == unifiedSamVec.sizeLocal() ));
  queso_require_msg(bRC, "invalid input data");

  ScalarSequence<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    unifiedSamVec[i] = data.unifiedMeanCltStd(m_vectorSpace.numOfProcsForStorage() == 1,
                                              0,
                                              numPos,
                                              unifiedMeanVec[i]);
  }

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::bmm(
  unsigned int initialPos,
  unsigned int batchLength,
  V&           bmmVec) const
{
  bool bRC = ((initialPos              <  this->subSequenceSize()            ) &&
              (batchLength             < (this->subSequenceSize()-initialPos)) &&
              (this->vectorSizeLocal() == bmmVec.sizeLocal()                 ));
  queso_require_msg(bRC, "invalid input data");

  ScalarSequence<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           this->subSequenceSize()-initialPos,
                           i,
                           data);
    bmmVec[i] = data.bmm(0,
                         batchLength);
  }

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::fftForward(
  unsigned int                        initialPos,
  unsigned int                        fftSize,
  unsigned int                        paramId,
  std::vector<std::complex<double> >& fftResult) const
{
  bool bRC = ((initialPos           <  this->subSequenceSize()) &&
              (paramId              <  this->vectorSizeLocal()) &&
              (0                    <  fftSize                ) &&
              ((initialPos+fftSize) <= this->subSequenceSize()) &&
              (fftSize              <  this->subSequenceSize()));
  queso_require_msg(bRC, "invalid input data");

  std::vector<double> rawData(fftSize,0.);
  this->extractRawData(initialPos,
                       1, // spacing
                       fftSize,
                       paramId,
                       rawData);

  m_fftObj->forward(rawData,fftSize,fftResult);

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::psd(
  unsigned int         initialPos,
  unsigned int         numBlocks,
  double               hopSizeRatio,
  unsigned int         paramId,
  std::vector<double>& psdResult) const
{
  bool bRC = ((initialPos < this->subSequenceSize()) &&
              (paramId    < this->vectorSizeLocal()));
  queso_require_msg(bRC, "invalid input data");

  ScalarSequence<double> data(m_env,0,"");

  this->extractScalarSeq(initialPos,
                         1, // spacing
                         this->subSequenceSize()-initialPos,
                         paramId,
                         data);
  data.psd(0,
           numBlocks,
           hopSizeRatio,
           psdResult);

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::psdAtZero(
  unsigned int initialPos,
  unsigned int numBlocks,
  double       hopSizeRatio,
  V&           psdVec) const
{
  bool bRC = ((initialPos         <  this->subSequenceSize()) &&
              (this->vectorSizeLocal() == psdVec.sizeLocal()));
  queso_require_msg(bRC, "invalid input data");

  ScalarSequence<double> data(m_env,0,"");
  std::vector<double> psdResult(0,0.); // size will be determined by 'data.psd()'

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           this->subSequenceSize()-initialPos,
                           i,
                           data);
    data.psd(0,
             numBlocks,
             hopSizeRatio,
             psdResult);
    psdVec[i] = psdResult[0];
    //*m_env.subDisplayFile() << "psdResult[0] = " << psdResult[0] << std::endl;
  }

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::geweke(
  unsigned int initialPos,
  double       ratioNa,
  double       ratioNb,
  V&           gewVec) const
{
  bool bRC = ((initialPos              <  this->subSequenceSize()) &&
              (this->vectorSizeLocal() == gewVec.sizeLocal()     ));
  queso_require_msg(bRC, "invalid input data");

  unsigned int numPos = this->subSequenceSize() - initialPos;
  ScalarSequence<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    gewVec[i] = data.geweke(0,
                            ratioNa,
                            ratioNb);
  }

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::meanStacc(
  unsigned int initialPos,
  V&           meanStaccVec) const
{
  bool bRC = ((initialPos              <  this->subSequenceSize() ) &&
              (this->vectorSizeLocal() == meanStaccVec.sizeLocal()));
  queso_require_msg(bRC, "invalid input data");

  unsigned int numPos = this->subSequenceSize() - initialPos;
  ScalarSequence<double> data(m_env,0,"");

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    meanStaccVec[i] = data.meanStacc(0);
  }

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::subCdfPercentageRange(
  unsigned int initialPos,
  unsigned int numPos,
  double       range, // \in [0,1]
  V&           lowerVec,
  V&           upperVec) const
{
  bool bRC = ((0                       <  numPos                 ) &&
              ((initialPos+numPos)     <= this->subSequenceSize()) &&
              (this->vectorSizeLocal() == lowerVec.sizeLocal()   ) &&
              (this->vectorSizeLocal() == upperVec.sizeLocal()   ));
  queso_require_msg(bRC, "invalid input data");

  unsigned int numParams = this->vectorSizeLocal();
  ScalarSequence<double> data(m_env,0,"");

  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    data.subCdfPercentageRange(0,
                               numPos,
                               range,
                               lowerVec[i],
                               upperVec[i]);
  }

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::unifiedCdfPercentageRange(
  unsigned int initialPos,
  unsigned int numPos,
  double       range, // \in [0,1]
  V&           lowerVec,
  V&           upperVec) const
{
  bool bRC = ((0                       <  numPos                 ) &&
              ((initialPos+numPos)     <= this->subSequenceSize()) &&
              (this->vectorSizeLocal() == lowerVec.sizeLocal()   ) &&
              (this->vectorSizeLocal() == upperVec.sizeLocal()   ));
  queso_require_msg(bRC, "invalid input data");

  unsigned int numParams = this->vectorSizeLocal();
  ScalarSequence<double> data(m_env,0,"");

  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    data.unifiedCdfPercentageRange(m_vectorSpace.numOfProcsForStorage() == 1,
                                   0,
                                   numPos,
                                   range,
                                   lowerVec[i],
                                   upperVec[i]);
  }

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::subCdfStacc(
  unsigned int     initialPos,
  std::vector<V*>& cdfStaccVecs,
  std::vector<V*>& cdfStaccVecsUp,
  std::vector<V*>& cdfStaccVecsLow,
  std::vector<V*>& sortedDataVecs) const
{
  bool bRC = (initialPos < this->subSequenceSize());
  queso_require_msg(bRC, "invalid input data");

  unsigned int numPos = this->subSequenceSize() - initialPos;
  unsigned int numEvals = numPos;
  for (unsigned int j = 0; j < numEvals; ++j) {
    cdfStaccVecs   [j] = new V(m_vectorSpace.zeroVector());
    cdfStaccVecsUp [j] = new V(m_vectorSpace.zeroVector());
    cdfStaccVecsLow[j] = new V(m_vectorSpace.zeroVector());
    sortedDataVecs [j] = new V(m_vectorSpace.zeroVector());
  }
  std::vector<double> cdfStaccs   (numEvals,0.);
  std::vector<double> cdfStaccsup (numEvals,0.);
  std::vector<double> cdfStaccslow(numEvals,0.);

  ScalarSequence<double> data      (m_env,0,"");
  ScalarSequence<double> sortedData(m_env,0,"");
  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    //std::cout << "x-data" << data<< std::endl;
    data.subSort(initialPos,sortedData);
    data.subCdfStacc(initialPos,
                     cdfStaccs,
                     cdfStaccsup,
                     cdfStaccslow,
                     sortedData);

    for (unsigned int j = 0; j < numEvals; ++j) {
      (*sortedDataVecs [j])[i] = sortedData  [j];
      (*cdfStaccVecs   [j])[i] = cdfStaccs   [j];
      (*cdfStaccVecsUp [j])[i] = cdfStaccsup [j];
      (*cdfStaccVecsLow[j])[i] = cdfStaccslow[j];
    }
  }

  return;
}
//---------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::subCdfStacc(
  unsigned int           initialPos,
  const std::vector<V*>& evalPositionsVecs,
  std::vector<V*>&       cdfStaccVecs) const
{
  bool bRC = ((initialPos               <  this->subSequenceSize() ) &&
              (0                        <  evalPositionsVecs.size()) &&
              (evalPositionsVecs.size() == cdfStaccVecs.size()     ));
  queso_require_msg(bRC, "invalid input data");

  unsigned int numPos = this->subSequenceSize() - initialPos;
  ScalarSequence<double> data(m_env,0,"");

  unsigned int numEvals = evalPositionsVecs.size();
  for (unsigned int j = 0; j < numEvals; ++j) {
    cdfStaccVecs[j] = new V(m_vectorSpace.zeroVector());
  }
  std::vector<double> evalPositions(numEvals,0.);
  std::vector<double> cdfStaccs    (numEvals,0.);

  unsigned int numParams = this->vectorSizeLocal();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);

    for (unsigned int j = 0; j < numEvals; ++j) {
      evalPositions[j] = (*evalPositionsVecs[j])[i];
    }

    data.subCdfStacc(0,
                     evalPositions,
                     cdfStaccs);

    for (unsigned int j = 0; j < numEvals; ++j) {
      (*cdfStaccVecs[j])[i] = cdfStaccs[j];
    }
  }

  return;
}
#endif // #ifdef QUESO_COMPUTES_EXTRA_POST_PROCESSING_STATISTICS

// --------------------------------------------------
// --------------------------------------------------
// --------------------------------------------------

#ifdef UQ_CODE_HAS_MONITORS
template <class V, class M>
void
SequenceOfVectors<V,M>::subMeanMonitorAlloc(unsigned int numberOfMonitorPositions)
{
  m_subMeanMonitorPosSeq = new ScalarSequence<double>(m_env,        numberOfMonitorPositions,(m_name+"_subMeanMonitorPosSeq").c_str());
  m_subMeanVecSeq        = new SequenceOfVectors<V,M>(m_vectorSpace,numberOfMonitorPositions,(m_name+"_subMeanVecSeq").c_str()       );
  m_subMeanCltStdSeq     = new SequenceOfVectors<V,M>(m_vectorSpace,numberOfMonitorPositions,(m_name+"_subMeanCltStdSeq").c_str()    );

  return;
}
// --------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::subMeanInter0MonitorAlloc(unsigned int numberOfMonitorPositions)
{
  m_subMeanInter0MonitorPosSeq = new ScalarSequence<double>(m_env,        numberOfMonitorPositions,(m_name+"_subMeanInter0MonitorPosSeq").c_str() );
  m_subMeanInter0Mean          = new SequenceOfVectors<V,M>(m_vectorSpace,numberOfMonitorPositions,(m_name+"_subMeanInter0MeanSeq").c_str()       );
  m_subMeanInter0Clt95         = new SequenceOfVectors<V,M>(m_vectorSpace,numberOfMonitorPositions,(m_name+"_subMeanInter0Clt95Seq").c_str()      );
  m_subMeanInter0Empirical90   = new SequenceOfVectors<V,M>(m_vectorSpace,numberOfMonitorPositions,(m_name+"_subMeanInter0Empirical90Seq").c_str());
  m_subMeanInter0Min           = new SequenceOfVectors<V,M>(m_vectorSpace,numberOfMonitorPositions,(m_name+"_subMeanInter0MinSeq").c_str()        );
  m_subMeanInter0Max           = new SequenceOfVectors<V,M>(m_vectorSpace,numberOfMonitorPositions,(m_name+"_subMeanInter0MaxSeq").c_str()        );

  return;
}
// --------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::unifiedMeanMonitorAlloc(unsigned int numberOfMonitorPositions)
{
  m_unifiedMeanMonitorPosSeq = new ScalarSequence<double>(m_env,        numberOfMonitorPositions,(m_name+"_unifiedMeanMonitorPosSeq").c_str());
  m_unifiedMeanVecSeq        = new SequenceOfVectors<V,M>(m_vectorSpace,numberOfMonitorPositions,(m_name+"_unifiedMeanVecSeq").c_str()       );
  m_unifiedMeanCltStdSeq     = new SequenceOfVectors<V,M>(m_vectorSpace,numberOfMonitorPositions,(m_name+"_unifiedMeanCltStdSeq").c_str()    );

  return;
}
// --------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::subMeanMonitorRun(unsigned int monitorPosition,
                                                 V&           subMeanVec,
                                                 V&           subMeanCltStd)
{
  this->subMeanExtra(0,
                     monitorPosition,
                     subMeanVec);

  this->subMeanCltStd(0,
                      monitorPosition,
                      subMeanVec,
                      subMeanCltStd);

  return;
}
// --------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::subMeanInter0MonitorRun(unsigned int monitorPosition,
                                                       V&           subMeanInter0Mean,
                                                       V&           subMeanInter0Clt95,
                                                       V&           subMeanInter0Empirical90,
                                                       V&           subMeanInter0Min,
                                                       V&           subMeanInter0Max)
{
  V subMeanVec(m_vectorSpace.zeroVector());
  this->subMeanExtra(0,
                     monitorPosition,
                     subMeanVec);

  subMeanVec.mpiAllReduce(RawValue_MPI_SUM,m_env.inter0Comm(),subMeanInter0Mean);
  subMeanInter0Mean /= ((double) m_env.inter0Comm().NumProc());

  V subMeanInter0CltVariance = subMeanVec-subMeanInter0Mean;
  subMeanInter0CltVariance *= subMeanInter0CltVariance;
  subMeanInter0CltVariance.mpiAllReduce(RawValue_MPI_SUM,m_env.inter0Comm(),subMeanInter0Clt95);
  subMeanInter0Clt95 /= ((double) (m_env.inter0Comm().NumProc()-1));
  subMeanInter0Clt95 /= ((double) (m_env.inter0Comm().NumProc()-1));
  subMeanInter0Clt95.cwSqrt();
  subMeanInter0Clt95 *= 3.;

  V subMeanInter0Quantile5(m_vectorSpace.zeroVector());
  subMeanVec.mpiAllQuantile(.05,m_env.inter0Comm(),subMeanInter0Quantile5);
  V subMeanInter0Quantile95(m_vectorSpace.zeroVector());
  subMeanVec.mpiAllQuantile(.95,m_env.inter0Comm(),subMeanInter0Quantile95);
  subMeanInter0Empirical90 = subMeanInter0Quantile95 - subMeanInter0Quantile5;

  subMeanVec.mpiAllReduce(RawValue_MPI_MIN,m_env.inter0Comm(),subMeanInter0Min);

  subMeanVec.mpiAllReduce(RawValue_MPI_MAX,m_env.inter0Comm(),subMeanInter0Max);

  return;
}
// --------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::unifiedMeanMonitorRun(unsigned int monitorPosition,
                                                     V&           unifiedMeanVec,
                                                     V&           unifiedMeanCltStd)
{
  this->unifiedMeanExtra(0,
                         monitorPosition,
                         unifiedMeanVec);

  this->unifiedMeanCltStd(0,
                          monitorPosition,
                          unifiedMeanVec,
                          unifiedMeanCltStd);
  return;
}
// --------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::subMeanMonitorStore(unsigned int i,
                                                   unsigned int monitorPosition,
                                                   const V&     subMeanVec,
                                                   const V&     subMeanCltStd)
{
  (*m_subMeanMonitorPosSeq)[i] = monitorPosition;
  m_subMeanVecSeq->setPositionValues(i,subMeanVec);
  m_subMeanCltStdSeq->setPositionValues(i,subMeanCltStd);

  return;
}
// --------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::subMeanInter0MonitorStore(unsigned int i,
                                                         unsigned int monitorPosition,
                                                         const V&     subMeanInter0Mean,
                                                         const V&     subMeanInter0Clt95,
                                                         const V&     subMeanInter0Empirical90,
                                                         const V&     subMeanInter0Min,
                                                         const V&     subMeanInter0Max)
{
  (*m_subMeanInter0MonitorPosSeq)[i] = monitorPosition;
  m_subMeanInter0Mean->setPositionValues(i,subMeanInter0Mean);
  m_subMeanInter0Clt95->setPositionValues(i,subMeanInter0Clt95);
  m_subMeanInter0Empirical90->setPositionValues(i,subMeanInter0Empirical90);
  m_subMeanInter0Min->setPositionValues(i,subMeanInter0Min);
  m_subMeanInter0Max->setPositionValues(i,subMeanInter0Max);

  return;
}
// --------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::unifiedMeanMonitorStore(unsigned int i,
                                                       unsigned int monitorPosition,
                                                       V&           unifiedMeanVec,
                                                       V&           unifiedMeanCltStd)
{
  (*m_unifiedMeanMonitorPosSeq)[i] = monitorPosition;
  m_unifiedMeanVecSeq->setPositionValues(i,unifiedMeanVec);
  m_unifiedMeanCltStdSeq->setPositionValues(i,unifiedMeanCltStd);

  return;
}
// --------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::subMeanMonitorWrite(std::ofstream& ofs)
{
  m_subMeanMonitorPosSeq->subWriteContents(0,m_subMeanMonitorPosSeq->subSequenceSize(),ofs,UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT); // Yes, always ".m"
  m_subMeanVecSeq->subWriteContents       (0,m_subMeanVecSeq->subSequenceSize(),       ofs,UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT); // Yes, always ".m"
  m_subMeanCltStdSeq->subWriteContents    (0,m_subMeanCltStdSeq->subSequenceSize(),    ofs,UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT); // Yes, always ".m"

  return;
}
// --------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::subMeanInter0MonitorWrite(std::ofstream& ofs)
{
  m_subMeanInter0MonitorPosSeq->subWriteContents(0,m_subMeanInter0MonitorPosSeq->subSequenceSize(),ofs,UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT); // Yes, always ".m"
  m_subMeanInter0Mean->subWriteContents         (0,m_subMeanInter0Mean->subSequenceSize(),         ofs,UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT); // Yes, always ".m"
  m_subMeanInter0Clt95->subWriteContents        (0,m_subMeanInter0Clt95->subSequenceSize(),        ofs,UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT); // Yes, always ".m"
  m_subMeanInter0Empirical90->subWriteContents  (0,m_subMeanInter0Empirical90->subSequenceSize(),  ofs,UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT); // Yes, always ".m"
  m_subMeanInter0Min->subWriteContents          (0,m_subMeanInter0Min->subSequenceSize(),          ofs,UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT); // Yes, always ".m"
  m_subMeanInter0Max->subWriteContents          (0,m_subMeanInter0Max->subSequenceSize(),          ofs,UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT); // Yes, always ".m"

  return;
}
// --------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::unifiedMeanMonitorWrite(std::ofstream& ofs)
{
  // std::set<unsigned int> tmpSet;
  // tmpSet.insert(0);
  m_unifiedMeanMonitorPosSeq->subWriteContents(0,m_unifiedMeanMonitorPosSeq->subSequenceSize(),ofs,UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT); // Yes, always ".m" // Yes, 'subWriteContents()'
  m_unifiedMeanVecSeq->subWriteContents       (0,m_unifiedMeanVecSeq->subSequenceSize(),       ofs,UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT); // Yes, always ".m" // Yes, 'subWriteContents()'
  m_unifiedMeanCltStdSeq->subWriteContents    (0,m_unifiedMeanCltStdSeq->subSequenceSize(),    ofs,UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT); // Yes, always ".m" // Yes, 'subWriteContents()'

  return;
}
// --------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::subMeanMonitorFree()
{
  delete m_subMeanMonitorPosSeq;
  m_subMeanMonitorPosSeq = NULL;
  delete m_subMeanVecSeq;
  m_subMeanVecSeq = NULL;
  delete m_subMeanCltStdSeq;
  m_subMeanCltStdSeq = NULL;

  return;
}
// --------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::subMeanInter0MonitorFree()
{
  delete m_subMeanInter0MonitorPosSeq;
  m_subMeanInter0MonitorPosSeq = NULL;
  delete m_subMeanInter0Mean;
  m_subMeanInter0Mean = NULL;
  delete m_subMeanInter0Clt95;
  m_subMeanInter0Clt95 = NULL;
  delete m_subMeanInter0Empirical90;
  m_subMeanInter0Empirical90 = NULL;
  delete m_subMeanInter0Min;
  m_subMeanInter0Min = NULL;
  delete m_subMeanInter0Max;
  m_subMeanInter0Max = NULL;

  return;
}
// --------------------------------------------------
template <class V, class M>
void
SequenceOfVectors<V,M>::unifiedMeanMonitorFree()
{
  delete m_unifiedMeanMonitorPosSeq;
  m_unifiedMeanMonitorPosSeq = NULL;
  delete m_unifiedMeanVecSeq;
  m_unifiedMeanVecSeq = NULL;
  delete m_unifiedMeanCltStdSeq;
  m_unifiedMeanCltStdSeq = NULL;

  return;
}
#endif // #ifdef UQ_CODE_HAS_MONITORS

}  // End namespace QUESO

template class QUESO::SequenceOfVectors<QUESO::GslVector, QUESO::GslMatrix>;
