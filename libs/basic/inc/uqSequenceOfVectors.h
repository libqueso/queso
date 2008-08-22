/* uq/libs/basic/inc/uqSequenceOfVectors.h
 *
 * Copyright (C) 2008 The PECOS Team, http://queso.ices.utexas.edu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __UQ_SEQUENCE_OF_VECTORS_H__
#define __UQ_SEQUENCE_OF_VECTORS_H__

#include <EpetraExt_DistArray.h>

typedef std::vector<double> ScalarSequenceType;

template <class V>
class uqSequenceOfVectorsClass
{
public:
  uqSequenceOfVectorsClass(unsigned int sequenceSize, const V& zeroVector);
 ~uqSequenceOfVectorsClass();

  const unsigned int sequenceSize      () const;
  const unsigned int vectorSize        () const;
        void         resizeSequence    (unsigned int newSequenceSize            );
        void         setPositionValues (unsigned int positionId, const V& vector);
        void         mean              (unsigned int initialPos,
                                        unsigned int numPos,
                                        V&           mean) const;
        void         sampleVariance    (unsigned int initialPos,
                                        unsigned int numPos,
                                        const V&     mean,
                                        V&           sampleVariance) const;
        void         populationVariance(unsigned int initialPos,
                                        unsigned int numPos,
                                        const V&     mean,
                                        V&           populVariance) const;
#if 0
  const V            operator[]                                 const;
        V            operator[]                                 const;
#endif
private:
  const uqEnvironmentClass&                 m_env;
  V*                                        m_zeroVector;
  EpetraExt::DistArray<ScalarSequenceType*> m_scalarSequences;
};

template <class V>
uqSequenceOfVectorsClass<V>::uqSequenceOfVectorsClass(
  unsigned int sequenceSize,
  const V&     zeroVector)
  :
  m_env            (zeroVector.env()  ),
  m_zeroVector     (new V(zeroVector) ),
  m_scalarSequences(zeroVector.map(),1)
{

  //if (m_env.rank() == 0) std::cout << "Entering uqSequenceOfVectorsClass<V>::constructor()"
  //                                 << std::endl;

  //if (m_env.rank() == 0) std::cout << "In uqSequenceOfVectorsClass<V>::constructor()"
  //                                 << "\n sequenceSize = "                 << sequenceSize
  //                                 << "\n m_scalarSequences.MyLength() = " << m_scalarSequences.MyLength()
  //                                 << std::endl;

  for (unsigned int i = 0; i < (unsigned int) m_scalarSequences.MyLength(); ++i) {
    m_scalarSequences(i,0) = new ScalarSequenceType(sequenceSize,0.);
  }

  //if (m_env.rank() == 0) std::cout << "Leaving uqSequenceOfVectorsClass<V>::constructor()"
  //                                 << std::endl;
}

template <class V>
uqSequenceOfVectorsClass<V>::~uqSequenceOfVectorsClass()
{
  for (unsigned int i = 0; i < (unsigned int) m_scalarSequences.MyLength(); ++i) {
    if (m_scalarSequences(i,0)) delete m_scalarSequences(i,0);
  }
  if (m_zeroVector) delete m_zeroVector;
}

template <class V>
const unsigned int
uqSequenceOfVectorsClass<V>::sequenceSize() const
{
  uqSequenceOfVectorsClass<V>* tmp = const_cast<uqSequenceOfVectorsClass<V>*>(this);

  return tmp->m_scalarSequences(0,0)->size();
}

template <class V>
const unsigned int
uqSequenceOfVectorsClass<V>::vectorSize() const
{
  return m_zeroVector->size();
}

template <class V>
void
uqSequenceOfVectorsClass<V>::resizeSequence(unsigned int newSequenceSize)
{
  for (unsigned int i = 0; i < (unsigned int) m_scalarSequences.MyLength(); ++i) {
    m_scalarSequences(i,0)->resize(newSequenceSize);
    std::vector<double>(*(m_scalarSequences(i,0))).swap(*(m_scalarSequences(i,0)));
  }

  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::setPositionValues(unsigned int positionId, const V& vector)
{
  for (unsigned int i = 0; i < (unsigned int) m_scalarSequences.MyLength(); ++i) {
    (*(m_scalarSequences(i,0)))[positionId] = vector[i];
  }

  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::mean(
  unsigned int initialPos,
  unsigned int numPos,
  V&           mean) const
{
  bool bRC = ((0                     <= initialPos             ) &&
              (0                     <  numPos                 ) &&
              ((initialPos+numPos-1) <= (this->sequenceSize()-1)));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqVectorSequenceMean<V>()",
                      "invalid initial position or number of positions");

  bRC = (this->vectorSize() == mean.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqVectorSequenceMean<V>()",
                      "incompatible sizes between mean vector and vectors in sequence");

  unsigned int loopSize      = numPos;
  unsigned int finalPosPlus1 = initialPos + loopSize;
  double doubleLoopSize = (double) loopSize;
  mean.cwSet(0.);

  uqSequenceOfVectorsClass<V>* tmp = const_cast<uqSequenceOfVectorsClass<V>*>(this);
  for (unsigned int i = 0; i < mean.size(); ++i) {
    const ScalarSequenceType& seq = *(tmp->m_scalarSequences(i,0));
    double                    result = 0.;
    for (unsigned int j = initialPos; j < finalPosPlus1; ++j) {
      result += seq[j];
    }
    mean[i] = result/doubleLoopSize;
  }

  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::sampleVariance(
  unsigned int initialPos,
  unsigned int numPos,
  const V&     mean,
  V&           sampleVariance) const
{
  bool bRC = ((0                     <= initialPos             ) &&
              (0                     <  numPos                 ) &&
              ((initialPos+numPos-1) <= (this->sequenceSize()-1)));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqVectorSequenceSampleVariance<V>()",
                      "invalid initial position or number of positions");

  bRC = (this->vectorSize() == sampleVariance.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqVectorSequenceSampleVariance<V>()",
                      "incompatible sizes between sampleVariance vector and vectors in sequence");

  bRC = (this->vectorSize() == mean.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqVectorSequenceSampleVariance<V>()",
                      "incompatible sizes between mean vector and vectors in sequence");

  unsigned int loopSize      = numPos;
  unsigned int finalPosPlus1 = initialPos + loopSize;
  double doubleLoopSize = (double) loopSize;
  sampleVariance.cwSet(0.);

  uqSequenceOfVectorsClass<V>* tmp = const_cast<uqSequenceOfVectorsClass<V>*>(this);
  for (unsigned int i = 0; i < sampleVariance.size(); ++i) {
    const ScalarSequenceType& seq = *(tmp->m_scalarSequences(i,0));
    double                    tmpMean = mean[i];
    double                    result = 0.;
    for (unsigned int j = initialPos; j < finalPosPlus1; ++j) {
      double diff = seq[j] - tmpMean;
      result += diff*diff;
    }
    sampleVariance[i] = result/(doubleLoopSize - 1.);
  }

  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::populationVariance(
  unsigned int initialPos,
  unsigned int numPos,
  const V&     mean,
  V&           populVariance) const
{
  bool bRC = ((0                     <= initialPos             ) &&
              (0                     <  numPos                 ) &&
              ((initialPos+numPos-1) <= (this->sequenceSize()-1)));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqVectorSequencePopulationVariance<V>()",
                      "invalid initial position or number of positions");

  bRC = (this->vectorSize() == populVariance.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqVectorSequencePopulationVariance<V>()",
                      "incompatible sizes between populVariance vector and vectors in sequence");

  bRC = (this->vectorSize() == mean.size());
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqVectorSequencePopulationVariance<V>()",
                      "incompatible sizes between mean vector and vectors in sequence");

  unsigned int loopSize      = numPos;
  unsigned int finalPosPlus1 = initialPos + loopSize;
  double doubleLoopSize = (double) loopSize;
  populVariance.cwSet(0.);

  uqSequenceOfVectorsClass<V>* tmp = const_cast<uqSequenceOfVectorsClass<V>*>(this);
  for (unsigned int i = 0; i < populVariance.size(); ++i) {
    const ScalarSequenceType& seq = *(tmp->m_scalarSequences(i,0));
    double                    tmpMean = mean[i];
    double                    result = 0.;
    for (unsigned int j = initialPos; j < finalPosPlus1; ++j) {
      double diff = seq[j] - tmpMean;
      result += diff*diff;
    }
    populVariance[i] = result/doubleLoopSize;
  }

  return;
}
#endif // __UQ_SEQUENCE_OF_VECTORS_H__

