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

#include <uqChain.h>
#define UQ_SEQ_VEC_USES_SCALAR_SEQ_CODE
#undef UQ_SEQ_VEC_USES_OPERATOR

template <class V>
class uqSequenceOfVectorsClass : public uqChainBaseClass<V>
{
public:
  typedef typename std::vector<const V*>::const_iterator seqVectorPositionConstIteratorTypedef;
  typedef typename std::vector<const V*>::iterator       seqVectorPositionIteratorTypedef;
  uqSequenceOfVectorsClass(unsigned int sequenceSize, const V& vectorExample);
 ~uqSequenceOfVectorsClass();

        unsigned int sequenceSize      () const;
        void         resizeSequence    (unsigned int newSequenceSize);
        void         resetValues       (unsigned int initialPos, unsigned int numPos);
        void         erasePositions    (unsigned int initialPos, unsigned int numPos);
#ifdef UQ_SEQ_VEC_USES_OPERATOR
	const V*     operator[]        (unsigned int posId) const;
	const V*&    operator[]        (unsigned int posId);
#endif
        void         getPositionValues (unsigned int posId,       V& vec) const;
        void         setPositionValues (unsigned int posId, const V& vec);
        void         setGaussian       (const gsl_rng* rng, const V& meanVec, const V& stdDevVec);
        void         setUniform        (const gsl_rng* rng, const V& aVec,    const V& bVec     );

        void         mean              (unsigned int              initialPos,
                                        unsigned int              numPos,
                                        V&                        meanVec) const;
        void         sampleVariance    (unsigned int              initialPos,
                                        unsigned int              numPos,
                                        const V&                  meanVec,
                                        V&                        samVec) const;
        void         populationVariance(unsigned int              initialPos,
                                        unsigned int              numPos,
                                        const V&                  meanVec,
                                        V&                        popVec) const;
        void         autoCovariance    (unsigned int              initialPos,
                                        unsigned int              numPos,
                                        const V&                  meanVec,
                                        unsigned int              lag,
                                        V&                        covVec) const;

        void         autoCorrViaDef    (unsigned int              initialPos,
                                        unsigned int              numPos,
                                        unsigned int              lag,
                                        V&                        corrVec) const;
        void         autoCorrViaFft    (unsigned int                     initialPos,
                                        unsigned int                     numPos,
                                        const std::vector<unsigned int>& lags,
                                        std::vector<V*>&                 corrVecs) const;
        void         autoCorrViaFft    (unsigned int              initialPos,
                                        unsigned int              numPos,
                                        unsigned int              numSum,
                                        V&                        autoCorrsSumVec) const;
        void         bmm               (unsigned int              initialPos,
                                        unsigned int              batchLength,
                                        V&                        bmmVec) const;
        void         fftForward        (unsigned int                        initialPos,
                                        unsigned int                        fftSize,
                                        unsigned int                        paramId,
                                        std::vector<std::complex<double> >& fftResult) const;
        //void         fftInverse        (unsigned int fftSize);
        void         psd               (unsigned int              initialPos,
                                        unsigned int              numBlocks,
                                        double                    hopSizeRatio,
                                        unsigned int              paramId,
                                        std::vector<double>&      psdResult) const;
        void         psdAtZero         (unsigned int              initialPos,
                                        unsigned int              numBlocks,
                                        double                    hopSizeRatio,
                                        V&                        psdVec) const;
        void         geweke            (unsigned int              initialPos,
                                        double                    ratioNa,
                                        double                    ratioNb,
                                        V&                        gewVec) const;
        void         minMax            (unsigned int              initialPos,
                                        V&                        minVec,
                                        V&                        maxVec) const;
        void         histogram         (unsigned int              initialPos,
                                        const V&                  minVec,
                                        const V&                  maxVec,
                                        std::vector<V*>&          centersForAllBins,
                                        std::vector<V*>&          binsForAllParams) const;
        void         interQuantileRange(unsigned int              initialPos,
                                        V&                        iqrVec) const;
        void         scalesForKDE      (unsigned int              initialPos,
                                        const V&                  iqrVec,
                                        V&                        scaleVec) const;
        void         gaussianKDE       (unsigned int              initialPos,
                                        const V&                  scaleVec,
                                        const std::vector<V*>&    evaluationPositions,
                                        std::vector<V*>&          densityValues) const;
        void         write             (const std::string&        name,
                                        std::ofstream&            ofs) const;
        void         select            (const std::vector<unsigned int>& idsOfUniquePositions);
        void         filter            (unsigned int              initialPos,
                                        unsigned int              spacing);

private:
        void         extractScalarSeq  (unsigned int                   initialPos,
                                        unsigned int                   spacing,
                                        unsigned int                   numPos,
                                        unsigned int                   paramId,
                                        uqScalarSequenceClass<double>& scalarSeq) const;
        void         extractRawData    (unsigned int                   initialPos,
                                        unsigned int                   spacing,
                                        unsigned int                   numPos,
                                        unsigned int                   paramId,
                                        std::vector<double>&           rawData) const;

  std::vector<const V*> m_seq;

  using uqChainBaseClass<V>::m_env;
  using uqChainBaseClass<V>::m_vectorExample;
  using uqChainBaseClass<V>::m_fftObj;
};

template <class V>
uqSequenceOfVectorsClass<V>::uqSequenceOfVectorsClass(
  unsigned int sequenceSize,
  const V&     vectorExample)
  :
  uqChainBaseClass<V>(sequenceSize,vectorExample),
  m_seq              (sequenceSize,NULL)
{

  //if (m_env.rank() == 0) std::cout << "Entering uqSequenceOfVectorsClass<V>::constructor()"
  //                                 << std::endl;

  //if (m_env.rank() == 0) std::cout << "Leaving uqSequenceOfVectorsClass<V>::constructor()"
  //                                 << std::endl;
}

template <class V>
uqSequenceOfVectorsClass<V>::~uqSequenceOfVectorsClass()
{
  for (unsigned int i = 0; i < (unsigned int) m_seq.size(); ++i) {
    if (m_seq[i]) delete m_seq[i];
  }
}

template <class V>
unsigned int
uqSequenceOfVectorsClass<V>::sequenceSize() const
{
  return m_seq.size();
}

template <class V>
void
uqSequenceOfVectorsClass<V>::resizeSequence(unsigned int newSequenceSize)
{
  if (newSequenceSize != this->sequenceSize()) {
    m_seq.resize(newSequenceSize,NULL);
    std::vector<const V*>(m_seq).swap(m_seq);
  }

 return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::resetValues(unsigned int initialPos, unsigned int numPos)
{
  bool bRC = ((initialPos          <  this->sequenceSize()) &&
              (0                   <  numPos              ) &&
              ((initialPos+numPos) <= this->sequenceSize()));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V>::resetValues()",
                      "invalid input data");

  for (unsigned int j = 0; j < numPos; ++j) {
    if (m_seq[initialPos+j] != NULL) {
      delete m_seq[initialPos+j];
      m_seq[initialPos+j] = NULL;
    }
  }

  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::erasePositions(unsigned int initialPos, unsigned int numPos)
{
  bool bRC = ((initialPos          <  this->sequenceSize()) &&
              (0                   <  numPos              ) &&
              ((initialPos+numPos) <= this->sequenceSize()));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V>::erasePositions()",
                      "invalid input data");

  for (unsigned int j = 0; j < numPos; ++j) {
    if (m_seq[initialPos+j] != NULL) delete m_seq[initialPos+j];
  }

  seqVectorPositionIteratorTypedef posIteratorBegin = m_seq.begin();
  if (initialPos < this->sequenceSize()) std::advance(posIteratorBegin,initialPos);
  else                                   posIteratorBegin = m_seq.end();

  unsigned int posEnd = initialPos + numPos - 1;
  seqVectorPositionIteratorTypedef posIteratorEnd = m_seq.begin();
  if (posEnd < this->sequenceSize()) std::advance(posIteratorEnd,posEnd);
  else                               posIteratorEnd = m_seq.end();

  unsigned int oldSequenceSize = this->sequenceSize();
  m_seq.erase(posIteratorBegin,posIteratorEnd);
  UQ_FATAL_TEST_MACRO((oldSequenceSize - numPos) != this->sequenceSize(),
                      m_env.rank(),
                      "uqSequenceOfVectors::erasePositions()",
                      "(oldSequenceSize - numPos) != this->sequenceSize()");

  return;
}

#ifdef UQ_SEQ_VEC_USES_OPERATOR
template <class V>
const V*
uqSequenceOfVectorsClass<V>::operator[](unsigned int posId) const
{
  UQ_FATAL_TEST_MACRO((posId >= this->sequenceSize()),
                      m_env.rank(),
                      "uqScalarSequences<V>::operator[] const",
                      "posId > sequenceSize()");

  return (const V*) (m_seq[posId]);
}

template <class V>
const V*&
uqSequenceOfVectorsClass<V>::operator[](unsigned int posId)
{
  UQ_FATAL_TEST_MACRO((posId >= this->sequenceSize()),
                      m_env.rank(),
                      "uqScalarSequences<V>::operator[] const",
                      "posId > sequenceSize()");

  return m_seq[posId];
}
#endif

template <class V>
void
uqSequenceOfVectorsClass<V>::getPositionValues(unsigned int posId, V& vec) const
{
  UQ_FATAL_TEST_MACRO((posId >= this->sequenceSize()),
                      m_env.rank(),
                      "uqScalarSequences<V>::getPositionValues()",
                      "posId > sequenceSize()");

  vec = *(const_cast<V*>(m_seq[posId]));

  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::setPositionValues(unsigned int posId, const V& vec)
{
  UQ_FATAL_TEST_MACRO((posId >= this->sequenceSize()),
                      m_env.rank(),
                      "uqScalarSequences<V>::setPositionValues()",
                      "posId > sequenceSize()");

  if (m_seq[posId] != NULL) delete m_seq[posId];
  m_seq[posId] = new V(vec);

  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::setGaussian(const gsl_rng* rng, const V& meanVec, const V& stdDevVec)
{
  V gaussianVector(m_vectorExample);
  for (unsigned int j = 0; j < this->sequenceSize(); ++j) {
    gaussianVector.cwSetGaussian(m_env.rng(),meanVec,stdDevVec);
    if (m_seq[j] != NULL) delete m_seq[j];
    m_seq[j] = new V(gaussianVector);
  }

  return;
}


template <class V>
void
uqSequenceOfVectorsClass<V>::setUniform(const gsl_rng* rng, const V& aVec, const V& bVec)
{
  V uniformVector(m_vectorExample);
  for (unsigned int j = 0; j < this->sequenceSize(); ++j) {
    uniformVector.cwSetUniform(m_env.rng(),aVec,bVec);
    if (m_seq[j] != NULL) delete m_seq[j];
    m_seq[j] = new V(uniformVector);
  }

  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::mean(
  unsigned int initialPos,
  unsigned int numPos,
  V&           meanVec) const
{
  bool bRC = ((initialPos          <  this->sequenceSize()) &&
              (0                   <  numPos              ) &&
              ((initialPos+numPos) <= this->sequenceSize()) &&
              (this->vectorSize()  == meanVec.size()      ));
  if ((bRC == false) && (m_env.rank() == 0)) {
    std::cout << "In uqSequenceOfVectorsClass<V>::mean()"
              << ", initialPos = "           << initialPos
              << ", this->sequenceSize() = " << this->sequenceSize()
              << ", numPos = "               << numPos
              << ", this->vectorSize() = "   << this->vectorSize()
              << ", meanVec.size() = "       << meanVec.size()
              << std::endl;
  }
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V>::mean()",
                      "invalid input data");

  uqScalarSequenceClass<double> data(m_env,0);

  unsigned int numParams = this->vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    meanVec[i] = data.mean(0,
                           numPos);
  }

  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::sampleVariance(
  unsigned int initialPos,
  unsigned int numPos,
  const V&     meanVec,
  V&           samVec) const
{
  bool bRC = ((initialPos          <  this->sequenceSize()) &&
              (0                   <  numPos              ) &&
              ((initialPos+numPos) <= this->sequenceSize()) &&
              (this->vectorSize()  == meanVec.size()      ) &&
              (this->vectorSize()  == samVec.size()       ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V>::sampleVariance()",
                      "invalid input data");

  uqScalarSequenceClass<double> data(m_env,0);

  unsigned int numParams = this->vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    samVec[i] = data.sampleVariance(0,
                                    numPos,
                                    meanVec[i]);
  }

  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::populationVariance(
  unsigned int initialPos,
  unsigned int numPos,
  const V&     meanVec,
  V&           popVec) const
{
  bool bRC = ((initialPos          <  this->sequenceSize()) &&
              (0                   <  numPos              ) &&
              ((initialPos+numPos) <= this->sequenceSize()) &&
              (this->vectorSize()  == meanVec.size()      ) &&
              (this->vectorSize()  == popVec.size()       ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V>::populationVariance()",
                      "invalid input data");

  uqScalarSequenceClass<double> data(m_env,0);

  unsigned int numParams = this->vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    popVec[i] = data.populationVariance(0,
                                        numPos,
                                        meanVec[i]);
  }

  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::autoCovariance(
  unsigned int initialPos,
  unsigned int numPos,
  const V&     meanVec,
  unsigned int lag,
  V&           covVec) const
{
  bool bRC = ((initialPos          <  this->sequenceSize()) &&
              (0                   <  numPos              ) &&
              ((initialPos+numPos) <= this->sequenceSize()) &&
              (this->vectorSize()  == meanVec.size()      ) &&
              (lag                 <  numPos              ) && // lag should not be too large
              (this->vectorSize()  == covVec.size()       ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V>::autoCovariance()",
                      "invalid input data");

  uqScalarSequenceClass<double> data(m_env,0);

  unsigned int numParams = this->vectorSize();
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

template <class V>
void
uqSequenceOfVectorsClass<V>::autoCorrViaDef(
  unsigned int initialPos,
  unsigned int numPos,
  unsigned int lag,
  V&           corrVec) const
{
  bool bRC = ((initialPos          <  this->sequenceSize()) &&
              (0                   <  numPos              ) &&
              ((initialPos+numPos) <= this->sequenceSize()) &&
              (lag                 <  numPos              ) && // lag should not be too large
              (this->vectorSize()  == corrVec.size()      ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V>::autoCorrViaDef()",
                      "invalid input data");

  uqScalarSequenceClass<double> data(m_env,0);

  unsigned int numParams = this->vectorSize();
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

template <class V>
void
uqSequenceOfVectorsClass<V>::autoCorrViaFft(
  unsigned int                     initialPos,
  unsigned int                     numPos,
  const std::vector<unsigned int>& lags,
  std::vector<V*>&                 corrVecs) const
{
  bool bRC = ((initialPos          <  this->sequenceSize()) &&
              (0                   <  numPos              ) &&
              ((initialPos+numPos) <= this->sequenceSize()) &&
              (0                   < lags.size()          ) &&
              (lags[lags.size()-1] <  numPos              )); // lag should not be too large
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V>::autoCorrViaFft()",
                      "invalid input data");

  for (unsigned int j = lags.size(); j < corrVecs.size(); ++j) {
    if (corrVecs[j] != NULL) delete corrVecs[j];
  }
  corrVecs.resize(lags.size(),NULL);
  for (unsigned int j = 0;           j < corrVecs.size(); ++j) {
    if (corrVecs[j] == NULL) corrVecs[j] = new V(m_vectorExample);
  }

  uqScalarSequenceClass<double> data(m_env,0);
  unsigned int maxLag = lags[lags.size()-1];
  std::vector<double> autoCorrs(maxLag+1,0.); // Yes, +1

  unsigned int numParams = this->vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);

    //if (m_env.rank() == 0) {
    //  std::cout << "In uqSequenceOfVectorsClass<V>::autoCorrViaFft()"
    //            << ": about to call data.autoCorrViaFft() for paramId = " << i
    //            << ", with numPos = "      << numPos
    //            << ", maxLag = "           << maxLag
    //            << ", autoCorrs.size() = " << autoCorrs.size()
    //            << std::endl;
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

template <class V>
void
uqSequenceOfVectorsClass<V>::autoCorrViaFft(
  unsigned int initialPos,
  unsigned int numPos,
  unsigned int numSum,
  V&           autoCorrsSumVec) const
{
  bool bRC = ((initialPos             <  this->sequenceSize()) &&
              (0                      <  numPos              ) &&
              ((initialPos+numPos)    <= this->sequenceSize()) &&
              (0                      <  numSum              ) &&
              (numSum                 <= numPos              ) &&
              (autoCorrsSumVec.size() == this->vectorSize()  ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V>::autoCorrViaFft(), for sum",
                      "invalid input data");

  uqScalarSequenceClass<double> data(m_env,0);

  unsigned int numParams = this->vectorSize();
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

template <class V>
void
uqSequenceOfVectorsClass<V>::bmm(
  unsigned int initialPos,
  unsigned int batchLength,
  V&           bmmVec) const
{
  bool bRC = ((initialPos          <  this->sequenceSize()            ) &&
              (batchLength         < (this->sequenceSize()-initialPos)) &&
              (this->vectorSize()  == bmmVec.size()                   ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V>::bmm()",
                      "invalid input data");

  uqScalarSequenceClass<double> data(m_env,0);

  unsigned int numParams = this->vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           this->sequenceSize()-initialPos,
                           i,
                           data);
    bmmVec[i] = data.bmm(0,
                         batchLength);
  }

  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::fftForward(
  unsigned int                        initialPos,
  unsigned int                        fftSize,
  unsigned int                        paramId,
  std::vector<std::complex<double> >& fftResult) const
{
  bool bRC = ((initialPos           <  this->sequenceSize()) &&
              (paramId              <  this->vectorSize()  ) &&
              (0                    <  fftSize             ) &&
              ((initialPos+fftSize) <= this->sequenceSize()) &&
              (fftSize              <  this->sequenceSize()));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V>::fftForward()",
                      "invalid input data");

  std::vector<double> rawData(fftSize,0.);
  this->extractRawData(initialPos,
                       1, // spacing
                       fftSize,
                       paramId,
                       rawData);

  m_fftObj->forward(rawData,fftSize,fftResult);

  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::psd(
  unsigned int         initialPos,
  unsigned int         numBlocks,
  double               hopSizeRatio,
  unsigned int         paramId,
  std::vector<double>& psdResult) const
{
  bool bRC = ((initialPos < this->sequenceSize()) &&
              (paramId    < this->vectorSize()  ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V>::psd()",
                      "invalid input data");

  uqScalarSequenceClass<double> data(m_env,0);

  this->extractScalarSeq(initialPos,
                         1, // spacing
                         this->sequenceSize()-initialPos,
                         paramId,
                         data);
  data.psd(0,
           numBlocks,
           hopSizeRatio,
           psdResult);

  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::psdAtZero(
  unsigned int initialPos,
  unsigned int numBlocks,
  double       hopSizeRatio,
  V&           psdVec) const
{
  bool bRC = ((initialPos         <  this->sequenceSize()) &&
              (this->vectorSize() == psdVec.size()       ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V>::psdAtZero()",
                      "invalid input data");

  uqScalarSequenceClass<double> data(m_env,0);
  std::vector<double> psdResult(0,0.); // size will be determined by 'data.psd()'

  unsigned int numParams = this->vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           this->sequenceSize()-initialPos,
                           i,
                           data);
    data.psd(0,
             numBlocks,
             hopSizeRatio,
             psdResult);
    psdVec[i] = psdResult[0];
    //std::cout << "psdResult[0] = " << psdResult[0] << std::endl;
  }

  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::geweke(
  unsigned int initialPos,
  double       ratioNa,
  double       ratioNb,
  V&           gewVec) const
{
  bool bRC = ((initialPos         <  this->sequenceSize()) &&
              (this->vectorSize() == gewVec.size()       ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V>::geweke()",
                      "invalid input data");

  unsigned int numPos = this->sequenceSize() - initialPos;
  uqScalarSequenceClass<double> data(m_env,0);

  unsigned int numParams = this->vectorSize();
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

template <class V>
void
uqSequenceOfVectorsClass<V>::minMax(
  unsigned int initialPos,
  V&           minVec,
  V&           maxVec) const
{
  bool bRC = ((initialPos         <  this->sequenceSize()) &&
              (this->vectorSize() == minVec.size()       ) &&
              (this->vectorSize() == maxVec.size()       ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V>::minMax()",
                      "invalid input data");

  unsigned int numPos = this->sequenceSize() - initialPos;
  unsigned int numParams = this->vectorSize();
  uqScalarSequenceClass<double> data(m_env,0);

  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    data.minMax(0,minVec[i],maxVec[i]);
  }

  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::histogram(
  unsigned int     initialPos,
  const V&         minVec,
  const V&         maxVec,
  std::vector<V*>& centersForAllBins,
  std::vector<V*>& binsForAllParams) const
{
  bool bRC = ((initialPos               <  this->sequenceSize()    ) &&
              (this->vectorSize()       == minVec.size()           ) &&
              (this->vectorSize()       == maxVec.size()           ) &&
              (0                        <  centersForAllBins.size()) &&
              (centersForAllBins.size() == binsForAllParams.size() ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V>::histogram()",
                      "invalid input data");

  for (unsigned int j = 0; j < binsForAllParams.size(); ++j) {
    centersForAllBins[j] = new V(m_vectorExample);
    binsForAllParams [j] = new V(m_vectorExample);
  }

  unsigned int dataSize = this->sequenceSize() - initialPos;
  unsigned int numParams = this->vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    uqScalarSequenceClass<double> data(m_env,dataSize);
    for (unsigned int j = 0; j < dataSize; ++j) {
      data[j] = (*(m_seq[initialPos+j]))[i];
    }

    std::vector<double      > centers(centersForAllBins.size(),0.);
    std::vector<unsigned int> bins   (binsForAllParams.size(), 0 );
    data.histogram(0,
                   minVec[i],
                   maxVec[i],
                   centers,
                   bins);

    for (unsigned int j = 0; j < bins.size(); ++j) {
      (*(centersForAllBins[j]))[i] = centers[j];
      (*(binsForAllParams [j]))[i] = (double) bins[j];
    }
  }

  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::interQuantileRange(
  unsigned int initialPos,
  V&           iqrVec) const
{
  bool bRC = ((initialPos         <  this->sequenceSize()) &&
              (this->vectorSize() == iqrVec.size()       ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V>::interQuantileRange()",
                      "invalid input data");

  unsigned int numPos = this->sequenceSize() - initialPos;
  uqScalarSequenceClass<double> data(m_env,0);

  unsigned int numParams = this->vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    iqrVec[i] = data.interQuantileRange(0);
  }

  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::scalesForKDE(
  unsigned int initialPos,
  const V&     iqrVec,
  V&           scaleVec) const
{
  bool bRC = ((initialPos         <  this->sequenceSize()) &&
              (this->vectorSize() == iqrVec.size()       ) &&
              (this->vectorSize() == scaleVec.size()     ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V>::scalesForKDE()",
                      "invalid input data");

  unsigned int numPos = this->sequenceSize() - initialPos;
  uqScalarSequenceClass<double> data(m_env,0);

  unsigned int numParams = this->vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);
    scaleVec[i] = data.scaleForKDE(0,
                                   iqrVec[i]);
  }

  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::gaussianKDE(
  unsigned int           initialPos,
  const V&               scaleVec,
  const std::vector<V*>& evaluationPositions,
  std::vector<V*>&       densityValues) const
{
  bool bRC = ((initialPos                 <  this->sequenceSize()      ) &&
              (this->vectorSize()         == scaleVec.size()           ) &&
              (0                          <  evaluationPositions.size()) &&
              (evaluationPositions.size() == densityValues.size()      ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      m_env.rank(),
                      "uqSequenceOfVectorsClass<V>::gaussianKDE()",
                      "invalid input data");

  unsigned int numPos = this->sequenceSize() - initialPos;
  uqScalarSequenceClass<double> data(m_env,0);

  unsigned int numEvals = evaluationPositions.size();
  for (unsigned int j = 0; j < numEvals; ++j) {
    densityValues[j] = new V(m_vectorExample);
  }
  std::vector<double> evalScalarPositions(numEvals,0.);
  std::vector<double> evaluatedDensities (numEvals,0.);

  unsigned int numParams = this->vectorSize();
  for (unsigned int i = 0; i < numParams; ++i) {
    this->extractScalarSeq(initialPos,
                           1, // spacing
                           numPos,
                           i,
                           data);

    for (unsigned int j = 0; j < numEvals; ++j) {
      evalScalarPositions[j] = (*evaluationPositions[j])[i];
    }

    data.gaussianKDE(0,
                     scaleVec[i],
                     evalScalarPositions,
                     evaluatedDensities);

    for (unsigned int j = 0; j < numEvals; ++j) {
      (*densityValues[j])[i] = evaluatedDensities[j];
    }
  }

  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::select(const std::vector<unsigned int>& idsOfUniquePositions)
{
  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::filter(
  unsigned int initialPos,
  unsigned int spacing)
{
  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::write(
  const std::string& chainName,
  std::ofstream&     ofs) const
{
  // Write chain
  ofs << chainName << " = zeros(" << this->sequenceSize()
      << ","                      << this->vectorSize()
      << ");"
      << std::endl;
  ofs << chainName << " = [";
  unsigned int chainSize = this->sequenceSize();
  for (unsigned int j = 0; j < chainSize; ++j) {
    ofs << *(m_seq[j])
        << std::endl;
  }
  ofs << "];\n";

  return;
}

template <class V>
void
uqSequenceOfVectorsClass<V>::extractScalarSeq(
  unsigned int                   initialPos,
  unsigned int                   spacing,
  unsigned int                   numPos,
  unsigned int                   paramId,
  uqScalarSequenceClass<double>& scalarSeq) const
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

template <class V>
void
uqSequenceOfVectorsClass<V>::extractRawData(
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
#endif // __UQ_SEQUENCE_OF_VECTORS_H__

