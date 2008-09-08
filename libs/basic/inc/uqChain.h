/* uq/libs/basic/inc/uqChain.h
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

#ifndef __UQ_CHAIN_H__
#define __UQ_CHAIN_H__

#include <uqScalarSequence.h>

template <class V>
class uqChainBaseClass
{
public:
           uqChainBaseClass(unsigned int sequenceSize, const V& vectorExample);
  virtual ~uqChainBaseClass();

  virtual const unsigned int sequenceSize  () const = 0;
          const unsigned int vectorSize    () const;
  virtual       void         resizeSequence(unsigned int newSequenceSize) = 0;
  virtual       void         resetValues   (unsigned int initialPos, unsigned int numPos) = 0;
  virtual       void         erasePositions(unsigned int initialPos, unsigned int numPos) = 0;
  virtual       void         getPositionValues(unsigned int posId,       V& vec) const = 0;
  virtual       void         setPositionValues(unsigned int posId, const V& vec) = 0;
  virtual       void         setGaussian   (gsl_rng* rng, const V& meanVec, const V& stdDevVec) = 0;

  virtual void         mean              (unsigned int              initialPos,
                                          unsigned int              numPos,
                                          V&                        meanVec) const = 0;
  virtual void         sampleVariance    (unsigned int              initialPos,
                                          unsigned int              numPos,
                                          const V&                  meanVec,
                                          V&                        samVec) const = 0;
  virtual void         populationVariance(unsigned int              initialPos,
                                          unsigned int              numPos,
                                          const V&                  meanVec,
                                          V&                        popVec) const = 0;
  virtual void         autoCovariance    (unsigned int              initialPos,
                                          unsigned int              numPos,
                                          const V&                  meanVec,
                                          unsigned int              lag,
                                          V&                        covVec) const = 0;

  virtual void         autoCorrViaDef    (unsigned int              initialPos,
                                          unsigned int              numPos,
                                          unsigned int              lag,
                                          V&                        corrVec) const = 0;
  virtual void         autoCorrViaFft    (unsigned int                     initialPos,
                                          unsigned int                     numPos,
                                          const std::vector<unsigned int>& lags,
                                          std::vector<V*>&                 corrVecs) const = 0;
  virtual void         bmm               (unsigned int              initialPos,
                                          unsigned int              batchLength,
                                          V&                        bmmVec) const = 0;
  virtual void         fftForward        (unsigned int                        initialPos,
                                          unsigned int                        fftSize,
                                          unsigned int                        paramId,
                                          std::vector<std::complex<double> >& fftResult) const = 0;
//virtual void         fftInverse        (unsigned int fftSize) = 0;
  virtual void         psd               (unsigned int              initialPos,
                                          unsigned int              numBlocks,
                                          double                    hopSizeRatio,
                                          unsigned int              paramId,
                                          std::vector<double>&      psdResult) const = 0;
  virtual void         psdAtZero         (unsigned int              initialPos,
                                          unsigned int              numBlocks,
                                          double                    hopSizeRatio,
                                          V&                        psdVec) const = 0;
  virtual void         geweke            (unsigned int              initialPos,
                                          double                    ratioNa,
                                          double                    ratioNb,
                                          V&                        gewVec) const = 0;
  virtual void         minMax            (unsigned int              initialPos,
                                          V&                        minVec,
                                          V&                        maxVec) const = 0;
  virtual void         histogram         (unsigned int              initialPos,
                                          const V&                  minVec,
                                          const V&                  maxVec,
                                          std::vector<V*>&          centersForAllBins,
                                          std::vector<V*>&          binsForAllParams) const = 0;
  virtual void         interQuantileRange(unsigned int              initialPos,
                                          V&                        iqrVec) const = 0;
  virtual void         scalesForKDE      (unsigned int              initialPos,
                                          const V&                  iqrVec,
                                          V&                        scaleVec) const = 0;
  virtual void         gaussianKDE       (unsigned int              initialPos,
                                          const std::vector<V*>&    evaluationPositions,
                                          const V&                  scaleVec,
                                          std::vector<V*>&          densityValues) const = 0;
  virtual void         write             (const std::string&        name,
                                          std::ofstream&            ofs) const = 0;
  virtual void         filter            (unsigned int              initialPos,
                                          unsigned int              spacing) = 0;

protected:
  virtual void         extractScalarSeq  (unsigned int                   initialPos,
                                          unsigned int                   spacing,
                                          unsigned int                   numPos,
                                          unsigned int                   paramId,
                                          uqScalarSequenceClass<double>& scalarSeq) const = 0;
  virtual void         extractRawData    (unsigned int                   initialPos,
                                          unsigned int                   spacing,
                                          unsigned int                   numPos,
                                          unsigned int                   paramId,
                                          std::vector<double>&           rawData) const = 0;

  const uqEnvironmentClass&   m_env;
  V                           m_vectorExample;
  mutable uqFftClass<double>* m_fftObj;
};

template <class V>
uqChainBaseClass<V>::uqChainBaseClass(
  unsigned int sequenceSize,
  const V&     vectorExample)
  :
  m_env          (vectorExample.env()),
  m_vectorExample(vectorExample),
  m_fftObj       (new uqFftClass<double>(m_env))
{
}

template <class V>
uqChainBaseClass<V>::~uqChainBaseClass()
{
  if (m_fftObj != NULL) delete m_fftObj;
}

template <class V>
const unsigned int
uqChainBaseClass<V>::vectorSize() const
{
  return m_vectorExample.size();
}

#endif // __UQ_CHAIN_H__
