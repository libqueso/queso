/* uq/libs/queso/inc/uqVectorRealizer.h
 *
 * Copyright (C) 2008 The QUESO Team, http://queso.ices.utexas.edu/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __UQ_REALIZER_H__
#define __UQ_REALIZER_H__

#include <uqVectorSequence.h>
#include <uqEnvironment.h>
#include <math.h>

//*****************************************************
// Classes to accomodate a probability density routine.
//*****************************************************

//*****************************************************
// Base class
//*****************************************************
template<class V, class M>
class uqBaseVectorRealizerClass {
public:
           uqBaseVectorRealizerClass(const char*                    prefix,
                                     const uqVectorSpaceClass<V,M>& imageSpace,
                                     unsigned int                   period,
                                     const V&                       imageMinValues,
                                     const V&                       imageMaxValues,
                                     const V&                       imageExpectedValues,
                                     const V&                       imageVarianceValues);
           uqBaseVectorRealizerClass(const char*                    prefix,
                                     const uqVectorSpaceClass<V,M>& imageSpace,
                                     unsigned int                   period,
                                     const V&                       imageMinValues,
                                     const V&                       imageMaxValues,
                                     const V&                       imageExpectedValues);
           uqBaseVectorRealizerClass(const char*                    prefix,
                                     const uqVectorSpaceClass<V,M>& imageSpace,
                                     unsigned int                   period,
                                     const V&                       imageMinValues,
                                     const V&                       imageMaxValues);
           uqBaseVectorRealizerClass(const char*                    prefix,
                                     const uqVectorSpaceClass<V,M>& imageSpace,
                                     unsigned int                   period);

  virtual ~uqBaseVectorRealizerClass();

  const   uqVectorSpaceClass<V,M>& imageSpace         ()              const;
          unsigned int             period             ()              const;
  virtual void                     realization        (V& nextValues) const = 0;

          bool                     outOfImageBounds   (const V& v)    const;
  const   V&                       imageMinValues     ()              const;
  const   V&                       imageMaxValues     ()              const;
  const   V&                       imageExpectedValues()              const;
  const   V&                       imageVarianceValues()              const;

protected:
  const uqBaseEnvironmentClass&      m_env;
        std::string              m_prefix;
  const uqVectorSpaceClass<V,M>& m_imageSpace;
        unsigned int             m_period;

        V*                       m_imageMinValues;
        V*                       m_imageMaxValues;
        V*                       m_imageExpectedValues;
        V*                       m_imageVarianceValues;
};

template<class V, class M>
uqBaseVectorRealizerClass<V,M>::uqBaseVectorRealizerClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& imageSpace,
  unsigned int                   period,
  const V&                       imageMinValues,
  const V&                       imageMaxValues,
  const V&                       imageExpectedValues,
  const V&                       imageVarianceValues)
  :
  m_env                (imageSpace.env()),
  m_prefix             ((std::string)(prefix)+"re_"),
  m_imageSpace         (imageSpace),
  m_period             (period),
  m_imageMinValues     (new V(imageMinValues     )),
  m_imageMaxValues     (new V(imageMaxValues     )),
  m_imageExpectedValues(new V(imageExpectedValues)),
  m_imageVarianceValues(new V(imageVarianceValues))
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqBaseVectorRealizerClass<V,M>::constructor() [1]"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqBaseVectorRealizerClass<V,M>::constructor() [1]"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V, class M>
uqBaseVectorRealizerClass<V,M>::uqBaseVectorRealizerClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& imageSpace,
  unsigned int                   period,
  const V&                       imageMinValues,
  const V&                       imageMaxValues,
  const V&                       imageExpectedValues)
  :
  m_env                (imageSpace.env()),
  m_prefix             ((std::string)(prefix)+"re_"),
  m_imageSpace         (imageSpace),
  m_period             (period),
  m_imageMinValues     (new V(imageMinValues)         ),
  m_imageMaxValues     (new V(imageMaxValues)         ),
  m_imageExpectedValues(new V(imageExpectedValues)    ),
  m_imageVarianceValues(imageSpace.newVector(INFINITY))
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqBaseVectorRealizerClass<V,M>::constructor() [2]"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqBaseVectorRealizerClass<V,M>::constructor() [2]"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V, class M>
uqBaseVectorRealizerClass<V,M>::uqBaseVectorRealizerClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& imageSpace,
  unsigned int                   period,
  const V&                       imageMinValues,
  const V&                       imageMaxValues)
  :
  m_env                (imageSpace.env()),
  m_prefix             ((std::string)(prefix)+"re_"),
  m_imageSpace         (imageSpace),
  m_period             (period),
  m_imageMinValues     (new V(imageMinValues)         ),
  m_imageMaxValues     (new V(imageMaxValues)         ),
  m_imageExpectedValues(imageSpace.newVector(      0.)),
  m_imageVarianceValues(imageSpace.newVector(INFINITY))
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqBaseVectorRealizerClass<V,M>::constructor() [3]"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqBaseVectorRealizerClass<V,M>::constructor() [3]"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V, class M>
uqBaseVectorRealizerClass<V,M>::uqBaseVectorRealizerClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& imageSpace,
  unsigned int                   period)
  :
  m_env                (imageSpace.env()),
  m_prefix             ((std::string)(prefix)+"re_"),
  m_imageSpace         (imageSpace),
  m_period             (period),
  m_imageMinValues     (imageSpace.newVector(-INFINITY)),
  m_imageMaxValues     (imageSpace.newVector( INFINITY)),
  m_imageExpectedValues(imageSpace.newVector(       0.)),
  m_imageVarianceValues(imageSpace.newVector( INFINITY))
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqBaseVectorRealizerClass<V,M>::constructor() [4]"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqBaseVectorRealizerClass<V,M>::constructor() [4]"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V, class M>
uqBaseVectorRealizerClass<V,M>::~uqBaseVectorRealizerClass()
{
}

template<class V, class M>
unsigned int
uqBaseVectorRealizerClass<V,M>::period() const
{
  return m_period;
}

template <class V, class M>
const V&
uqBaseVectorRealizerClass<V,M>::imageMinValues() const
{
  return *m_imageMinValues;
}

template <class V, class M>
const V&
uqBaseVectorRealizerClass<V,M>::imageMaxValues() const
{
  return *m_imageMaxValues;
}

template <class V, class M>
const V&
uqBaseVectorRealizerClass<V,M>::imageExpectedValues() const
{
  return *m_imageExpectedValues;
}

template <class V, class M>
const V&
uqBaseVectorRealizerClass<V,M>::imageVarianceValues() const
{
  return *m_imageVarianceValues;
}

template <class V, class M>
bool
uqBaseVectorRealizerClass<V,M>::outOfImageBounds(const V& v) const
{
  return (v.atLeastOneComponentSmallerThan(this->imageMinValues()) ||
          v.atLeastOneComponentBiggerThan (this->imageMaxValues()));
}

template<class V, class M>
const uqVectorSpaceClass<V,M>&
uqBaseVectorRealizerClass<V,M>::imageSpace() const
{
  return m_imageSpace;
}

//*****************************************************
// Generic class
//*****************************************************
template<class V, class M>
class uqGenericVectorRealizerClass : public uqBaseVectorRealizerClass<V,M> {
public:
  uqGenericVectorRealizerClass(const char*                    prefix,
                               const uqVectorSpaceClass<V,M>& imageSpace,
                               unsigned int                   period,
                               double (*routinePtr)(const void* routineDataPtr, V& nextParamValues),
                               const void* routineDataPtr);
 ~uqGenericVectorRealizerClass();

  void realization(V& nextValues) const;

private:
  double (*m_routinePtr)(const void* routineDataPtr, V& nextParamValues);
  const void* m_routineDataPtr;

  using uqBaseVectorRealizerClass<V,M>::m_env;
  using uqBaseVectorRealizerClass<V,M>::m_prefix;
  using uqBaseVectorRealizerClass<V,M>::m_imageSpace;
  using uqBaseVectorRealizerClass<V,M>::m_period;
};

template<class V, class M>
uqGenericVectorRealizerClass<V,M>::uqGenericVectorRealizerClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& imageSpace,
  unsigned int                   period,
  double (*routinePtr)(const void* routineDataPtr, V& nextParamValues),
  const void* routineDataPtr)
  :
  uqBaseVectorRealizerClass<V,M>(((std::string)(prefix)+"gen").c_str(),imageSpace,period),
  m_routinePtr    (routinePtr),
  m_routineDataPtr(routineDataPtr)
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqGenericVectorRealizerClass<V,M>::constructor()"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqGenericVectorRealizerClass<V,M>::constructor()"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V, class M>
uqGenericVectorRealizerClass<V,M>::~uqGenericVectorRealizerClass()
{
}

template<class V, class M>
void
uqGenericVectorRealizerClass<V,M>::realization(V& nextValues) const
{
  m_routinePtr(m_routineDataPtr,nextValues);
  return;
}

//*****************************************************
// Gaussian class
//*****************************************************
template<class V, class M>
class uqGaussianVectorRealizerClass : public uqBaseVectorRealizerClass<V,M> {
public:
  uqGaussianVectorRealizerClass(const char* prefix,
				const uqVectorSpaceClass<V,M>& imageSpace,
				const V& expectedValues, // vector of mean values
				const M& lowerCholCovMatrix); // lower triangular matrix resulting from Cholesky decomposition of the covariance matrix

  ~uqGaussianVectorRealizerClass();

  void realization(V& nextValues) const;
  void updateExpectedValues(const V& newExpectedValues);
  void updateLowerCholCovMatrix(const V& newLowerCholCovMatrix);
    
private:
  M* m_lowerCholCovMatrix;

  using uqBaseVectorRealizerClass<V,M>::m_env;
  using uqBaseVectorRealizerClass<V,M>::m_prefix;
  using uqBaseVectorRealizerClass<V,M>::m_imageSpace;
  using uqBaseVectorRealizerClass<V,M>::m_period;
  using uqBaseVectorRealizerClass<V,M>::m_imageExpectedValues;
};

template<class V, class M>
uqGaussianVectorRealizerClass<V,M>::uqGaussianVectorRealizerClass(const char* prefix,
								  const uqVectorSpaceClass<V,M>& imageSpace,
								  const V& expectedValues,
								  const M& lowerCholCovMatrix)
  :
  uqBaseVectorRealizerClass<V,M>( ((std::string)(prefix)+"gau").c_str(), imageSpace, 0 ),
  m_lowerCholCovMatrix(new M(lowerCholCovMatrix))
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqGaussianVectorRealizerClass<V,M>::constructor()"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  *m_imageExpectedValues = expectedValues;

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqGaussianVectorRealizerClass<V,M>::constructor()"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}
								  
template<class V, class M>
uqGaussianVectorRealizerClass<V,M>::~uqGaussianVectorRealizerClass()
{
  delete m_lowerCholCovMatrix;
}

template<class V, class M>
void
uqGaussianVectorRealizerClass<V,M>::realization(V& nextValues) const
{
  V iidGaussianVector(m_imageSpace.zeroVector());
  iidGaussianVector.cwSetGaussian(m_env.rng(), 0.0, 1.0);

  nextValues = (*m_imageExpectedValues) + (*m_lowerCholCovMatrix)*iidGaussianVector;
  return;
}

template<class V, class M>
void
uqGaussianVectorRealizerClass<V,M>::updateExpectedValues(const V& newExpectedValues)
{
  // delete old expected values (alloced at construction or last call to this function)
  delete m_imageExpectedValues;
  m_imageExpectedValues = new V(newExpectedValues);
  return;
}

template<class V, class M>
void
uqGaussianVectorRealizerClass<V,M>::updateLowerCholCovMatrix(const V& newLowerCholCovMatrix)
{
  // delete old expected values (alloced at construction or last call to this function)
  delete m_lowerCholCovMatrix;
  m_lowerCholCovMatrix = new V(newLowerCholCovMatrix);
  return;
}

//*****************************************************
// Sequential class
//*****************************************************
template<class V, class M>
class uqSequentialVectorRealizerClass : public uqBaseVectorRealizerClass<V,M> {
public:
  uqSequentialVectorRealizerClass(const char*                           prefix,
                                  const uqBaseVectorSequenceClass<V,M>& chain);
 ~uqSequentialVectorRealizerClass();

  void realization(V& nextValues) const;

private:
  const uqBaseVectorSequenceClass<V,M>& m_chain;
  mutable unsigned int                  m_currentChainPos;

  using uqBaseVectorRealizerClass<V,M>::m_env;
  using uqBaseVectorRealizerClass<V,M>::m_prefix;
  using uqBaseVectorRealizerClass<V,M>::m_imageSpace;
  using uqBaseVectorRealizerClass<V,M>::m_period;
};

template<class V, class M>
uqSequentialVectorRealizerClass<V,M>::uqSequentialVectorRealizerClass(
  const char*                           prefix,
  const uqBaseVectorSequenceClass<V,M>& chain)
  :
  uqBaseVectorRealizerClass<V,M>(((std::string)(prefix)+"seq").c_str(),chain.vectorSpace(),chain.sequenceSize(),chain.minValues(),chain.maxValues(),chain.meanValues(),chain.sampleVarianceValues()),
  m_chain          (chain),
  m_currentChainPos(0)
{
  if ((m_env.verbosity() >= 0) && (m_env.rank() == 0)) {
    std::cout << "In uqSequentialVectorRealizerClass<V,M>::constructor()"
              << ": m_chain.sequenceSize() = " << m_chain.sequenceSize()
              << std::endl;
  }
}

template<class V, class M>
uqSequentialVectorRealizerClass<V,M>::~uqSequentialVectorRealizerClass()
{
}

template<class V, class M>
void
uqSequentialVectorRealizerClass<V,M>::realization(V& nextParamValues) const
{
  m_chain.getPositionValues(m_currentChainPos++,nextParamValues);
  if (m_currentChainPos >= m_period) m_currentChainPos = 0;

  return;
}
#endif // __UQ_REALIZER_H__
