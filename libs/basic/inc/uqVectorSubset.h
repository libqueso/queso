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

#ifndef __UQ_VECTOR_SUBSET_H__
#define __UQ_VECTOR_SUBSET_H__

#include <uqVectorSpace.h>

//*****************************************************
// Base class
//*****************************************************
template <class V, class M>
class uqVectorSubsetClass : public uqVectorSetClass<V,M>
{
public:
           uqVectorSubsetClass();
           uqVectorSubsetClass(const char*                    prefix,
                               const uqVectorSpaceClass<V,M>& vectorSpace,
                               double                         volume);
  virtual ~uqVectorSubsetClass();

           const uqVectorSpaceClass<V,M>& vectorSpace()                 const;
  virtual        bool                     contains   (const V& vec)     const = 0;
  virtual        void                     print      (std::ostream& os) const;

protected:
  using uqVectorSetClass<V,M>::m_env;
  using uqVectorSetClass<V,M>::m_prefix;

  const uqVectorSpaceClass<V,M>* m_vectorSpace;
};

template <class V, class M>
uqVectorSubsetClass<V,M>::uqVectorSubsetClass()
  :
  uqVectorSetClass<V,M>(),
  m_vectorSpace        (NULL)
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.fullRank(),
                      "uqVectorSubsetClass<V,M>::constructor(), default",
                      "should not be used by user");
}

template <class V, class M>
uqVectorSubsetClass<V,M>::uqVectorSubsetClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& vectorSpace,
  double                         volume)
  :
  uqVectorSetClass<V,M>(vectorSpace.env(),prefix,volume),
  m_vectorSpace        (&vectorSpace)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqVectorSubsetClass<V,M>::constructor()"
              << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqVectorSubsetClass<V,M>::constructor()"
              << std::endl;
  }
}

template <class V, class M>
uqVectorSubsetClass<V,M>::~uqVectorSubsetClass()
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqVectorSubsetClass<V,M>::destructor()"
                           << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqVectorSubsetClass<V,M>::destructor()"
                           << std::endl;
  }
}

template <class V, class M>
const uqVectorSpaceClass<V,M>&
uqVectorSubsetClass<V,M>::vectorSpace() const
{
  return *m_vectorSpace;
}

template <class V, class M>
void
uqVectorSubsetClass<V,M>::print(std::ostream& os) const
{
  return;
}

//*****************************************************
// Box class
//*****************************************************
template<class V, class M>
class uqBoxSubsetClass : public uqVectorSubsetClass<V,M> {
public:
  uqBoxSubsetClass(const char*                    prefix,
                   const uqVectorSpaceClass<V,M>& vectorSpace,
                   const V&                       minValues,
                   const V&                       maxValues);
 ~uqBoxSubsetClass();

        bool contains (const V& vec)     const;
  const V&   minValues()                 const;
  const V&   maxValues()                 const;
        void print    (std::ostream& os) const;

protected:
  using uqVectorSetClass   <V,M>::m_env;
  using uqVectorSetClass   <V,M>::m_prefix;
  using uqVectorSetClass   <V,M>::m_volume;
  using uqVectorSubsetClass<V,M>::m_vectorSpace;

  V m_minValues;
  V m_maxValues;
};

template<class V, class M>
uqBoxSubsetClass<V,M>::uqBoxSubsetClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& vectorSpace,
  const V&                       minValues,
  const V&                       maxValues)
  :
  uqVectorSubsetClass<V,M>(prefix,vectorSpace,0.),
  m_minValues(minValues),
  m_maxValues(maxValues)
{
  UQ_FATAL_TEST_MACRO(minValues.sizeLocal() != maxValues.sizeLocal(),
                      m_env.fullRank(),
                      "uqBoxSubsetClass<V,M>::uqBoxSubsetClass()",
                      "vectors 'minValues' and 'maxValues' should have the same size");
  UQ_FATAL_TEST_MACRO(minValues.sizeLocal() != vectorSpace.dimLocal(),
                      m_env.fullRank(),
                      "uqBoxSubsetClass<V,M>::uqBoxSubsetClass()",
                      "sizes of vectors 'minValues' and 'maxValues' should be equal to dimension of the vector space");
  for (unsigned int i = 0; i < m_vectorSpace->dimLocal(); ++i) {
    UQ_FATAL_TEST_MACRO(minValues[i] > maxValues[i],
                        m_env.fullRank(),
                        "uqBoxSubsetClass<V,M>::uqBoxSubsetClass()",
                        "it should happen minValue <= maxValue for all dimensions");
  }

  m_volume = 1.;
  for (unsigned int i = 0; i < m_vectorSpace->dimLocal(); ++i) {
    m_volume *= (m_maxValues[i] - m_minValues[i]);
  }
}

template<class V, class M>
uqBoxSubsetClass<V,M>::~uqBoxSubsetClass()
{
}

template<class V, class M>
bool
uqBoxSubsetClass<V,M>::contains(const V& vec) const
{
  //bool result = true;

  //for (unsigned int i = 0; (i < m_vectorSpace->dimLocal()) && (result == true); ++i) {
  //  result = (m_maxValues[i] <= vec[i]) && (vec[i] <= m_minValues[i]);
  //}

  //return result;

  return (!vec.atLeastOneComponentSmallerThan(m_minValues) &&
          !vec.atLeastOneComponentBiggerThan (m_maxValues));
}

template<class V, class M>
const V&
uqBoxSubsetClass<V,M>::minValues() const
{
  return m_minValues;
}

template<class V, class M>
const V&
uqBoxSubsetClass<V,M>::maxValues() const
{
  return m_maxValues;
}

template <class V, class M>
void
uqBoxSubsetClass<V,M>::print(std::ostream& os) const
{
  os << "In uqBoxSubsetClass<V,M>::print()"
     << ": m_minValues = " << m_minValues
     << ", m_maxValues = " << m_maxValues
     << ", m_volume = "    << m_volume
     << std::endl;

  return;
}

//*****************************************************
// Intersection class
//*****************************************************
template<class V, class M>
class uqIntersectionSubsetClass : public uqVectorSubsetClass<V,M> {
public:
  uqIntersectionSubsetClass(const char*                    prefix,
                            const uqVectorSpaceClass<V,M>& vectorSpace,
                                  double                   volume,
                            const uqVectorSetClass<V,M>&   set1,
                            const uqVectorSetClass<V,M>&   set2);
 ~uqIntersectionSubsetClass();

  bool contains(const V& vec)     const;
  void print   (std::ostream& os) const;

protected:
  using uqVectorSetClass   <V,M>::m_env;
  using uqVectorSetClass   <V,M>::m_prefix;
  using uqVectorSetClass   <V,M>::m_volume;
  using uqVectorSubsetClass<V,M>::m_vectorSpace;

  const uqVectorSetClass<V,M>& m_set1;
  const uqVectorSetClass<V,M>& m_set2;
};

template<class V, class M>
uqIntersectionSubsetClass<V,M>::uqIntersectionSubsetClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& vectorSpace,
        double                   volume,
  const uqVectorSetClass<V,M>&   set1,
  const uqVectorSetClass<V,M>&   set2)
  :
  uqVectorSubsetClass<V,M>(prefix,vectorSpace,volume),
  m_set1                  (set1),
  m_set2                  (set2)
{
}

template<class V, class M>
uqIntersectionSubsetClass<V,M>::~uqIntersectionSubsetClass()
{
}

template<class V, class M>
bool
uqIntersectionSubsetClass<V,M>::contains(const V& vec) const
{
  return (m_set1.contains(vec) && m_set2.contains(vec));
}

template <class V, class M>
void
uqIntersectionSubsetClass<V,M>::print(std::ostream& os) const
{
  os << "In uqIntersectionSubsetClass<V,M>::print()"
     << ": m_set1 = " << m_set1
     << ", m_set2 = " << m_set2
     << std::endl;

  return;
}

//*****************************************************
// Concatenation class
//*****************************************************
template<class V, class M>
class uqConcatenationSubsetClass : public uqVectorSubsetClass<V,M> {
public:
  uqConcatenationSubsetClass(const char*                    prefix,
                             const uqVectorSpaceClass<V,M>& vectorSpace,
                             const uqVectorSetClass<V,M>&   set1,
                             const uqVectorSetClass<V,M>&   set2);
 ~uqConcatenationSubsetClass();

  bool contains(const V& vec)     const;
  void print   (std::ostream& os) const;

protected:
  using uqVectorSetClass   <V,M>::m_env;
  using uqVectorSetClass   <V,M>::m_prefix;
  using uqVectorSetClass   <V,M>::m_volume;
  using uqVectorSubsetClass<V,M>::m_vectorSpace;

  const uqVectorSetClass<V,M>& m_set1;
  const uqVectorSetClass<V,M>& m_set2;
};

template<class V, class M>
uqConcatenationSubsetClass<V,M>::uqConcatenationSubsetClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& vectorSpace,
  const uqVectorSetClass<V,M>&   set1,
  const uqVectorSetClass<V,M>&   set2)
  :
  uqVectorSubsetClass<V,M>(prefix,vectorSpace,set1.volume()*set2.volume()),
  m_set1                  (set1),
  m_set2                  (set2)
{
}

template<class V, class M>
uqConcatenationSubsetClass<V,M>::~uqConcatenationSubsetClass()
{
}

template<class V, class M>
bool
uqConcatenationSubsetClass<V,M>::contains(const V& vec) const
{
  V v1(m_set1.vectorSpace().zeroVector());
  V v2(m_set2.vectorSpace().zeroVector());

  UQ_FATAL_TEST_MACRO((vec.sizeLocal() != (v1.sizeLocal()+v2.sizeLocal())),
                      m_env.fullRank(),
                      "uqConcatenationSubsetClass<V,M>::contains()",
                      "incompatible vector sizes");

  vec.cwExtract(             0,v1);
  vec.cwExtract(v1.sizeLocal(),v2);

  return (m_set1.contains(v1) && m_set2.contains(v2));
}

template <class V, class M>
void
uqConcatenationSubsetClass<V,M>::print(std::ostream& os) const
{
  os << "In uqConcatenationSubsetClass<V,M>::print()"
     << ": m_set1 = " << m_set1
     << ", m_set2 = " << m_set2
     << std::endl;

  return;
}

//*****************************************************
// Discrete class
//*****************************************************
template<class V, class M>
class uqDiscreteSubsetClass : public uqVectorSubsetClass<V,M> {
public:
  uqDiscreteSubsetClass(const char*                    prefix,
                        const uqVectorSpaceClass<V,M>& vectorSpace,
                        const std::vector<V*>&         elements);
 ~uqDiscreteSubsetClass();

        bool contains (const V& vec)     const;
        void print    (std::ostream& os) const;

protected:
  using uqVectorSetClass   <V,M>::m_env;
  using uqVectorSetClass   <V,M>::m_prefix;
  using uqVectorSetClass   <V,M>::m_volume;
  using uqVectorSubsetClass<V,M>::m_vectorSpace;

  std::vector<V*> m_elements;
};

template<class V, class M>
uqDiscreteSubsetClass<V,M>::uqDiscreteSubsetClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& vectorSpace,
  const std::vector<V*>&         elements)
  :
  uqVectorSubsetClass<V,M>(prefix,vectorSpace,0.),
  m_elements(elements.size(),NULL)
{
  m_volume = 0.;
  UQ_FATAL_TEST_MACRO(true,
                      m_env.fullRank(),
                      "uqDiscreteSubsetClass<V,M>::contains()",
                      "incomplete code");
}

template<class V, class M>
uqDiscreteSubsetClass<V,M>::~uqDiscreteSubsetClass()
{
}

template<class V, class M>
bool
uqDiscreteSubsetClass<V,M>::contains(const V& vec) const
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.fullRank(),
                      "uqDiscreteSubsetClass<V,M>::contains()",
                      "incomplete code");

  return false;
}

template <class V, class M>
void
uqDiscreteSubsetClass<V,M>::print(std::ostream& os) const
{
  os << "In uqBoxSubsetClass<V,M>::print()"
     << ": nothing to print"
     << std::endl;

  return;
}
#endif // __UQ_VECTOR_SUBSET_H__

