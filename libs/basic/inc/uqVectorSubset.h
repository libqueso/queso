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
                               const uqVectorSpaceClass<V,M>& vectorSpace);
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
                      m_env.rank(),
                      "uqVectorSubsetClass<V,M>::constructor(), default",
                      "should not be used by user");
}

template <class V, class M>
uqVectorSubsetClass<V,M>::uqVectorSubsetClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& vectorSpace)
  :
  uqVectorSetClass<V,M>(vectorSpace.env(),prefix,0.),
  m_vectorSpace        (&vectorSpace)
{
  if ((m_env.subScreenFile()) && (m_env.verbosity() >= 5)) {
    *m_env.subScreenFile() << "Entering uqVectorSubsetClass<V,M>::constructor()"
              << std::endl;
  }

  if ((m_env.subScreenFile()) && (m_env.verbosity() >= 5)) {
    *m_env.subScreenFile() << "Leaving uqVectorSubsetClass<V,M>::constructor()"
              << std::endl;
  }
}

template <class V, class M>
uqVectorSubsetClass<V,M>::~uqVectorSubsetClass()
{
  if ((m_env.subScreenFile()) && (m_env.verbosity() >= 5)) {
    *m_env.subScreenFile() << "Entering uqVectorSubsetClass<V,M>::destructor()"
                           << std::endl;
  }

  if ((m_env.subScreenFile()) && (m_env.verbosity() >= 5)) {
    *m_env.subScreenFile() << "Leaving uqVectorSubsetClass<V,M>::destructor()"
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
  uqVectorSubsetClass<V,M>(prefix,vectorSpace),
  m_minValues(minValues),
  m_maxValues(maxValues)
{
  m_volume = 1.;
  for (unsigned int i = 0; i < m_vectorSpace->dim(); ++i) {
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

  //for (unsigned int i = 0; (i < m_vectorSpace->dim()) && (result == true); ++i) {
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
  uqVectorSubsetClass<V,M>(prefix,vectorSpace),
  m_set1                  (set1),
  m_set2                  (set2)
{
  m_volume = volume;
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
  return;
}
#endif // __UQ_VECTOR_SUBSET_H__

