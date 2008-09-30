/* uq/libs/mcmc/inc/uqVectorRV.h
 *
 * Copyright (C) 2008 The PECOS Team, http://www.ices.utexas.edu/centers/pecos
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

#ifndef __UQ_VECTOR_RV_H__
#define __UQ_VECTOR_RV_H__

#include <uqScalarRV.h>
#include <uqVectorSpace.h>
#include <uqVectorProbDensity.h>
#include <uqVectorRealizer.h>
#include <uqSequenceOfVectors.h>
#include <uqArrayOfSequences.h>

//*****************************************************
// Base class
//*****************************************************
template<class V, class M>
class uqBaseVectorRVClass {
public:
  uqBaseVectorRVClass(const char*                    prefix,
                      const uqVectorSpaceClass<V,M>& imageSpace);
  virtual ~uqBaseVectorRVClass();

  const   uqEnvironmentClass&                env           () const;
  const   uqVectorSpaceClass          <V,M>& imageSpace    () const;
  const   uqBaseVectorProbDensityClass<V,M>& probDensity   () const;
  const   uqBaseVectorRealizerClass   <V,M>& realizer      () const;
          uqBaseVectorSequenceClass<V>&      chain         ();

  virtual void                               setProbDensity(const uqBaseVectorProbDensityClass<V,M>& probDensity) = 0;
  virtual void                               setRealizer   (const uqBaseVectorRealizerClass   <V,M>& realizer   ) = 0;

  virtual void                               print         (std::ostream& os) const;

protected:
  const   uqEnvironmentClass&                m_env;
          std::string                        m_prefix;
  const   uqVectorSpaceClass          <V,M>& m_imageSpace;
  const   uqBaseVectorProbDensityClass<V,M>* m_probDensity;
  const   uqBaseVectorRealizerClass   <V,M>* m_realizer;

          bool                               m_chainUse2;
          uqSequenceOfVectorsClass<V>        m_chain1;
          uqArrayOfSequencesClass<V>         m_chain2;
};

template<class V, class M>
uqBaseVectorRVClass<V,M>::uqBaseVectorRVClass(
  const char*                              prefix,
  const uqVectorSpaceClass          <V,M>& imageSpace)
  //const uqBaseVectorProbDensityClass<V,M>* probDensity,
  //const uqBaseVectorRealizerClass   <V,M>* realizer)
  :
  m_env            (imageSpace.env()),
  m_prefix         ((std::string)(prefix)+"rv_"),
  m_imageSpace     (imageSpace),
  m_probDensity    (NULL),//(probDensity),
  m_realizer       (NULL),//(realizer)
  m_chainUse2      (false),
  m_chain1         (0,m_imageSpace.zeroVector()),
  m_chain2         (0,m_imageSpace.zeroVector())
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqBaseVectorRVClass<V,M>::constructor()"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqBaseVectorRVClass<V,M>::constructor()"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V, class M>
uqBaseVectorRVClass<V,M>::~uqBaseVectorRVClass()
{
  m_chain1.clear();
  m_chain2.clear();
}

template<class V, class M>
const uqVectorSpaceClass<V,M>&
uqBaseVectorRVClass<V,M>::imageSpace() const
{
  return m_imageSpace;
}

template<class V, class M>
const uqBaseVectorProbDensityClass<V,M>&
uqBaseVectorRVClass<V,M>::probDensity() const
{
  UQ_FATAL_TEST_MACRO(m_probDensity == NULL,
                      m_env.rank(),
                      "uqBaseVectorRVClass<V,M>::probDensity()",
                      "m_probDensity is NULL");

  return *m_probDensity;
}

template<class V, class M>
const uqBaseVectorRealizerClass<V,M>&
uqBaseVectorRVClass<V,M>::realizer() const
{
  UQ_FATAL_TEST_MACRO(m_realizer == NULL,
                      m_env.rank(),
                      "uqBaseVectorRVClass<V,M>::realizer()",
                      "m_realizer is NULL");

  return *m_realizer;
}

template<class V,class M>
uqBaseVectorSequenceClass<V>&
uqBaseVectorRVClass<V,M>::chain()
{
  if (m_chainUse2) return m_chain2;
  return m_chain1;
}

template <class V, class M>
const uqEnvironmentClass&
uqBaseVectorRVClass<V,M>::env() const
{
  return m_env;
}

template <class V, class M>
void
uqBaseVectorRVClass<V,M>::print(std::ostream& os) const
{
#if 0
  os << "\nComponents are:"
     << std::endl;
  for (unsigned int i = 0; i < m_components.size(); ++i) {
    os << i << " ";
    if (m_components[i]) {
      os << *(m_components[i]);
    }
    else {
      os << "NULL";
    }
    os << std::endl;
  }
#endif
  return;
}

template<class V, class M>
std::ostream& operator<<(std::ostream& os, const uqBaseVectorRVClass<V,M>& obj)
{
  obj.print(os);

  return os;
}

//*****************************************************
// Generic class
//*****************************************************
template<class V, class M>
class uqGenericVectorRVClass : public uqBaseVectorRVClass<V,M> {
public:
  uqGenericVectorRVClass(const char*                              prefix,
                         const uqVectorSpaceClass          <V,M>& imageSpace,
                         const uqBaseVectorProbDensityClass<V,M>* probDensity,
                         const uqBaseVectorRealizerClass   <V,M>* realizer);
  virtual ~uqGenericVectorRVClass();

          void setProbDensity(const uqBaseVectorProbDensityClass<V,M>& probDensity);
          void setRealizer   (const uqBaseVectorRealizerClass   <V,M>& realizer   );

private:
  using uqBaseVectorRVClass<V,M>::m_env;
  using uqBaseVectorRVClass<V,M>::m_prefix;
  using uqBaseVectorRVClass<V,M>::m_imageSpace;
  using uqBaseVectorRVClass<V,M>::m_probDensity;
  using uqBaseVectorRVClass<V,M>::m_realizer;
};

template<class V, class M>
uqGenericVectorRVClass<V,M>::uqGenericVectorRVClass(
  const char*                              prefix,
  const uqVectorSpaceClass          <V,M>& imageSpace,
  const uqBaseVectorProbDensityClass<V,M>* probDensity,
  const uqBaseVectorRealizerClass   <V,M>* realizer)
  :
  uqBaseVectorRVClass<V,M>(prefix,imageSpace)
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqGenericVectorRVClass<V,M>::constructor()"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  m_probDensity = probDensity;
  m_realizer    = realizer;

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqGenericVectorRVClass<V,M>::constructor()"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V, class M>
uqGenericVectorRVClass<V,M>::~uqGenericVectorRVClass()
{
}

template<class V, class M>
void
uqGenericVectorRVClass<V,M>::setProbDensity(const uqBaseVectorProbDensityClass<V,M>& probDensity)
{
  m_probDensity = &probDensity;
  return;
}

template<class V, class M>
void
uqGenericVectorRVClass<V,M>::setRealizer(const uqBaseVectorRealizerClass<V,M>& realizer)
{
  m_realizer = &realizer;
  return;
}

//*****************************************************
// Gaussian class
//*****************************************************
template<class V, class M>
class uqGaussianVectorRVClass : public uqBaseVectorRVClass<V,M> {
public:
  uqGaussianVectorRVClass(const char*                    prefix,
                          const uqVectorSpaceClass<V,M>& imageSpace,
                          const V&                       imageMinValues,
                          const V&                       imageMaxValues,
                          const V&                       imageExpectValues,
                          const V&                       imageStdDevValues);
  uqGaussianVectorRVClass(const char*                    prefix,
                          const uqVectorSpaceClass<V,M>& imageSpace,
                          const V&                       imageMinValues,
                          const V&                       imageMaxValues,
                          const V&                       imageExpectValues,
                          const M&                       covMatrix);
  virtual ~uqGaussianVectorRVClass();

          void setProbDensity(const uqBaseVectorProbDensityClass<V,M>& probDensity);
          void setRealizer   (const uqBaseVectorRealizerClass   <V,M>& realizer   );

private:
  using uqBaseVectorRVClass<V,M>::m_env;
  using uqBaseVectorRVClass<V,M>::m_prefix;
  using uqBaseVectorRVClass<V,M>::m_imageSpace;
  using uqBaseVectorRVClass<V,M>::m_probDensity;
  using uqBaseVectorRVClass<V,M>::m_realizer;
};

template<class V, class M>
uqGaussianVectorRVClass<V,M>::uqGaussianVectorRVClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& imageSpace,
  const V&                       imageMinValues,
  const V&                       imageMaxValues,
  const V&                       imageExpectValues,
  const V&                       imageStdDevValues)
  :
  uqBaseVectorRVClass<V,M>(prefix,imageSpace)
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqGaussianVectorRVClass<V,M>::constructor() [1]"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  m_probDensity = new uqGaussianVectorProbDensityClass<V,M>(m_prefix.c_str(),
                                                            m_imageSpace,
                                                            imageMinValues,
                                                            imageMaxValues,
                                                            imageExpectValues,
                                                            imageStdDevValues);
  m_realizer    = NULL; // FIX ME: complete code

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqGaussianVectorRVClass<V,M>::constructor() [1]"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V, class M>
uqGaussianVectorRVClass<V,M>::uqGaussianVectorRVClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& imageSpace,
  const V&                       imageMinValues,
  const V&                       imageMaxValues,
  const V&                       imageExpectValues,
  const M&                       covMatrix)
  :
  uqBaseVectorRVClass<V,M>(prefix,imageSpace)
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqGaussianVectorRVClass<V,M>::constructor() [2]"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  m_probDensity = new uqGaussianVectorProbDensityClass<V,M>(m_prefix.c_str(),
                                                            m_imageSpace,
                                                            imageMinValues,
                                                            imageMaxValues,
                                                            imageExpectValues,
                                                            covMatrix);
  m_realizer    = NULL; // FIX ME: complete code

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqGaussianVectorRVClass<V,M>::constructor() [2]"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V, class M>
uqGaussianVectorRVClass<V,M>::~uqGaussianVectorRVClass()
{
}

template<class V, class M>
void
uqGaussianVectorRVClass<V,M>::setProbDensity(const uqBaseVectorProbDensityClass<V,M>& probDensity)
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.rank(),
                      "uqGaussianVectorRVClass<V,M>::setProbDensity()",
                      "it does not make sense to call such routine for this class");
  return;
}

template<class V, class M>
void
uqGaussianVectorRVClass<V,M>::setRealizer(const uqBaseVectorRealizerClass<V,M>& realizer)
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.rank(),
                      "uqGaussianVectorRVClass<V,M>::setRealizer()",
                      "it does not make sense to call such routine for this class");
  return;
}

#endif // __UQ_VECTOR_RV_H__
