/* uq/libs/mcmc/inc/uqVectorRealizer.h
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
                                     const uqVectorSpaceClass<V,M>& domainSpace,
                                     unsigned int                   period);

           uqBaseVectorRealizerClass(const uqBaseVectorSequenceClass<V>* chain);
  virtual ~uqBaseVectorRealizerClass();

          unsigned int period     ()              const;
  virtual void         realization(V& nextValues) const = 0;

protected:
  const uqEnvironmentClass&      m_env;
        std::string              m_prefix;
  const uqVectorSpaceClass<V,M>& m_domainSpace;
        unsigned int             m_period;

  const uqBaseVectorSequenceClass<V>* m_chain;
  mutable unsigned int                m_currentChainPos;
};

template<class V, class M>
uqBaseVectorRealizerClass<V,M>::uqBaseVectorRealizerClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& domainSpace,
  unsigned int                   period)
  :
  m_env        (domainSpace.env()),
  m_prefix     ((std::string)(prefix)+"re_"),
  m_domainSpace(domainSpace),
  m_period     (period)
{
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

//*****************************************************
// Generic class
//*****************************************************
template<class V, class M>
class uqGenericVectorRealizerClass : public uqBaseVectorRealizerClass<V,M> {
public:
  uqGenericVectorRealizerClass(const char*                    prefix,
                               const uqVectorSpaceClass<V,M>& domainSpace,
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
  using uqBaseVectorRealizerClass<V,M>::m_domainSpace;
  using uqBaseVectorRealizerClass<V,M>::m_period;
};

template<class V, class M>
uqGenericVectorRealizerClass<V,M>::uqGenericVectorRealizerClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& domainSpace,
  unsigned int                   period,
  double (*routinePtr)(const void* routineDataPtr, V& nextParamValues),
  const void* routineDataPtr)
  :
  uqBaseVectorRealizerClass<V,M>(((std::string)(prefix)+"gen").c_str(),domainSpace,period),
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
// Sequential class
//*****************************************************
template<class V, class M>
class uqSequentialVectorRealizerClass : public uqBaseVectorRealizerClass<V,M> {
public:
  uqSequentialVectorRealizerClass(const char*                         prefix,
                                  const uqVectorSpaceClass<V,M>&      domainSpace,
                                  const uqBaseVectorSequenceClass<V>& chain);
 ~uqSequentialVectorRealizerClass();

  void realization(V& nextValues) const;

private:
  const uqBaseVectorSequenceClass<V>& m_chain;
  mutable unsigned int                m_currentChainPos;

  using uqBaseVectorRealizerClass<V,M>::m_env;
  using uqBaseVectorRealizerClass<V,M>::m_prefix;
  using uqBaseVectorRealizerClass<V,M>::m_domainSpace;
  using uqBaseVectorRealizerClass<V,M>::m_period;
};

template<class V, class M>
uqSequentialVectorRealizerClass<V,M>::uqSequentialVectorRealizerClass(
  const char*                         prefix,
  const uqVectorSpaceClass<V,M>&      domainSpace,
  const uqBaseVectorSequenceClass<V>& chain)
  :
  uqBaseVectorRealizerClass<V,M>(((std::string)(prefix)+"seq").c_str(),domainSpace,chain.sequenceSize()),
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
