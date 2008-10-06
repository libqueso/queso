/* uq/libs/mcmc/inc/uqVectorMdf.h
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

#ifndef __UQ_VECTOR_MARGINAL_DENSITY_FUNCTION_H__
#define __UQ_VECTOR_MARGINAL_DENSITY_FUNCTION_H__

#include <uqArrayOfOneDGrids.h>
#include <uqArrayOfOneDTables.h>
#include <uqEnvironment.h>
#include <math.h>

//*****************************************************
// Classes to accomodate a marginal density function
//*****************************************************

//*****************************************************
// Base class
//*****************************************************
template<class V, class M>
class uqBaseVectorMdfClass {
public:
           uqBaseVectorMdfClass(const char*                    prefix,
                                const uqVectorSpaceClass<V,M>& domainSpace,
                                const V&                       domainMinValues,
                                const V&                       domainMaxValues);
           uqBaseVectorMdfClass(const char*                    prefix,
                                const uqVectorSpaceClass<V,M>& domainSpace);
  virtual ~uqBaseVectorMdfClass();

  const   uqVectorSpaceClass<V,M>& domainSpace      ()                const;
  virtual void                     values           (const V& paramValues,
                                                           V& mdfVec) const = 0;
          bool                     outOfDomainBounds(const V& v)      const;
  const   V&                       domainMinValues  ()                const;
  const   V&                       domainMaxValues  ()                const;
  virtual void                     printContents    (std::ostream& os) const = 0;

protected:

  const   uqEnvironmentClass&      m_env;
          std::string              m_prefix;
  const   uqVectorSpaceClass<V,M>& m_domainSpace;

          V*                       m_domainMinValues;
          V*                       m_domainMaxValues;
};

template<class V, class M>
uqBaseVectorMdfClass<V,M>::uqBaseVectorMdfClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& domainSpace,
  const V&                       domainMinValues,
  const V&                       domainMaxValues)
  :
  m_env            (domainSpace.env()),
  m_prefix         ((std::string)(prefix)+"mdf_"),
  m_domainSpace    (domainSpace),
  m_domainMinValues(new V(domainMinValues)),
  m_domainMaxValues(new V(domainMaxValues))
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqBaseVectorMdfClass<V,M>::constructor() [3]"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqBaseVectorMdfClass<V,M>::constructor() [3]"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V, class M>
uqBaseVectorMdfClass<V,M>::uqBaseVectorMdfClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& domainSpace)
  :
  m_env            (domainSpace.env()),
  m_prefix         ((std::string)(prefix)+"pd_"),
  m_domainSpace    (domainSpace),
  m_domainMinValues(domainSpace.newVector(-INFINITY)),
  m_domainMaxValues(domainSpace.newVector( INFINITY))
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqBaseVectorMdfClass<V,M>::constructor() [4]"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqBaseVectorMdfClass<V,M>::constructor() [4]"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V, class M>
uqBaseVectorMdfClass<V,M>::~uqBaseVectorMdfClass()
{
  delete m_domainMinValues;
  delete m_domainMaxValues;
}

template <class V, class M>
const V&
uqBaseVectorMdfClass<V,M>::domainMinValues() const
{
  return *m_domainMinValues;
}

template <class V, class M>
const V&
uqBaseVectorMdfClass<V,M>::domainMaxValues() const
{
  return *m_domainMaxValues;
}

template <class V, class M>
bool
uqBaseVectorMdfClass<V,M>::outOfDomainBounds(const V& v) const
{
  return (v.atLeastOneComponentSmallerThan(this->domainMinValues()) ||
          v.atLeastOneComponentBiggerThan (this->domainMaxValues()));
}

template<class V, class M>
const uqVectorSpaceClass<V,M>&
uqBaseVectorMdfClass<V,M>::domainSpace() const
{
  return m_domainSpace;
}

//*****************************************************
// Generic probability distibution function class
//*****************************************************
template<class V, class M>
class uqGenericVectorMdfClass : public uqBaseVectorMdfClass<V,M> {
public:
  uqGenericVectorMdfClass(const char*                    prefix,
                          const uqVectorSpaceClass<V,M>& domainSpace,
                          const V&                       domainMinValues,
                          const V&                       domainMaxValues,
                          double (*routinePtr)(const V& paramValues, const void* routineDataPtr, V& mdfVec),
                          const void* routineDataPtr);
  uqGenericVectorMdfClass(const char*                    prefix,
                          const uqVectorSpaceClass<V,M>& domainSpace,
                          double (*routinePtr)(const V& paramValues, const void* routineDataPtr, V& mdfVec),
                          const void* routineDataPtr);
 ~uqGenericVectorMdfClass();

  void values       (const V& paramValues, V& mdfVec) const;
  void printContents(std::ostream& os) const;

protected:
  double (*m_routinePtr)(const V& paramValues, const void* routineDataPtr, V& mdfVec);
  const void* m_routineDataPtr;

  using uqBaseVectorMdfClass<V,M>::m_env;
  using uqBaseVectorMdfClass<V,M>::m_prefix;
  using uqBaseVectorMdfClass<V,M>::m_domainSpace;
  using uqBaseVectorMdfClass<V,M>::m_domainMinValues;
  using uqBaseVectorMdfClass<V,M>::m_domainMaxValues;
};

template<class V, class M>
uqGenericVectorMdfClass<V,M>::uqGenericVectorMdfClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& domainSpace,
  const V&                       domainMinValues,
  const V&                       domainMaxValues,
  double (*routinePtr)(const V& paramValues, const void* routineDataPtr, V& mdfVec),
  const void* routineDataPtr)
  :
  uqBaseVectorMdfClass<V,M>(((std::string)(prefix)+"gen").c_str(),domainSpace,domainMinValues,domainMaxValues),
  m_routinePtr    (routinePtr),
  m_routineDataPtr(routineDataPtr)
{
}

template<class V, class M>
uqGenericVectorMdfClass<V,M>::uqGenericVectorMdfClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& domainSpace,
  double (*routinePtr)(const V& paramValues, const void* routineDataPtr, V& mdfVec),
  const void* routineDataPtr)
  :
  uqBaseVectorMdfClass<V,M>(((std::string)(prefix)+"gen").c_str(),domainSpace),
  m_routinePtr    (routinePtr),
  m_routineDataPtr(routineDataPtr)
{
}

template<class V, class M>
uqGenericVectorMdfClass<V,M>::~uqGenericVectorMdfClass()
{
}

template<class V, class M>
void
uqGenericVectorMdfClass<V,M>::values(
  const V& paramValues,
        V& mdfVec) const
{
  m_routinePtr(paramValues, m_routineDataPtr, mdfVec);
  return;
}

template <class V, class M>
void
uqGenericVectorMdfClass<V,M>::printContents(std::ostream& os) const
{
  return;
}

//*****************************************************
// Gaussian probability distribution function class
//*****************************************************
template<class V, class M>
class uqGaussianVectorMdfClass : public uqBaseVectorMdfClass<V,M> {
public:
  uqGaussianVectorMdfClass(const char*                    prefix,
                           const uqVectorSpaceClass<V,M>& domainSpace,
                           const V&                       domainMinValues,
                           const V&                       domainMaxValues,
                           const V&                       domainExpectedValues,
                           const V&                       domainVarianceValues);
  uqGaussianVectorMdfClass(const char*                    prefix,
                           const uqVectorSpaceClass<V,M>& domainSpace,
                           const V&                       domainMinValues,
                           const V&                       domainMaxValues,
                           const V&                       domainExpectedValues,
                           const M&                       covMatrix);
 ~uqGaussianVectorMdfClass();

  void values       (const V& paramValues, V& mdfVec) const;
  void printContents(std::ostream& os) const;

protected:
  const M*                         m_covMatrix;

  using uqBaseVectorMdfClass<V,M>::m_env;
  using uqBaseVectorMdfClass<V,M>::m_prefix;
  using uqBaseVectorMdfClass<V,M>::m_domainSpace;
  using uqBaseVectorMdfClass<V,M>::m_domainMinValues;
  using uqBaseVectorMdfClass<V,M>::m_domainMaxValues;

  void commonConstructor();
};

template<class V,class M>
uqGaussianVectorMdfClass<V,M>::uqGaussianVectorMdfClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& domainSpace,
  const V&                       domainMinValues,
  const V&                       domainMaxValues,
  const V&                       domainExpectedValues,
  const V&                       domainVarianceValues)
  :
  uqBaseVectorMdfClass<V,M>(((std::string)(prefix)+"gau").c_str(),domainSpace,domainMinValues,domainMaxValues),
  m_covMatrix              (m_domainSpace.newDiagMatrix(domainVarianceValues*domainVarianceValues))
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqGaussianVectorMdfClass<V,M>::constructor() [1]"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  commonConstructor();

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqGaussianVectorMdfClass<V,M>::constructor() [1]"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V,class M>
uqGaussianVectorMdfClass<V,M>::uqGaussianVectorMdfClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& domainSpace,
  const V&                       domainMinValues,
  const V&                       domainMaxValues,
  const V&                       domainExpectedValues,
  const M&                       covMatrix)
  :
  uqBaseVectorMdfClass<V,M>(((std::string)(prefix)+"gau").c_str(),domainSpace,domainMinValues,domainMaxValues),
  m_covMatrix                      (new M(covMatrix))
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqGaussianVectorMdfClass<V,M>::constructor() [2]"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  commonConstructor();

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqGaussianVectorMdfClass<V,M>::constructor() [2]"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V,class M>
void
uqGaussianVectorMdfClass<V,M>::commonConstructor()
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.rank(),
                      "uqGaussianVectorMdfClass<V,M>::commonConstructor()",
                      "incomplete code");
  return;
}

template<class V,class M>
uqGaussianVectorMdfClass<V,M>::~uqGaussianVectorMdfClass()
{
  delete m_covMatrix;
}

template<class V, class M>
void
uqGaussianVectorMdfClass<V,M>::values(
  const V& paramValues,
        V& mdfVec) const
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.rank(),
                      "uqGaussianVectorMdfClass<V,M>::mdfVec()",
                      "incomplete code");
  return;
}

template <class V, class M>
void
uqGaussianVectorMdfClass<V,M>::printContents(std::ostream& os) const
{
  return;
}

//*****************************************************
// Sampled probability distribution function class
//*****************************************************
template<class V, class M>
class uqSampledVectorMdfClass : public uqBaseVectorMdfClass<V,M> {
public:
  uqSampledVectorMdfClass(const char*                          prefix,
                          const uqArrayOfOneDGridsClass <V,M>& oneDGrids,
                          const uqArrayOfOneDTablesClass<V,M>& mdfValues);
 ~uqSampledVectorMdfClass();

  void values       (const V& paramValues, V& mdfVec) const;
  void printContents(std::ostream& os) const;

protected:
  using uqBaseVectorMdfClass<V,M>::m_env;
  using uqBaseVectorMdfClass<V,M>::m_prefix;
  using uqBaseVectorMdfClass<V,M>::m_domainSpace;
  using uqBaseVectorMdfClass<V,M>::m_domainMinValues;
  using uqBaseVectorMdfClass<V,M>::m_domainMaxValues;

  const uqArrayOfOneDGridsClass <V,M>& m_oneDGrids;
  const uqArrayOfOneDTablesClass<V,M>& m_mdfValues;
};

template<class V,class M>
uqSampledVectorMdfClass<V,M>::uqSampledVectorMdfClass(
  const char*                          prefix,
  const uqArrayOfOneDGridsClass <V,M>& oneDGrids,
  const uqArrayOfOneDTablesClass<V,M>& mdfValues)
  :
  uqBaseVectorMdfClass<V,M>(((std::string)(prefix)+"sam").c_str(),oneDGrids.rowSpace(),oneDGrids.minPositions(),oneDGrids.maxPositions()),
  m_oneDGrids(oneDGrids),
  m_mdfValues(mdfValues)
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqSampledVectorMdfClass<V,M>::constructor()"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqSampledVectorMdfClass<V,M>::constructor()"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V,class M>
uqSampledVectorMdfClass<V,M>::~uqSampledVectorMdfClass()
{
}

template<class V, class M>
void
uqSampledVectorMdfClass<V,M>::values(
  const V& paramValues,
        V& mdfVec) const
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.rank(),
                      "uqSampledVectorMdfClass<V,M>::mdfVec()",
                      "incomplete code");
  return;
}

template <class V, class M>
void
uqSampledVectorMdfClass<V,M>::printContents(std::ostream& os) const
{
  // Print values *of* grid points
  os << m_oneDGrids;

  // Print *mdf* values *at* grid points
  os << m_mdfValues;

  return;
}

#endif // __UQ_VECTOR_MARGINAL_DENSITY_FUNCTION_H__
