/* uq/libs/mcmc/inc/uqVectorCdf.h
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

#ifndef __UQ_VECTOR_CUMULATIVE_DISTRIBUTION_FUNCTION_H__
#define __UQ_VECTOR_CUMULATIVE_DISTRIBUTION_FUNCTION_H__

#include <uqArrayOfOneDUniformGrids.h>
#include <uqArrayOfScalarSets.h>
#include <uqEnvironment.h>
#include <math.h>

//*****************************************************
// Classes to accomodate a cumulative distribution function
//*****************************************************

//*****************************************************
// Base class
//*****************************************************
template<class V, class M>
class uqBaseVectorCdfClass {
public:
           uqBaseVectorCdfClass(const char*                    prefix,
                                const uqVectorSpaceClass<V,M>& domainSpace,
                                const V&                       domainMinValues,
                                const V&                       domainMaxValues);
           uqBaseVectorCdfClass(const char*                    prefix,
                                const uqVectorSpaceClass<V,M>& domainSpace);
  virtual ~uqBaseVectorCdfClass();

  const   uqVectorSpaceClass<V,M>& domainSpace      ()                const;
  virtual void                     values           (const V& paramValues,
                                                           V& cdfVec) const = 0;
          bool                     outOfDomainBounds(const V& v)      const;
  const   V&                       domainMinValues  ()                const;
  const   V&                       domainMaxValues  ()                const;
  virtual void                     printContents    (std::ofstream& ofs) const = 0;

protected:

  const   uqEnvironmentClass&      m_env;
          std::string              m_prefix;
  const   uqVectorSpaceClass<V,M>& m_domainSpace;

          V*                       m_domainMinValues;
          V*                       m_domainMaxValues;
};

template<class V, class M>
uqBaseVectorCdfClass<V,M>::uqBaseVectorCdfClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& domainSpace,
  const V&                       domainMinValues,
  const V&                       domainMaxValues)
  :
  m_env            (domainSpace.env()),
  m_prefix         ((std::string)(prefix)+"cdf_"),
  m_domainSpace    (domainSpace),
  m_domainMinValues(new V(domainMinValues)),
  m_domainMaxValues(new V(domainMaxValues))
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqBaseVectorCdfClass<V,M>::constructor() [3]"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqBaseVectorCdfClass<V,M>::constructor() [3]"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V, class M>
uqBaseVectorCdfClass<V,M>::uqBaseVectorCdfClass(
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
    std::cout << "Entering uqBaseVectorCdfClass<V,M>::constructor() [4]"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqBaseVectorCdfClass<V,M>::constructor() [4]"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V, class M>
uqBaseVectorCdfClass<V,M>::~uqBaseVectorCdfClass()
{
  delete m_domainMinValues;
  delete m_domainMaxValues;
}

template <class V, class M>
const V&
uqBaseVectorCdfClass<V,M>::domainMinValues() const
{
  return *m_domainMinValues;
}

template <class V, class M>
const V&
uqBaseVectorCdfClass<V,M>::domainMaxValues() const
{
  return *m_domainMaxValues;
}

template <class V, class M>
bool
uqBaseVectorCdfClass<V,M>::outOfDomainBounds(const V& v) const
{
  return (v.atLeastOneComponentSmallerThan(this->domainMinValues()) ||
          v.atLeastOneComponentBiggerThan (this->domainMaxValues()));
}

template<class V, class M>
const uqVectorSpaceClass<V,M>&
uqBaseVectorCdfClass<V,M>::domainSpace() const
{
  return m_domainSpace;
}

//*****************************************************
// Generic probability distibution function class
//*****************************************************
template<class V, class M>
class uqGenericVectorCdfClass : public uqBaseVectorCdfClass<V,M> {
public:
  uqGenericVectorCdfClass(const char*                    prefix,
                          const uqVectorSpaceClass<V,M>& domainSpace,
                          const V&                       domainMinValues,
                          const V&                       domainMaxValues,
                          double (*routinePtr)(const V& paramValues, const void* routineDataPtr, V& cdfVec),
                          const void* routineDataPtr);
  uqGenericVectorCdfClass(const char*                    prefix,
                          const uqVectorSpaceClass<V,M>& domainSpace,
                          double (*routinePtr)(const V& paramValues, const void* routineDataPtr, V& cdfVec),
                          const void* routineDataPtr);
 ~uqGenericVectorCdfClass();

  void values       (const V& paramValues, V& cdfVec) const;
  void printContents(std::ofstream& ofs) const;

protected:
  double (*m_routinePtr)(const V& paramValues, const void* routineDataPtr, V& cdfVec);
  const void* m_routineDataPtr;

  using uqBaseVectorCdfClass<V,M>::m_env;
  using uqBaseVectorCdfClass<V,M>::m_prefix;
  using uqBaseVectorCdfClass<V,M>::m_domainSpace;
  using uqBaseVectorCdfClass<V,M>::m_domainMinValues;
  using uqBaseVectorCdfClass<V,M>::m_domainMaxValues;
};

template<class V, class M>
uqGenericVectorCdfClass<V,M>::uqGenericVectorCdfClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& domainSpace,
  const V&                       domainMinValues,
  const V&                       domainMaxValues,
  double (*routinePtr)(const V& paramValues, const void* routineDataPtr, V& cdfVec),
  const void* routineDataPtr)
  :
  uqBaseVectorCdfClass<V,M>(((std::string)(prefix)+"gen").c_str(),domainSpace,domainMinValues,domainMaxValues),
  m_routinePtr    (routinePtr),
  m_routineDataPtr(routineDataPtr)
{
}

template<class V, class M>
uqGenericVectorCdfClass<V,M>::uqGenericVectorCdfClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& domainSpace,
  double (*routinePtr)(const V& paramValues, const void* routineDataPtr, V& cdfVec),
  const void* routineDataPtr)
  :
  uqBaseVectorCdfClass<V,M>(((std::string)(prefix)+"gen").c_str(),domainSpace),
  m_routinePtr    (routinePtr),
  m_routineDataPtr(routineDataPtr)
{
}

template<class V, class M>
uqGenericVectorCdfClass<V,M>::~uqGenericVectorCdfClass()
{
}

template<class V, class M>
void
uqGenericVectorCdfClass<V,M>::values(
  const V& paramValues,
        V& cdfVec) const
{
  m_routinePtr(paramValues, m_routineDataPtr, cdfVec);
  return;
}

template <class V, class M>
void
uqGenericVectorCdfClass<V,M>::printContents(std::ofstream& ofs) const
{
  return;
}

//*****************************************************
// Gaussian probability distribution function class
//*****************************************************
template<class V, class M>
class uqGaussianVectorCdfClass : public uqBaseVectorCdfClass<V,M> {
public:
  uqGaussianVectorCdfClass(const char*                    prefix,
                           const uqVectorSpaceClass<V,M>& domainSpace,
                           const V&                       domainMinValues,
                           const V&                       domainMaxValues,
                           const V&                       domainExpectedValues,
                           const V&                       domainVarianceValues);
  uqGaussianVectorCdfClass(const char*                    prefix,
                           const uqVectorSpaceClass<V,M>& domainSpace,
                           const V&                       domainMinValues,
                           const V&                       domainMaxValues,
                           const V&                       domainExpectedValues,
                           const M&                       covMatrix);
 ~uqGaussianVectorCdfClass();

  void values       (const V& paramValues, V& cdfVec) const;
  void printContents(std::ofstream& ofs) const;

protected:
  const M*                         m_covMatrix;

  using uqBaseVectorCdfClass<V,M>::m_env;
  using uqBaseVectorCdfClass<V,M>::m_prefix;
  using uqBaseVectorCdfClass<V,M>::m_domainSpace;
  using uqBaseVectorCdfClass<V,M>::m_domainMinValues;
  using uqBaseVectorCdfClass<V,M>::m_domainMaxValues;

  void commonConstructor();
};

template<class V,class M>
uqGaussianVectorCdfClass<V,M>::uqGaussianVectorCdfClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& domainSpace,
  const V&                       domainMinValues,
  const V&                       domainMaxValues,
  const V&                       domainExpectedValues,
  const V&                       domainVarianceValues)
  :
  uqBaseVectorCdfClass<V,M>(((std::string)(prefix)+"gau").c_str(),domainSpace,domainMinValues,domainMaxValues),
  m_covMatrix              (m_domainSpace.newDiagMatrix(domainVarianceValues*domainVarianceValues))
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqGaussianVectorCdfClass<V,M>::constructor() [1]"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  commonConstructor();

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqGaussianVectorCdfClass<V,M>::constructor() [1]"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V,class M>
uqGaussianVectorCdfClass<V,M>::uqGaussianVectorCdfClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& domainSpace,
  const V&                       domainMinValues,
  const V&                       domainMaxValues,
  const V&                       domainExpectedValues,
  const M&                       covMatrix)
  :
  uqBaseVectorCdfClass<V,M>(((std::string)(prefix)+"gau").c_str(),domainSpace,domainMinValues,domainMaxValues),
  m_covMatrix                      (new M(covMatrix))
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqGaussianVectorCdfClass<V,M>::constructor() [2]"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  commonConstructor();

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqGaussianVectorCdfClass<V,M>::constructor() [2]"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V,class M>
void
uqGaussianVectorCdfClass<V,M>::commonConstructor()
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.rank(),
                      "uqGaussianVectorCdfClass<V,M>::commonConstructor()",
                      "incomplete code");
  return;
}

template<class V,class M>
uqGaussianVectorCdfClass<V,M>::~uqGaussianVectorCdfClass()
{
  delete m_covMatrix;
}

template<class V, class M>
void
uqGaussianVectorCdfClass<V,M>::values(
  const V& paramValues,
        V& cdfVec) const
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.rank(),
                      "uqGaussianVectorCdfClass<V,M>::cdfVec()",
                      "incomplete code");
  return;
}

template <class V, class M>
void
uqGaussianVectorCdfClass<V,M>::printContents(std::ofstream& ofs) const
{
  return;
}

//*****************************************************
// Sampled probability distribution function class
//*****************************************************
template<class V, class M>
class uqSampledVectorCdfClass : public uqBaseVectorCdfClass<V,M> {
public:
  uqSampledVectorCdfClass(const char*                                prefix,
                          const uqArrayOfOneDUniformGridsClass<V,M>& oneDGrids,
                          const uqArrayOfScalarSetsClass      <V,M>& cdfValues);
 ~uqSampledVectorCdfClass();

  void values       (const V& paramValues, V& cdfVec) const;
  void printContents(std::ofstream& ofs) const;

protected:
  using uqBaseVectorCdfClass<V,M>::m_env;
  using uqBaseVectorCdfClass<V,M>::m_prefix;
  using uqBaseVectorCdfClass<V,M>::m_domainSpace;
  using uqBaseVectorCdfClass<V,M>::m_domainMinValues;
  using uqBaseVectorCdfClass<V,M>::m_domainMaxValues;

  const uqArrayOfOneDUniformGridsClass<V,M>& m_oneDGrids;
  const uqArrayOfScalarSetsClass      <V,M>& m_cdfValues;
};

template<class V,class M>
uqSampledVectorCdfClass<V,M>::uqSampledVectorCdfClass(
  const char*                                prefix,
  const uqArrayOfOneDUniformGridsClass<V,M>& oneDGrids,
  const uqArrayOfScalarSetsClass      <V,M>& cdfValues)
  :
  uqBaseVectorCdfClass<V,M>(((std::string)(prefix)+"sam").c_str(),oneDGrids.rowSpace(),oneDGrids.minValues(),oneDGrids.maxValues()),
  m_oneDGrids(oneDGrids),
  m_cdfValues(cdfValues)
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqSampledVectorCdfClass<V,M>::constructor()"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqSampledVectorCdfClass<V,M>::constructor()"
              << ": prefix = " << m_prefix
              << std::endl;
  }
}

template<class V,class M>
uqSampledVectorCdfClass<V,M>::~uqSampledVectorCdfClass()
{
}

template<class V, class M>
void
uqSampledVectorCdfClass<V,M>::values(
  const V& paramValues,
        V& cdfVec) const
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.rank(),
                      "uqSampledVectorCdfClass<V,M>::cdfVec()",
                      "incomplete code");
  return;
}

template <class V, class M>
void
uqSampledVectorCdfClass<V,M>::printContents(std::ofstream& ofs) const
{
  for (unsigned int i = 0; i < m_domainSpace.dim(); ++i) {
    // Print values *of* grid points
    std::vector<double> aGrid;
    m_oneDGrids.getAGrid(i,aGrid);
    ofs << m_prefix << i << "_grid = zeros(" << aGrid.size()
        << ","                               << 1
        << ");"
        << std::endl;
    ofs << m_prefix << i << "_grid = [";
    for (unsigned int j = 0; j < aGrid.size(); ++j) {
      ofs << aGrid[j] << " ";
    }
    ofs << "];"
        << std::endl;

    // Print *cdf* values *at* grid points
    const std::vector<double>& aSetOfCdfValues = m_cdfValues.scalarSet(i);
    ofs << m_prefix << i << "_values = zeros(" << aSetOfCdfValues.size()
        << ","                                 << 1
        << ");"
        << std::endl;
    ofs << m_prefix << i << "_values = [";
    for (unsigned int j = 0; j < aSetOfCdfValues.size(); ++j) {
      ofs << aSetOfCdfValues[j] << " ";
    }
    ofs << "];"
        << std::endl;
  }

  return;
}

#endif // __UQ_VECTOR_CUMULATIVE_DISTRIBUTION_FUNCTION_H__
