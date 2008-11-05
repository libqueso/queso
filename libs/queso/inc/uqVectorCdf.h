/* uq/libs/queso/inc/uqVectorCdf.h
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

#ifndef __UQ_VECTOR_CUMULATIVE_DISTRIBUTION_FUNCTION_H__
#define __UQ_VECTOR_CUMULATIVE_DISTRIBUTION_FUNCTION_H__

#include <uqArrayOfOneDGrids.h>
#include <uqArrayOfOneDTables.h>
#include <uqScalarCdf.h>
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

  const   uqVectorSpaceClass<V,M>&            domainSpace    ()                                const;
  virtual void                                values         (const V& paramValues, V& cdfVec) const = 0;
  virtual const uqBaseScalarCdfClass<double>& cdf            (unsigned int rowId)              const = 0;
  const   V&                                  domainMinValues()                                const;
  const   V&                                  domainMaxValues()                                const;
  virtual void                                print          (std::ostream& os)                const = 0;

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
  m_prefix         ((std::string)(prefix)+"cdf_"),
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

template<class V, class M>
const uqVectorSpaceClass<V,M>&
uqBaseVectorCdfClass<V,M>::domainSpace() const
{
  return m_domainSpace;
}

template <class V, class M>
  std::ostream& operator<< (std::ostream& os, const uqBaseVectorCdfClass<V,M>& obj)
{
  obj.print(os);
  return os;
}

//*****************************************************
// Generic cumulative distibution function class
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

  void values(const V& paramValues, V& cdfVec) const;
  void print (std::ostream& os) const;

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
  uqBaseVectorCdfClass<V,M>(prefix,domainSpace,domainMinValues,domainMaxValues),
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
  uqBaseVectorCdfClass<V,M>(prefix,domainSpace),
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
uqGenericVectorCdfClass<V,M>::print(std::ostream& os) const
{
  return;
}

//*****************************************************
// Gaussian cumulative distribution function class
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

  void values(const V& paramValues, V& cdfVec) const;
  void print (std::ostream& os) const;

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
  uqBaseVectorCdfClass<V,M>(prefix,domainSpace,domainMinValues,domainMaxValues),
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
  uqBaseVectorCdfClass<V,M>(prefix,domainSpace,domainMinValues,domainMaxValues),
  m_covMatrix              (new M(covMatrix))
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
uqGaussianVectorCdfClass<V,M>::print(std::ostream& os) const
{
  return;
}

//*****************************************************
// Sampled cumulative distribution function class
//*****************************************************
template<class V, class M>
class uqSampledVectorCdfClass : public uqBaseVectorCdfClass<V,M> {
public:
  uqSampledVectorCdfClass(const char*                          prefix,
                          const uqArrayOfOneDGridsClass <V,M>& oneDGrids,
                          const uqArrayOfOneDTablesClass<V,M>& cdfValues);
 ~uqSampledVectorCdfClass();

        void                          values(const V& paramValues, V& cdfVec) const;
  const uqBaseScalarCdfClass<double>& cdf   (unsigned int rowId)              const;
        void                          print (std::ostream& os)                const;

protected:
  using uqBaseVectorCdfClass<V,M>::m_env;
  using uqBaseVectorCdfClass<V,M>::m_prefix;
  using uqBaseVectorCdfClass<V,M>::m_domainSpace;
  using uqBaseVectorCdfClass<V,M>::m_domainMinValues;
  using uqBaseVectorCdfClass<V,M>::m_domainMaxValues;

  EpetraExt::DistArray<uqSampledScalarCdfClass<double>*> m_cdfs;
};

template<class V,class M>
uqSampledVectorCdfClass<V,M>::uqSampledVectorCdfClass(
  const char*                          prefix,
  const uqArrayOfOneDGridsClass <V,M>& oneDGrids,
  const uqArrayOfOneDTablesClass<V,M>& cdfValues)
  :
  uqBaseVectorCdfClass<V,M>(prefix,oneDGrids.rowSpace(),oneDGrids.minPositions(),oneDGrids.maxPositions()),
  m_cdfs(m_domainSpace.map(),1)
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqSampledVectorCdfClass<V,M>::constructor()"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  char strI[65];
  for (unsigned int i = 0; i < (unsigned int) m_cdfs.MyLength(); ++i) {
    sprintf(strI,"%d_",i);
    m_cdfs(i,0) = new uqSampledScalarCdfClass<double>(m_env,
                                                      ((std::string)(m_prefix)+strI).c_str(),
                                                      oneDGrids.grid(i),
                                                      cdfValues.oneDTable(i));
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
  for (unsigned int i = 0; i < (unsigned int) m_cdfs.MyLength(); ++i) {
    if (m_cdfs(i,0)) delete m_cdfs(i,0);
  }
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

template<class V, class M>
const uqBaseScalarCdfClass<double>&
uqSampledVectorCdfClass<V,M>::cdf(unsigned int rowId) const
{
  UQ_FATAL_TEST_MACRO(rowId >= m_domainSpace.dim(),
                      m_env.rank(),
                      "uqSampledVectorCdfClass<T>::cdf()",
                      "rowId is out of range");

  uqSampledVectorCdfClass<V,M>* tmp = const_cast<uqSampledVectorCdfClass<V,M>*>(this);
  return *(tmp->m_cdfs(rowId,0));
  
}


template <class V, class M>
void
uqSampledVectorCdfClass<V,M>::print(std::ostream& os) const
{
  uqSampledVectorCdfClass<V,M>* tmp = const_cast<uqSampledVectorCdfClass<V,M>*>(this);
  for (unsigned int i = 0; i < (unsigned int) m_cdfs.MyLength(); ++i) {
    os << (*tmp->m_cdfs(i,0))
       << std::endl;
  }

  return;
}

template <class V, class M>
void
horizontalDistances(const uqBaseVectorCdfClass<V,M>& cdf1,
                    const uqBaseVectorCdfClass<V,M>& cdf2,
                    const V& epsilonVec,
                    V&       distances)
{
  for (unsigned int i = 0; i < cdf1.domainSpace().dim(); ++i) {
    distances[i] = horizontalDistance(cdf1.cdf(i),
                                      cdf2.cdf(i),
                                      epsilonVec[i]);
  }

  return;
}
#endif // __UQ_VECTOR_CUMULATIVE_DISTRIBUTION_FUNCTION_H__
