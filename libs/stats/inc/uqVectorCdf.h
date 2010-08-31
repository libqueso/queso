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
           uqBaseVectorCdfClass(const char*                  prefix,
                                const uqVectorSetClass<V,M>& pdfSupport);
  virtual ~uqBaseVectorCdfClass();

          const uqVectorSetClass<V,M>&        pdfSupport      ()                                const;
  virtual void                                values          (const V& paramValues, V& cdfVec) const = 0;
  virtual const uqBaseScalarCdfClass<double>& cdf             (unsigned int rowId)              const = 0;
  virtual void                                print           (std::ostream& os)                const = 0;
  virtual void                                subWriteContents(const std::string&            varNamePrefix,
                                                               const std::string&            fileName,
                                                               const std::string&            fileType,
                                                               const std::set<unsigned int>& allowedSubEnvIds) const;

protected:

  const   uqBaseEnvironmentClass& m_env;
          std::string             m_prefix;
  const   uqVectorSetClass<V,M>&  m_pdfSupport;
};

template<class V, class M>
uqBaseVectorCdfClass<V,M>::uqBaseVectorCdfClass(
  const char*                  prefix,
  const uqVectorSetClass<V,M>& pdfSupport)
  :
  m_env       (pdfSupport.env()),
  m_prefix    ((std::string)(prefix)+"Cdf_"),
  m_pdfSupport(pdfSupport)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqBaseVectorCdfClass<V,M>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqBaseVectorCdfClass<V,M>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }
}

template<class V, class M>
uqBaseVectorCdfClass<V,M>::~uqBaseVectorCdfClass()
{
}

template<class V, class M>
const uqVectorSetClass<V,M>&
uqBaseVectorCdfClass<V,M>::pdfSupport() const
{
  return m_pdfSupport;
}

template<class V, class M>
void
uqBaseVectorCdfClass<V,M>::subWriteContents(
  const std::string&            varNamePrefix,
  const std::string&            fileName,
  const std::string&            fileType,
  const std::set<unsigned int>& allowedSubEnvIds) const
{
  std::cerr << "WARNING: uqBaseVectorCdfClass<V,M>::subWriteContents() being used..."
            << std::endl;
  return;
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
  uqGenericVectorCdfClass(const char*                  prefix,
                          const uqVectorSetClass<V,M>& pdfSupport,
                          double (*routinePtr)(const V& paramValues, const void* routineDataPtr, V& cdfVec),
                          const void*                  routineDataPtr);
 ~uqGenericVectorCdfClass();

  void values(const V& paramValues, V& cdfVec) const;
  void print (std::ostream& os) const;

protected:
  double (*m_routinePtr)(const V& paramValues, const void* routineDataPtr, V& cdfVec);
  const void* m_routineDataPtr;

  using uqBaseVectorCdfClass<V,M>::m_env;
  using uqBaseVectorCdfClass<V,M>::m_prefix;
  using uqBaseVectorCdfClass<V,M>::m_pdfSupport;
};

template<class V, class M>
uqGenericVectorCdfClass<V,M>::uqGenericVectorCdfClass(
  const char*                    prefix,
  const uqVectorSetClass<V,M>& pdfSupport,
  double (*routinePtr)(const V& paramValues, const void* routineDataPtr, V& cdfVec),
  const void* routineDataPtr)
  :
  uqBaseVectorCdfClass<V,M>(prefix,pdfSupport),
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
  uqGaussianVectorCdfClass(const char*                  prefix,
                           const uqVectorSetClass<V,M>& pdfSupport,
                           const V&                     domainExpectedValues,
                           const V&                     domainVarianceValues);
  uqGaussianVectorCdfClass(const char*                  prefix,
                           const uqVectorSetClass<V,M>& pdfSupport,
                           const V&                     domainExpectedValues,
                           const M&                     covMatrix);
 ~uqGaussianVectorCdfClass();

  void values(const V& paramValues, V& cdfVec) const;
  void print (std::ostream& os) const;

protected:
  const M*                         m_covMatrix;

  using uqBaseVectorCdfClass<V,M>::m_env;
  using uqBaseVectorCdfClass<V,M>::m_prefix;
  using uqBaseVectorCdfClass<V,M>::m_pdfSupport;

  void commonConstructor();
};

template<class V,class M>
uqGaussianVectorCdfClass<V,M>::uqGaussianVectorCdfClass(
  const char*                    prefix,
  const uqVectorSetClass<V,M>& pdfSupport,
  const V&                       domainExpectedValues,
  const V&                       domainVarianceValues)
  :
  uqBaseVectorCdfClass<V,M>(prefix,pdfSupport),
  m_covMatrix              (m_pdfSupport.newDiagMatrix(domainVarianceValues*domainVarianceValues))
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqGaussianVectorCdfClass<V,M>::constructor() [1]"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  commonConstructor();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqGaussianVectorCdfClass<V,M>::constructor() [1]"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }
}

template<class V,class M>
uqGaussianVectorCdfClass<V,M>::uqGaussianVectorCdfClass(
  const char*                    prefix,
  const uqVectorSetClass<V,M>& pdfSupport,
  const V&                       domainExpectedValues,
  const M&                       covMatrix)
  :
  uqBaseVectorCdfClass<V,M>(prefix,pdfSupport),
  m_covMatrix              (new M(covMatrix))
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqGaussianVectorCdfClass<V,M>::constructor() [2]"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  commonConstructor();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqGaussianVectorCdfClass<V,M>::constructor() [2]"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }
}

template<class V,class M>
void
uqGaussianVectorCdfClass<V,M>::commonConstructor()
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.fullRank(),
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
                      m_env.fullRank(),
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
        void                          subWriteContents(const std::string&            varNamePrefix,
                                                       const std::string&            fileName,
                                                       const std::string&            fileType,
                                                       const std::set<unsigned int>& allowedSubEnvIds) const;

protected:
  using uqBaseVectorCdfClass<V,M>::m_env;
  using uqBaseVectorCdfClass<V,M>::m_prefix;
  using uqBaseVectorCdfClass<V,M>::m_pdfSupport;

  EpetraExt::DistArray<uqSampledScalarCdfClass<double>*> m_cdfs;
};

template<class V,class M>
uqSampledVectorCdfClass<V,M>::uqSampledVectorCdfClass(
  const char*                          prefix,
  const uqArrayOfOneDGridsClass <V,M>& oneDGrids,
  const uqArrayOfOneDTablesClass<V,M>& cdfValues)
  :
  uqBaseVectorCdfClass<V,M>(prefix,oneDGrids.rowSpace()),
  m_cdfs(m_pdfSupport.vectorSpace().map(),1)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqSampledVectorCdfClass<V,M>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  char strI[65];
  for (unsigned int i = 0; i < (unsigned int) m_cdfs.MyLength(); ++i) {
    sprintf(strI,"%u_",i);
    m_cdfs(i,0) = new uqSampledScalarCdfClass<double>(m_env,
                                                      ((std::string)(m_prefix)+strI).c_str(),
                                                      oneDGrids.grid(i),
                                                      cdfValues.oneDTable(i));
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqSampledVectorCdfClass<V,M>::constructor()"
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
                      m_env.fullRank(),
                      "uqSampledVectorCdfClass<V,M>::cdfVec()",
                      "incomplete code");
  return;
}

template<class V, class M>
const uqBaseScalarCdfClass<double>&
uqSampledVectorCdfClass<V,M>::cdf(unsigned int rowId) const
{
  UQ_FATAL_TEST_MACRO(rowId >= m_pdfSupport.vectorSpace().dimLocal(),
                      m_env.fullRank(),
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

template<class V, class M>
void
uqSampledVectorCdfClass<V,M>::subWriteContents(
  const std::string&            varNamePrefix,
  const std::string&            fileName,
  const std::string&            fileType,
  const std::set<unsigned int>& allowedSubEnvIds) const
{
  UQ_FATAL_TEST_MACRO(m_env.subRank() < 0,
                      m_env.fullRank(),
                      "uqSampledVectorCdfClass<V,M>::subWriteContents()",
                      "unexpected subRank");

  uqSampledVectorCdfClass<V,M>* tmp = const_cast<uqSampledVectorCdfClass<V,M>*>(this);
  char compId[16+1];
  for (unsigned int i = 0; i < (unsigned int) m_cdfs.MyLength(); ++i) {
    sprintf(compId,"%d",i);
    tmp->m_cdfs(i,0)->subWriteContents(varNamePrefix+compId,fileName,fileType,allowedSubEnvIds);
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
  for (unsigned int i = 0; i < cdf1.pdfSupport().vectorSpace().dimLocal(); ++i) {
    distances[i] = horizontalDistance(cdf1.cdf(i),
                                      cdf2.cdf(i),
                                      epsilonVec[i]);
  }

  return;
}
#endif // __UQ_VECTOR_CUMULATIVE_DISTRIBUTION_FUNCTION_H__
