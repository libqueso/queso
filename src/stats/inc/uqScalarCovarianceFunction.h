//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, 
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-
// 
// $Id$
//
//--------------------------------------------------------------------------

#ifndef __UQ_SCALAR_COVARIANCE_FUNCTION_H__
#define __UQ_SCALAR_COVARIANCE_FUNCTION_H__

#include <uqVectorSet.h>
#include <uqEnvironment.h>
#include <cmath>

template<class V,class M>
class uqBaseScalarCovarianceFunctionClass {
public:
           uqBaseScalarCovarianceFunctionClass(const char*                  prefix,
                                               const uqVectorSetClass<V,M>& basicDomainSet);
  virtual ~uqBaseScalarCovarianceFunctionClass();

          const uqVectorSetClass<V,M>& basicDomainSet   ()                                               const;
  virtual       double                 value            (const V& domainVector1, const V& domainVector2) const = 0;

protected:
  const uqBaseEnvironmentClass& m_env;
        std::string             m_prefix;
  const uqVectorSetClass<V,M>&  m_basicDomainSet;
};

template<class V, class M>
uqBaseScalarCovarianceFunctionClass<V,M>::uqBaseScalarCovarianceFunctionClass(
  const char*                  prefix,
  const uqVectorSetClass<V,M>& basicDomainSet)
  :
  m_env           (basicDomainSet.env()),
  m_prefix        ((std::string)(prefix)+"cov_func_"),
  m_basicDomainSet(basicDomainSet)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqBaseScalarCovarianceFunctionClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqBaseScalarCovarianceFunctionClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}

template<class V, class M>
uqBaseScalarCovarianceFunctionClass<V,M>::~uqBaseScalarCovarianceFunctionClass()
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqBaseScalarCovarianceFunctionClass<V,M>::destructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqBaseScalarCovarianceFunctionClass<V,M>::destructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}

template<class V, class M>
const uqVectorSetClass<V,M>&
uqBaseScalarCovarianceFunctionClass<V,M>::basicDomainSet() const
{
  return m_basicDomainSet;
}

//*****************************************************
// Exponential class
//*****************************************************
template<class V,class M>
class uqExponentialScalarCovarianceFunctionClass : public uqBaseScalarCovarianceFunctionClass<V,M> {
public:
           uqExponentialScalarCovarianceFunctionClass(const char*                  prefix,
                                                      const uqVectorSetClass<V,M>& basicDomainSet,
                                                      double                       sigma,
                                                      double                       a);
  virtual ~uqExponentialScalarCovarianceFunctionClass();

           double value(const V& domainVector1, const V& domainVector2) const;

protected:
  using uqBaseScalarCovarianceFunctionClass<V,M>::m_env;
  using uqBaseScalarCovarianceFunctionClass<V,M>::m_prefix;
  using uqBaseScalarCovarianceFunctionClass<V,M>::m_basicDomainSet;

  double m_sigma;
  double m_a;
};

template<class V,class M>
uqExponentialScalarCovarianceFunctionClass<V,M>::uqExponentialScalarCovarianceFunctionClass(
  const char*                  prefix,
  const uqVectorSetClass<V,M>& basicDomainSet,
  double                       sigma,
  double                       a)
  : 
  uqBaseScalarCovarianceFunctionClass<V,M>(prefix,basicDomainSet),
  m_sigma(sigma),
  m_a    (a)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqExponentialScalarCovarianceFunctionClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqExponentialScalarCovarianceFunctionClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}

template<class V,class M>
uqExponentialScalarCovarianceFunctionClass<V,M>::~uqExponentialScalarCovarianceFunctionClass()
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqExponentialScalarCovarianceFunctionClass<V,M>::destructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqExponentialScalarCovarianceFunctionClass<V,M>::destructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}

template<class V,class M>
double
uqExponentialScalarCovarianceFunctionClass<V,M>::value(const V& domainVector1, const V& domainVector2) const
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqExponentialScalarCovarianceFunctionClass<V,M>::value()"
                          << std::endl;
  }

  double result = 0.;

  double exponent = -(domainVector1 - domainVector2).norm2Sq()/(m_sigma*m_sigma);

  result = m_a*std::exp(exponent);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqExponentialScalarCovarianceFunctionClass<V,M>::value()"
                          << std::endl;
  }

  return result;
}

//*****************************************************
// Generic class
//*****************************************************
template<class V,class M>
class uqGenericScalarCovarianceFunctionClass : public uqBaseScalarCovarianceFunctionClass<V,M> {
public:
           uqGenericScalarCovarianceFunctionClass(const char*                  prefix,
                                                  const uqVectorSetClass<V,M>& domainSet,
                                                  double (*covRoutinePtr)(const V& positionVector1, const V& positionVector2, const void* routineDataPtr),
                                                  const void*                  routinesDataPtr);
  virtual ~uqGenericScalarCovarianceFunctionClass();

           double value(const V& positionVector1, const V& positionVector2) const;

protected:
  using uqBaseScalarCovarianceFunctionClass<V,M>::m_env;
  using uqBaseScalarCovarianceFunctionClass<V,M>::m_prefix;
  using uqBaseScalarCovarianceFunctionClass<V,M>::m_basicDomainSet;

  double (*m_covRoutinePtr)(const V& positionVector1, const V& positionVector2, const void* routineDataPtr);
  const void* m_routineDataPtr;
};

template<class V,class M>
uqGenericScalarCovarianceFunctionClass<V,M>::uqGenericScalarCovarianceFunctionClass(
  const char*                  prefix,
  const uqVectorSetClass<V,M>& domainSet,
  double (*covRoutinePtr)(const V& positionVector1, const V& positionVector2, const void* routineDataPtr),
  const void*                  routinesDataPtr)
  : 
  uqBaseScalarCovarianceFunctionClass<V,M>(prefix,domainSet),
  m_covRoutinePtr                         (covRoutinePtr),
  m_routineDataPtr                        (routinesDataPtr)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqGenericScalarCovarianceFunctionClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqGenericScalarCovarianceFunctionClass<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}

template<class V,class M>
uqGenericScalarCovarianceFunctionClass<V,M>::~uqGenericScalarCovarianceFunctionClass()
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqGenericScalarCovarianceFunctionClass<V,M>::destructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqGenericScalarCovarianceFunctionClass<V,M>::destructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}

template<class V,class M>
double
uqGenericScalarCovarianceFunctionClass<V,M>::value(const V& positionVector1, const V& positionVector2) const
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqGenericScalarCovarianceFunctionClass<V,M>::value()"
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO(m_covRoutinePtr == NULL,
                      m_env.worldRank(),
                      "uqGenericScalarCovarianceFunctionClass<V,M>::value()",
                      "m_covRoutinePtr = NULL");

  double result = m_covRoutinePtr(positionVector1, positionVector2, m_routineDataPtr);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqGenericScalarCovarianceFunctionClass<V,M>::value()"
                            << std::endl;
  }

  return result;
}

#endif // __UQ_SCALAR_COVARIANCE_FUNCTION_H__ 
