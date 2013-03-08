//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
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

#ifndef __UQ_MATRIX_COVARIANCE_FUNCTION_H__
#define __UQ_MATRIX_COVARIANCE_FUNCTION_H__

#include <uqVectorSet.h>
#include <uqEnvironment.h>
#include <cmath>

template<class P_V, class P_M, class Q_V, class Q_M>
class uqBaseMatrixCovarianceFunctionClass {
public:
           uqBaseMatrixCovarianceFunctionClass(const char*                      prefix,
                                               const uqVectorSetClass<P_V,P_M>& basicDomainSet,
                                               const uqVectorSetClass<Q_V,Q_M>& imageSet);
  virtual ~uqBaseMatrixCovarianceFunctionClass();

          const uqVectorSetClass<P_V,P_M>& basicDomainSet()                                                                     const;
  virtual       void                       covMatrix     (const P_V& domainVector1, const P_V& domainVector2, Q_M& imageMatrix) const = 0;

protected:
  const uqBaseEnvironmentClass&    m_env;
        std::string                m_prefix;
  const uqVectorSetClass<P_V,P_M>& m_basicDomainSet;
  const uqVectorSetClass<Q_V,Q_M>& m_imageSet;
};

template<class P_V, class P_M, class Q_V, class Q_M>
uqBaseMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::uqBaseMatrixCovarianceFunctionClass(
  const char*                      prefix,
  const uqVectorSetClass<P_V,P_M>& basicDomainSet,
  const uqVectorSetClass<Q_V,Q_M>& imageSet)
  :
  m_env           (basicDomainSet.env()),
  m_prefix        ((std::string)(prefix)+"cov_func_"),
  m_basicDomainSet(basicDomainSet),
  m_imageSet      (imageSet)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqBaseMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqBaseMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}

template<class P_V, class P_M, class Q_V, class Q_M>
uqBaseMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::~uqBaseMatrixCovarianceFunctionClass()
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqBaseMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::destructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqBaseMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::destructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}

template<class P_V, class P_M, class Q_V, class Q_M>
const uqVectorSetClass<P_V,P_M>&
uqBaseMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::basicDomainSet() const
{
  return m_basicDomainSet;
}

//*****************************************************
// Exponential class
//*****************************************************
template<class P_V, class P_M, class Q_V, class Q_M>
class uqExponentialMatrixCovarianceFunctionClass : public uqBaseMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M> {
public:
           uqExponentialMatrixCovarianceFunctionClass(const char*                      prefix,
                                                      const uqVectorSetClass<P_V,P_M>& basicDomainSet,
                                                      const uqVectorSetClass<Q_V,Q_M>& imageSet,
                                                      const Q_M&                       sigmas,
                                                      const Q_M&                       as);
  virtual ~uqExponentialMatrixCovarianceFunctionClass();

           void covMatrix(const P_V& domainVector1, const P_V& domainVector2, Q_M& imageMatrix) const;

protected:
  using uqBaseMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::m_env;
  using uqBaseMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::m_prefix;
  using uqBaseMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::m_basicDomainSet;
  using uqBaseMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::m_imageSet;

  Q_M* m_sigmas;
  Q_M* m_as;
};

template<class P_V, class P_M, class Q_V, class Q_M>
uqExponentialMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::uqExponentialMatrixCovarianceFunctionClass(
  const char*                      prefix,
  const uqVectorSetClass<P_V,P_M>& basicDomainSet,
  const uqVectorSetClass<Q_V,Q_M>& imageSet,
  const Q_M&                       sigmas,
  const Q_M&                       as)
  : 
  uqBaseMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>(prefix,basicDomainSet,imageSet),
  m_sigmas(NULL),
  m_as    (NULL)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqExponentialMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  m_sigmas = new Q_M(sigmas);
  m_as     = new Q_M(as);

  unsigned int matrixOrder = m_imageSet.vectorSpace().dimLocal();

  UQ_FATAL_TEST_MACRO(m_sigmas->numRowsLocal() != matrixOrder,
                      m_env.worldRank(),
                      "uqExponentialMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::Constructor()",
                      "m_sigmas has invalid number of rows");

  UQ_FATAL_TEST_MACRO(m_sigmas->numCols() != matrixOrder,
                      m_env.worldRank(),
                      "uqExponentialMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::Constructor()",
                      "m_sigmas has invalid number of columns");

  UQ_FATAL_TEST_MACRO(m_as->numRowsLocal() != matrixOrder,
                      m_env.worldRank(),
                      "uqExponentialMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::Constructor()",
                      "m_as has invalid number of rows");

  UQ_FATAL_TEST_MACRO(m_as->numCols() != matrixOrder,
                      m_env.worldRank(),
                      "uqExponentialMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::Constructor()",
                      "m_as has invalid number of columns");

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqExponentialMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}

template<class P_V, class P_M, class Q_V, class Q_M>
uqExponentialMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::~uqExponentialMatrixCovarianceFunctionClass()
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqExponentialMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::destructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  delete m_as;
  delete m_sigmas;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqExponentialMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::destructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}

template<class P_V, class P_M, class Q_V, class Q_M>
void
uqExponentialMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::covMatrix(const P_V& domainVector1, const P_V& domainVector2, Q_M& imageMatrix) const
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqExponentialMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::covMatrix()"
                          << std::endl;
  }

  unsigned int matrixOrder = m_imageSet.vectorSpace().dimLocal();

  UQ_FATAL_TEST_MACRO(imageMatrix.numRowsLocal() != matrixOrder,
                      m_env.worldRank(),
                      "uqExponentialMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::covMatrix()",
                      "imageMatrix has invalid number of rows");

  UQ_FATAL_TEST_MACRO(imageMatrix.numCols() != matrixOrder,
                      m_env.worldRank(),
                      "uqExponentialMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::covMatrix()",
                      "imageMatrix has invalid number of columns");

  double tmpSq = -(domainVector1 - domainVector2).norm2Sq();

  for (unsigned int i = 0; i < matrixOrder; ++i) {
    for (unsigned int j = 0; j < matrixOrder; ++j) {
      double tmp = tmpSq/( (*m_sigmas)(i,j) * (*m_sigmas)(i,j) );
      imageMatrix(i,j) = (*m_as)(i,j) * std::exp(tmp);
    }
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqExponentialMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::covMatrix()"
                          << std::endl;
  }

  return;
}

//*****************************************************
// Generic class
//*****************************************************
template<class P_V, class P_M, class Q_V, class Q_M>
class uqGenericMatrixCovarianceFunctionClass : public uqBaseMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M> {
public:
           uqGenericMatrixCovarianceFunctionClass(const char*                      prefix,
                                                  const uqVectorSetClass<P_V,P_M>& basicDomainSet,
                                                  const uqVectorSetClass<Q_V,Q_M>& imageSet,
                                                  void (*covRoutinePtr)(const P_V& positionVector1, const P_V& positionVector2, const void* routineDataPtr, Q_M& imageMatrix),
                                                  const void*                      routinesDataPtr);
  virtual ~uqGenericMatrixCovarianceFunctionClass();

           void covMatrix(const P_V& positionVector1, const P_V& positionVector2, Q_M& imageMatrix) const;

protected:
  using uqBaseMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::m_env;
  using uqBaseMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::m_prefix;
  using uqBaseMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::m_basicDomainSet;
  using uqBaseMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::m_imageSet;

  void (*m_covRoutinePtr)(const P_V& positionVector1, const P_V& positionVector2, const void* routineDataPtr, Q_M& imageMatrix);
  const void* m_routineDataPtr;
};

template<class P_V, class P_M, class Q_V, class Q_M>
uqGenericMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::uqGenericMatrixCovarianceFunctionClass(
  const char*                      prefix,
  const uqVectorSetClass<P_V,P_M>& basicDomainSet,
  const uqVectorSetClass<Q_V,Q_M>& imageSet,
  void (*covRoutinePtr)(const P_V& positionVector1, const P_V& positionVector2, const void* routineDataPtr, Q_M& imageMatrix),
  const void*                      routinesDataPtr)
  : 
  uqBaseMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>(prefix,basicDomainSet,imageSet),
  m_covRoutinePtr                                     (covRoutinePtr),
  m_routineDataPtr                                    (routinesDataPtr)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqGenericMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqGenericMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}

template<class P_V, class P_M, class Q_V, class Q_M>
uqGenericMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::~uqGenericMatrixCovarianceFunctionClass()
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqGenericMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::destructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqGenericMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::destructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}

template<class P_V, class P_M, class Q_V, class Q_M>
void
uqGenericMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::covMatrix(const P_V& positionVector1, const P_V& positionVector2, Q_M& imageMatrix) const
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqGenericMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::covMatrix()"
                            << std::endl;
  }

  unsigned int matrixOrder = m_imageSet.vectorSpace().dimLocal();

  UQ_FATAL_TEST_MACRO(imageMatrix.numRowsLocal() != matrixOrder,
                      m_env.worldRank(),
                      "uqGenericMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::covMatrix()",
                      "imageMatrix has invalid number of rows");

  UQ_FATAL_TEST_MACRO(imageMatrix.numCols() != matrixOrder,
                      m_env.worldRank(),
                      "uqGenericMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::covMatrix()",
                      "imageMatrix has invalid number of columns");

  UQ_FATAL_TEST_MACRO(m_covRoutinePtr == NULL,
                      m_env.worldRank(),
                      "uqGenericMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::covMatrix()",
                      "m_covRoutinePtr = NULL");

  m_covRoutinePtr(positionVector1, positionVector2, m_routineDataPtr, imageMatrix);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqGenericMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::covMatrix()"
                            << std::endl;
  }

  return;
}

#endif // __UQ_MATRIX_COVARIANCE_FUNCTION_H__ 
