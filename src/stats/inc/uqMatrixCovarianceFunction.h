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

namespace QUESO {

/*! \file uqMatrixCovarianceFunction.h
 * \brief Classes to accommodate covariance matrix of random vector functions.
 * 
 * \class BaseMatrixCovarianceFunctionClass
 * \brief A templated (base) class to accommodate covariance matrix of (random) vector functions.
 *
 * This class allows the mathematical definition of a multivariate covariance function, i.e.
 * a covariance matrix of random vector functions. 
 * Sometimes the covariance matrix of a multivariate random variable is not known but has 
 * to be estimated. Estimation of covariance matrices then deals with the question of how 
 * to approximate the actual covariance matrix on the basis of a sample from the multivariate 
 * distribution. */
 
 /* The covariance for two random variates \f$ X \f$ and \f$ Y \f$, 
 * each with sample size \f$ N \f$, is defined by the expectation value:
 * \f$ cov (X,Y) = <(X-\mu_X)(Y-\mu_Y)> = < X Y > - \mu_X \mu_Y \f$ 
 * where \f$ \mu_X \f$ and \f$ \mu_Y \f$ are the respective means, which can be written out 
 * explicitly as \f[ cov (X,Y) = \sum_{i=1}^{N} \frac{(x_i - \bar{x})(y_i - \bar{y})}{N}\f] */

template<class P_V, class P_M, class Q_V, class Q_M>
class BaseMatrixCovarianceFunctionClass {
public:
    //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor.
  /*! Instantiates an object of the class  given a prefix, the domain set and the image set.*/
  BaseMatrixCovarianceFunctionClass(const char*                      prefix,
				      const VectorSetClass<P_V,P_M>& basicDomainSet,
				      const VectorSetClass<Q_V,Q_M>& imageSet);
  
  //! Virtual destructor
  virtual ~BaseMatrixCovarianceFunctionClass();
  //@}
  
  //! @name Math methods
  //@{
  //! Domain set; access to private attribute m_basicDomainSet.
  const VectorSetClass<P_V,P_M>& basicDomainSet()             const;
  
  //! Calculates the covariance matrix. See template specialization.
  virtual void                     covMatrix     (const P_V& domainVector1, const P_V& domainVector2, Q_M& imageMatrix) const = 0;
  //@}
  
protected:
  const BaseEnvironmentClass&    m_env;
        std::string                m_prefix;
  const VectorSetClass<P_V,P_M>& m_basicDomainSet;
  const VectorSetClass<Q_V,Q_M>& m_imageSet;
};
// Default constructor -----------------------------
template<class P_V, class P_M, class Q_V, class Q_M>
BaseMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::BaseMatrixCovarianceFunctionClass(
  const char*                      prefix,
  const VectorSetClass<P_V,P_M>& basicDomainSet,
  const VectorSetClass<Q_V,Q_M>& imageSet)
  :
  m_env           (basicDomainSet.env()),
  m_prefix        ((std::string)(prefix)+"cov_func_"),
  m_basicDomainSet(basicDomainSet),
  m_imageSet      (imageSet)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering BaseMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving BaseMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor ---------------------------------------
template<class P_V, class P_M, class Q_V, class Q_M>
BaseMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::~BaseMatrixCovarianceFunctionClass()
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering BaseMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::destructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving BaseMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::destructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Math methods -------------------------------------
template<class P_V, class P_M, class Q_V, class Q_M>
const VectorSetClass<P_V,P_M>&
BaseMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::basicDomainSet() const
{
  return m_basicDomainSet;
}

//*****************************************************
// Exponential class
//*****************************************************
/*!
 * \class ExponentialMatrixCovarianceFunctionClass
 * \brief A class for exponential covariance matrices.
 *  
 * This class implements squared exponential covariance matrices of the form:
 * \f[ cov = a \exp{(-d^2/\sigma^2)}\f], where \f$ d=d(x,y) \f$ is the distance between two vectors, 
 * \f$ \sigma^2 \f$ is the variance matrix and \f$ a \f$ is the length scale ().*/
 
template<class P_V, class P_M, class Q_V, class Q_M>
class ExponentialMatrixCovarianceFunctionClass : public BaseMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M> {
public:
  //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor.
  /*! Instantiates an object of the class given a prefix, the domain and image sets, the variances scale factors.*/
  ExponentialMatrixCovarianceFunctionClass(const char*                      prefix,
					     const VectorSetClass<P_V,P_M>& basicDomainSet,
					     const VectorSetClass<Q_V,Q_M>& imageSet,
					     const Q_M&                       sigmas,
					     const Q_M&                       as);
  
  //! Virtual destructor
  virtual ~ExponentialMatrixCovarianceFunctionClass();
  //@}
  
  //!@ \name Math methods
  //@{
  //! Calculates the covariance matrix, given two parameter domains.
  void covMatrix(const P_V& domainVector1, const P_V& domainVector2, Q_M& imageMatrix) const;
  //@}
  
protected:
  using BaseMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::m_env;
  using BaseMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::m_prefix;
  using BaseMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::m_basicDomainSet;
  using BaseMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::m_imageSet;

  Q_M* m_sigmas;
  Q_M* m_as;
};
// Default constructor -----------------------------
template<class P_V, class P_M, class Q_V, class Q_M>
ExponentialMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::ExponentialMatrixCovarianceFunctionClass(
  const char*                      prefix,
  const VectorSetClass<P_V,P_M>& basicDomainSet,
  const VectorSetClass<Q_V,Q_M>& imageSet,
  const Q_M&                       sigmas,
  const Q_M&                       as)
  : 
  BaseMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>(prefix,basicDomainSet,imageSet),
  m_sigmas(NULL),
  m_as    (NULL)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering ExponentialMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  m_sigmas = new Q_M(sigmas);
  m_as     = new Q_M(as);

  unsigned int matrixOrder = m_imageSet.vectorSpace().dimLocal();

  UQ_FATAL_TEST_MACRO(m_sigmas->numRowsLocal() != matrixOrder,
                      m_env.worldRank(),
                      "ExponentialMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::Constructor()",
                      "m_sigmas has invalid number of rows");

  UQ_FATAL_TEST_MACRO(m_sigmas->numCols() != matrixOrder,
                      m_env.worldRank(),
                      "ExponentialMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::Constructor()",
                      "m_sigmas has invalid number of columns");

  UQ_FATAL_TEST_MACRO(m_as->numRowsLocal() != matrixOrder,
                      m_env.worldRank(),
                      "ExponentialMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::Constructor()",
                      "m_as has invalid number of rows");

  UQ_FATAL_TEST_MACRO(m_as->numCols() != matrixOrder,
                      m_env.worldRank(),
                      "ExponentialMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::Constructor()",
                      "m_as has invalid number of columns");

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving ExponentialMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor ---------------------------------------
template<class P_V, class P_M, class Q_V, class Q_M>
ExponentialMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::~ExponentialMatrixCovarianceFunctionClass()
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering ExponentialMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::destructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  delete m_as;
  delete m_sigmas;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving ExponentialMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::destructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Math methods -------------------------------------
template<class P_V, class P_M, class Q_V, class Q_M>
void
ExponentialMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::covMatrix(const P_V& domainVector1, const P_V& domainVector2, Q_M& imageMatrix) const
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering ExponentialMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::covMatrix()"
                          << std::endl;
  }

  unsigned int matrixOrder = m_imageSet.vectorSpace().dimLocal();

  UQ_FATAL_TEST_MACRO(imageMatrix.numRowsLocal() != matrixOrder,
                      m_env.worldRank(),
                      "ExponentialMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::covMatrix()",
                      "imageMatrix has invalid number of rows");

  UQ_FATAL_TEST_MACRO(imageMatrix.numCols() != matrixOrder,
                      m_env.worldRank(),
                      "ExponentialMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::covMatrix()",
                      "imageMatrix has invalid number of columns");

  double tmpSq = -(domainVector1 - domainVector2).norm2Sq();

  for (unsigned int i = 0; i < matrixOrder; ++i) {
    for (unsigned int j = 0; j < matrixOrder; ++j) {
      double tmp = tmpSq/( (*m_sigmas)(i,j) * (*m_sigmas)(i,j) );
      imageMatrix(i,j) = (*m_as)(i,j) * std::exp(tmp);
    }
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving ExponentialMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::covMatrix()"
                          << std::endl;
  }

  return;
}

//*****************************************************
// Generic class
//*****************************************************
/*! 
 * \class GenericMatrixCovarianceFunctionClass
 * \brief A class for generic covariance matrices.
 *  
 * This class implements a generic covariance matrices by calling a routine (via pointer).*/

template<class P_V, class P_M, class Q_V, class Q_M>
class GenericMatrixCovarianceFunctionClass : public BaseMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M> {
public:
  //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor.
  /*! Instantiates an object of the class given a prefix, the domain set, the pointer to the routine. */
  GenericMatrixCovarianceFunctionClass(const char*                      prefix,
					 const VectorSetClass<P_V,P_M>& basicDomainSet,
					 const VectorSetClass<Q_V,Q_M>& imageSet,
					 void (*covRoutinePtr)(const P_V& positionVector1, const P_V& positionVector2, const void* routineDataPtr, Q_M& imageMatrix),
					 const void*                      routinesDataPtr);
  //! Virtual destructor
  virtual ~GenericMatrixCovarianceFunctionClass();
  //@}
  
  //! @name Math methods
  //@{
  //! Calculates the value of the generic covariance matrix.
  /*! This function accesses the routine that calculates the covariance function. */
  void covMatrix(const P_V& positionVector1, const P_V& positionVector2, Q_M& imageMatrix) const;

protected:
  using BaseMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::m_env;
  using BaseMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::m_prefix;
  using BaseMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::m_basicDomainSet;
  using BaseMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::m_imageSet;

  void (*m_covRoutinePtr)(const P_V& positionVector1, const P_V& positionVector2, const void* routineDataPtr, Q_M& imageMatrix);
  const void* m_routineDataPtr;
};
// Default constructor -----------------------------

template<class P_V, class P_M, class Q_V, class Q_M>
GenericMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::GenericMatrixCovarianceFunctionClass(
  const char*                      prefix,
  const VectorSetClass<P_V,P_M>& basicDomainSet,
  const VectorSetClass<Q_V,Q_M>& imageSet,
  void (*covRoutinePtr)(const P_V& positionVector1, const P_V& positionVector2, const void* routineDataPtr, Q_M& imageMatrix),
  const void*                      routinesDataPtr)
  : 
  BaseMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>(prefix,basicDomainSet,imageSet),
  m_covRoutinePtr                                     (covRoutinePtr),
  m_routineDataPtr                                    (routinesDataPtr)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering GenericMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving GenericMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor ---------------------------------------
template<class P_V, class P_M, class Q_V, class Q_M>
GenericMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::~GenericMatrixCovarianceFunctionClass()
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering GenericMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::destructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving GenericMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::destructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Math methods -------------------------------------
template<class P_V, class P_M, class Q_V, class Q_M>
void
GenericMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::covMatrix(const P_V& positionVector1, const P_V& positionVector2, Q_M& imageMatrix) const
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering GenericMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::covMatrix()"
                            << std::endl;
  }

  unsigned int matrixOrder = m_imageSet.vectorSpace().dimLocal();

  UQ_FATAL_TEST_MACRO(imageMatrix.numRowsLocal() != matrixOrder,
                      m_env.worldRank(),
                      "GenericMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::covMatrix()",
                      "imageMatrix has invalid number of rows");

  UQ_FATAL_TEST_MACRO(imageMatrix.numCols() != matrixOrder,
                      m_env.worldRank(),
                      "GenericMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::covMatrix()",
                      "imageMatrix has invalid number of columns");

  UQ_FATAL_TEST_MACRO(m_covRoutinePtr == NULL,
                      m_env.worldRank(),
                      "GenericMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::covMatrix()",
                      "m_covRoutinePtr = NULL");

  m_covRoutinePtr(positionVector1, positionVector2, m_routineDataPtr, imageMatrix);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving GenericMatrixCovarianceFunctionClass<P_V,P_M,Q_V,Q_M>::covMatrix()"
                            << std::endl;
  }

  return;
}

}  // End namespace QUESO

#endif // __UQ_MATRIX_COVARIANCE_FUNCTION_H__ 
