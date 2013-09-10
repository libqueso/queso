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

#ifndef __UQ_SCALAR_COVARIANCE_FUNCTION_H__
#define __UQ_SCALAR_COVARIANCE_FUNCTION_H__

#include <queso/VectorSet.h>
#include <queso/Environment.h>
#include <cmath>

namespace QUESO {

//*****************************************************
// Base class
//*****************************************************
/*! \file uqScalarCovarianceFunction.h
 * \brief Classes to accommodate covariance of scalar functions (random variables).
 * 
 * \class BaseScalarCovarianceFunction
 * \brief A templated (base) class to accommodate scalar covariance functions (of random variables).
 *
 * This class allows the mathematical definition of the covariance function of a random variable
 * Covariance provides a measure of the strength of the correlation between two or more sets of 
 * random variates. 
 * The covariance for two random variates \f$ X \f$ and \f$ Y \f$, each with sample size 
 * \f$ N \f$, is defined by the expectation value:
 * \f$ cov (X,Y) = <(X-\mu_X)(Y-\mu_Y)> = < X Y > - \mu_X \mu_Y \f$ 
 * where \f$ \mu_X \f$ and \f$ \mu_Y \f$ are the respective means, which can be written out 
 * explicitly as \f[ cov (X,Y) = \sum_{i=1}^{N} \frac{(x_i - \bar{x})(y_i - \bar{y})}{N}\f] */

template<class V,class M>
class BaseScalarCovarianceFunction {
public:
  //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor.
  /*! Instantiates an object of the class  given a prefix and the domain set.*/
  BaseScalarCovarianceFunction(const char*                  prefix,
                                      const VectorSet<V,M>& basicDomainSet);
  
  //! Virtual destructor.
  virtual ~BaseScalarCovarianceFunction();
//@}
  
  //! @name Math methods
  //@{
  //! Domain set; access to private attribute m_basicDomainSet.
  const VectorSet<V,M>& basicDomainSet   ()       const;
  
  //! The value of the covariance function. See template specialization.
  virtual       double                 value    (const V& domainVector1, const V& domainVector2) const = 0;
  //@}

protected:
  const BaseEnvironment& m_env;
        std::string             m_prefix;
  const VectorSet<V,M>&  m_basicDomainSet;
};
// Default constructor -----------------------------
template<class V, class M>
BaseScalarCovarianceFunction<V,M>::BaseScalarCovarianceFunction(
  const char*                  prefix,
  const VectorSet<V,M>& basicDomainSet)
  :
  m_env           (basicDomainSet.env()),
  m_prefix        ((std::string)(prefix)+"cov_func_"),
  m_basicDomainSet(basicDomainSet)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering BaseScalarCovarianceFunction<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving BaseScalarCovarianceFunction<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor ---------------------------------------
template<class V, class M>
BaseScalarCovarianceFunction<V,M>::~BaseScalarCovarianceFunction()
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering BaseScalarCovarianceFunction<V,M>::destructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving BaseScalarCovarianceFunction<V,M>::destructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Math methods -------------------------------------
template<class V, class M>
const VectorSet<V,M>&
BaseScalarCovarianceFunction<V,M>::basicDomainSet() const
{
  return m_basicDomainSet;
}

//*****************************************************
// Exponential class
//*****************************************************
/*! 
 * \class ExponentialScalarCovarianceFunction
 * \brief A class for exponential covariances.
 *  
 * This class implements squared exponential covariance functions of the form:
 * \f[ cov = a \exp{(-d^2/\sigma^2)}\f], where \f$ d=d(x,y) \f$ is the distance between two points, 
 * \f$ \sigma^2 \f$ is the variance and \f$ a \f$ is the length scale. 
 * This is a stationary covariance function with smooth sample paths. Exponential covariance 
 * functions are largely employed in Gaussian processes.  */
 
template<class V,class M>
class ExponentialScalarCovarianceFunction : public BaseScalarCovarianceFunction<V,M> {
public:
  //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor.
  /*! Instantiates an object of the class given a prefix, the domain set, the variance and a scale factor. */
  ExponentialScalarCovarianceFunction(const char*                  prefix,
					     const VectorSet<V,M>& basicDomainSet,
					     double                       sigma,
					     double                       a);
  
  //! Virtual destructor.
  virtual ~ExponentialScalarCovarianceFunction();
  
  //! @name Math methods
  //@{
  //! Calculates the value of the exponential covariance function.
  /*! The value of the exponential covariance function is: \f$ cov= a exp (-d^2/sigma)\f$, with 
   * \f$ d= \sqrt{(domainVector1 - domainVector1)^2} \f$*/
  double value(const V& domainVector1, const V& domainVector2) const;
  //@}
  
protected:
  using BaseScalarCovarianceFunction<V,M>::m_env;
  using BaseScalarCovarianceFunction<V,M>::m_prefix;
  using BaseScalarCovarianceFunction<V,M>::m_basicDomainSet;

  double m_sigma;
  double m_a;
};
// Default constructor -----------------------------
template<class V,class M>
ExponentialScalarCovarianceFunction<V,M>::ExponentialScalarCovarianceFunction(
  const char*                  prefix,
  const VectorSet<V,M>& basicDomainSet,
  double                       sigma,
  double                       a)
  : 
  BaseScalarCovarianceFunction<V,M>(prefix,basicDomainSet),
  m_sigma(sigma),
  m_a    (a)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering ExponentialScalarCovarianceFunction<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving ExponentialScalarCovarianceFunction<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor ---------------------------------------
template<class V,class M>
ExponentialScalarCovarianceFunction<V,M>::~ExponentialScalarCovarianceFunction()
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering ExponentialScalarCovarianceFunction<V,M>::destructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving ExponentialScalarCovarianceFunction<V,M>::destructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Math methods -------------------------------------
template<class V,class M>
double
ExponentialScalarCovarianceFunction<V,M>::value(const V& domainVector1, const V& domainVector2) const
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering ExponentialScalarCovarianceFunction<V,M>::value()"
                          << std::endl;
  }

  double result = 0.;

  double exponent = -(domainVector1 - domainVector2).norm2Sq()/(m_sigma*m_sigma);

  result = m_a*std::exp(exponent);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving ExponentialScalarCovarianceFunction<V,M>::value()"
                          << std::endl;
  }

  return result;
}

//*****************************************************
// Generic class
//*****************************************************
/*! 
 * \class GenericScalarCovarianceFunction
 * \brief A class for generic covariances.
 *  
 * This class implements a generic covariance functions, by calling a routine (via pointer).*/

template<class V,class M>
class GenericScalarCovarianceFunction : public BaseScalarCovarianceFunction<V,M> {
public:
  //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor.
  /*! Instantiates an object of the class given a prefix, the domain set, the pointer to the routine. */
  GenericScalarCovarianceFunction(const char*                  prefix,
					 const VectorSet<V,M>& domainSet,
					 double (*covRoutinePtr)(const V& positionVector1, const V& positionVector2, const void* routineDataPtr),
					 const void*                  routinesDataPtr);
  
  //! Virtual destructor
  virtual ~GenericScalarCovarianceFunction();
  //@}
    //! @name Math methods
  //@{
  //! Calculates the value of the generic covariance function.
  /*! This function accesses the routine that calculates the covariance function. */
  double value(const V& positionVector1, const V& positionVector2) const;

protected:
  using BaseScalarCovarianceFunction<V,M>::m_env;
  using BaseScalarCovarianceFunction<V,M>::m_prefix;
  using BaseScalarCovarianceFunction<V,M>::m_basicDomainSet;

  double (*m_covRoutinePtr)(const V& positionVector1, const V& positionVector2, const void* routineDataPtr);
  const void* m_routineDataPtr;
};
// Default constructor -----------------------------
template<class V,class M>
GenericScalarCovarianceFunction<V,M>::GenericScalarCovarianceFunction(
  const char*                  prefix,
  const VectorSet<V,M>& domainSet,
  double (*covRoutinePtr)(const V& positionVector1, const V& positionVector2, const void* routineDataPtr),
  const void*                  routinesDataPtr)
  : 
  BaseScalarCovarianceFunction<V,M>(prefix,domainSet),
  m_covRoutinePtr                         (covRoutinePtr),
  m_routineDataPtr                        (routinesDataPtr)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering GenericScalarCovarianceFunction<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving GenericScalarCovarianceFunction<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor ---------------------------------------
template<class V,class M>
GenericScalarCovarianceFunction<V,M>::~GenericScalarCovarianceFunction()
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering GenericScalarCovarianceFunction<V,M>::destructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving GenericScalarCovarianceFunction<V,M>::destructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Math methods -------------------------------------
template<class V,class M>
double
GenericScalarCovarianceFunction<V,M>::value(const V& positionVector1, const V& positionVector2) const
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering GenericScalarCovarianceFunction<V,M>::value()"
                            << std::endl;
  }

  UQ_FATAL_TEST_MACRO(m_covRoutinePtr == NULL,
                      m_env.worldRank(),
                      "GenericScalarCovarianceFunction<V,M>::value()",
                      "m_covRoutinePtr = NULL");

  double result = m_covRoutinePtr(positionVector1, positionVector2, m_routineDataPtr);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving GenericScalarCovarianceFunction<V,M>::value()"
                            << std::endl;
  }

  return result;
}

}  // End namespace QUESO

#endif // __UQ_SCALAR_COVARIANCE_FUNCTION_H__ 
