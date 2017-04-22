//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2017 The PECOS Development Team
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

#ifndef UQ_EXP_SCALAR_COVARIANCE_FUNCTION_H
#define UQ_EXP_SCALAR_COVARIANCE_FUNCTION_H

#include <queso/ScalarCovarianceFunction.h>
#include <queso/VectorSet.h>
#include <queso/Environment.h>
#include <cmath>

namespace QUESO {

class GslVector;
class GslMatrix;

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

template <class V = GslVector, class M = GslMatrix>
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

}  // End namespace QUESO

#endif // UQ_EXP_SCALAR_COVARIANCE_FUNCTION_H
