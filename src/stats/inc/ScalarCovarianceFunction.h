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

#ifndef UQ_SCALAR_COVARIANCE_FUNCTION_H
#define UQ_SCALAR_COVARIANCE_FUNCTION_H

#include <queso/VectorSet.h>
#include <queso/Environment.h>
#include <cmath>

namespace QUESO {

class GslVector;
class GslMatrix;

/*! \file ScalarCovarianceFunction.h
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

template <class V = GslVector, class M = GslMatrix>
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

}  // End namespace QUESO

#endif // UQ_SCALAR_COVARIANCE_FUNCTION_H
