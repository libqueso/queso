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

#ifndef UQ_MATRIX_COVARIANCE_FUNCTION_H
#define UQ_MATRIX_COVARIANCE_FUNCTION_H

#include <queso/VectorSet.h>
#include <queso/Environment.h>
#include <cmath>

namespace QUESO {

class GslVector;
class GslMatrix;

/*! \file MatrixCovarianceFunction.h
 * \brief Classes to accommodate covariance matrix of random vector functions.
 *
 * \class BaseMatrixCovarianceFunction
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

template <class P_V = GslVector, class P_M = GslMatrix, class Q_V = GslVector, class Q_M = GslMatrix>
class BaseMatrixCovarianceFunction {
public:
    //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  /*! Instantiates an object of the class  given a prefix, the domain set and the image set.*/
  BaseMatrixCovarianceFunction(const char*                      prefix,
				      const VectorSet<P_V,P_M>& basicDomainSet,
				      const VectorSet<Q_V,Q_M>& imageSet);

  //! Virtual destructor
  virtual ~BaseMatrixCovarianceFunction();
  //@}

  //! @name Math methods
  //@{
  //! Domain set; access to private attribute m_basicDomainSet.
  const VectorSet<P_V,P_M>& basicDomainSet()             const;

  //! Calculates the covariance matrix. See template specialization.
  virtual void                     covMatrix     (const P_V& domainVector1, const P_V& domainVector2, Q_M& imageMatrix) const = 0;
  //@}

protected:
  const BaseEnvironment&    m_env;
        std::string                m_prefix;
  const VectorSet<P_V,P_M>& m_basicDomainSet;
  const VectorSet<Q_V,Q_M>& m_imageSet;
};

}  // End namespace QUESO

#endif // UQ_MATRIX_COVARIANCE_FUNCTION_H
