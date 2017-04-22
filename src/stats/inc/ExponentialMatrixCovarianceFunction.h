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

#ifndef UQ_EXP_MATRIX_COVARIANCE_FUNCTION_H
#define UQ_EXP_MATRIX_COVARIANCE_FUNCTION_H

#include <queso/MatrixCovarianceFunction.h>
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
 * \class ExponentialMatrixCovarianceFunction
 * \brief A class for exponential covariance matrices.
 *
 * This class implements squared exponential covariance matrices of the form:
 * \f[ cov = a \exp{(-d^2/\sigma^2)}\f], where \f$ d=d(x,y) \f$ is the distance between two vectors,
 * \f$ \sigma^2 \f$ is the variance matrix and \f$ a \f$ is the length scale ().*/

template <class P_V = GslVector, class P_M = GslMatrix, class Q_V = GslVector, class Q_M = GslMatrix>
class ExponentialMatrixCovarianceFunction : public BaseMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M> {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  /*! Instantiates an object of the class given a prefix, the domain and image sets, the variances scale factors.*/
  ExponentialMatrixCovarianceFunction(const char*                      prefix,
               const VectorSet<P_V,P_M>& basicDomainSet,
               const VectorSet<Q_V,Q_M>& imageSet,
               const Q_M&                       sigmas,
               const Q_M&                       as);

  //! Virtual destructor
  virtual ~ExponentialMatrixCovarianceFunction();
  //@}

  //!@ \name Math methods
  //@{
  //! Calculates the covariance matrix, given two parameter domains.
  void covMatrix(const P_V& domainVector1, const P_V& domainVector2, Q_M& imageMatrix) const;
  //@}

protected:
  using BaseMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>::m_env;
  using BaseMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>::m_prefix;
  using BaseMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>::m_basicDomainSet;
  using BaseMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>::m_imageSet;

  Q_M* m_sigmas;
  Q_M* m_as;
};

}  // End namespace QUESO

#endif // UQ_EXP_MATRIX_COVARIANCE_FUNCTION_H
