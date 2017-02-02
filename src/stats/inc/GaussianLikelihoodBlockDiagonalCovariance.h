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

#ifndef UQ_GAUSSIAN_LIKELIHOOD_BLOCK_DIAG_COV_H
#define UQ_GAUSSIAN_LIKELIHOOD_BLOCK_DIAG_COV_H

#include <vector>
#include <queso/GslBlockMatrix.h>
#include <queso/LikelihoodBase.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*!
 * \file GaussianLikelihoodBlockDiagonalCovariance.h
 *
 * \class GaussianLikelihoodBlockDiagonalCovariance
 * \brief A class representing a Gaussian likelihood with block-diagonal covariance matrix
 */

template <class V = GslVector, class M = GslMatrix>
class GaussianLikelihoodBlockDiagonalCovariance : public LikelihoodBase<V, M> {
public:
  //! @name Constructor/Destructor methods.
  //@{
  //! Default constructor.
  /*!
   * Instantiates a Gaussian likelihood function, given a prefix, its domain, a
   * vector of observations and a block diagonal covariance matrix.
   * The diagonal covariance matrix is of type \c GslBlockMatrix.  Each block
   * in the block diagonal matrix is an object of type \c GslMatrix.
   *
   * Furthermore, each block comes with a multiplicative coefficient which
   * defaults to 1.0.
   */
  GaussianLikelihoodBlockDiagonalCovariance(const char * prefix,
      const VectorSet<V, M> & domainSet, const V & observations,
      const GslBlockMatrix & covariance);

  //! Destructor
  virtual ~GaussianLikelihoodBlockDiagonalCovariance();
  //@}

  //! Get (non-const) multiplicative coefficient for block \c i
  double & blockCoefficient(unsigned int i);

  //! Get (const) multiplicative coefficient for block \c i
  const double & getBlockCoefficient(unsigned int i) const;

  //! Logarithm of the value of the scalar function.
  virtual double lnValue(const V & domainVector) const;

private:
  std::vector<double> m_covarianceCoefficients;
  const GslBlockMatrix & m_covariance;
};

}  // End namespace QUESO

#endif  // UQ_GAUSSIAN_LIKELIHOOD_BLOCK_DIAG_COV_H
