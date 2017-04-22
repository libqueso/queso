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

#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/VectorSet.h>
#include <queso/VectorSpace.h>
#include <queso/GaussianLikelihoodFullCovarianceRandomCoefficient.h>

namespace QUESO {

template<class V, class M>
GaussianLikelihoodFullCovarianceRandomCoefficient<V, M>::GaussianLikelihoodFullCovarianceRandomCoefficient(
    const char * prefix, const VectorSet<V, M> & domainSet,
    const V & observations, const M & covariance)
  : LikelihoodBase<V, M>(prefix, domainSet, observations),
    m_covariance(covariance)
{
  if (covariance.numRowsLocal() != observations.sizeLocal()) {
    queso_error_msg("Covariance matrix not same size as observation vector");
  }
}

template<class V, class M>
GaussianLikelihoodFullCovarianceRandomCoefficient<V, M>::~GaussianLikelihoodFullCovarianceRandomCoefficient()
{
}

template<class V, class M>
double
GaussianLikelihoodFullCovarianceRandomCoefficient<V, M>::lnValue(const V & domainVector) const
{
  V modelOutput(this->m_observations, 0, 0);  // At least it's not a copy
  V weightedMisfit(this->m_observations, 0, 0);  // At least it's not a copy

  this->evaluateModel(domainVector, modelOutput);

  // Compute misfit G(x) - y
  modelOutput -= this->m_observations;

  // Solve \Sigma u = G(x) - y for u
  this->m_covariance.invertMultiply(modelOutput, weightedMisfit);

  // Compute (G(x) - y)^T \Sigma^{-1} (G(x) - y)
  modelOutput *= weightedMisfit;

  // This is square of 2-norm
  double norm2_squared = modelOutput.sumOfComponents();

  // Get the determinant of the covariance matrix |\Sigma|
  double deter_cov = this->m_covariance.determinant();

  deter_cov = std::sqrt(deter_cov);

  // Set the right hyperparameter coefficient
  // The last element of domainVector is the multiplicative coefficient of the
  // covariance matrix
  double cov_coeff = domainVector[domainVector.sizeLocal()-1];
  cov_coeff = std::pow(std::sqrt(cov_coeff), this->m_observations.sizeLocal());

  return -0.5 * norm2_squared / cov_coeff - std::log(cov_coeff * deter_cov);
}

}  // End namespace QUESO

template class QUESO::GaussianLikelihoodFullCovarianceRandomCoefficient<QUESO::GslVector, QUESO::GslMatrix>;
