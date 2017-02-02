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

#include <cmath>

#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/VectorSet.h>
#include <queso/GaussianLikelihoodFullCovariance.h>

namespace QUESO {

template<class V, class M>
GaussianLikelihoodFullCovariance<V, M>::GaussianLikelihoodFullCovariance(
    const char * prefix, const VectorSet<V, M> & domainSet,
    const V & observations, const M & covariance, double covarianceCoefficient)
  : LikelihoodBase<V, M>(prefix, domainSet, observations),
    m_covarianceCoefficient(covarianceCoefficient),
    m_covariance(covariance)
{
  if (covariance.numRowsLocal() != observations.sizeLocal()) {
    queso_error_msg("Covariance matrix not same size as observation vector");
  }
}

template<class V, class M>
GaussianLikelihoodFullCovariance<V, M>::~GaussianLikelihoodFullCovariance()
{
}

template<class V, class M>
double
GaussianLikelihoodFullCovariance<V, M>::lnValue(const V & domainVector) const
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

  return -0.5 * norm2_squared / (this->m_covarianceCoefficient);
}

}  // End namespace QUESO

template class QUESO::GaussianLikelihoodFullCovariance<QUESO::GslVector, QUESO::GslMatrix>;
