//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2015 The PECOS Development Team
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
  : BaseGaussianLikelihood<V, M>(prefix, domainSet, observations),
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
GaussianLikelihoodFullCovarianceRandomCoefficient<V, M>::actualValue(
    const V & domainVector, const V * domainDirection, V * gradVector,
    M * hessianMatrix, V * hessianEffect) const
{
  return std::exp(this->lnValue(domainVector, domainDirection, gradVector,
        hessianMatrix, hessianEffect));
}

template<class V, class M>
double
GaussianLikelihoodFullCovarianceRandomCoefficient<V, M>::lnValue(
    const V & domainVector, const V * domainDirection, V * gradVector,
    M * hessianMatrix, V * hessianEffect) const
{
  V modelOutput(this->m_observations, 0, 0);  // At least it's not a copy
  V weightedMisfit(this->m_observations, 0, 0);  // At least it's not a copy

  this->evaluateModel(domainVector, domainDirection, modelOutput, gradVector,
      hessianMatrix, hessianEffect);

  // Compute misfit G(x) - y
  modelOutput -= this->m_observations;

  // Solve \Sigma u = G(x) - y for u
  this->m_covariance.invertMultiply(modelOutput, weightedMisfit);

  // Compute (G(x) - y)^T \Sigma^{-1} (G(x) - y)
  modelOutput *= weightedMisfit;

  // This is square of 2-norm
  double norm2_squared = modelOutput.sumOfComponents();

  // The last element of domainVector is the multiplicative coefficient of the
  // covariance matrix
  return -0.5 * norm2_squared / (domainVector[domainVector.sizeLocal()-1]);
}

}  // End namespace QUESO

template class QUESO::GaussianLikelihoodFullCovarianceRandomCoefficient<QUESO::GslVector, QUESO::GslMatrix>;
