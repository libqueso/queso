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
#include <queso/GaussianLikelihoodBlockDiagonalCovarianceRandomCoefficients.h>

namespace QUESO {

template<class V, class M>
GaussianLikelihoodBlockDiagonalCovarianceRandomCoefficients<V, M>::GaussianLikelihoodBlockDiagonalCovarianceRandomCoefficients(
    const char * prefix, const VectorSet<V, M> & domainSet,
    const V & observations, const GslBlockMatrix & covariance)
  : LikelihoodBase<V, M>(prefix, domainSet, observations),
    m_covariance(covariance)
{
  unsigned int totalDim = 0;

  for (unsigned int i = 0; i < this->m_covariance.numBlocks(); i++) {
    totalDim += this->m_covariance.getBlock(i).numRowsLocal();
  }

  if (totalDim != observations.sizeLocal()) {
    queso_error_msg("Covariance matrix not same dimension as observation vector");
  }
}

template<class V, class M>
GaussianLikelihoodBlockDiagonalCovarianceRandomCoefficients<V, M>::~GaussianLikelihoodBlockDiagonalCovarianceRandomCoefficients()
{
}

template<class V, class M>
double
GaussianLikelihoodBlockDiagonalCovarianceRandomCoefficients<V, M>::lnValue(const V & domainVector) const
{
  V modelOutput(this->m_observations, 0, 0);  // At least it's not a copy
  V weightedMisfit(this->m_observations, 0, 0);  // At least it's not a copy

  this->evaluateModel(domainVector, modelOutput);

  // Compute misfit G(x) - y
  modelOutput -= this->m_observations;

  // Solve \Sigma u = G(x) - y for u
  this->m_covariance.invertMultiply(modelOutput, weightedMisfit);

  // Deal with the multiplicative coefficients for each of the blocks
  unsigned int numBlocks = this->m_covariance.numBlocks();
  unsigned int offset = 0;

  // For each block...
  double cov_norm_factor = 0.0;
  for (unsigned int i = 0; i < this->m_covariance.numBlocks(); i++) {

    // ...find the right hyperparameter
    unsigned int index = domainVector.sizeLocal() + (i - numBlocks);
    double coefficient = domainVector[index];

    // ...divide the appropriate parts of the solution by the coefficient
    unsigned int blockDim = this->m_covariance.getBlock(i).numRowsLocal();
    for (unsigned int j = 0; j < blockDim; j++) {
      // 'coefficient' is a variance, so we divide by it
      modelOutput[offset+j] /= coefficient;
    }

    // Keep track of the part of the covariance matrix that appears in the
    // normalising constant because of the hyperparameter
    double cov_determinant = this->m_covariance.getBlock(i).determinant();
    cov_determinant = std::sqrt(cov_determinant);

    coefficient = std::sqrt(coefficient);
    cov_norm_factor += std::log(std::pow(coefficient, blockDim) * cov_determinant);

    offset += blockDim;
  }

  // Compute (G(x) - y)^T \Sigma^{-1} (G(x) - y)
  modelOutput *= weightedMisfit;

  double norm2_squared = modelOutput.sumOfComponents();  // This is square of 2-norm

  return -0.5 * norm2_squared - cov_norm_factor;
}

}  // End namespace QUESO

template class QUESO::GaussianLikelihoodBlockDiagonalCovarianceRandomCoefficients<QUESO::GslVector, QUESO::GslMatrix>;
