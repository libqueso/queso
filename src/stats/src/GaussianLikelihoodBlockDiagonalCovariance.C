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
#include <queso/GaussianLikelihoodBlockDiagonalCovariance.h>

namespace QUESO {

template<class V, class M>
GaussianLikelihoodBlockDiagonalCovariance<V, M>::GaussianLikelihoodBlockDiagonalCovariance(
    const char * prefix, const VectorSet<V, M> & domainSet,
    const V & observations, const GslBlockMatrix & covariance)
  : LikelihoodBase<V, M>(prefix, domainSet, observations),
    m_covarianceCoefficients(covariance.numBlocks(), 1.0),
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
GaussianLikelihoodBlockDiagonalCovariance<V, M>::~GaussianLikelihoodBlockDiagonalCovariance()
{
}

template<class V, class M>
double &
GaussianLikelihoodBlockDiagonalCovariance<V, M>::blockCoefficient(
    unsigned int i)
{
  return this->m_covarianceCoefficients[i];
}

template<class V, class M>
const double &
GaussianLikelihoodBlockDiagonalCovariance<V, M>::getBlockCoefficient(
    unsigned int i) const
{
  return this->m_covarianceCoefficients[i];
}

template<class V, class M>
double
GaussianLikelihoodBlockDiagonalCovariance<V, M>::lnValue(const V & domainVector) const
{
  V modelOutput(this->m_observations, 0, 0);  // At least it's not a copy
  V weightedMisfit(this->m_observations, 0, 0);  // At least it's not a copy

  this->evaluateModel(domainVector, modelOutput);

  // Compute misfit G(x) - y
  modelOutput -= this->m_observations;

  // Solve \Sigma u = G(x) - y for u
  this->m_covariance.invertMultiply(modelOutput, weightedMisfit);

  // Deal with the multiplicative coefficients for each of the blocks
  unsigned int offset = 0;

  // For each block...
  for (unsigned int i = 0; i < this->m_covariance.numBlocks(); i++) {
    // ...divide the appropriate parts of the solution by the coefficient
    unsigned int blockDim = this->m_covariance.getBlock(i).numRowsLocal();
    for (unsigned int j = 0; j < blockDim; j++) {
      // coefficient is a variance, so we divide by it
      modelOutput[offset+j] /= this->m_covarianceCoefficients[i];
    }
    offset += blockDim;
  }

  // Compute (G(x) - y)^T \Sigma^{-1} (G(x) - y)
  modelOutput *= weightedMisfit;

  double norm2_squared = modelOutput.sumOfComponents();  // This is square of 2-norm

  return -0.5 * norm2_squared;
}

}  // End namespace QUESO

template class QUESO::GaussianLikelihoodBlockDiagonalCovariance<QUESO::GslVector, QUESO::GslMatrix>;
