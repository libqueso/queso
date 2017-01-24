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
#include <queso/GaussianLikelihoodDiagonalCovariance.h>

namespace QUESO {

template<class V, class M>
GaussianLikelihoodDiagonalCovariance<V, M>::GaussianLikelihoodDiagonalCovariance(
    const char * prefix, const VectorSet<V, M> & domainSet,
    const V & observations, const V & covariance)
  : LikelihoodBase<V, M>(prefix, domainSet, observations),
    m_covariance(covariance)
{
  if (covariance.sizeLocal() != observations.sizeLocal()) {
    queso_error_msg("Covariance matrix not same size as observation vector");
  }
}

template<class V, class M>
GaussianLikelihoodDiagonalCovariance<V, M>::~GaussianLikelihoodDiagonalCovariance()
{
}

template<class V, class M>
double
GaussianLikelihoodDiagonalCovariance<V, M>::lnValue(const V & domainVector) const
{
  V modelOutput(this->m_observations, 0, 0);  // At least it's not a copy

  this->evaluateModel(domainVector, modelOutput);

  modelOutput -= this->m_observations;  // Compute misfit
  modelOutput *= modelOutput;
  modelOutput /= this->m_covariance;  // Multiply by inverse covriance matrix

  double norm2_squared = modelOutput.sumOfComponents();  // This is square of 2-norm

  return -0.5 * norm2_squared;
}

}  // End namespace QUESO

template class QUESO::GaussianLikelihoodDiagonalCovariance<QUESO::GslVector, QUESO::GslMatrix>;
