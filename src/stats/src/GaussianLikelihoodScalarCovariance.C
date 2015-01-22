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
#include <queso/GaussianLikelihoodScalarCovariance.h>

namespace QUESO {

template<class V, class M>
GaussianLikelihoodScalarCovariance<V, M>::GaussianLikelihoodScalarCovariance(
    const char * prefix, const VectorSet<V, M> & domainSet,
    const std::vector<double> & observations, double variance)
  : BaseGaussianLikelihood<V, M>(prefix, domainSet, observations),
    m_variance(variance)
{
}

template<class V, class M>
GaussianLikelihoodScalarCovariance<V, M>::~GaussianLikelihoodScalarCovariance()
{
}

template<class V, class M>
double
GaussianLikelihoodScalarCovariance<V, M>::actualValue(const V & domainVector,
    const V * domainDirection, V * gradVector, M * hessianMatrix,
    V * hessianEffect) const
{
  return std::exp(this->lnValue(domainVector, domainDirection, gradVector,
        hessianMatrix, hessianEffect));
}

template<class V, class M>
double
GaussianLikelihoodScalarCovariance<V, M>::lnValue(const V & domainVector,
    const V * domainDirection, V * gradVector, M * hessianMatrix,
    V * hessianEffect) const
{
  double misfit = 0.0;

  for (unsigned int i = 0; i < this->m_observations.size(); i++) {
    double diff = this->m_modelOutput[i] - this->m_observations[i];
    misfit += diff * diff;
  }

  return -0.5 * misfit / m_variance;
}

}  // End namespace QUESO

template class QUESO::GaussianLikelihoodScalarCovariance<QUESO::GslVector, QUESO::GslMatrix>;
