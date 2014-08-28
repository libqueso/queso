//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
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

#include <queso/InvLogitGaussianVectorRealizer.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

template<class V, class M>
InvLogitGaussianVectorRealizer<V, M>::InvLogitGaussianVectorRealizer(
    const char * prefix,
    const BoxSubset<V, M> & unifiedImageBoxSubset,
    const V & lawExpVector,
    const M & lowerCholLawCovMatrix)
  : BaseVectorRealizer<V, M>(((std::string)(prefix)+"invlogit_gau").c_str(),
      unifiedImageBoxSubset, std::numeric_limits<unsigned int>::max()),
    m_gaussianRealizer(((std::string)(prefix)+"invlogit_gau").c_str(),
        unifiedImageBoxSubset, lawExpVector, lowerCholLawCovMatrix)
{
}

template<class V, class M>
InvLogitGaussianVectorRealizer<V, M>::InvLogitGaussianVectorRealizer(
    const char * prefix,
    const BoxSubset<V, M> & unifiedImageBoxSubset,
    const V & lawExpVector,
    const M & matU,
    const V & vecSsqrt,
    const M & matVt)
  : BaseVectorRealizer<V, M>(((std::string)(prefix)+"invlogit_gau").c_str(),
      unifiedImageBoxSubset, std::numeric_limits<unsigned int>::max()),
    m_gaussianRealizer(((std::string)(prefix)+"invlogit_gau").c_str(),
        unifiedImageBoxSubset, lawExpVector, matU, vecSsqrt, matVt)
{
}

template<class V, class M>
InvLogitGaussianVectorRealizer<V, M>::~InvLogitGaussianVectorRealizer()
{
}

template <class V, class M>
const V &
InvLogitGaussianVectorRealizer<V, M>::unifiedLawExpVector() const
{
  return this->m_gaussianRealizer.unifiedLawExpVector();
}

template <class V, class M>
const V &
InvLogitGaussianVectorRealizer<V, M>::unifiedLawVarVector() const
{
  return this->m_gaussianRealizer.unifiedLawVarVector();
}

template<class V, class M>
void
InvLogitGaussianVectorRealizer<V, M>::realization(V & nextValues) const
{
  this->m_gaussianRealizer.realization(nextValues);

  for (unsigned int i = 0; i < nextValues.sizeLocal(); i++) {
    double temp = std::exp(nextValues[i]);
    nextValues[i] = temp / (1.0 + temp);
  }
}

template<class V, class M>
void
InvLogitGaussianVectorRealizer<V, M>::updateLawExpVector(
    const V & newLawExpVector)
{
  this->m_gaussianRealizer.updateLawExpVector(newLawExpVector);
}

template<class V, class M>
void
InvLogitGaussianVectorRealizer<V, M>::updateLowerCholLawCovMatrix(
    const M & newLowerCholLawCovMatrix)
{
  this->m_gaussianRealizer.updateLowerCholLawCovMatrix(newLowerCholLawCovMatrix);
}

template<class V, class M>
void
InvLogitGaussianVectorRealizer<V, M>::updateLowerCholLawCovMatrix(
  const M & matU,
  const V & vecSsqrt,
  const M & matVt)
{
  this->m_gaussianRealizer.updateLowerCholLawCovMatrix(matU, vecSsqrt, matVt);
}

}  // End namespace QUESO

template class QUESO::InvLogitGaussianVectorRealizer<QUESO::GslVector, QUESO::GslMatrix>;
