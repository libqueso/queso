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

#include <queso/InvLogitGaussianJointPdf.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Constructor -------------------------------------
template<class V,class M>
InvLogitGaussianJointPdf<V,M>::InvLogitGaussianJointPdf(
  const char*                  prefix,
  const BoxSubset<V,M>& domainBoxSubset,
  const V&                     lawExpVector,
  const V&                     lawVarVector)
  :
  BaseJointPdf<V,M>(((std::string)(prefix)+"invlogit_gau").c_str(),
      domainBoxSubset),
  m_gaussianPdf(((std::string)(prefix)+"invlogit_gau").c_str(),
      domainBoxSubset, lawExpVector, lawVarVector)
{
}


template<class V,class M>
InvLogitGaussianJointPdf<V,M>::InvLogitGaussianJointPdf(
  const char*                  prefix,
  const BoxSubset<V,M>& domainBoxSubset,
  const V&                     lawExpVector,
  const M&                     lawCovMatrix)
  :
  BaseJointPdf<V,M>(((std::string)(prefix)+"invlogit_gau").c_str(),
      domainBoxSubset),
  m_gaussianPdf(((std::string)(prefix)+"invlogit_gau").c_str(),
      domainBoxSubset, lawExpVector, lawCovMatrix)
{
}

template<class V,class M>
InvLogitGaussianJointPdf<V,M>::~InvLogitGaussianJointPdf()
{
}

template <class V, class M>
const V&
InvLogitGaussianJointPdf<V,M>::lawExpVector() const
{
  return this->m_gaussianPdf.lawExpVector();
}

template <class V, class M>
const V&
InvLogitGaussianJointPdf<V,M>::lawVarVector() const
{
  return this->m_gaussianPdf.lawVarVector();
}

template<class V, class M>
double
InvLogitGaussianJointPdf<V,M>::actualValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  double returnValue = 0.;

  if (this->m_domainSet.contains(domainVector) == false) { // prudenci 2011-Oct-04
    returnValue = 0.;
  }
  else {
    returnValue = std::exp(this->lnValue(domainVector,domainDirection,gradVector,hessianMatrix,hessianEffect));
  }

  return returnValue;
}

template<class V, class M>
double
InvLogitGaussianJointPdf<V,M>::lnValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  double returnValue;
  V transformedDomainVector(domainVector);

  double jacobian = 0.0;
  for (unsigned int i = 0; i < domainVector.sizeLocal(); i++) {
    transformedDomainVector[i] = std::log(domainVector[i] /
        (1 - domainVector[i]));
    jacobian -= std::log(domainVector[i]) + std::log(1.0 - domainVector[i]);
  }

  returnValue = this->m_gaussianPdf.lnValue(domainVector, domainDirection,
      gradVector, hessianMatrix, hessianEffect) + jacobian;

  return returnValue;
}

template<class V, class M>
double
InvLogitGaussianJointPdf<V,M>::computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const
{
  return this->m_gaussianPdf.computeLogOfNormalizationFactor(numSamples,
      updateFactorInternally);
}

template<class V, class M>
void
InvLogitGaussianJointPdf<V,M>::updateLawExpVector(const V& newLawExpVector)
{
  this->m_gaussianPdf.updateLawExpVector(newLawExpVector);
}

template<class V, class M>
void
InvLogitGaussianJointPdf<V,M>::updateLawCovMatrix(const M& newLawCovMatrix)
{
  this->m_gaussianPdf.updateLawCovMatrix(newLawCovMatrix);
}

template<class V, class M>
const M&
InvLogitGaussianJointPdf<V,M>::lawCovMatrix() const
{
  return this->m_gaussianPdf.lawCovMatrix();
}

}  // End namespace QUESO

template class QUESO::InvLogitGaussianJointPdf<QUESO::GslVector, QUESO::GslMatrix>;
