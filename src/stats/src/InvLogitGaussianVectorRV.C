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

#include <queso/InvLogitGaussianVectorRV.h>
#include <queso/InvLogitGaussianVectorRealizer.h>
#include <queso/InvLogitGaussianJointPdf.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Constructor---------------------------------------
template<class V, class M>
InvLogitGaussianVectorRV<V,M>::InvLogitGaussianVectorRV(
    const char * prefix,
    const BoxSubset<V, M> & imageBoxSubset,
    const V & lawExpVector,
    const V & lawVarVector)
  : BaseVectorRV<V, M>(((std::string)(prefix)+"invlogit_gau").c_str(),
      imageBoxSubset)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering InvLogitGaussianVectorRV<V,M>::constructor() [1]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  queso_require_greater_msg(lawVarVector.getMinValue(), 0.0, "Covariance matrix is not symmetric positive definite.");

  m_pdf = new InvLogitGaussianJointPdf<V,M>(m_prefix.c_str(),
      dynamic_cast<const BoxSubset<V, M> & >(m_imageSet), lawExpVector,
      lawVarVector);

  V cholDiag(lawVarVector);
  cholDiag.cwSqrt();
  M lowerCholLawCovMatrix(cholDiag);
  lowerCholLawCovMatrix.zeroUpper(false);

  m_realizer = new InvLogitGaussianVectorRealizer<V,M>(m_prefix.c_str(),
      dynamic_cast<const BoxSubset<V, M> & >(m_imageSet), lawExpVector,
      lowerCholLawCovMatrix);

  m_subCdf     = NULL; // FIX ME: complete code
  m_unifiedCdf = NULL; // FIX ME: complete code
  m_mdf        = NULL; // FIX ME: complete code

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving InvLogitGaussianVectorRV<V,M>::constructor() [1]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}

template<class V, class M>
InvLogitGaussianVectorRV<V, M>::InvLogitGaussianVectorRV(
    const char * prefix,
    const BoxSubset<V, M> & imageBoxSubset,
    const V & lawExpVector,
    const M & lawCovMatrix)
  : BaseVectorRV<V, M>(((std::string)(prefix)+"invlogit_gau").c_str(),
      imageBoxSubset)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering InvLogitGaussianVectorRV<V,M>::constructor() [2]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  m_pdf = new InvLogitGaussianJointPdf<V, M>(m_prefix.c_str(),
      dynamic_cast<const BoxSubset<V, M> & >(m_imageSet),
      lawExpVector, lawCovMatrix);

  M lowerCholLawCovMatrix(lawCovMatrix);
  int iRC = lowerCholLawCovMatrix.chol();
  lowerCholLawCovMatrix.zeroUpper(false);
  if (iRC) {
    std::cerr << "In InvLogitGaussianVectorRV<V,M>::constructor() [2]: chol failed, will use svd\n";
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In InvLogitGaussianVectorRV<V,M>::constructor() [2]: chol failed; will use svd; lawCovMatrix contents are\n";
      *m_env.subDisplayFile() << lawCovMatrix; // FIX ME: might demand parallelism
      *m_env.subDisplayFile() << std::endl;
    }
    M matU (lawCovMatrix);
    M matVt(m_imageSet.vectorSpace().zeroVector());
    V vecS (m_imageSet.vectorSpace().zeroVector());
    iRC = lawCovMatrix.svd(matU,vecS,matVt);
    queso_require_msg(!(iRC), "Cholesky decomposition of covariance matrix failed.");

    vecS.cwSqrt();
    m_realizer = new InvLogitGaussianVectorRealizer<V,M>(m_prefix.c_str(),
        dynamic_cast<const BoxSubset<V, M> & >(m_imageSet), lawExpVector, matU,
        vecS, matVt);
  }
  else {
    m_realizer = new InvLogitGaussianVectorRealizer<V, M>(m_prefix.c_str(),
        dynamic_cast<const BoxSubset<V, M> & >(m_imageSet), lawExpVector, lowerCholLawCovMatrix);
  }

  m_subCdf     = NULL; // FIX ME: complete code
  m_unifiedCdf = NULL; // FIX ME: complete code
  m_mdf        = NULL; // FIX ME: complete code

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving InvLogitGaussianVectorRV<V,M>::constructor() [2]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}

template<class V, class M>
InvLogitGaussianVectorRV<V, M>::~InvLogitGaussianVectorRV()
{
  delete m_mdf;
  delete m_unifiedCdf;
  delete m_subCdf;
  delete m_realizer;
  delete m_pdf;
}

template<class V, class M>
void
InvLogitGaussianVectorRV<V, M>::updateLawExpVector(const V & newLawExpVector)
{
  // We are sure that m_pdf (and m_realizer, etc) point to associated Gaussian
  // classes, so all is well
  (dynamic_cast<InvLogitGaussianJointPdf<V, M> * >(m_pdf))->updateLawExpVector(
      newLawExpVector);
  (dynamic_cast<InvLogitGaussianVectorRealizer<V, M> * >(m_realizer))->
    updateLawExpVector(newLawExpVector);
  return;
}

template<class V, class M>
void
InvLogitGaussianVectorRV<V, M>::updateLawCovMatrix(const M & newLawCovMatrix)
{
  // We are sure that m_pdf (and m_realizer, etc) point to associated Gaussian
  // classes, so all is well
  (dynamic_cast<InvLogitGaussianJointPdf<V,M> * >(m_pdf))->updateLawCovMatrix(
      newLawCovMatrix);

  M newLowerCholLawCovMatrix(newLawCovMatrix);
  int iRC = newLowerCholLawCovMatrix.chol();
  newLowerCholLawCovMatrix.zeroUpper(false);
  if (iRC) {
    std::cerr << "In InvLogitGaussianVectorRV<V,M>::updateLawCovMatrix(): chol failed, will use svd\n";
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In InvLogitGaussianVectorRV<V,M>::updateLawCovMatrix(): chol failed; will use svd; newLawCovMatrix contents are\n";
      *m_env.subDisplayFile() << newLawCovMatrix; // FIX ME: might demand parallelism
      *m_env.subDisplayFile() << std::endl;
    }
    M matU (newLawCovMatrix);
    M matVt(m_imageSet.vectorSpace().zeroVector());
    V vecS (m_imageSet.vectorSpace().zeroVector());
    iRC = newLawCovMatrix.svd(matU,vecS,matVt);
    queso_require_msg(!(iRC), "Cholesky decomposition of covariance matrix failed.");

    vecS.cwSqrt();
    (dynamic_cast<InvLogitGaussianVectorRealizer<V, M> * >(m_realizer))->
      updateLowerCholLawCovMatrix(matU, vecS, matVt);
  }
  else {
    (dynamic_cast<InvLogitGaussianVectorRealizer<V, M> * >(m_realizer))->
      updateLowerCholLawCovMatrix(newLowerCholLawCovMatrix);
  }
  return;
}

template <class V, class M>
void
InvLogitGaussianVectorRV<V, M>::print(std::ostream & os) const
{
  os << "InvLogitGaussianVectorRV<V,M>::print() says, 'Please implement me.'" << std::endl;
  return;
}

}  // End namespace QUESO

template class QUESO::InvLogitGaussianVectorRV<QUESO::GslVector,QUESO::GslMatrix>;
