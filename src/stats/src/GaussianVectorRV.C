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

#include <queso/GaussianVectorRV.h>
#include <queso/GaussianVectorRealizer.h>
#include <queso/GaussianJointPdf.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Constructor---------------------------------------
template<class V, class M>
GaussianVectorRV<V,M>::GaussianVectorRV(
  const char*                  prefix,
  const VectorSet<V,M>& imageSet,
  const V&                     lawExpVector,
  const V&                     lawVarVector)
  :
  BaseVectorRV<V,M>(((std::string)(prefix)+"gau").c_str(),imageSet)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering GaussianVectorRV<V,M>::constructor() [1]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  queso_require_greater_msg(lawVarVector.getMinValue(), 0.0, "Covariance matrix is not symmetric positive definite.");

  m_pdf = new GaussianJointPdf<V,M>(m_prefix.c_str(),
                                            m_imageSet,
                                            lawExpVector,
                                            lawVarVector);

  V cholDiag(lawVarVector);
  cholDiag.cwSqrt();
  M lowerCholLawCovMatrix(cholDiag);
  lowerCholLawCovMatrix.zeroUpper(false);

  m_realizer = new GaussianVectorRealizer<V,M>(m_prefix.c_str(),
                                                      m_imageSet,
                                                      lawExpVector,
                                                      lowerCholLawCovMatrix);

  m_subCdf     = NULL; // FIX ME: complete code
  m_unifiedCdf = NULL; // FIX ME: complete code
  m_mdf        = NULL; // FIX ME: complete code

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving GaussianVectorRV<V,M>::constructor() [1]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Constructor---------------------------------------
template<class V, class M>
GaussianVectorRV<V,M>::GaussianVectorRV(
  const char*                  prefix,
  const VectorSet<V,M>& imageSet,
  const V&                     lawExpVector,
  const M&                     lawCovMatrix)
  :
  BaseVectorRV<V,M>(((std::string)(prefix)+"gau").c_str(),imageSet)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering GaussianVectorRV<V,M>::constructor() [2]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  m_pdf = new GaussianJointPdf<V,M>(m_prefix.c_str(),
                                           m_imageSet,
                                           lawExpVector,
                                           lawCovMatrix);

  M lowerCholLawCovMatrix(lawCovMatrix);
  int iRC = lowerCholLawCovMatrix.chol();
  lowerCholLawCovMatrix.zeroUpper(false);
  if (iRC) {
    std::cerr << "In GaussianVectorRV<V,M>::constructor() [2]: chol failed, will use svd\n";
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GaussianVectorRV<V,M>::constructor() [2]: chol failed; will use svd; lawCovMatrix contents are\n";
      *m_env.subDisplayFile() << lawCovMatrix; // FIX ME: might demand parallelism
      *m_env.subDisplayFile() << std::endl;
    }
    M matU (lawCovMatrix);
    M matVt(m_imageSet.vectorSpace().zeroVector());
    V vecS (m_imageSet.vectorSpace().zeroVector());
    iRC = lawCovMatrix.svd(matU,vecS,matVt);
    queso_require_msg(!(iRC), "Cholesky decomposition of covariance matrix failed.");

    vecS.cwSqrt();
    m_realizer = new GaussianVectorRealizer<V,M>(m_prefix.c_str(),
                                                        m_imageSet,
                                                        lawExpVector,
                                                        matU,
                                                        vecS, // already square rooted
                                                        matVt);
  }
  else {
    m_realizer = new GaussianVectorRealizer<V,M>(m_prefix.c_str(),
                                                        m_imageSet,
                                                        lawExpVector,
                                                        lowerCholLawCovMatrix);
  }

  m_subCdf     = NULL; // FIX ME: complete code
  m_unifiedCdf = NULL; // FIX ME: complete code
  m_mdf        = NULL; // FIX ME: complete code

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving GaussianVectorRV<V,M>::constructor() [2]"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor ---------------------------------------
template<class V, class M>
GaussianVectorRV<V,M>::~GaussianVectorRV()
{
  delete m_mdf;
  delete m_unifiedCdf;
  delete m_subCdf;
  delete m_realizer;
  delete m_pdf;
}
// Statistical methods-------------------------------
template<class V, class M>
void
GaussianVectorRV<V,M>::updateLawExpVector(const V& newLawExpVector)
{
  // We are sure that m_pdf (and m_realizer, etc) point to associated Gaussian classes, so all is well
  ( dynamic_cast< GaussianJointPdf      <V,M>* >(m_pdf     ) )->updateLawExpVector(newLawExpVector);
  ( dynamic_cast< GaussianVectorRealizer<V,M>* >(m_realizer) )->updateLawExpVector(newLawExpVector);
  return;
}
//---------------------------------------------------
template<class V, class M>
void
GaussianVectorRV<V,M>::updateLawCovMatrix(const M& newLawCovMatrix)
{
  // We are sure that m_pdf (and m_realizer, etc) point to associated Gaussian classes, so all is well
  ( dynamic_cast< GaussianJointPdf<V,M>* >(m_pdf) )->updateLawCovMatrix(newLawCovMatrix);

  M newLowerCholLawCovMatrix(newLawCovMatrix);
  int iRC = newLowerCholLawCovMatrix.chol();
  newLowerCholLawCovMatrix.zeroUpper(false);
  if (iRC) {
    std::cerr << "In GaussianVectorRV<V,M>::updateLawCovMatrix(): chol failed, will use svd\n";
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GaussianVectorRV<V,M>::updateLawCovMatrix(): chol failed; will use svd; newLawCovMatrix contents are\n";
      *m_env.subDisplayFile() << newLawCovMatrix; // FIX ME: might demand parallelism
      *m_env.subDisplayFile() << std::endl;
    }
    M matU (newLawCovMatrix);
    M matVt(m_imageSet.vectorSpace().zeroVector());
    V vecS (m_imageSet.vectorSpace().zeroVector());
    iRC = newLawCovMatrix.svd(matU,vecS,matVt);
    queso_require_msg(!(iRC), "Cholesky decomposition of covariance matrix failed.");

    vecS.cwSqrt();
    ( dynamic_cast< GaussianVectorRealizer<V,M>* >(m_realizer) )->updateLowerCholLawCovMatrix(matU,
                                                                                                     vecS, // already square rooted
                                                                                                     matVt);
  }
  else {
    ( dynamic_cast< GaussianVectorRealizer<V,M>* >(m_realizer) )->updateLowerCholLawCovMatrix(newLowerCholLawCovMatrix);
  }
  return;
}
// I/O methods---------------------------------------
template <class V, class M>
void
GaussianVectorRV<V,M>::print(std::ostream& os) const
{
  os << "GaussianVectorRV<V,M>::print() says, 'Please implement me.'" << std::endl;
  return;
}


//---------------------------------------------------
// Method declared outside class definition ---------
//---------------------------------------------------
template<class V, class M>
void
ComputeConditionalGaussianVectorRV(
  const V& muVec1,
  const V& muVec2,
  const M& sigmaMat11,
  const M& sigmaMat12,
  const M& sigmaMat21,
  const M& sigmaMat22,
  const V& sampleVec2,
        V& muVec1_cond_on_2,
        M& sigmaMat11_cond_on_2)
{
  const BaseEnvironment& env = muVec1.env();
  unsigned int dim1 = muVec1.sizeLocal();
  unsigned int dim2 = muVec2.sizeLocal();

  queso_require_msg(!((sigmaMat11.numRowsLocal() != dim1) || (sigmaMat11.numCols() != dim1)), "invalid sigmaMat11");

  queso_require_msg(!((sigmaMat12.numRowsLocal() != dim1) || (sigmaMat12.numCols() != dim2)), "invalid sigmaMat12");

  queso_require_msg(!((sigmaMat21.numRowsLocal() != dim2) || (sigmaMat21.numCols() != dim1)), "invalid sigmaMat21");

  queso_require_msg(!((sigmaMat22.numRowsLocal() != dim2) || (sigmaMat22.numCols() != dim2)), "invalid sigmaMat22");

  // Check transpose operation
  M mat_tt(sigmaMat12);
  mat_tt.cwSet(0.);
  mat_tt.fillWithTranspose(0,0,sigmaMat21,true,true);
  double auxNorm = (mat_tt - sigmaMat12).normFrob();
  if (auxNorm >= 1.e-12) {
    if (env.subDisplayFile()) {
      *env.subDisplayFile() << "In ComputeConditionalGaussianVectorRV()"
                            << ": WARNING, ||sigmaMat21^T - sigmaMat12||_2 = " << auxNorm
                            << std::endl;
    }
  }
  queso_require_less_msg(auxNorm, 1.e-12, "sigmaMat12 and sigmaMat21 are not transpose of each other");

  queso_require_equal_to_msg(sampleVec2.sizeLocal(), dim2, "invalid sampleVec2");

  queso_require_equal_to_msg(muVec1_cond_on_2.sizeLocal(), dim1, "invalid muVec1_cond_on_2");

  queso_require_msg(!((sigmaMat11_cond_on_2.numRowsLocal() != dim1) || (sigmaMat11_cond_on_2.numCols() != dim1)), "invalid sigmaMat11_cond_on_2");

  muVec1_cond_on_2     = muVec1     + sigmaMat12 * sigmaMat22.invertMultiply(sampleVec2 - muVec2);
  sigmaMat11_cond_on_2 = sigmaMat11 - sigmaMat12 * sigmaMat22.invertMultiply(sigmaMat21);

  return;
}

}  // End namespace QUESO

template class QUESO::GaussianVectorRV<QUESO::GslVector,QUESO::GslMatrix>;

template void QUESO::ComputeConditionalGaussianVectorRV<QUESO::GslVector, QUESO::GslMatrix>(QUESO::GslVector const&, QUESO::GslVector const&, QUESO::GslMatrix const&, QUESO::GslMatrix const&, QUESO::GslMatrix const&, QUESO::GslMatrix const&, QUESO::GslVector const&, QUESO::GslVector&, QUESO::GslMatrix&);
