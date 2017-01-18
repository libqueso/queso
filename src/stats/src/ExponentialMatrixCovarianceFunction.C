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

#include <queso/ExponentialMatrixCovarianceFunction.h>

namespace QUESO {

// Default constructor -----------------------------
template<class P_V, class P_M, class Q_V, class Q_M>
ExponentialMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>::ExponentialMatrixCovarianceFunction(
  const char*                      prefix,
  const VectorSet<P_V,P_M>& basicDomainSet,
  const VectorSet<Q_V,Q_M>& imageSet,
  const Q_M&                       sigmas,
  const Q_M&                       as)
  :
  BaseMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>(prefix,basicDomainSet,imageSet),
  m_sigmas(NULL),
  m_as    (NULL)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering ExponentialMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  m_sigmas = new Q_M(sigmas);
  m_as     = new Q_M(as);

  unsigned int matrixOrder = m_imageSet.vectorSpace().dimLocal();

  queso_require_equal_to_msg(m_sigmas->numRowsLocal(), matrixOrder, "m_sigmas has invalid number of rows");

  queso_require_equal_to_msg(m_sigmas->numCols(), matrixOrder, "m_sigmas has invalid number of columns");

  queso_require_equal_to_msg(m_as->numRowsLocal(), matrixOrder, "m_as has invalid number of rows");

  queso_require_equal_to_msg(m_as->numCols(), matrixOrder, "m_as has invalid number of columns");

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving ExponentialMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor ---------------------------------------
template<class P_V, class P_M, class Q_V, class Q_M>
ExponentialMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>::~ExponentialMatrixCovarianceFunction()
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering ExponentialMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>::destructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  delete m_as;
  delete m_sigmas;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving ExponentialMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>::destructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Math methods -------------------------------------
template<class P_V, class P_M, class Q_V, class Q_M>
void
ExponentialMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>::covMatrix(const P_V& domainVector1, const P_V& domainVector2, Q_M& imageMatrix) const
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering ExponentialMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>::covMatrix()"
                          << std::endl;
  }

  unsigned int matrixOrder = m_imageSet.vectorSpace().dimLocal();

  queso_require_equal_to_msg(imageMatrix.numRowsLocal(), matrixOrder, "imageMatrix has invalid number of rows");

  queso_require_equal_to_msg(imageMatrix.numCols(), matrixOrder, "imageMatrix has invalid number of columns");

  double tmpSq = -(domainVector1 - domainVector2).norm2Sq();

  for (unsigned int i = 0; i < matrixOrder; ++i) {
    for (unsigned int j = 0; j < matrixOrder; ++j) {
      double tmp = tmpSq/( (*m_sigmas)(i,j) * (*m_sigmas)(i,j) );
      imageMatrix(i,j) = (*m_as)(i,j) * std::exp(tmp);
    }
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving ExponentialMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>::covMatrix()"
                          << std::endl;
  }

  return;
}

}  // End namespace QUESO
