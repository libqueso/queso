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

#include <queso/MatrixCovarianceFunction.h>

namespace QUESO {

// Default constructor -----------------------------
template<class P_V, class P_M, class Q_V, class Q_M>
BaseMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>::BaseMatrixCovarianceFunction(
  const char*                      prefix,
  const VectorSet<P_V,P_M>& basicDomainSet,
  const VectorSet<Q_V,Q_M>& imageSet)
  :
  m_env           (basicDomainSet.env()),
  m_prefix        ((std::string)(prefix)+"cov_func_"),
  m_basicDomainSet(basicDomainSet),
  m_imageSet      (imageSet)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering BaseMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving BaseMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor ---------------------------------------
template<class P_V, class P_M, class Q_V, class Q_M>
BaseMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>::~BaseMatrixCovarianceFunction()
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering BaseMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>::destructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving BaseMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>::destructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Math methods -------------------------------------
template<class P_V, class P_M, class Q_V, class Q_M>
const VectorSet<P_V,P_M>&
BaseMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>::basicDomainSet() const
{
  return m_basicDomainSet;
}

}  // End namespace QUESO
