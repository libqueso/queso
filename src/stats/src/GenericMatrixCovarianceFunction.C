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

#include <queso/GenericMatrixCovarianceFunction.h>

namespace QUESO {

// Default constructor -----------------------------

template<class P_V, class P_M, class Q_V, class Q_M>
GenericMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>::GenericMatrixCovarianceFunction(
  const char*                      prefix,
  const VectorSet<P_V,P_M>& basicDomainSet,
  const VectorSet<Q_V,Q_M>& imageSet,
  void (*covRoutinePtr)(const P_V& positionVector1, const P_V& positionVector2, const void* routineDataPtr, Q_M& imageMatrix),
  const void*                      routinesDataPtr)
  :
  BaseMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>(prefix,basicDomainSet,imageSet),
  m_covRoutinePtr                                     (covRoutinePtr),
  m_routineDataPtr                                    (routinesDataPtr)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering GenericMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving GenericMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor ---------------------------------------
template<class P_V, class P_M, class Q_V, class Q_M>
GenericMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>::~GenericMatrixCovarianceFunction()
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering GenericMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>::destructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving GenericMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>::destructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Math methods -------------------------------------
template<class P_V, class P_M, class Q_V, class Q_M>
void
GenericMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>::covMatrix(const P_V& positionVector1, const P_V& positionVector2, Q_M& imageMatrix) const
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering GenericMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>::covMatrix()"
                            << std::endl;
  }

  unsigned int matrixOrder = m_imageSet.vectorSpace().dimLocal();

  queso_require_equal_to_msg(imageMatrix.numRowsLocal(), matrixOrder, "imageMatrix has invalid number of rows");

  queso_require_equal_to_msg(imageMatrix.numCols(), matrixOrder, "imageMatrix has invalid number of columns");

  queso_require_msg(m_covRoutinePtr, "m_covRoutinePtr = NULL");

  m_covRoutinePtr(positionVector1, positionVector2, m_routineDataPtr, imageMatrix);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving GenericMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>::covMatrix()"
                            << std::endl;
  }

  return;
}

}  // End namespace QUESO
