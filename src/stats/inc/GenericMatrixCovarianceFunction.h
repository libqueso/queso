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

#ifndef UQ_GENERIC_MATRIX_COVARIANCE_FUNCTION_H
#define UQ_GENERIC_MATRIX_COVARIANCE_FUNCTION_H

#include <queso/MatrixCovarianceFunction.h>
#include <queso/VectorSet.h>
#include <queso/Environment.h>
#include <cmath>

namespace QUESO {

class GslVector;
class GslMatrix;

//*****************************************************
// Generic class
//*****************************************************
/*!
 * \class GenericMatrixCovarianceFunction
 * \brief A class for generic covariance matrices.
 *
 * This class implements a generic covariance matrices by calling a routine (via pointer).*/

template <class P_V = GslVector, class P_M = GslMatrix, class Q_V = GslVector, class Q_M = GslMatrix>
class GenericMatrixCovarianceFunction : public BaseMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M> {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  /*! Instantiates an object of the class given a prefix, the domain set, the pointer to the routine. */
  GenericMatrixCovarianceFunction(const char*                      prefix,
           const VectorSet<P_V,P_M>& basicDomainSet,
           const VectorSet<Q_V,Q_M>& imageSet,
           void (*covRoutinePtr)(const P_V& positionVector1, const P_V& positionVector2, const void* routineDataPtr, Q_M& imageMatrix),
           const void*                      routinesDataPtr);
  //! Virtual destructor
  virtual ~GenericMatrixCovarianceFunction();
  //@}

  //! @name Math methods
  //@{
  //! Calculates the value of the generic covariance matrix.
  /*! This function accesses the routine that calculates the covariance function. */
  void covMatrix(const P_V& positionVector1, const P_V& positionVector2, Q_M& imageMatrix) const;

protected:
  using BaseMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>::m_env;
  using BaseMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>::m_prefix;
  using BaseMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>::m_basicDomainSet;
  using BaseMatrixCovarianceFunction<P_V,P_M,Q_V,Q_M>::m_imageSet;

  void (*m_covRoutinePtr)(const P_V& positionVector1, const P_V& positionVector2, const void* routineDataPtr, Q_M& imageMatrix);
  const void* m_routineDataPtr;
};

}  // End namespace QUESO

#endif // UQ_GENERIC_MATRIX_COVARIANCE_FUNCTION_H
