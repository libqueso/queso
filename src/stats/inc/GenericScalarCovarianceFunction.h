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

#ifndef UQ_GENERIC_SCALAR_COVARIANCE_FUNCTION_H
#define UQ_GENERIC_SCALAR_COVARIANCE_FUNCTION_H

#include <queso/ScalarCovarianceFunction.h>
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
 * \class GenericScalarCovarianceFunction
 * \brief A class for generic covariances.
 *
 * This class implements a generic covariance functions, by calling a routine (via pointer).*/

template <class V = GslVector, class M = GslMatrix>
class GenericScalarCovarianceFunction : public BaseScalarCovarianceFunction<V,M> {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  /*! Instantiates an object of the class given a prefix, the domain set, the pointer to the routine. */
  GenericScalarCovarianceFunction(const char*                  prefix,
           const VectorSet<V,M>& domainSet,
           double (*covRoutinePtr)(const V& positionVector1, const V& positionVector2, const void* routineDataPtr),
           const void*                  routinesDataPtr);

  //! Virtual destructor
  virtual ~GenericScalarCovarianceFunction();
  //@}
    //! @name Math methods
  //@{
  //! Calculates the value of the generic covariance function.
  /*! This function accesses the routine that calculates the covariance function. */
  double value(const V& positionVector1, const V& positionVector2) const;

protected:
  using BaseScalarCovarianceFunction<V,M>::m_env;
  using BaseScalarCovarianceFunction<V,M>::m_prefix;
  using BaseScalarCovarianceFunction<V,M>::m_basicDomainSet;

  double (*m_covRoutinePtr)(const V& positionVector1, const V& positionVector2, const void* routineDataPtr);
  const void* m_routineDataPtr;
};

}  // End namespace QUESO

#endif // UQ_GENERIC_SCALAR_COVARIANCE_FUNCTION_H
