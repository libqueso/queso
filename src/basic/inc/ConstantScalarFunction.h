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

#ifndef UQ_CONSTANT_SCALAR_FUNCTION_H
#define UQ_CONSTANT_SCALAR_FUNCTION_H

#include <queso/ScalarFunction.h>
#include <queso/VectorSet.h>
#include <queso/VectorSubset.h>
#include <queso/Environment.h>
#include <queso/Defines.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*!\class ConstantScalarFunction
 * \brief A class for handling scalar functions which image is a constant (real number).
 *
 * This class allows the mathematical definition of a scalar function which image set
 * is a constant (real number). */

template <class V = GslVector, class M = GslMatrix>
class ConstantScalarFunction : public BaseScalarFunction<V,M> {
public:
    //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  /*! Instantiates an object of the class, i.e. a scalar function, given a prefix, its domain and constant-valued image.*/
  ConstantScalarFunction(const char*                  prefix,
                                const VectorSet<V,M>& domainSet,
                                double                       constantValue);
  //! Virtual destructor
  virtual ~ConstantScalarFunction();

  //! @name Mathematical method
  //@{
  //! Calculates the actual value of this scalar function.
  double actualValue      (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;

  //! Calculates the logarithm of the value of this scalar function (which is zero).
  double lnValue          (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  //@}
protected:
  using BaseScalarFunction<V,M>::m_env;
  using BaseScalarFunction<V,M>::m_prefix;
  using BaseScalarFunction<V,M>::m_domainSet;

  //! Constant value is the image set of this scalar function.
  double m_constantValue;
};

}  // End namespace QUESO

#endif // UQ_CONSTANT_SCALAR_FUNCTION_H
