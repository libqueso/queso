//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2015 The PECOS Development Team
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

#ifndef UQ_SCALAR_FUNCTION_H
#define UQ_SCALAR_FUNCTION_H

#include <queso/VectorSet.h>
#include <queso/VectorSubset.h>
#include <queso/Environment.h>
#include <queso/Defines.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*! \file ScalarFunction.h
 * \brief Set of classes for handling vector functions.
 *
 * \class BaseScalarFunction
 * \brief A templated (base) class for handling scalar functions.
 *
 * This class allows the mathematical definition of a scalar function such as:
 * \f$ f: B \subset R \rightarrow R \f$. A function of one or more variables
 * has always one-dimensional range.  PDFs (marginal, joint) and CDFs are examples
 * of scalar functions.
 */

template <class V = GslVector, class M = GslMatrix>
class BaseScalarFunction {
public:
  //! @name Constructor/Destructor methods.
  //@{
  //! Default constructor.
  /*!
   * Instantiates an object of the class, i.e. a scalar function, given a
   * prefix and its domain.
   */
  BaseScalarFunction(const char * prefix, const VectorSet<V, M> & domainSet);

  //! Destructor
  virtual ~BaseScalarFunction();
  //@}

  //! @name Mathematical methods.
  //@{
  //! Access to the protected attribute \c m_domainSet: domain set of the scalar function.
  const VectorSet<V, M> & domainSet() const;

  //! Actual value of the scalar function.
  virtual double actualValue(const V & domainVector, const V * domainDirection,
      V * gradVector, M * hessianMatrix, V * hessianEffect) const = 0;

  //! Logarithm of the value of the scalar function.
  virtual double lnValue(const V & domainVector, const V * domainDirection,
      V * gradVector, M * hessianMatrix, V * hessianEffect) const = 0;
  //@}
protected:
  const BaseEnvironment & m_env;
  std::string m_prefix;

  //! Domain set of the scalar function.
  const VectorSet<V, M> & m_domainSet;
};

}  // End namespace QUESO

#endif // UQ_SCALAR_FUNCTION_H
