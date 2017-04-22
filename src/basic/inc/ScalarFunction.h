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

#ifndef UQ_SCALAR_FUNCTION_H
#define UQ_SCALAR_FUNCTION_H

#include <queso/Environment.h>
#include <queso/VectorSet.h>
#include <queso/VectorSubset.h>
#include <queso/ScopedPtr.h>

namespace QUESO {

class GslVector;
class GslMatrix;
class BoostInputOptionsParser;

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
  //@}

  //! @name Evaluation methods
  //@{
  //! Logarithm of the value of the scalar function.  Deprecated.
  /*!
   * Pointers will be NULL if derivative information is not required by QUESO.
   *
   * Default implementation throws an exception.
   */
  virtual double lnValue(const V & domainVector,
                         const V * domainDirection,
                         V * gradVector,
                         M * hessianMatrix,
                         V * hessianEffect) const;

  //! Returns the logarithm of the function at \c domainVector
  /*!
   * Default implementation calls the deprecated method.
   *
   * QUESO calls this method when it needs to evaluate the function but doesn't
   * need derivative information about that function.
   */
  virtual double lnValue(const V & domainVector) const;

  //! Returns the logarithm of the function and its gradient at \c domainVector
  /*!
   * Default implementation calls above method successively to fill up
   * \c gradVector with a finite difference approximation.
   *
   * Note, gradVector should be filled with the gradient of the logarithm of
   * the function, not the gradient of the function.
   *
   * QUESO calls this method when it needs to evaluate the function, needs
   * first order derivative information, but doesn't need second order
   * derivative information.
   */
  virtual double lnValue(const V & domainVector, V & gradVector) const;

  //! Returns the logarithm of the function \c domainVector, the gradient of
  //! the logarithm at \c domainVector, and the effect of the hessian of the
  //! logarithm in the direction at \c domainDirection
  /*!
   * Default implementation throws an exception because we don't expect the
   * user to tolerate a finite difference approximation of the hessian.
   *
   * Note, gradVector should be filled with the gradient of the logarithm of
   * the function, not the gradient of the function.
   *
   * The 'hessian' referred to here should be the hessian of the logarithm
   * of the function at \c domainVector.
   *
   * QUESO calls this method when it needs to evaluate the function, needs
   * first order derivative information, and also needs second-order deriative
   * information
   */
  virtual double lnValue(const V & domainVector,
                         V & gradVector,
                         const V & domainDirection,
                         V & hessianEffect) const;

  //! Actual value of the scalar function.
  virtual double actualValue(const V & domainVector, const V * domainDirection,
      V * gradVector, M * hessianMatrix, V * hessianEffect) const = 0;
  //@}

  //! Sets the step size for finite differencing gradients
  /*!
   * If the function is multi-dimensional then the same finite difference step
   * size is used in every direction.
   */
  void setFiniteDifferenceStepSize(double fdStepSize);

  //! Sets the step size for the i-th component of the finite differencing
  //! vector
  /*!
   * i is a zero-based index.
   *
   * If the function is one-dimensional, the only allowed value for i is 0.
   */
  void setFiniteDifferenceStepSize(unsigned int i, double fdStepSize);

protected:
  const BaseEnvironment & m_env;
  std::string m_prefix;

  //! Domain set of the scalar function.
  const VectorSet<V, M> & m_domainSet;

private:
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  //! Input parser
  ScopedPtr<BoostInputOptionsParser>::Type m_parser;
#endif

  //! Finite different step size
  std::vector<double> m_fdStepSize;
};

}  // End namespace QUESO

#endif // UQ_SCALAR_FUNCTION_H
