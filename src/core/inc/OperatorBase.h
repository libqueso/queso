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

#ifndef QUESO_OPERATOR_BASE_H
#define QUESO_OPERATOR_BASE_H

#include <string>
#include <vector>
#include <queso/SharedPtr.h>
#include <queso/FunctionBase.h>

namespace QUESO {

/*!
 * \file OperatorBase.h
 * \brief Abstract base class for operator objects
 *
 * \class OperatorBase
 * \brief Abstract base class for operator objects.  Operators are assumed to be symmetric and positive-definite.
 */

class FunctionOperatorBuilder;

class OperatorBase {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Constructor
  OperatorBase();

  //! Destructor
  virtual ~OperatorBase();
  //@}

  //! Return eigenvalue \c i.
  /*!
   * You can store them however you want, but having some kind of order to them
   * is useful for \c InfiniteDimensionalMeasure
   */
  virtual double get_eigenvalue(unsigned int i) const = 0;

  //! Return the reciprocal of eigenvalue \c i.
  virtual double get_inverted_eigenvalue(unsigned int i) const = 0;

  //! Return the number of converged eigenpairs
  virtual unsigned int get_num_converged() const = 0;

  //! Given coefficients \c xi, computes the Karhunen-Loeve transform
  /*!
   *  This transform goes from coefficient space to physical space using
   *  \c this as the precision operator:
   *  \sum_k \xi_k / pow(\lambda_k, \c alpha / 2.0) \phi_k(x)
   *  where the lambda are eigenvalues of the precision operator, \c this, and
   *  the \phi(x) are eigenfunctions of the precision operator, \c this
   */
  virtual SharedPtr<FunctionBase>::Type
  inverse_kl_transform(std::vector<double>& xi, double alpha) const = 0;
};

}  // End namespace QUESO

#endif // QUESO_OPERATOR_BASE_H
