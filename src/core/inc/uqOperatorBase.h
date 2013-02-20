//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012 The PECOS Development Team
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
// 
// $Id: $
//
//--------------------------------------------------------------------------

#ifndef __QUESO_OPERATOR_BASE__
#define __QUESO_OPERATOR_BASE__

/*!
 * \file uqOperatorBase.h
 * \brief Abstract base class for operator objects. Operators are assumed to
 *        be symmetric and positive-definite, so the eigenvalues are real and
 *        positive.
 */

#include <string>
#include <vector>
#include <memory>
#include <uqFunctionBase.h>

class uqFunctionOperatorBuilder;

class uqOperatorBase {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Construct with a \c builder
  /*!
   * A \c builder object is just one that a FEM library backend can use to set
   * up various options. Polynomial type, polynomial order, and the number of
   * eigenpairs to request are good examples.
   */
  uqOperatorBase(const uqFunctionOperatorBuilder & builder);

  //! Destructor
  virtual ~uqOperatorBase();
  //@}

  //! Return eigenvalue \c i.
  /*!
   * You can store them however you want, but having some kind of order to them
   * is useful for \c uqInfiniteDimensionalMeasure
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
  virtual std::auto_ptr<uqFunctionBase>
  inverse_kl_transform(std::vector<double>& xi, double alpha) const = 0;

protected:
  //! A reference to the builder object
  const uqFunctionOperatorBuilder & builder;
};

#endif // __QUESO_OPERATOR_BASE__
