//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
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

#ifndef UQ_GSL_OPTIMIZER_H
#define UQ_GSL_OPTIMIZER_H

#include <queso/Optimizer.h>

namespace QUESO {

/*!
 * \file GslOptimizer.h
 * \brief Class for handling optimization of scalar functions
 *
 * \class GslOptimizer
 * \brief A base class for handling optimisation of scalar functions
 *
 * WRITE DOCS HERE
 */

class Vector;
class GslVector;
class GslMatrix;

template <class V, class M>
class BaseScalarFunction;

class GslOptimizer : public BaseOptimizer {
public:
  //! Constructs an object that will maximize a scalar function
  /*!
   * The function \c objectiveFunction is the function that will be maximized.
   */
  GslOptimizer(
      const BaseScalarFunction<GslVector, GslMatrix> & objectiveFunction);

  //! Destructor
  virtual ~GslOptimizer();

  //! Minimize the objective function, starting at \c m_initialPoint
  /*!
   * m_initialPoint is handled in the derived class
   */
  virtual void minimize();

  //! Returns the objective function
  const BaseScalarFunction<GslVector, GslMatrix> & objectiveFunction() const;

private:
  const BaseScalarFunction<GslVector, GslMatrix> & m_objectiveFunction;
};

}  // End namespace QUESO

#endif // UQ_GSL_OPTIMIZER_H
