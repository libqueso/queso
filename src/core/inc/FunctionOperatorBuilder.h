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

#ifndef QUESO_FUNCTIONOPERATORBUILDER_BASE_H
#define QUESO_FUNCTIONOPERATORBUILDER_BASE_H

#include <string>

namespace QUESO {

/*!
 * \file FunctionOperatorBuilder.h
 * \brief Helper class for function and operator objects. This class is meant
 *        to hold common FEM library backend options in a library-agnostic
 *        fashion
 */

class FunctionOperatorBuilder {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.  Does nothing.
  FunctionOperatorBuilder();

  //! Destructor.  Does nothing.
  ~FunctionOperatorBuilder();
  //@}

  //! String to store the polynomial family to use. Default is "LAGRANGE".
  std::string family;

  //! String to store the polynomial order to use. Default is "FIRST".
  std::string order;

  //! Number of eigenpairs to request when building an operator. Default is 0.
  unsigned int num_req_eigenpairs;
};

}  // End namespace QUESO

#endif // QUESO_FUNCTIONOPERATORBUILDER_BASE_H
