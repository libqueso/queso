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

#ifndef __QUESO_FUNCTION_BASE__
#define __QUESO_FUNCTION_BASE__

#include <string>

/*!
 * \file uqFunctionBase.h
 * \brief Abstract base class for function objects
 */

class uqFunctionBase {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor. Zero everywhere.
  uqFunctionBase();

  //   libmesh object
  //   read from a file

  //! Destructor
  virtual ~uqFunctionBase();
  //@}

  //! Save the current function to an Exodus file
  virtual void save_function(const std::string & filename) const = 0;

protected:
  // Number of degrees of freedom
  unsigned int ndofs;
};

#endif // __QUESO_FUNCTION_BASE__
