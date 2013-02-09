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
 * \brief Abstract base class for operator objects
 */

#include <string>

class uqOperatorBase {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor
  uqOperatorBase();

  //! Construct using the mesh from file \c filename
  uqOperatorBase(const std::string& filename);

  //! Destructor
  ~uqOperatorBase();
  //@}
};

#endif // __QUESO_OPERATOR_BASE__
