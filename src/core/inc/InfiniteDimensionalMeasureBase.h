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

#ifndef __QUESO_INFINITEDIMENSIONALMEASURE_BASE__
#define __QUESO_INFINITEDIMENSIONALMEASURE_BASE__

#include <boost/shared_ptr.hpp>
#include <queso/FunctionBase.h>

/*!
 * \file InfiniteDimensionalMeasureBase.h
 * \brief Abstract base class for infinite dimensional measures
 */

namespace QUESO {

class FunctionBase;

class InfiniteDimensionalMeasureBase {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor
  InfiniteDimensionalMeasureBase();

  //! Destructor
  virtual ~InfiniteDimensionalMeasureBase();
  //@}

  //! Draw from the measure, and store the result in the public member varaible. This
  //! updates the public memeber variable current draw
  virtual boost::shared_ptr<FunctionBase> draw() const = 0;

  //! Return the coefficient \c i of the KL expansion of the current draw
  /*!
   * You need to make a draw before you call this
   */
  virtual double get_kl_coefficient(unsigned int i) const = 0;
};

}  // End namespace QUESO

#endif // __QUESO_INFINITEDIMENSIONALMEASURE_BASE__
