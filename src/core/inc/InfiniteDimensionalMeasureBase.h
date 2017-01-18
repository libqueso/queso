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

#ifndef QUESO_INFINITEDIMENSIONALMEASURE_BASE_H
#define QUESO_INFINITEDIMENSIONALMEASURE_BASE_H

#include <queso/SharedPtr.h>
#include <queso/FunctionBase.h>

namespace QUESO {

/*!
 * \file InfiniteDimensionalMeasureBase.h
 * \brief Abstract base class for infinite dimensional measures
 *
 * \class InfiniteDimensionalMeasureBase
 * \brief Abstract base class for infinite dimensional measures
 */

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

  //! Draw from the measure, and then return a shared pointer to the draw
  virtual SharedPtr<FunctionBase>::Type draw() = 0;

  //! Return coefficient \c i of the KL expansion of the current draw.  Must be called after draw()
  virtual double get_kl_coefficient(unsigned int i) const = 0;
};

}  // End namespace QUESO

#endif // QUESO_INFINITEDIMENSIONALMEASURE_BASE_H
