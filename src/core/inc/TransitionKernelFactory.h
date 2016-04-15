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

#ifndef QUESO_TK_FACTORY_H
#define QUESO_TK_FACTORY_H

#include <queso/Factory.h>

namespace QUESO
{

/**
 * TransitionKernelFactory class defintion.
 */
template <class Base>
class TransitionKernelFactory : public Factory<Base>
{
public:
  /**
   * Constructor. Takes the name to be mapped.
   */
  TransitionKernelFactory(const std::string & name)
    : Factory<Base>(name)
  {}

  /**
   * Destructor. (Empty.)
   */
  virtual ~TransitionKernelFactory() {}

  /**
   * Create a Base class.  Force this to be implemented
   * later.
   */
  virtual typename SharedPtr<Base>::Type create() = 0;
};

} // namespace QUESO

#endif // QUESO_TK_FACTORY_H
