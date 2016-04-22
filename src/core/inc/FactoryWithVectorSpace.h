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

#ifndef QUESO_FACTORYWITHVECTORSPACE_H
#define QUESO_FACTORYWITHVECTORSPACE_H

#include <queso/Factory.h>
#include <queso/VectorSpace.h>

namespace QUESO
{

/**
 * FactoryWithVectorSpace class defintion.
 */
template <class Base>
class FactoryWithVectorSpace : public Factory<Base>
{
public:
  FactoryWithVectorSpace(const std::string & name)
    : Factory<Base>(name)
  {}

  virtual ~FactoryWithVectorSpace() {}

  static void set_vectorspace(const VectorSpace<GslVector, GslMatrix> & v)
  {
    m_vectorSpace = &v;
  }

protected:
  // We are not taking ownership of this, so we don't want it to be deleted.
  static const VectorSpace<GslVector, GslMatrix> * m_vectorSpace;
};

} // namespace QUESO

#endif // QUESO_FACTORYWITHVECTORSPACE_H
