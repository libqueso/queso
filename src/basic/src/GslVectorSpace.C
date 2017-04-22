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

#include <queso/VectorSpace.h>
#include <queso/GslMatrix.h>

namespace QUESO {

template <>
Map*
VectorSpace<GslVector, GslMatrix>::newMap()
{
  return new Map(m_dimGlobal,0,m_env.selfComm());
}

template<>
GslVector*
VectorSpace<GslVector,GslMatrix>::newVector() const
{
  return new GslVector(m_env,*m_map);
}

template<>
GslVector*
VectorSpace<GslVector,GslMatrix>::newVector(double value) const
{
  return new GslVector(m_env,*m_map,value);
}

template<>
GslMatrix*
VectorSpace<GslVector,GslMatrix>::newMatrix() const
{
  return new GslMatrix(m_env,*m_map,this->dimGlobal());
}

template<>
GslMatrix*
VectorSpace<GslVector,GslMatrix>::newDiagMatrix(double diagValue) const
{
  return new GslMatrix(m_env,*m_map,diagValue);
}

}  // End namespace QUESO
