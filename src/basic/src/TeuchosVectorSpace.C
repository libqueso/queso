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

#ifdef QUESO_HAS_TRILINOS

#include <queso/TeuchosMatrix.h>

namespace QUESO {

template <>
Map*
VectorSpace<TeuchosVector, TeuchosMatrix>::newMap()
{
  return new Map(m_dimGlobal,0,m_env.selfComm());
}

template<>
TeuchosVector*
VectorSpace<TeuchosVector,TeuchosMatrix>::newVector() const
{
  return new TeuchosVector(m_env,*m_map);
}

template<>
TeuchosVector*
VectorSpace<TeuchosVector,TeuchosMatrix>::newVector(double value) const
{
  return new TeuchosVector(m_env,*m_map,value);
}

template<>
TeuchosMatrix*
VectorSpace<TeuchosVector,TeuchosMatrix>::newMatrix() const
{
  return new TeuchosMatrix(m_env,*m_map,this->dimGlobal());
}

template<>
TeuchosMatrix*
VectorSpace<TeuchosVector,TeuchosMatrix>::newDiagMatrix(double diagValue) const
{
  return new TeuchosMatrix(m_env,*m_map,diagValue);
}

}  // End namespace QUESO

#endif // ifdef QUESO_HAS_TRILINOS
