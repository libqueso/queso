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
// 
// $Id$
//
//--------------------------------------------------------------------------

#include <uqVectorSpace.h>

#ifdef QUESO_HAS_TRILINOS

#include <uqTeuchosMatrix.h>

namespace QUESO {

template <>
MapClass*
VectorSpaceClass<TeuchosVectorClass, TeuchosMatrixClass>::newMap()
{
  return new MapClass(m_dimGlobal,0,m_env.selfComm());
}

template<>
TeuchosVectorClass*
VectorSpaceClass<TeuchosVectorClass,TeuchosMatrixClass>::newVector() const
{
  return new TeuchosVectorClass(m_env,*m_map);
}

template<>
TeuchosVectorClass*
VectorSpaceClass<TeuchosVectorClass,TeuchosMatrixClass>::newVector(double value) const
{
  return new TeuchosVectorClass(m_env,*m_map,value);
}

template<>
TeuchosMatrixClass*
VectorSpaceClass<TeuchosVectorClass,TeuchosMatrixClass>::newMatrix() const
{
  return new TeuchosMatrixClass(m_env,*m_map,this->dimGlobal());
}

template<>
TeuchosMatrixClass*
VectorSpaceClass<TeuchosVectorClass,TeuchosMatrixClass>::newDiagMatrix(double diagValue) const
{
  return new TeuchosMatrixClass(m_env,*m_map,diagValue);
}

}  // End namespace QUESO

#endif // ifdef QUESO_HAS_TRILINOS
