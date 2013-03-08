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
// $Id:$
//
//--------------------------------------------------------------------------

#include <uqVectorSpace.h>
#include <uqGslMatrix.h>

template <>
uqMapClass*
uqVectorSpaceClass<uqGslVectorClass, uqGslMatrixClass>::newMap()
{
  return new uqMapClass(m_dimGlobal,0,m_env.selfComm());
}

template<>
uqGslVectorClass*
uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass>::newVector() const
{
  return new uqGslVectorClass(m_env,*m_map);
}

template<>
uqGslVectorClass*
uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass>::newVector(double value) const
{
  return new uqGslVectorClass(m_env,*m_map,value);
}

template<>
uqGslMatrixClass*
uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass>::newMatrix() const
{
  return new uqGslMatrixClass(m_env,*m_map,this->dimGlobal());
}

template<>
uqGslMatrixClass*
uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass>::newDiagMatrix(double diagValue) const
{
  return new uqGslMatrixClass(m_env,*m_map,diagValue);
}
