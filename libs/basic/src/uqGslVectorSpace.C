/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of the QUESO Library (Quantification of Uncertainty
 * for Estimation, Simulation and Optimization).
 *
 * QUESO is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QUESO is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QUESO. If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * $Id$
 *
 * Brief description of this file: 
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#include <uqVectorSpace.h>
#include <uqGslMatrix.h>

template <>
Epetra_Map*
uqVectorSpaceClass<class uqGslVectorClass, class uqGslMatrixClass>::newMap()
{
  return new Epetra_Map(m_dim,0,m_env.selfComm());
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
  return new uqGslMatrixClass(m_env,*m_map,this->dim());
}

template<>
uqGslMatrixClass*
uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass>::newDiagMatrix(double diagValue) const
{
  return new uqGslMatrixClass(m_env,*m_map,diagValue);
}
