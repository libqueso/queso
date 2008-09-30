/* uq/libs/mcmc/src/uqGslVectorSpace.C
 *
 * Copyright (C) 2008 The PECOS Team, http://queso.ices.utexas.edu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <uqVectorSpace.h>
#include <uqGslMatrix.h>

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
