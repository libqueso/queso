/* uq/libs/queso/src/uqTrilinosFinDimLinearSpace.C
 *
 * Copyright (C) 2008 The QUESO Team, http://queso.ices.utexas.edu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#if 0
#include <uqFinDimLinearSpace.h>
#include <uqTrilinosMatrix.h>

template<>
uqTrilinosVectorClass*
uqFinDimLinearSpaceClass<uqTrilinosVectorClass,uqTrilinosMatrixClass>::newVector() const
{
  return new uqTrilinosVectorClass(m_env,this->map());
}

template<>
uqTrilinosMatrixClass*
uqFinDimLinearSpaceClass<uqTrilinosVectorClass,uqTrilinosMatrixClass>::newMatrix() const
{
  return new uqTrilinosMatrixClass(m_env,this->map(),this->dim());
}

template<>
uqTrilinosMatrixClass*
uqFinDimLinearSpaceClass<uqTrilinosVectorClass,uqTrilinosMatrixClass>::newDiagMatrix(double diagValue) const
{
  uqTrilinosVectorClass v(m_env,this->map());
  v.cwSet(diagValue);
  return new uqTrilinosMatrixClass(v);
}
#endif
