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
#include <uqTrilinosMatrix.h>

template<>
uqTrilinosVectorClass*
uqVectorSpaceClass<uqTrilinosVectorClass,uqTrilinosMatrixClass>::newVector() const
{
  return new uqTrilinosVectorClass(m_env,*m_map);
}

template<>
uqTrilinosVectorClass*
uqVectorSpaceClass<uqTrilinosVectorClass,uqTrilinosMatrixClass>::newVector(double value) const
{
  return new uqTrilinosVectorClass(m_env,*m_map,value);
}

template<>
uqTrilinosMatrixClass*
uqVectorSpaceClass<uqTrilinosVectorClass,uqTrilinosMatrixClass>::newMatrix() const
{
  return new uqTrilinosMatrixClass(m_env,*m_map,this->dim());
}

template<>
uqTrilinosMatrixClass*
uqVectorSpaceClass<uqTrilinosVectorClass,uqTrilinosMatrixClass>::newDiagMatrix(double diagValue) const
{
  return new uqTrilinosMatrixClass(m_env,*m_map,diagValue);
}
