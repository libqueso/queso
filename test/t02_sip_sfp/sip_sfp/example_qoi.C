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

#include <example_qoi.h>

void
qoiRoutine(
  const uqGslVectorClass&                        paramValues,
  const uqGslVectorClass*                        paramDirection,
  const void*                                    functionDataPtr,
        uqGslVectorClass&                        qoiValues,
        EpetraExt::DistArray<uqGslVectorClass*>* gradVectors,
        EpetraExt::DistArray<uqGslMatrixClass*>* hessianMatrices,
        EpetraExt::DistArray<uqGslVectorClass*>* hessianEffects)
{
  // Logic just to avoid warnings from INTEL compiler
  const uqGslVectorClass* aux1 = paramDirection;
  if (aux1) {};
  EpetraExt::DistArray<uqGslVectorClass*>* aux2 = gradVectors;
  if (aux2) {};
  aux2 = hessianEffects;
  EpetraExt::DistArray<uqGslMatrixClass*>* aux3 = hessianMatrices;
  if (aux3) {};

  // Actual code
  double coef1 = ((qoiRoutine_DataType *) functionDataPtr)->coef1;
  double coef2 = ((qoiRoutine_DataType *) functionDataPtr)->coef2;

  qoiValues[0] = (coef1*paramValues[0] + coef2*paramValues[1]);

  return;
}
