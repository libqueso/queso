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

#ifndef __UQ_NORMAL_LIKELIHOOD_H__
#define __UQ_NORMAL_LIKELIHOOD_H__

//********************************************************
// The likelihood routine: provided by user and called by QUESO
//********************************************************
template<class P_V,class P_M>
struct
likelihoodRoutine_DataType
{
  const P_V* paramMeans;
  const P_M* matrix;
  bool       applyMatrixInvert;
};

template<class P_V,class P_M>
double
likelihoodRoutine(
  const P_V&  paramValues,
  const P_V*  paramDirection,
  const void* functionDataPtr,
  P_V*        gradVector,
  P_M*        hessianMatrix,
  P_V*        hessianEffect)
{
  double result = 0.;

  const P_V& paramMeans        = *((likelihoodRoutine_DataType<P_V,P_M> *) functionDataPtr)->paramMeans;
  const P_M& matrix            = *((likelihoodRoutine_DataType<P_V,P_M> *) functionDataPtr)->matrix;
  bool       applyMatrixInvert =  ((likelihoodRoutine_DataType<P_V,P_M> *) functionDataPtr)->applyMatrixInvert;

  P_V diffVec(paramValues - paramMeans);
  if (applyMatrixInvert) {
    result = scalarProduct(diffVec, matrix.invertMultiply(diffVec));
  }
  else {
    result = scalarProduct(diffVec, matrix * diffVec);
  }

  return result;
}
#endif // __UQ_NORMAL_LIKELIHOOD_H__
