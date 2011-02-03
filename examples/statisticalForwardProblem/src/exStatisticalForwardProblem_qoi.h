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

#ifndef __EX_STATISTICAL_FORWARD_PROBLEM_QOI_H__
#define __EX_STATISTICAL_FORWARD_PROBLEM_QOI_H__

#include <uqEnvironment.h>
#include <uqDistArray.h>

//********************************************************
// The qoi routine: provided by user and called by QUESO
//********************************************************
template<class P_V,class P_M, class Q_V, class Q_M>
struct
qoiRoutine_DataType
{
  double p1MultiplicativeFactor;
  double p1ExponentFactor;
  double p2MultiplicativeFactor;
  double p2ExponentFactor;
};

template<class P_V,class P_M, class Q_V, class Q_M>
void
qoiRoutine(
  const P_V&                                   paramValues,
  const P_V*                                   paramDirection,
  const void*                                  functionDataPtr,
        Q_V&                                   qoiValues,
        typename uqDistArrayClass<P_V*>::type* gradVectors,
        typename uqDistArrayClass<P_M*>::type* hessianMatrices,
        typename uqDistArrayClass<P_V*>::type* hessianEffects)
{
  //double a1 = ((qoiRoutine_DataType<P_V,P_M,Q_V,Q_M> *) functionDataPtr)->p1MultiplicativeFactor;
  //double e1 = ((qoiRoutine_DataType<P_V,P_M,Q_V,Q_M> *) functionDataPtr)->p1ExponentFactor;
  //double a2 = ((qoiRoutine_DataType<P_V,P_M,Q_V,Q_M> *) functionDataPtr)->p2MultiplicativeFactor;
  //double e2 = ((qoiRoutine_DataType<P_V,P_M,Q_V,Q_M> *) functionDataPtr)->p2ExponentFactor;

  double p1 = paramValues[0];
  double p2 = paramValues[1];

  //qoiValues[0] = a1*pow(p1,e1) + a2*pow(p2,e2);
  qoiValues[0] = p1+p2;

  //std::cout << "Leaving qoiRoutine()"
  //          << ": a1 = " << a1
  //          << ", e1 = " << e1
  //          << ", a2 = " << a2
  //          << ", e2 = " << e2
  //          << ", p1 = " << p1
  //          << ", p2 = " << p2
  //          << ", qoi0 = " << qoiValues[0]
  //          << std::endl;

  return;
}

#endif // __EX_STATISTICAL_FORWARD_PROBLEM_QOI_H__
