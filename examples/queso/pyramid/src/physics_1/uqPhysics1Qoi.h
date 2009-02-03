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

#ifndef __UQ_TGA_QOI_H__
#define __UQ_TGA_QOI_H__

#include <uqDefines.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>

// The (user defined) data class that carries the data needed by the (user defined) qoi routine
template<class P_V,class P_M,class Q_V, class Q_M>
struct
uqPhysics1QoiInfoStruct
{
  uqPhysics1QoiInfoStruct(const uqVectorSpaceClass<P_V,P_M>& paramSpace,
                          double                             initialTemp,
                          double                             beta,
                          double                             criticalW,
                          double                             criticalTime);
 ~uqPhysics1QoiInfoStruct();

  const uqVectorSpaceClass<P_V,P_M>& m_paramSpace;
  uqBase1D1DFunctionClass*           m_temperatureFunctionObj;
  double                             m_criticalW;
  double                             m_criticalTime;
};

template<class P_V,class P_M,class Q_V,class Q_M>
uqPhysics1QoiInfoStruct<P_V,P_M,Q_V,Q_M>::uqPhysics1QoiInfoStruct(
  const uqVectorSpaceClass<P_V,P_M>& paramSpace,
  double                             initialTemp,
  double                             beta,
  double                             criticalW,
  double                             criticalTime)
  :
  m_paramSpace            (paramSpace),
  m_temperatureFunctionObj(new uqLinear1D1DFunctionClass(-INFINITY,INFINITY,0.,initialTemp,beta)),
  m_criticalW             (criticalW),
  m_criticalTime          (criticalTime)
{
}

template<class P_V,class P_M,class Q_V,class Q_M>
uqPhysics1QoiInfoStruct<P_V,P_M,Q_V,Q_M>::~uqPhysics1QoiInfoStruct()
{
  delete m_temperatureFunctionObj;
}

// The actual (user defined) qoi routine
template<class P_V,class P_M,class Q_V,class Q_M>
void uqPhysics1QoiRoutine(const P_V& paramValues, const void* functionDataPtr, Q_V& qoiValues)
{
  const uqPhysics1QoiInfoStruct<P_V,P_M,Q_V,Q_M>& info = *((const uqPhysics1QoiInfoStruct<P_V,P_M,Q_V,Q_M> *) functionDataPtr);

  qoiValues[0] = 0.*info.m_criticalW;

  return;
}

#endif // __UQ_TGA_QOI_H__
