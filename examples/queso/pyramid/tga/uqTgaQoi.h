/* uq/examples/queso/pyramid/uqTgaQoi.h
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

#ifndef __UQ_TGA_QOI_H__
#define __UQ_TGA_QOI_H__

#include <uqTgaDefines.h>
#include <uqTgaW.h>
#include <uqDefines.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>

//********************************************************
// QoI function object for both forward problems of the validation cycle.
// A QoI function object is provided by user and is called by the UQ library.
// This QoI function object consists of data and routine.
//********************************************************
// The (user defined) data class that carries the data needed by the (user defined) qoi routine
template<class P_V,class P_M,class Q_V, class Q_M>
struct
uqTgaQoiInfoStruct
{
  uqTgaQoiInfoStruct(const uqVectorSpaceClass<P_V,P_M>& paramSpace,
                     double                             initialTemp,
                     double                             beta,
                     double                             criticalW,
                     double                             criticalTime);
 ~uqTgaQoiInfoStruct();

  const uqVectorSpaceClass<P_V,P_M>& m_paramSpace;
  uqBase1D1DFunctionClass*           m_temperatureFunctionObj;
  double                             m_criticalW;
  double                             m_criticalTime;
  uqTgaWClass<P_V,P_M>*              m_wObj;
};

template<class P_V,class P_M,class Q_V,class Q_M>
uqTgaQoiInfoStruct<P_V,P_M,Q_V,Q_M>::uqTgaQoiInfoStruct(
  const uqVectorSpaceClass<P_V,P_M>& paramSpace,
  double                             initialTemp,
  double                             beta,
  double                             criticalW,
  double                             criticalTime)
  :
  m_paramSpace            (paramSpace),
  m_temperatureFunctionObj(new uqLinear1D1DFunctionClass(-INFINITY,INFINITY,0.,initialTemp,beta)),
  m_criticalW             (criticalW),
  m_criticalTime          (criticalTime),
  m_wObj                  (new uqTgaWClass<P_V,P_M>(m_paramSpace,*m_temperatureFunctionObj))
{
}

template<class P_V,class P_M,class Q_V,class Q_M>
uqTgaQoiInfoStruct<P_V,P_M,Q_V,Q_M>::~uqTgaQoiInfoStruct()
{
  delete m_wObj;
  delete m_temperatureFunctionObj;
}

// The actual (user defined) qoi routine
template<class P_V,class P_M,class Q_V,class Q_M>
void uqTgaQoiRoutine(const P_V& paramValues, const void* functionDataPtr, Q_V& qoiValues)
{
  const uqTgaQoiInfoStruct<P_V,P_M,Q_V,Q_M>& info = *((const uqTgaQoiInfoStruct<P_V,P_M,Q_V,Q_M> *) functionDataPtr);

  uqTgaWClass<P_V,P_M>& wObj = *(info.m_wObj);
                       
#ifdef QUESO_TGA_USES_OLD_COMPATIBLE_CODE
  wObj.computeUsingTemp(paramValues,
                        info.m_criticalTime*info.m_temperatureFunctionObj->deriv(0.), // Should add initialTemp --> COMPATIBILITY WITH OLD VERSION
                        NULL,  // referenceW
                        NULL);
#else
  wObj.compute(paramValues,
               info.m_criticalTime,
               0.,
               false, // computeGradAlso
               NULL,  // referenceW
               NULL,  // weightFunction
               NULL,  // misfitValue
               NULL); // diffFunction
#endif

  unsigned int tmpSize = wObj.ws().size();
  qoiValues[0] = wObj.ws()[tmpSize-1]; // QoI = mass fraction remaining at critical time

  return;
}

#endif // __UQ_TGA_QOI_H__
