/* uq/libs/mcmc/inc/uqQoIFunction.h
 *
 * Copyright (C) 2008 The PECOS Team, http://www.ices.utexas.edu/centers/pecos
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

#ifndef __UQ_QOI_FUNCTION_H__
#define __UQ_QOI_FUNCTION_H__

#include <uqEnvironment.h>
#include <uqDefines.h>

//*****************************************************
// Class to accomodate the qoi function.
// User must provide a qoi() routine.
//*****************************************************

//*****************************************************
// Base class
//*****************************************************
template<class P_V,class P_M,class Q_V,class Q_M>
class uqQoIFunction_BaseClass {
public:
  uqQoIFunction_BaseClass(void (*routinePtr)(const P_V& paramValues, const void* functionDataPtr, Q_V& qoiValues),
                          const void* functionDataPtr);
  virtual ~uqQoIFunction_BaseClass();
  virtual void computeQoIs(const P_V& paramValues, Q_V& qoiValues) const;

protected:
  void (*m_routinePtr)(const P_V& paramValues, const void* functionDataPtr, Q_V& qoiValues);
  const void* m_routineDataPtr;
};

template<class P_V,class P_M,class Q_V,class Q_M>
uqQoIFunction_BaseClass<P_V,P_M,Q_V,Q_M>::uqQoIFunction_BaseClass(
  void (*routinePtr)(const P_V& paramValues, const void* functionDataPtr, Q_V& qoiValues),
  const void* functionDataPtr)
  :
  m_routinePtr(routinePtr),
  m_routineDataPtr(functionDataPtr)
{
}

template<class P_V,class P_M,class Q_V,class Q_M>
uqQoIFunction_BaseClass<P_V,P_M,Q_V,Q_M>::~uqQoIFunction_BaseClass()
{
}

template<class P_V,class P_M,class Q_V,class Q_M>
void
uqQoIFunction_BaseClass<P_V,P_M,Q_V,Q_M>::computeQoIs(const P_V& paramValues, Q_V& qoiValues) const
{
  //UQ_FATAL_TEST_MACRO(false,
  //                    paramValues.env().rank(),
  //                    "uqQoIFunction_BaseClass<P_V,P_M,Q_V,Q_M>::computeQoIs()",
  //                    "this method should not be called in the case of this class");

  m_routinePtr(paramValues, m_routineDataPtr, qoiValues);

  return;
}

#endif // __UQ_QOI_FUNCTION_H__
