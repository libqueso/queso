/* uq/libs/mcmc/inc/uqVectorFunction.h
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

#ifndef __UQ_VECTOR_FUNCTION_H__
#define __UQ_VECTOR_FUNCTION_H__

#include <uqEnvironment.h>
#include <uqDefines.h>

//*****************************************************
// Base class
//*****************************************************
template<class P_V,class P_M,class Q_V,class Q_M>
class uqVectorFunctionClass {
public:
  uqVectorFunctionClass(void (*routinePtr)(const P_V& domainVector, const void* functionDataPtr, Q_V& imageVector),
                        const void* functionDataPtr);
  virtual ~uqVectorFunctionClass();
  virtual void compute(const P_V& domainVector, Q_V& imageVector) const;

protected:
  void (*m_routinePtr)(const P_V& domainVector, const void* functionDataPtr, Q_V& imageVector);
  const void* m_routineDataPtr;
};

template<class P_V,class P_M,class Q_V,class Q_M>
uqVectorFunctionClass<P_V,P_M,Q_V,Q_M>::uqVectorFunctionClass(
  void (*routinePtr)(const P_V& domainVector, const void* functionDataPtr, Q_V& imageVector),
  const void* functionDataPtr)
  :
  m_routinePtr(routinePtr),
  m_routineDataPtr(functionDataPtr)
{
}

template<class P_V,class P_M,class Q_V,class Q_M>
uqVectorFunctionClass<P_V,P_M,Q_V,Q_M>::~uqVectorFunctionClass()
{
}

template<class P_V,class P_M,class Q_V,class Q_M>
void
uqVectorFunctionClass<P_V,P_M,Q_V,Q_M>::compute(const P_V& domainVector, Q_V& imageVector) const
{
  //UQ_FATAL_TEST_MACRO(false,
  //                    domainVector.env().rank(),
  //                    "uqVectorFunctionClass<P_V,P_M,Q_V,Q_M>::compute()",
  //                    "this method should not be called in the case of this class");

  m_routinePtr(domainVector, m_routineDataPtr, imageVector);

  return;
}

#endif // __UQ_VECTOR_FUNCTION_H__
