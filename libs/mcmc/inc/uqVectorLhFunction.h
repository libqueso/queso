/* uq/libs/mcmc/inc/uqVectorLhFunction.h
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

#ifndef __UQ_VECTOR_LH_FUNCTION_H__
#define __UQ_VECTOR_LH_FUNCTION_H__

#include <uqEnvironment.h>
#include <uqDefines.h>

//*****************************************************
// Class to accomodate the likehood function.
// User must provide a likelihood() routine.
//*****************************************************

//*****************************************************
// Base class
//*****************************************************
template<class P_V,class P_M,class L_V,class L_M>
class uqVectorLhFunction_BaseClass {
public:
  uqVectorLhFunction_BaseClass(void (*routinePtr)(const P_V& paramValues, const void* functionDataPtr, L_V& likelihoodValues),
                                  const void* functionDataPtr);
  virtual ~uqVectorLhFunction_BaseClass();
  virtual void computeActualLikelihoods  (const P_V& paramValues, L_V& likelihoodValues) const = 0;
  virtual void computeMinus2LnLikelihoods(const P_V& paramValues, L_V& likelihoodValues) const = 0;
  virtual void computeMisfits            (const P_V& paramValues, L_V& likelihoodValues) const = 0;

protected:
  void (*m_routinePtr)(const P_V& paramValues, const void* functionDataPtr, L_V& likelihoodValues);
  const void* m_routineDataPtr;
};

template<class P_V,class P_M,class L_V,class L_M>
uqVectorLhFunction_BaseClass<P_V,P_M,L_V,L_M>::uqVectorLhFunction_BaseClass(
  void (*routinePtr)(const P_V& paramValues, const void* functionDataPtr, L_V& likelihoodValues),
  const void* functionDataPtr)
  :
  m_routinePtr(routinePtr),
  m_routineDataPtr(functionDataPtr)
{
}

template<class P_V,class P_M,class L_V,class L_M>
uqVectorLhFunction_BaseClass<P_V,P_M,L_V,L_M>::~uqVectorLhFunction_BaseClass()
{
}

//*****************************************************
// Misfit class
//*****************************************************
template<class P_V,class P_M,class L_V,class L_M>
class uqMisfitVectorLhFunction_Class : public uqVectorLhFunction_BaseClass<P_V,P_M,L_V,L_M> {
public:
  uqMisfitVectorLhFunction_Class(void (*routinePtr)(const P_V& paramValues, const void* functionDataPtr, L_V& likelihoodValues),
                                    const void* functionDataPtr);
 ~uqMisfitVectorLhFunction_Class();
  void computeActualLikelihoods  (const P_V& paramValues, L_V& likelihoodValues) const;
  void computeMinus2LnLikelihoods(const P_V& paramValues, L_V& likelihoodValues) const;
  void computeMisfits            (const P_V& paramValues, L_V& likelihoodValues) const;

  using uqVectorLhFunction_BaseClass<P_V,P_M,L_V,L_M>::m_routinePtr;
  using uqVectorLhFunction_BaseClass<P_V,P_M,L_V,L_M>::m_routineDataPtr;
};

template<class P_V,class P_M,class L_V,class L_M>
uqMisfitVectorLhFunction_Class<P_V,P_M,L_V,L_M>::uqMisfitVectorLhFunction_Class(
  void (*routinePtr)(const P_V& paramValues, const void* functionDataPtr, L_V& likelihoodValues),
  const void* functionDataPtr)
  :
  uqVectorLhFunction_BaseClass<P_V,P_M,L_V,L_M>(routinePtr,functionDataPtr)
{
}

template<class P_V,class P_M,class L_V,class L_M>
uqMisfitVectorLhFunction_Class<P_V,P_M,L_V,L_M>::~uqMisfitVectorLhFunction_Class()
{
}

template<class P_V,class P_M,class L_V,class L_M>
void
uqMisfitVectorLhFunction_Class<P_V,P_M,L_V,L_M>::computeActualLikelihoods(const P_V& paramValues, L_V& likelihoodValues) const
{
  UQ_FATAL_TEST_MACRO(false,
                      paramValues.env().rank(),
                      "uqMisfitVectorLhFunction_Class<P_V,P_M,L_V,L_M>::computeActualLikelihoods()",
                      "this method should not be called in the case of this class");
  return;
}

template<class P_V,class P_M,class L_V,class L_M>
void
uqMisfitVectorLhFunction_Class<P_V,P_M,L_V,L_M>::computeMinus2LnLikelihoods(const P_V& paramValues, L_V& likelihoodValues) const
{
  UQ_FATAL_TEST_MACRO(false,
                      paramValues.env().rank(),
                      "uqMisfitVectorLhFunction_Class<P_V,P_M,L_V,L_M>::computeMinus2LnLikelihoods()",
                      "this method should not be called in the case of this class");
  return;
}

template<class P_V,class P_M,class L_V,class L_M>
void
uqMisfitVectorLhFunction_Class<P_V,P_M,L_V,L_M>::computeMisfits(const P_V& paramValues, L_V& likelihoodValues) const
{
  return m_routinePtr(paramValues, m_routineDataPtr, likelihoodValues);
}

//*****************************************************
// Complete class
//*****************************************************
template<class P_V,class P_M,class L_V,class L_M>
class uqCompleteVectorLhFunction_Class : public uqVectorLhFunction_BaseClass<P_V,P_M,L_V,L_M> {
public:
  uqCompleteVectorLhFunction_Class(void (*routinePtr)(const P_V& paramValues, const void* functionDataPtr, L_V& likelihoodValues),
                                      const void* functionDataPtr,
                                      bool routineComputesMinus2LogOfLikelihood);
 ~uqCompleteVectorLhFunction_Class();
  void computeActualLikelihoods  (const P_V& paramValues, L_V& likelihoodValues) const;
  void computeMinus2LnLikelihoods(const P_V& paramValues, L_V& likelihoodValues) const;
  void computeMisfits            (const P_V& paramValues, L_V& likelihoodValues) const;

  using uqVectorLhFunction_BaseClass<P_V,P_M,L_V,L_M>::m_routinePtr;
  using uqVectorLhFunction_BaseClass<P_V,P_M,L_V,L_M>::m_routineDataPtr;

private:
  bool m_routineComputesMinus2LogOfLikelihood;
};

template<class P_V,class P_M,class L_V,class L_M>
uqCompleteVectorLhFunction_Class<P_V,P_M,L_V,L_M>::uqCompleteVectorLhFunction_Class(
  void (*routinePtr)(const P_V& paramValues, const void* functionDataPtr, L_V& likelihoodValues),
  const void* functionDataPtr,
  bool        routineComputesMinus2LogOfLikelihood)
  :
  uqVectorLhFunction_BaseClass<P_V,P_M,L_V,L_M>(routinePtr,functionDataPtr),
  m_routineComputesMinus2LogOfLikelihood(routineComputesMinus2LogOfLikelihood)
{
}

template<class P_V,class P_M,class L_V,class L_M>
uqCompleteVectorLhFunction_Class<P_V,P_M,L_V,L_M>::~uqCompleteVectorLhFunction_Class()
{
}

template<class P_V,class P_M,class L_V,class L_M>
void
uqCompleteVectorLhFunction_Class<P_V,P_M,L_V,L_M>::computeActualLikelihoods(const P_V& paramValues, L_V& likelihoodValues) const
{
  m_routinePtr(paramValues, m_routineDataPtr, likelihoodValues);
  if (m_routineComputesMinus2LogOfLikelihood) {
    for (unsigned int i = 0; i < likelihoodValues.size(); ++i) {
      likelihoodValues[i] = exp(-.5*likelihoodValues[i]);
    }
  }

  return;
}

template<class P_V,class P_M,class L_V,class L_M>
void
uqCompleteVectorLhFunction_Class<P_V,P_M,L_V,L_M>::computeMinus2LnLikelihoods(const P_V& paramValues, L_V& likelihoodValues) const
{
  m_routinePtr(paramValues, m_routineDataPtr, likelihoodValues);
  if (m_routineComputesMinus2LogOfLikelihood == false) {
    for (unsigned int i = 0; i < likelihoodValues.size(); ++i) {
      likelihoodValues[i] = -2.*log(likelihoodValues[i]);
    }
  }

  return;
}

template<class P_V,class P_M,class L_V,class L_M>
void
uqCompleteVectorLhFunction_Class<P_V,P_M,L_V,L_M>::computeMisfits(const P_V& paramValues, L_V& likelihoodValues) const
{
  UQ_FATAL_TEST_MACRO(false,
                      paramValues.env().rank(),
                      "uqCompleteVectorLhFunction_Class<P_V,P_M,L_V,L_M>::computeMisfits()",
                      "this method should not be called in the case of this class");
  return;
}
#endif // __UQ_VECTOR_LH_FUNCTION_H__
