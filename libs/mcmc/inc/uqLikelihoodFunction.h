/* uq/libs/mcmc/inc/uqLikelihoodFunction.h
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

#ifndef __UQ_LIKELIHOOD_FUNCTION_H__
#define __UQ_LIKELIHOOD_FUNCTION_H__

#include <uqEnvironment.h>
#include <uqDefines.h>

//*****************************************************
// Class to accomodate the likehood function.
// User must provide a likelihood() routine.
//*****************************************************

//*****************************************************
// Base class
//*****************************************************
template<class V, class M>
class uq_LikelihoodFunction_BaseClass {
public:
  uq_LikelihoodFunction_BaseClass(void (*functionPtr)(const V& paramValues, const void* functionDataPtr, V& likelihoodValues),
                                 const void* functionDataPtr);
  virtual ~uq_LikelihoodFunction_BaseClass();
  virtual void computeActualLikelihoods  (const V& paramValues, V& likelihoodValues) const = 0;
  virtual void computeMinus2LnLikelihoods(const V& paramValues, V& likelihoodValues) const = 0;
  virtual void computeMisfits            (const V& paramValues, V& likelihoodValues) const = 0;

protected:
  void (*m_functionPtr)(const V& paramValues, const void* functionDataPtr, V& likelihoodValues);
  const void* m_functionDataPtr;
};

template<class V, class M>
uq_LikelihoodFunction_BaseClass<V,M>::uq_LikelihoodFunction_BaseClass(
  void (*functionPtr)(const V& paramValues, const void* functionDataPtr, V& likelihoodValues),
  const void* functionDataPtr)
  :
  m_functionPtr(functionPtr),
  m_functionDataPtr(functionDataPtr)
{
}

template<class V, class M>
uq_LikelihoodFunction_BaseClass<V,M>::~uq_LikelihoodFunction_BaseClass()
{
}

//*****************************************************
// Misfit class
//*****************************************************
template<class V, class M>
class uq_MisfitLikelihoodFunction_Class : public uq_LikelihoodFunction_BaseClass<V,M> {
public:
  uq_MisfitLikelihoodFunction_Class(void (*functionPtr)(const V& paramValues, const void* functionDataPtr, V& likelihoodValues),
                                    const void* functionDataPtr);
 ~uq_MisfitLikelihoodFunction_Class();
  void computeActualLikelihoods  (const V& paramValues, V& likelihoodValues) const;
  void computeMinus2LnLikelihoods(const V& paramValues, V& likelihoodValues) const;
  void computeMisfits            (const V& paramValues, V& likelihoodValues) const;

  using uq_LikelihoodFunction_BaseClass<V,M>::m_functionPtr;
  using uq_LikelihoodFunction_BaseClass<V,M>::m_functionDataPtr;
};

template<class V, class M>
uq_MisfitLikelihoodFunction_Class<V,M>::uq_MisfitLikelihoodFunction_Class(
  void (*functionPtr)(const V& paramValues, const void* functionDataPtr, V& likelihoodValues),
  const void* functionDataPtr)
  :
  uq_LikelihoodFunction_BaseClass<V,M>(functionPtr,functionDataPtr)
{
}

template<class V, class M>
uq_MisfitLikelihoodFunction_Class<V,M>::~uq_MisfitLikelihoodFunction_Class()
{
}

template<class V, class M>
void
uq_MisfitLikelihoodFunction_Class<V,M>::computeActualLikelihoods(const V& paramValues, V& likelihoodValues) const
{
  UQ_FATAL_TEST_MACRO(false,
                      paramValues.env().rank(),
                      "uq_MisfitLikelihoodFunction_Class<V,M>::computeActualLikelihoods()",
                      "this method should not be called in the case of this class");
  return;
}

template<class V, class M>
void
uq_MisfitLikelihoodFunction_Class<V,M>::computeMinus2LnLikelihoods(const V& paramValues, V& likelihoodValues) const
{
  UQ_FATAL_TEST_MACRO(false,
                      paramValues.env().rank(),
                      "uq_MisfitLikelihoodFunction_Class<V,M>::computeMinus2LnLikelihoods()",
                      "this method should not be called in the case of this class");
  return;
}

template<class V, class M>
void
uq_MisfitLikelihoodFunction_Class<V,M>::computeMisfits(const V& paramValues, V& likelihoodValues) const
{
  return m_functionPtr(paramValues, m_functionDataPtr, likelihoodValues);
}
#endif // __UQ_LIKELIHOOD_FUNCTION_H__
