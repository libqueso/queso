/* uq/libs/mcmc/inc/uqProbDensity.h
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

#ifndef __UQ_PROB_DENSITY_H__
#define __UQ_PROB_DENSITY_H__

#include <uqScalarLhFunction.h>
#include <uqEnvironment.h>
#include <math.h>

//*****************************************************
// Classes to accomodate a probability density routine.
//*****************************************************

//*****************************************************
// Base class
//*****************************************************
template<class V, class M>
class uqProbDensity_BaseClass {
public:
  uqProbDensity_BaseClass(double (*routinePtr)(const V& paramValues, const void* routineDataPtr),
                          const void* routineDataPtr);
  uqProbDensity_BaseClass(const uqProbDensity_BaseClass     <V,M>* priorDensity,
                          const uqScalarLhFunction_BaseClass<V,M>* likelihoodFunction); 
  virtual ~uqProbDensity_BaseClass();
  virtual double actualDensity  (const V& paramValues) const = 0;
  virtual double minus2LnDensity(const V& paramValues) const = 0;

protected:
  double (*m_routinePtr)(const V& paramValues, const void* routineDataPtr);
  const void* m_routineDataPtr;

  const uqProbDensity_BaseClass     <V,M>* m_priorDensity;
  const uqScalarLhFunction_BaseClass<V,M>* m_likelihoodFunction;
};

template<class V, class M>
uqProbDensity_BaseClass<V,M>::uqProbDensity_BaseClass(
  double (*routinePtr)(const V& paramValues, const void* routineDataPtr),
  const void* routineDataPtr)
  :
  m_routinePtr        (routinePtr),
  m_routineDataPtr    (routineDataPtr),
  m_priorDensity      (NULL),
  m_likelihoodFunction(NULL)
{
}

template<class V,class M>
uqProbDensity_BaseClass<V,M>::uqProbDensity_BaseClass(
  const uqProbDensity_BaseClass     <V,M>* priorDensity,
  const uqScalarLhFunction_BaseClass<V,M>* likelihoodFunction)
  :
  m_routinePtr        (NULL),
  m_routineDataPtr    (NULL),
  m_priorDensity      (priorDensity),
  m_likelihoodFunction(likelihoodFunction)
{
}

template<class V, class M>
uqProbDensity_BaseClass<V,M>::~uqProbDensity_BaseClass()
{
}

//*****************************************************
// M2l class
//*****************************************************
template<class V, class M>
class uqM2lProbDensity_Class : public uqProbDensity_BaseClass<V,M> {
public:
  uqM2lProbDensity_Class(double (*routinePtr)(const V& paramValues, const void* routineDataPtr),
                         const void* routineDataPtr);
  uqM2lProbDensity_Class(const uqProbDensity_BaseClass     <V,M>* priorDensity,
                         const uqScalarLhFunction_BaseClass<V,M>* likelihoodFunction); 
 ~uqM2lProbDensity_Class();

  double actualDensity  (const V& paramValues) const;
  double minus2LnDensity(const V& paramValues) const;

  using uqProbDensity_BaseClass<V,M>::m_routinePtr;
  using uqProbDensity_BaseClass<V,M>::m_routineDataPtr;

  using uqProbDensity_BaseClass<V,M>::m_priorDensity;
  using uqProbDensity_BaseClass<V,M>::m_likelihoodFunction;
};

template<class V, class M>
uqM2lProbDensity_Class<V,M>::uqM2lProbDensity_Class(
  double (*routinePtr)(const V& paramValues, const void* routineDataPtr),
  const void* routineDataPtr)
  :
  uqProbDensity_BaseClass<V,M>(routinePtr,routineDataPtr)
{
}

template<class V,class M>
uqM2lProbDensity_Class<V,M>::uqM2lProbDensity_Class(
  const uqProbDensity_BaseClass     <V,M>* priorDensity,
  const uqScalarLhFunction_BaseClass<V,M>* likelihoodFunction)
  :
  uqProbDensity_BaseClass<V,M>(priorDensity,likelihoodFunction)
{
}

template<class V, class M>
uqM2lProbDensity_Class<V,M>::~uqM2lProbDensity_Class()
{
}

template<class V, class M>
double
uqM2lProbDensity_Class<V,M>::minus2LnDensity(const V& paramValues) const
{
  double value = 0.;

  if (m_routinePtr) {
    value = m_routinePtr(paramValues, m_routineDataPtr);
  }
  else {
    value = m_priorDensity->minus2LnDensity(paramValues);
    value += m_likelihoodFunction->minus2LnLikelihood(paramValues);
  }

  return value;
}

template<class V, class M>
double
uqM2lProbDensity_Class<V,M>::actualDensity(const V& paramValues) const
{
  double value = 0.;

  if (m_routinePtr) {
    value = m_routinePtr(paramValues, m_routineDataPtr);
    value = exp(-.5*value);
  }
  else {
    value = m_priorDensity->actualDensity(paramValues);
    value *= m_likelihoodFunction->actualLikelihood(paramValues);
  }

  return value;
}
#endif // __UQ_PROB_DENSITY_H__
