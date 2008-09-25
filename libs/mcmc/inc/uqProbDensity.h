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
           uqProbDensity_BaseClass();
  virtual ~uqProbDensity_BaseClass();

  virtual double actualDensity  (const V& paramValues) const = 0;
  virtual double minus2LnDensity(const V& paramValues) const = 0;
};

template<class V, class M>
uqProbDensity_BaseClass<V,M>::uqProbDensity_BaseClass()
{
}

template<class V, class M>
uqProbDensity_BaseClass<V,M>::~uqProbDensity_BaseClass()
{
}

//*****************************************************
// Routine probability density class
//*****************************************************
template<class V, class M>
class uqRoutineProbDensity_Class : public uqProbDensity_BaseClass<V,M> {
public:
  uqRoutineProbDensity_Class(double (*routinePtr)(const V& paramValues, const void* routineDataPtr),
                             const void* routineDataPtr,
                             bool routineComputesMinus2LogOfDensity);
 ~uqRoutineProbDensity_Class();

  double actualDensity  (const V& paramValues) const;
  double minus2LnDensity(const V& paramValues) const;

protected:
  double (*m_routinePtr)(const V& paramValues, const void* routineDataPtr);
  const void* m_routineDataPtr;

  bool m_routineComputesMinus2LogOfDensity;

};

template<class V, class M>
uqRoutineProbDensity_Class<V,M>::uqRoutineProbDensity_Class(
  double (*routinePtr)(const V& paramValues, const void* routineDataPtr),
  const void* routineDataPtr,
  bool        routineComputesMinus2LogOfDensity)
  :
  uqProbDensity_BaseClass<V,M>(),
  m_routinePtr                       (routinePtr),
  m_routineDataPtr                   (routineDataPtr),
  m_routineComputesMinus2LogOfDensity(routineComputesMinus2LogOfDensity)
{
}

template<class V, class M>
uqRoutineProbDensity_Class<V,M>::~uqRoutineProbDensity_Class()
{
}

template<class V, class M>
double
uqRoutineProbDensity_Class<V,M>::minus2LnDensity(const V& paramValues) const
{
  double value = m_routinePtr(paramValues, m_routineDataPtr);
  if (m_routineComputesMinus2LogOfDensity == false) {
    value = -2.*log(value);
  }

  return value;
}

template<class V, class M>
double
uqRoutineProbDensity_Class<V,M>::actualDensity(const V& paramValues) const
{
  double value = m_routinePtr(paramValues, m_routineDataPtr);
  if (m_routineComputesMinus2LogOfDensity) {
    value = exp(-.5*value);
  }

  return value;
}

//*****************************************************
// Bayesian probability density class
//*****************************************************
template<class V, class M>
class uqBayesianProbDensity_Class : public uqProbDensity_BaseClass<V,M> {
public:
  uqBayesianProbDensity_Class(const uqProbDensity_BaseClass<V,M>* priorDensity,
                              const uqProbDensity_BaseClass<V,M>* likelihoodFunction); 
 ~uqBayesianProbDensity_Class();

  double actualDensity  (const V& paramValues) const;
  double minus2LnDensity(const V& paramValues) const;

protected:
  const uqProbDensity_BaseClass<V,M>* m_priorDensity;
  const uqProbDensity_BaseClass<V,M>* m_likelihoodFunction;
};

template<class V,class M>
uqBayesianProbDensity_Class<V,M>::uqBayesianProbDensity_Class(
  const uqProbDensity_BaseClass<V,M>* priorDensity,
  const uqProbDensity_BaseClass<V,M>* likelihoodFunction)
  :
  uqProbDensity_BaseClass<V,M>(),
  m_priorDensity      (priorDensity),
  m_likelihoodFunction(likelihoodFunction)
{
}

template<class V,class M>
uqBayesianProbDensity_Class<V,M>::~uqBayesianProbDensity_Class()
{
}

template<class V, class M>
double
uqBayesianProbDensity_Class<V,M>::minus2LnDensity(const V& paramValues) const
{
  double value = 0.;

  value = m_priorDensity->minus2LnDensity(paramValues);
  value += m_likelihoodFunction->minus2LnDensity(paramValues);

  return value;
}

template<class V, class M>
double
uqBayesianProbDensity_Class<V,M>::actualDensity(const V& paramValues) const
{
  double value = 0.;

  value = m_priorDensity->actualDensity(paramValues);
  value *= m_likelihoodFunction->actualDensity(paramValues);

  return value;
}

#endif // __UQ_PROB_DENSITY_H__
