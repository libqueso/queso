/* uq/libs/mcmc/inc/uqVectorProbDensity.h
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

#ifndef __UQ_VECTOR_PROB_DENSITY_H__
#define __UQ_VECTOR_PROB_DENSITY_H__

#include <uqEnvironment.h>
#include <math.h>

//*****************************************************
// Classes to accomodate a probability density routine.
//*****************************************************

//*****************************************************
// Base class
//*****************************************************
template<class V, class M>
class uqBaseVectorProbDensityClass {
public:
           uqBaseVectorProbDensityClass();
  virtual ~uqBaseVectorProbDensityClass();

  virtual double actualDensity  (const V& paramValues) const = 0;
  virtual double minus2LnDensity(const V& paramValues) const = 0;
};

template<class V, class M>
uqBaseVectorProbDensityClass<V,M>::uqBaseVectorProbDensityClass()
{
}

template<class V, class M>
uqBaseVectorProbDensityClass<V,M>::~uqBaseVectorProbDensityClass()
{
}

//*****************************************************
// Routine probability density class
//*****************************************************
template<class V, class M>
class uqRoutineVectorProbDensityClass : public uqBaseVectorProbDensityClass<V,M> {
public:
  uqRoutineVectorProbDensityClass(double (*routinePtr)(const V& paramValues, const void* routineDataPtr),
                                  const void* routineDataPtr,
                                  bool routineComputesMinus2LogOfDensity);
 ~uqRoutineVectorProbDensityClass();

  double actualDensity  (const V& paramValues) const;
  double minus2LnDensity(const V& paramValues) const;

protected:
  double (*m_routinePtr)(const V& paramValues, const void* routineDataPtr);
  const void* m_routineDataPtr;

  bool m_routineComputesMinus2LogOfDensity;

};

template<class V, class M>
uqRoutineVectorProbDensityClass<V,M>::uqRoutineVectorProbDensityClass(
  double (*routinePtr)(const V& paramValues, const void* routineDataPtr),
  const void* routineDataPtr,
  bool        routineComputesMinus2LogOfDensity)
  :
  uqBaseVectorProbDensityClass<V,M>(),
  m_routinePtr                       (routinePtr),
  m_routineDataPtr                   (routineDataPtr),
  m_routineComputesMinus2LogOfDensity(routineComputesMinus2LogOfDensity)
{
}

template<class V, class M>
uqRoutineVectorProbDensityClass<V,M>::~uqRoutineVectorProbDensityClass()
{
}

template<class V, class M>
double
uqRoutineVectorProbDensityClass<V,M>::minus2LnDensity(const V& paramValues) const
{
  double value = m_routinePtr(paramValues, m_routineDataPtr);
  if (m_routineComputesMinus2LogOfDensity == false) {
    value = -2.*log(value);
  }

  return value;
}

template<class V, class M>
double
uqRoutineVectorProbDensityClass<V,M>::actualDensity(const V& paramValues) const
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
class uqBayesianVectorProbDensityClass : public uqBaseVectorProbDensityClass<V,M> {
public:
  uqBayesianVectorProbDensityClass(const uqBaseVectorProbDensityClass<V,M>* priorDensity,
                                   const uqBaseVectorProbDensityClass<V,M>* likelihoodFunction); 
 ~uqBayesianVectorProbDensityClass();

  double actualDensity  (const V& paramValues) const;
  double minus2LnDensity(const V& paramValues) const;

protected:
  const uqBaseVectorProbDensityClass<V,M>* m_priorDensity;
  const uqBaseVectorProbDensityClass<V,M>* m_likelihoodFunction;
};

template<class V,class M>
uqBayesianVectorProbDensityClass<V,M>::uqBayesianVectorProbDensityClass(
  const uqBaseVectorProbDensityClass<V,M>* priorDensity,
  const uqBaseVectorProbDensityClass<V,M>* likelihoodFunction)
  :
  uqBaseVectorProbDensityClass<V,M>(),
  m_priorDensity      (priorDensity),
  m_likelihoodFunction(likelihoodFunction)
{
}

template<class V,class M>
uqBayesianVectorProbDensityClass<V,M>::~uqBayesianVectorProbDensityClass()
{
}

template<class V, class M>
double
uqBayesianVectorProbDensityClass<V,M>::minus2LnDensity(const V& paramValues) const
{
  double value = 0.;

  double value1 = m_priorDensity->minus2LnDensity(paramValues);
  double value2 = m_likelihoodFunction->minus2LnDensity(paramValues);

  //if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
  //  std::cout << "In uqBayesianVectorProbDensityClass<P_V,P_M>::minus2LnDensity()"
  //            << ", -2ln(prior) = " << value1
  //            << ", -2ln(like) = "  << value2
  //            << std::endl;
  //}

  return value;
}

template<class V, class M>
double
uqBayesianVectorProbDensityClass<V,M>::actualDensity(const V& paramValues) const
{
  double value = 0.;

  value = m_priorDensity->actualDensity(paramValues);
  value *= m_likelihoodFunction->actualDensity(paramValues);

  return value;
}

#endif // __UQ_VECTOR_PROB_DENSITY_H__
