/* uq/libs/mcmc/inc/uqBayesProbDensity.h
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

#ifndef __UQ_BAYES_PROB_DENSITY_H__
#define __UQ_BAYES_PROB_DENSITY_H__

#include <uqProbDensity.h>
#include <uqScalarLhFunction.h>

//*****************************************************
// Base class
//*****************************************************
template<class P_V,class P_M>
class uqBayesProbDensity_BaseClass {
public:
  uqBayesProbDensity_BaseClass(const uqProbDensity_BaseClass     <P_V,P_M>* priorDensity,
                               const uqScalarLhFunction_BaseClass<P_V,P_M>* likelihoodFunction); 
  virtual ~uqBayesProbDensity_BaseClass();
  virtual double actualDensity  (const P_V& paramValues) const = 0;
  virtual double minus2LnDensity(const P_V& paramValues) const = 0;

protected:
  const uqProbDensity_BaseClass     <P_V,P_M>* m_priorDensity;
  const uqScalarLhFunction_BaseClass<P_V,P_M>* m_likelihoodFunction;
};

template<class P_V,class P_M>
uqBayesProbDensity_BaseClass<P_V,P_M>::uqBayesProbDensity_BaseClass(
  const uqProbDensity_BaseClass     <P_V,P_M>* priorDensity,
  const uqScalarLhFunction_BaseClass<P_V,P_M>* likelihoodFunction)
  :
  m_priorDensity      (priorDensity),
  m_likelihoodFunction(likelihoodFunction)
{
}

template<class P_V,class P_M>
uqBayesProbDensity_BaseClass<P_V,P_M>::~uqBayesProbDensity_BaseClass()
{
}

//*****************************************************
// M2l class
//*****************************************************
template<class P_V,class P_M>
class uqM2lBayesProbDensity_Class : public uqBayesProbDensity_BaseClass<P_V,P_M> {
public:
  uqM2lBayesProbDensity_Class(const uqProbDensity_BaseClass     <P_V,P_M>* priorDensity,
                              const uqScalarLhFunction_BaseClass<P_V,P_M>* likelihoodFunction); 
 ~uqM2lBayesProbDensity_Class();

  double actualDensity  (const P_V& paramValues) const;
  double minus2LnDensity(const P_V& paramValues) const;

  using uqBayesProbDensity_BaseClass<P_V,P_M>::m_priorDensity;
  using uqBayesProbDensity_BaseClass<P_V,P_M>::m_likelihoodFunction;
};

template<class P_V,class P_M>
uqM2lBayesProbDensity_Class<P_V,P_M>::uqM2lBayesProbDensity_Class(
  const uqProbDensity_BaseClass     <P_V,P_M>* priorDensity,
  const uqScalarLhFunction_BaseClass<P_V,P_M>* likelihoodFunction)
  :
  uqBayesProbDensity_BaseClass<P_V,P_M>(priorDensity,likelihoodFunction)
{
}

template<class P_V,class P_M>
uqM2lBayesProbDensity_Class<P_V,P_M>::~uqM2lBayesProbDensity_Class()
{
}

template<class P_V,class P_M>
double
uqM2lBayesProbDensity_Class<P_V,P_M>::minus2LnDensity(const P_V& paramValues) const
{
  double value = m_priorDensity->minus2LnDensity(paramValues);
  value += m_likelihoodFunction->minus2LnLikelihood(paramValues);

  return value;
}

template<class P_V,class P_M>
double
uqM2lBayesProbDensity_Class<P_V,P_M>::actualDensity(const P_V& paramValues) const
{
  double value = m_priorDensity->actualDensity(paramValues);
  value *= m_likelihoodFunction->actualLikelihood(paramValues);

  return value;
}
#endif // __UQ_PROB_DENSITY_H__
