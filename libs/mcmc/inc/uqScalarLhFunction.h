/* uq/libs/mcmc/inc/uqScalarLhFunction.h
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

#ifndef __UQ_SCALAR_LIKELIHOOD_FUNCTION_H__
#define __UQ_SCALAR_LIKELIHOOD_FUNCTION_H__

#include <uqEnvironment.h>
#include <uqDefines.h>

//*****************************************************
// Class to accomodate the likehood function.
// User must provide a likelihood() routine.
//*****************************************************

//*****************************************************
// Base class
//*****************************************************
template<class P_V,class P_M>
class uqScalarLhFunction_BaseClass {
public:
  uqScalarLhFunction_BaseClass(double (*routinePtr)(const P_V& paramValues, const void* functionDataPtr),
                               const void* functionDataPtr);
  virtual ~uqScalarLhFunction_BaseClass();
  virtual double actualLikelihood  (const P_V& paramValues) const = 0;
  virtual double minus2LnLikelihood(const P_V& paramValues) const = 0;
  virtual double misfit            (const P_V& paramValues) const = 0;

protected:
  double (*m_routinePtr)(const P_V& paramValues, const void* functionDataPtr);
  const void* m_routineDataPtr;
};

template<class P_V,class P_M>
uqScalarLhFunction_BaseClass<P_V,P_M>::uqScalarLhFunction_BaseClass(
  double (*routinePtr)(const P_V& paramValues, const void* functionDataPtr),
  const void* functionDataPtr)
  :
  m_routinePtr    (routinePtr),
  m_routineDataPtr(functionDataPtr)
{
}

template<class P_V,class P_M>
uqScalarLhFunction_BaseClass<P_V,P_M>::~uqScalarLhFunction_BaseClass()
{
}

//*****************************************************
// Misfit class
//*****************************************************
template<class P_V,class P_M>
class uqMisfitScalarLhFunction_Class : public uqScalarLhFunction_BaseClass<P_V,P_M> {
public:
  uqMisfitScalarLhFunction_Class(double (*routinePtr)(const P_V& paramValues, const void* functionDataPtr),
                                 const void* functionDataPtr);
 ~uqMisfitScalarLhFunction_Class();
  double actualLikelihood  (const P_V& paramValues) const;
  double minus2LnLikelihood(const P_V& paramValues) const;
  double misfit            (const P_V& paramValues) const;

  using uqScalarLhFunction_BaseClass<P_V,P_M>::m_routinePtr;
  using uqScalarLhFunction_BaseClass<P_V,P_M>::m_routineDataPtr;
};

template<class P_V,class P_M>
uqMisfitScalarLhFunction_Class<P_V,P_M>::uqMisfitScalarLhFunction_Class(
  double (*routinePtr)(const P_V& paramValues, const void* functionDataPtr),
  const void* functionDataPtr)
  :
  uqScalarLhFunction_BaseClass<P_V,P_M>(routinePtr,functionDataPtr)
{
}

template<class P_V,class P_M>
uqMisfitScalarLhFunction_Class<P_V,P_M>::~uqMisfitScalarLhFunction_Class()
{
}

template<class P_V,class P_M>
double
uqMisfitScalarLhFunction_Class<P_V,P_M>::actualLikelihood(const P_V& paramValues) const
{
  UQ_FATAL_TEST_MACRO(false,
                      paramValues.env().rank(),
                      "uqMisfitScalarLhFunction_Class<P_V,P_M>::actualLikelihood()",
                      "this method should not be called in the case of this class");
  return 0.;
}

template<class P_V,class P_M>
double
uqMisfitScalarLhFunction_Class<P_V,P_M>::minus2LnLikelihood(const P_V& paramValues) const
{
  UQ_FATAL_TEST_MACRO(false,
                      paramValues.env().rank(),
                      "uqMisfitScalarLhFunction_Class<P_V,P_M>::minus2LnLikelihood()",
                      "this method should not be called in the case of this class");
  return 0.;
}

template<class P_V,class P_M>
double
uqMisfitScalarLhFunction_Class<P_V,P_M>::misfit(const P_V& paramValues) const
{
  return m_routinePtr(paramValues, m_routineDataPtr);
}

//*****************************************************
// Complete class
//*****************************************************
template<class P_V,class P_M>
class uqCompleteScalarLhFunction_Class : public uqScalarLhFunction_BaseClass<P_V,P_M> {
public:
  uqCompleteScalarLhFunction_Class(double (*routinePtr)(const P_V& paramValues, const void* functionDataPtr),
                                   const void* functionDataPtr,
                                   bool routineComputesMinus2LogOfLikelihood);
 ~uqCompleteScalarLhFunction_Class();
  double actualLikelihood  (const P_V& paramValues) const;
  double minus2LnLikelihood(const P_V& paramValues) const;
  double misfit            (const P_V& paramValues) const;

  using uqScalarLhFunction_BaseClass<P_V,P_M>::m_routinePtr;
  using uqScalarLhFunction_BaseClass<P_V,P_M>::m_routineDataPtr;

private:
  bool m_routineComputesMinus2LogOfLikelihood;
};

template<class P_V,class P_M>
uqCompleteScalarLhFunction_Class<P_V,P_M>::uqCompleteScalarLhFunction_Class(
  double (*routinePtr)(const P_V& paramValues, const void* functionDataPtr),
  const void* functionDataPtr,
  bool        routineComputesMinus2LogOfLikelihood)
  :
  uqScalarLhFunction_BaseClass<P_V,P_M>(routinePtr,functionDataPtr),
  m_routineComputesMinus2LogOfLikelihood(routineComputesMinus2LogOfLikelihood)
{
}

template<class P_V,class P_M>
uqCompleteScalarLhFunction_Class<P_V,P_M>::~uqCompleteScalarLhFunction_Class()
{
}

template<class P_V,class P_M>
double
uqCompleteScalarLhFunction_Class<P_V,P_M>::actualLikelihood(const P_V& paramValues) const
{
  double value = m_routinePtr(paramValues, m_routineDataPtr);
  if (m_routineComputesMinus2LogOfLikelihood) {
    value = exp(-.5*value);
  }

  return value;
}

template<class P_V,class P_M>
double
uqCompleteScalarLhFunction_Class<P_V,P_M>::minus2LnLikelihood(const P_V& paramValues) const
{
  double value = m_routinePtr(paramValues, m_routineDataPtr);
  if (m_routineComputesMinus2LogOfLikelihood == false) {
    value = -2.*log(value);
  }

  return value;
}

template<class P_V,class P_M>
double
uqCompleteScalarLhFunction_Class<P_V,P_M>::misfit(const P_V& paramValues) const
{
  UQ_FATAL_TEST_MACRO(false,
                      paramValues.env().rank(),
                      "uqCompleteScalarLhFunction_Class<P_V,P_M>::misfit()",
                      "this method should not be called in the case of this class");
  return 0.;
}
#endif // __UQ_SCALAR_LIKELIHOOD_FUNCTION_H__
