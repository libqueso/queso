/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of the QUESO Library (Quantification of Uncertainty
 * for Estimation, Simulation and Optimization).
 *
 * QUESO is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QUESO is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QUESO. If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * $Id$
 *
 * Brief description of this file: 
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __UQ_SCALAR_FUNCTION_H__
#define __UQ_SCALAR_FUNCTION_H__

#include <uqVectorSet.h>
#include <uqVectorSubset.h>
#include <uqEnvironment.h>
#include <uqDefines.h>

//*****************************************************
// Base class
//*****************************************************
template<class V,class M>
class uqBaseScalarFunctionClass {
public:
           uqBaseScalarFunctionClass(const char*                  prefix,
                                     const uqVectorSetClass<V,M>& domainSet);
  virtual ~uqBaseScalarFunctionClass();

          const uqVectorSetClass<V,M>& domainSet        ()                                                                                                   const;
  virtual       double                 actualValue      (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const = 0;
  virtual       double                 minus2LnValue    (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const = 0;

protected:
  const uqBaseEnvironmentClass& m_env;
        std::string             m_prefix;
  const uqVectorSetClass<V,M>&  m_domainSet;
};

template<class V,class M>
uqBaseScalarFunctionClass<V,M>::uqBaseScalarFunctionClass(
  const char*                  prefix,
  const uqVectorSetClass<V,M>& domainSet)
  :
  m_env      (domainSet.env()),
  m_prefix   ((std::string)(prefix)+"func_"),
  m_domainSet(domainSet)
{
}

template<class V,class M>
uqBaseScalarFunctionClass<V,M>::~uqBaseScalarFunctionClass()
{
}

template<class V,class M>
const uqVectorSetClass<V,M>&
uqBaseScalarFunctionClass<V,M>::domainSet() const
{
  return m_domainSet;
}

//*****************************************************
// Generic class
//*****************************************************
template<class V,class M>
class uqGenericScalarFunctionClass : public uqBaseScalarFunctionClass<V,M> {
public:
  uqGenericScalarFunctionClass(const char*                  prefix,
                               const uqVectorSetClass<V,M>& domainSet,
                               double (*valueRoutinePtr)(const V& domainVector, const V* domainDirection, const void* routinesDataPtr, V* gradVector, M* hessianMatrix, V* hessianEffect),
                               const void* routinesDataPtr,
                               bool routinesAreForMinus2Ln);
  virtual ~uqGenericScalarFunctionClass();

  double actualValue      (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  double minus2LnValue    (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;

protected:
  using uqBaseScalarFunctionClass<V,M>::m_env;
  using uqBaseScalarFunctionClass<V,M>::m_prefix;
  using uqBaseScalarFunctionClass<V,M>::m_domainSet;

  double (*m_valueRoutinePtr)(const V& domainVector, const V* domainDirection, const void* routinesDataPtr, V* gradVector, M* hessianMatrix, V* hessianEffect);
  const void* m_routinesDataPtr;
  bool m_routinesAreForMinus2Ln;
};

template<class V,class M>
uqGenericScalarFunctionClass<V,M>::uqGenericScalarFunctionClass(
  const char*                  prefix,
  const uqVectorSetClass<V,M>& domainSet,
  double (*valueRoutinePtr)(const V& domainVector, const V* domainDirection, const void* routinesDataPtr, V* gradVector, M* hessianMatrix, V* hessianEffect),
  const void* routinesDataPtr,
  bool routinesAreForMinus2Ln)
  :
  uqBaseScalarFunctionClass<V,M>(((std::string)(prefix)+"gen").c_str(), domainSet),
  m_valueRoutinePtr             (valueRoutinePtr),
  m_routinesDataPtr             (routinesDataPtr),
  m_routinesAreForMinus2Ln      (routinesAreForMinus2Ln)
{
}

template<class V,class M>
uqGenericScalarFunctionClass<V,M>::~uqGenericScalarFunctionClass()
{
}

template<class V,class M>
double
uqGenericScalarFunctionClass<V,M>::actualValue(const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const
{
  UQ_FATAL_TEST_MACRO(m_valueRoutinePtr == NULL,
                      m_env.fullRank(),
                      "uqGenericScalarFunctionClass<V,M>::actualValue()",
                      "m_valueRoutinePtr = NULL");

  double value = m_valueRoutinePtr(domainVector, domainDirection, m_routinesDataPtr, gradVector, hessianMatrix, hessianEffect);
  if (m_routinesAreForMinus2Ln) {
    value = std::exp(-.5*value);
    UQ_FATAL_TEST_MACRO(true,
                        m_env.fullRank(),
                        "uqGenericScalarFunctionClass<V,M>::gradOfActual()",
                        "INCOMPLETE CODE");
  }
  return value;
}

template<class V,class M>
double
uqGenericScalarFunctionClass<V,M>::minus2LnValue(const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const
{
  UQ_FATAL_TEST_MACRO(m_valueRoutinePtr == NULL,
                      m_env.fullRank(),
                      "uqGenericScalarFunctionClass<V,M>::minus2LnValue()",
                      "m_valueRoutinePtr = NULL");

  double value = m_valueRoutinePtr(domainVector, domainDirection, m_routinesDataPtr, gradVector, hessianMatrix, hessianEffect);
  if (m_routinesAreForMinus2Ln == false) {
    value = -2.*log(value);
    UQ_FATAL_TEST_MACRO(true,
                        m_env.fullRank(),
                        "uqGenericScalarFunctionClass<V,M>::gradOfMinus2Ln()",
                        "INCOMPLETE CODE");
  }
  return value;
}
#endif // __UQ_SCALAR_FUNCTION_H__
