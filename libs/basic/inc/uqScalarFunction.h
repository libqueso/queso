/* uq/libs/queso/inc/uqScalarFunction.h
 *
 * Copyright (C) 2008 The QUESO Team, http://queso.ices.utexas.edu/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __UQ_SCALAR_FUNCTION_H__
#define __UQ_SCALAR_FUNCTION_H__

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

          const uqVectorSetClass<V,M>& domainSet        ()                                        const;
  virtual       double                 actualValue      (const V& domainVector)                   const = 0;
  virtual       double                 minus2LnValue    (const V& domainVector)                   const = 0;
  virtual       void                   gradOfActual     (const V& domainVector, V& gradVector)    const = 0;
  virtual       void                   gradOfMinus2Ln   (const V& domainVector, V& gradVector)    const = 0;
  virtual       void                   hessianOfActual  (const V& domainVector, M& hessianMatrix) const = 0;
  virtual       void                   hessianOfMinus2Ln(const V& domainVector, M& hessianMatrix) const = 0;

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
                               double (*valueRoutinePtr  )(const V& domainVector, const void* routinesDataPtr),
                               void   (*gradRoutinePtr   )(const V& domainVector, const void* routinesDataPtr, V& gradVector),
                               void   (*hessianRoutinePtr)(const V& domainVector, const void* routinesDataPtr, M& hessianMatrix),
                               const void* routinesDataPtr,
                               bool routinesAreForMinus2Ln);
  virtual ~uqGenericScalarFunctionClass();

  double actualValue      (const V& domainVector)                   const;
  double minus2LnValue    (const V& domainVector)                   const;
  void   gradOfActual     (const V& domainVector, V& gradVector)    const;
  void   gradOfMinus2Ln   (const V& domainVector, V& gradVector)    const;
  void   hessianOfActual  (const V& domainVector, M& hessianMatrix) const;
  void   hessianOfMinus2Ln(const V& domainVector, M& hessianMatrix) const;

protected:
  using uqBaseScalarFunctionClass<V,M>::m_env;
  using uqBaseScalarFunctionClass<V,M>::m_prefix;
  using uqBaseScalarFunctionClass<V,M>::m_domainSet;

  double (*m_valueRoutinePtr  )(const V& domainVector, const void* routinesDataPtr);
  void   (*m_gradRoutinePtr   )(const V& domainVector, const void* routinesDataPtr, V& gradVector);
  void   (*m_hessianRoutinePtr)(const V& domainVector, const void* routinesDataPtr, M& hessianMatrix);
  const void* m_routinesDataPtr;
  bool m_routinesAreForMinus2Ln;
};

template<class V,class M>
uqGenericScalarFunctionClass<V,M>::uqGenericScalarFunctionClass(
  const char*                  prefix,
  const uqVectorSetClass<V,M>& domainSet,
  double (*valueRoutinePtr)(const V& domainVector, const void* routinesDataPtr),
  void   (*gradRoutinePtr   )(const V& domainVector, const void* routinesDataPtr, V& gradVector),
  void   (*hessianRoutinePtr)(const V& domainVector, const void* routinesDataPtr, M& hessianMatrix),
  const void* routinesDataPtr,
  bool routinesAreForMinus2Ln)
  :
  uqBaseScalarFunctionClass<V,M>(((std::string)(prefix)+"gen").c_str(), domainSet),
  m_valueRoutinePtr             (valueRoutinePtr),
  m_gradRoutinePtr              (gradRoutinePtr),
  m_hessianRoutinePtr           (hessianRoutinePtr),
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
uqGenericScalarFunctionClass<V,M>::actualValue(const V& domainVector) const
{
  double value = m_valueRoutinePtr(domainVector, m_routinesDataPtr);
  if (m_routinesAreForMinus2Ln) {
    value = exp(-.5*value);
  }
  return value;
}

template<class V,class M>
double
uqGenericScalarFunctionClass<V,M>::minus2LnValue(const V& domainVector) const
{
  double value = m_valueRoutinePtr(domainVector, m_routinesDataPtr);
  if (m_routinesAreForMinus2Ln == false) {
    value = -2.*log(value);
  }
  return value;
}

template<class V,class M>
void
uqGenericScalarFunctionClass<V,M>::gradOfActual(const V& domainVector, V& gradVector) const
{
  return;
}

template<class V,class M>
void
uqGenericScalarFunctionClass<V,M>::gradOfMinus2Ln(const V& domainVector, V& gradVector) const
{
  return;
}

template<class V,class M>
void
uqGenericScalarFunctionClass<V,M>::hessianOfActual(const V& domainVector, M& hessianMatrix) const
{
  return;
}

template<class V,class M>
void
uqGenericScalarFunctionClass<V,M>::hessianOfMinus2Ln(const V& domainVector, M& hessianMatrix) const
{
  return;
}
#endif // __UQ_SCALAR_FUNCTION_H__
