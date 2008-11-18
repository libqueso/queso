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
template<class P_V,class P_M,class Q_V,class Q_M>
class uqBaseScalarFunctionClass {
public:
           uqBaseScalarFunctionClass(const char*                      prefix,
                                     const uqVectorSetClass<P_V,P_M>& domainSet,
                                     const uqVectorSetClass<Q_V,Q_M>& imageSet);
  virtual ~uqBaseScalarFunctionClass();

          const uqVectorSetClass<P_V,P_M>& domainSet()                        const;
          const uqVectorSetClass<Q_V,Q_M>& imageSet ()                        const;
  virtual       double                     value    (const P_V& domainVector) const = 0;

protected:
  const uqBaseEnvironmentClass&    m_env;
        std::string                m_prefix;
  const uqVectorSetClass<P_V,P_M>& m_domainSet;
  const uqVectorSetClass<Q_V,Q_M>& m_imageSet;
};

template<class P_V,class P_M,class Q_V,class Q_M>
uqBaseScalarFunctionClass<P_V,P_M,Q_V,Q_M>::uqBaseScalarFunctionClass(
  const char*                      prefix,
  const uqVectorSetClass<P_V,P_M>& domainSet,
  const uqVectorSetClass<Q_V,Q_M>& imageSet)
  :
  m_env      (domainSet.env()),
  m_prefix   ((std::string)(prefix)+"func_"),
  m_domainSet(domainSet),
  m_imageSet (imageSet)
{
}

template<class P_V,class P_M,class Q_V,class Q_M>
uqBaseScalarFunctionClass<P_V,P_M,Q_V,Q_M>::~uqBaseScalarFunctionClass()
{
}

template<class P_V,class P_M,class Q_V,class Q_M>
const uqVectorSetClass<P_V,P_M>&
uqBaseScalarFunctionClass<P_V,P_M,Q_V,Q_M>::domainSet() const
{
  return m_domainSet;
}

template<class P_V,class P_M,class Q_V,class Q_M>
const uqVectorSetClass<Q_V,Q_M>&
uqBaseScalarFunctionClass<P_V,P_M,Q_V,Q_M>::imageSet() const
{
  return m_imageSet;
}

//*****************************************************
// Generic class
//*****************************************************
template<class P_V,class P_M,class Q_V,class Q_M>
class uqGenericScalarFunctionClass : public uqBaseScalarFunctionClass<P_V,P_M,Q_V,Q_M> {
public:
  uqGenericScalarFunctionClass(const char*                      prefix,
                               const uqVectorSetClass<P_V,P_M>& domainSet,
                               const uqVectorSetClass<Q_V,Q_M>& imageSet,
                               double (*routinePtr)(const P_V& domainVector, const void* functionDataPtr),
                               const void* functionDataPtr);
  virtual ~uqGenericScalarFunctionClass();

  double value(const P_V& domainVector) const;

protected:
  double (*m_routinePtr)(const P_V& domainVector, const void* functionDataPtr);
  const void* m_routineDataPtr;

  using uqBaseScalarFunctionClass<P_V,P_M,Q_V,Q_M>::m_env;
  using uqBaseScalarFunctionClass<P_V,P_M,Q_V,Q_M>::m_prefix;
  using uqBaseScalarFunctionClass<P_V,P_M,Q_V,Q_M>::m_domainSet;
  using uqBaseScalarFunctionClass<P_V,P_M,Q_V,Q_M>::m_imageSet;
};

template<class P_V,class P_M,class Q_V,class Q_M>
uqGenericScalarFunctionClass<P_V,P_M,Q_V,Q_M>::uqGenericScalarFunctionClass(
  const char*                      prefix,
  const uqVectorSetClass<P_V,P_M>& domainSet,
  const uqVectorSetClass<Q_V,Q_M>& imageSet,
  double (*routinePtr)(const P_V& domainVector, const void* functionDataPtr),
  const void* functionDataPtr)
  :
  uqBaseScalarFunctionClass<P_V,P_M,Q_V,Q_M>(((std::string)(prefix)+"gen").c_str(),
                                             domainSet,
                                             imageSet),
  m_routinePtr    (routinePtr),
  m_routineDataPtr(functionDataPtr)
{
}

template<class P_V,class P_M,class Q_V,class Q_M>
uqGenericScalarFunctionClass<P_V,P_M,Q_V,Q_M>::~uqGenericScalarFunctionClass()
{
}

template<class P_V,class P_M,class Q_V,class Q_M>
double
uqGenericScalarFunctionClass<P_V,P_M,Q_V,Q_M>::value(const P_V& domainVector) const
{
  return m_routinePtr(domainVector, m_routineDataPtr);
}

#endif // __UQ_SCALAR_FUNCTION_H__
