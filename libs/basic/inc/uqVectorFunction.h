/* uq/libs/queso/inc/uqVectorFunction.h
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

#ifndef __UQ_VECTOR_FUNCTION_H__
#define __UQ_VECTOR_FUNCTION_H__

#include <uqEnvironment.h>
#include <uqDefines.h>

//*****************************************************
// Base class
//*****************************************************
template<class P_V,class P_M,class Q_V,class Q_M>
class uqBaseVectorFunctionClass {
public:
           uqBaseVectorFunctionClass(const char*                        prefix,
                                     const uqVectorSpaceClass<P_V,P_M>& domainSpace,
                                     const uqVectorSpaceClass<Q_V,Q_M>& imageSpace);
  virtual ~uqBaseVectorFunctionClass();

          const uqVectorSpaceClass<P_V,P_M>& domainSpace()                                          const;
          const uqVectorSpaceClass<Q_V,Q_M>& imageSpace ()                                          const;
  virtual       void                         compute    (const P_V& domainVector, Q_V& imageVector) const = 0;

protected:
  const uqEnvironmentClass&          m_env;
        std::string                  m_prefix;
  const uqVectorSpaceClass<P_V,P_M>& m_domainSpace;
  const uqVectorSpaceClass<Q_V,Q_M>& m_imageSpace;
};

template<class P_V,class P_M,class Q_V,class Q_M>
uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M>::uqBaseVectorFunctionClass(
  const char*                        prefix,
  const uqVectorSpaceClass<P_V,P_M>& domainSpace,
  const uqVectorSpaceClass<Q_V,Q_M>& imageSpace)
  :
  m_env           (domainSpace.env()),
  m_prefix        ((std::string)(prefix)+"func_"),
  m_domainSpace   (domainSpace),
  m_imageSpace    (imageSpace)
{
}

template<class P_V,class P_M,class Q_V,class Q_M>
uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M>::~uqBaseVectorFunctionClass()
{
}

template<class P_V,class P_M,class Q_V,class Q_M>
const uqVectorSpaceClass<P_V,P_M>&
uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M>::domainSpace() const
{
  return m_domainSpace;
}

template<class P_V,class P_M,class Q_V,class Q_M>
const uqVectorSpaceClass<Q_V,Q_M>&
uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M>::imageSpace() const
{
  return m_imageSpace;
}

//*****************************************************
// Generic class
//*****************************************************
template<class P_V,class P_M,class Q_V,class Q_M>
class uqGenericVectorFunctionClass : public uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M> {
public:
  uqGenericVectorFunctionClass(const char*                        prefix,
                               const uqVectorSpaceClass<P_V,P_M>& domainSpace,
                               const uqVectorSpaceClass<Q_V,Q_M>& imageSpace,
                               void (*routinePtr)(const P_V& domainVector, const void* functionDataPtr, Q_V& imageVector),
                               const void* functionDataPtr);
  virtual ~uqGenericVectorFunctionClass();

  void compute(const P_V& domainVector, Q_V& imageVector) const;

protected:
  void (*m_routinePtr)(const P_V& domainVector, const void* functionDataPtr, Q_V& imageVector);
  const void* m_routineDataPtr;

  using uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M>::m_env;
  using uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M>::m_prefix;
  using uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M>::m_domainSpace;
  using uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M>::m_imageSpace;
};

template<class P_V,class P_M,class Q_V,class Q_M>
uqGenericVectorFunctionClass<P_V,P_M,Q_V,Q_M>::uqGenericVectorFunctionClass(
  const char*                        prefix,
  const uqVectorSpaceClass<P_V,P_M>& domainSpace,
  const uqVectorSpaceClass<Q_V,Q_M>& imageSpace,
  void (*routinePtr)(const P_V& domainVector, const void* functionDataPtr, Q_V& imageVector),
  const void* functionDataPtr)
  :
  uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M>(((std::string)(prefix)+"gen").c_str(),
                                             domainSpace,
                                             imageSpace),
  m_routinePtr    (routinePtr),
  m_routineDataPtr(functionDataPtr)
{
}

template<class P_V,class P_M,class Q_V,class Q_M>
uqGenericVectorFunctionClass<P_V,P_M,Q_V,Q_M>::~uqGenericVectorFunctionClass()
{
}

template<class P_V,class P_M,class Q_V,class Q_M>
void
uqGenericVectorFunctionClass<P_V,P_M,Q_V,Q_M>::compute(const P_V& domainVector, Q_V& imageVector) const
{
  //UQ_FATAL_TEST_MACRO(false,
  //                    domainVector.env().rank(),
  //                    "uqGenericVectorFunctionClass<P_V,P_M,Q_V,Q_M>::compute()",
  //                    "this method should not be called in the case of this class");

  m_routinePtr(domainVector, m_routineDataPtr, imageVector);

  return;
}

#endif // __UQ_VECTOR_FUNCTION_H__
