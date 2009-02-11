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

#ifndef __UQ_VECTOR_FUNCTION_H__
#define __UQ_VECTOR_FUNCTION_H__

#include <uqVectorSet.h>
#include <uqEnvironment.h>
#include <uqDefines.h>
#include <EpetraExt_DistArray.h>

//*****************************************************
// Base class
//*****************************************************
template<class P_V,class P_M,class Q_V,class Q_M>
class uqBaseVectorFunctionClass {
public:
           uqBaseVectorFunctionClass(const char*                      prefix,
                                     const uqVectorSetClass<P_V,P_M>& domainSet,
                                     const uqVectorSetClass<Q_V,Q_M>& imageSet);
  virtual ~uqBaseVectorFunctionClass();

          const uqVectorSetClass<P_V,P_M>& domainSet() const;
          const uqVectorSetClass<Q_V,Q_M>& imageSet () const;
  virtual       void                       compute  (const P_V&                        domainVector,
                                                     const P_V*                        domainDirection,
                                                           Q_V&                        imageVector,
                                                           EpetraExt::DistArray<P_V*>* gradVectors,     // Yes, 'P_V'
                                                           EpetraExt::DistArray<P_M*>* hessianMatrices, // Yes, 'P_M'
                                                           EpetraExt::DistArray<P_V*>* hessianEffects) const = 0;

protected:
  const uqBaseEnvironmentClass&    m_env;
        std::string                m_prefix;
  const uqVectorSetClass<P_V,P_M>& m_domainSet;
  const uqVectorSetClass<Q_V,Q_M>& m_imageSet;
};

template<class P_V,class P_M,class Q_V,class Q_M>
uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M>::uqBaseVectorFunctionClass(
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
uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M>::~uqBaseVectorFunctionClass()
{
}

template<class P_V,class P_M,class Q_V,class Q_M>
const uqVectorSetClass<P_V,P_M>&
uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M>::domainSet() const
{
  return m_domainSet;
}

template<class P_V,class P_M,class Q_V,class Q_M>
const uqVectorSetClass<Q_V,Q_M>&
uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M>::imageSet() const
{
  return m_imageSet;
}

//*****************************************************
// Generic class
//*****************************************************
template<class P_V,class P_M,class Q_V,class Q_M>
class uqGenericVectorFunctionClass : public uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M> {
public:
  uqGenericVectorFunctionClass(const char*                      prefix,
                               const uqVectorSetClass<P_V,P_M>& domainSet,
                               const uqVectorSetClass<Q_V,Q_M>& imageSet,
                               void (*routinePtr)(const P_V&                        domainVector,
                                                  const P_V*                        domainDirection,
                                                  const void*                       functionDataPtr,
                                                        Q_V&                        imageVector,
                                                        EpetraExt::DistArray<P_V*>* gradVectors,
                                                        EpetraExt::DistArray<P_M*>* hessianMatrices,
                                                        EpetraExt::DistArray<P_V*>* hessianEffects),
                               const void* functionDataPtr);
  virtual ~uqGenericVectorFunctionClass();

  void compute  (const P_V&                        domainVector,
                 const P_V*                        domainDirection,
                       Q_V&                        imageVector,
                       EpetraExt::DistArray<P_V*>* gradVectors,     // Yes, 'P_V'
                       EpetraExt::DistArray<P_M*>* hessianMatrices, // Yes, 'P_M'
                       EpetraExt::DistArray<P_V*>* hessianEffects) const;

protected:
  void (*m_routinePtr)(const P_V&                        domainVector,
                       const P_V*                        domainDirection,
                       const void*                       functionDataPtr,
                             Q_V&                        imageVector,
                             EpetraExt::DistArray<P_V*>* gradVectors,
                             EpetraExt::DistArray<P_M*>* hessianMatrices,
                             EpetraExt::DistArray<P_V*>* hessianEffects);
  const void* m_routineDataPtr;

  using uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M>::m_env;
  using uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M>::m_prefix;
  using uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M>::m_domainSet;
  using uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M>::m_imageSet;
};

template<class P_V,class P_M,class Q_V,class Q_M>
uqGenericVectorFunctionClass<P_V,P_M,Q_V,Q_M>::uqGenericVectorFunctionClass(
  const char*                      prefix,
  const uqVectorSetClass<P_V,P_M>& domainSet,
  const uqVectorSetClass<Q_V,Q_M>& imageSet,
  void (*routinePtr)(const P_V&                        domainVector,
                     const P_V*                        domainDirection,
                     const void*                       functionDataPtr,
                           Q_V&                        imageVector,
                           EpetraExt::DistArray<P_V*>* gradVectors,
                           EpetraExt::DistArray<P_M*>* hessianMatrices,
                           EpetraExt::DistArray<P_V*>* hessianEffects),
  const void* functionDataPtr)
  :
  uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M>(((std::string)(prefix)+"gen").c_str(),
                                             domainSet,
                                             imageSet),
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
uqGenericVectorFunctionClass<P_V,P_M,Q_V,Q_M>::compute(
  const P_V&                        domainVector,
  const P_V*                        domainDirection,
        Q_V&                        imageVector,
        EpetraExt::DistArray<P_V*>* gradVectors,     // Yes, 'P_V'
        EpetraExt::DistArray<P_M*>* hessianMatrices, // Yes, 'P_M'
        EpetraExt::DistArray<P_V*>* hessianEffects) const
{
  //UQ_FATAL_TEST_MACRO(false,
  //                    domainVector.env().rank(),
  //                    "uqGenericVectorFunctionClass<P_V,P_M,Q_V,Q_M>::compute()",
  //                    "this method should not be called in the case of this class");

  m_routinePtr(domainVector, domainDirection, m_routineDataPtr, imageVector, gradVectors, hessianMatrices, hessianEffects);

  return;
}

#endif // __UQ_VECTOR_FUNCTION_H__
