//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2017 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-

#include <queso/GenericVectorFunction.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Default constructor -----------------------------
template<class P_V,class P_M,class Q_V,class Q_M>
GenericVectorFunction<P_V,P_M,Q_V,Q_M>::GenericVectorFunction(
  const char*                      prefix,
  const VectorSet<P_V,P_M>& domainSet,
  const VectorSet<Q_V,Q_M>& imageSet,
  void (*routinePtr)(const P_V&                    domainVector,
                     const P_V*                    domainDirection,
                     const void*                   functionDataPtr,
                           Q_V&                    imageVector,
                           DistArray<P_V*>* gradVectors,
                           DistArray<P_M*>* hessianMatrices,
                           DistArray<P_V*>* hessianEffects),
  const void* functionDataPtr)
  :
  BaseVectorFunction<P_V,P_M,Q_V,Q_M>(((std::string)(prefix)+"gen").c_str(),
                                             domainSet,
                                             imageSet),
  m_routinePtr    (routinePtr),
  m_routineDataPtr(functionDataPtr)
{
}

// Destructor ---------------------------------------
template<class P_V,class P_M,class Q_V,class Q_M>
GenericVectorFunction<P_V,P_M,Q_V,Q_M>::~GenericVectorFunction()
{
}

// Math methods -------------------------------------
template<class P_V,class P_M,class Q_V,class Q_M>
void
GenericVectorFunction<P_V,P_M,Q_V,Q_M>::compute(
  const P_V&                    domainVector,
  const P_V*                    domainDirection,
        Q_V&                    imageVector,
        DistArray<P_V*>* gradVectors,     // Yes, 'P_V'
        DistArray<P_M*>* hessianMatrices, // Yes, 'P_M'
        DistArray<P_V*>* hessianEffects) const
{

  //                    domainVector.env().worldRank(),
  //                    "GenericVectorFunction<P_V,P_M,Q_V,Q_M>::compute()",
  //                    "this method should not be called in the case of this class");

  m_routinePtr(domainVector, domainDirection, m_routineDataPtr, imageVector, gradVectors, hessianMatrices, hessianEffects);

  return;
}

}  // End namespace QUESO

template class QUESO::GenericVectorFunction<QUESO::GslVector, QUESO::GslMatrix, QUESO::GslVector, QUESO::GslMatrix>;
