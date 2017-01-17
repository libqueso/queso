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

#include <queso/ConstantVectorFunction.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Default constructor -----------------------------
template<class P_V,class P_M,class Q_V,class Q_M>
ConstantVectorFunction<P_V,P_M,Q_V,Q_M>::ConstantVectorFunction(
  const char*                      prefix,
  const VectorSet<P_V,P_M>& domainSet,
  const VectorSet<Q_V,Q_M>& imageSet,
  const Q_V&                       constantImageVector)
  :
  BaseVectorFunction<P_V,P_M,Q_V,Q_M>(((std::string)(prefix)+"gen").c_str(),
                                             domainSet,
                                             imageSet),
  m_constantImageVector(NULL)
{
  m_constantImageVector = new Q_V(constantImageVector);
}

// Destructor ---------------------------------------
template<class P_V,class P_M,class Q_V,class Q_M>
ConstantVectorFunction<P_V,P_M,Q_V,Q_M>::~ConstantVectorFunction()
{
  delete m_constantImageVector;
}

// Math methods -------------------------------------
template<class P_V,class P_M,class Q_V,class Q_M>
void
ConstantVectorFunction<P_V,P_M,Q_V,Q_M>::compute(
  const P_V&                    domainVector,
  const P_V*                    domainDirection,
        Q_V&                    imageVector,
        DistArray<P_V*>* gradVectors,     // Yes, 'P_V'
        DistArray<P_M*>* hessianMatrices, // Yes, 'P_M'
        DistArray<P_V*>* hessianEffects) const
{
  queso_require_msg(m_constantImageVector, "m_constantImageVector is NULL");

  imageVector = *m_constantImageVector;

  return;
}

}  // End namespace QUESO

template class QUESO::ConstantVectorFunction<QUESO::GslVector, QUESO::GslMatrix, QUESO::GslVector, QUESO::GslMatrix>;
