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

#include <queso/VectorFunction.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Default constructor -----------------------------
template<class P_V,class P_M,class Q_V,class Q_M>
BaseVectorFunction<P_V,P_M,Q_V,Q_M>::BaseVectorFunction(
  const char*                      prefix,
  const VectorSet<P_V,P_M>& domainSet,
  const VectorSet<Q_V,Q_M>& imageSet)
  :
  m_env      (domainSet.env()),
  m_prefix   ((std::string)(prefix)+"func_"),
  m_domainSet(domainSet),
  m_imageSet (imageSet)
{
}

// Destructor ---------------------------------------
template<class P_V,class P_M,class Q_V,class Q_M>
BaseVectorFunction<P_V,P_M,Q_V,Q_M>::~BaseVectorFunction()
{
}

// Math methods -------------------------------------
template<class P_V,class P_M,class Q_V,class Q_M>
const VectorSet<P_V,P_M>&
BaseVectorFunction<P_V,P_M,Q_V,Q_M>::domainSet() const
{
  return m_domainSet;
}


template<class P_V,class P_M,class Q_V,class Q_M>
const VectorSet<Q_V,Q_M>&
BaseVectorFunction<P_V,P_M,Q_V,Q_M>::imageSet() const
{
  return m_imageSet;
}

}  // End namespace QUESO

template class QUESO::BaseVectorFunction<QUESO::GslVector, QUESO::GslMatrix, QUESO::GslVector, QUESO::GslMatrix>;
