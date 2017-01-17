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

#include <queso/VectorSet.h>
#include <queso/VectorSubset.h>
#include <queso/Environment.h>
#include <queso/Defines.h>
#include <queso/ConstantScalarFunction.h>

namespace QUESO {

// Default constructor
template<class V,class M>
ConstantScalarFunction<V,M>::ConstantScalarFunction(const char* prefix,
    const VectorSet<V,M>& domainSet,
    double constantValue)
  : BaseScalarFunction<V,M>(((std::string)(prefix)+"gen").c_str(), domainSet),
    m_constantValue(constantValue)
{
}

// Destructor
template<class V,class M>
ConstantScalarFunction<V,M>::~ConstantScalarFunction()
{
}

// Math methods
template<class V,class M>
double ConstantScalarFunction<V,M>::actualValue(const V& domainVector,
    const V* domainDirection,
    V* gradVector,
    M* hessianMatrix,
    V* hessianEffect) const
{
  return m_constantValue;
}


template<class V,class M>
double ConstantScalarFunction<V,M>::lnValue(const V& domainVector,
    const V* domainDirection,
    V* gradVector,
    M* hessianMatrix,
    V* hessianEffect) const
{
  return 0.;
}

}  // End namespace QUESO
