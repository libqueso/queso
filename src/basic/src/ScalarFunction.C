//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2015 The PECOS Development Team
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
#include <queso/ScalarFunction.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Default constructor
template<class V, class M>
BaseScalarFunction<V, M>::BaseScalarFunction(const char * prefix,
    const VectorSet<V, M> & domainSet)
  : m_env(domainSet.env()),
    m_prefix((std::string)(prefix) + "func_"),
    m_domainSet(domainSet)
{
}

// Destructor
template<class V, class M>
BaseScalarFunction<V, M>::~BaseScalarFunction()
{
}

// Math methods
template<class V, class M>
const VectorSet<V, M> & BaseScalarFunction<V, M>::domainSet() const
{
  return m_domainSet;
}

}  // End namespace QUESO

template class QUESO::BaseScalarFunction<QUESO::GslVector, QUESO::GslMatrix>;
