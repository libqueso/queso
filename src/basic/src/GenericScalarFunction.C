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

#include <cmath>

#include <queso/Defines.h>
#include <queso/VectorSet.h>
#include <queso/VectorSubset.h>
#include <queso/Environment.h>
#include <queso/GenericScalarFunction.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Default constructor
template<class V,class M>
GenericScalarFunction<V,M>::GenericScalarFunction(const char* prefix,
    const VectorSet<V,M>& domainSet,
    double (*valueRoutinePtr)(const V& domainVector, const V* domainDirection, const void* routinesDataPtr, V* gradVector, M* hessianMatrix, V* hessianEffect),
    const void* routinesDataPtr,
    bool routineIsForLn)
  : BaseScalarFunction<V,M>(((std::string)(prefix)+"gen").c_str(), domainSet),
    m_valueRoutinePtr             (valueRoutinePtr),
    m_routinesDataPtr             (routinesDataPtr),
    m_routineIsForLn              (routineIsForLn)
{
}

// Destructor
template<class V,class M>
GenericScalarFunction<V,M>::~GenericScalarFunction()
{
}

// Math methods
template<class V,class M>
double GenericScalarFunction<V,M>::actualValue(const V& domainVector,
    const V* domainDirection,
    V* gradVector,
    M* hessianMatrix,
    V* hessianEffect) const
{
  queso_require_msg(m_valueRoutinePtr, "m_valueRoutinePtr = NULL");

  double value = m_valueRoutinePtr(domainVector, domainDirection, m_routinesDataPtr, gradVector, hessianMatrix, hessianEffect);
  if (m_routineIsForLn) {
    value = std::exp(value);

    if ((domainDirection != NULL) ||
        (gradVector      != NULL) ||
        (hessianMatrix   != NULL) ||
        (hessianEffect   != NULL))
      queso_not_implemented();
  }
  return value;
}


template<class V,class M>
double GenericScalarFunction<V,M>::lnValue(const V& domainVector,
    const V* domainDirection,
    V* gradVector,
    M* hessianMatrix,
    V* hessianEffect) const
{
  queso_require_msg(m_valueRoutinePtr, "m_valueRoutinePtr = NULL");

  double value = m_valueRoutinePtr(domainVector, domainDirection, m_routinesDataPtr, gradVector, hessianMatrix, hessianEffect);
  if (m_routineIsForLn == false) {
    value = log(value);

    if ((domainDirection != NULL) ||
        (gradVector      != NULL) ||
        (hessianMatrix   != NULL) ||
        (hessianEffect   != NULL))
      queso_not_implemented();
  }
  return value;
}

}  // End namespace QUESO

template class QUESO::GenericScalarFunction<QUESO::GslVector, QUESO::GslMatrix>;
