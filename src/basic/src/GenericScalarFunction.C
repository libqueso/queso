//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
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
// 
// $Id$
//
//--------------------------------------------------------------------------

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
  UQ_FATAL_TEST_MACRO(m_valueRoutinePtr == NULL,
                      m_env.worldRank(),
                      "GenericScalarFunction<V,M>::actualValue()",
                      "m_valueRoutinePtr = NULL");

  double value = m_valueRoutinePtr(domainVector, domainDirection, m_routinesDataPtr, gradVector, hessianMatrix, hessianEffect);
  if (m_routineIsForLn) {
#ifdef QUESO_EXPECTS_LN_LIKELIHOOD_INSTEAD_OF_MINUS_2_LN
    value = std::exp(value);
#else
    value = std::exp(-.5*value);
#endif
    UQ_FATAL_TEST_MACRO((domainDirection != NULL) ||
                        (gradVector      != NULL) ||
                        (hessianMatrix   != NULL) ||
                        (hessianEffect   != NULL),
                        m_env.worldRank(),
                        "GenericScalarFunction<V,M>::gradOfActual()",
                        "INCOMPLETE CODE");
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
  UQ_FATAL_TEST_MACRO(m_valueRoutinePtr == NULL,
                      m_env.worldRank(),
                      "GenericScalarFunction<V,M>::lnValue()",
                      "m_valueRoutinePtr = NULL");

  double value = m_valueRoutinePtr(domainVector, domainDirection, m_routinesDataPtr, gradVector, hessianMatrix, hessianEffect);
  if (m_routineIsForLn == false) {
#ifdef QUESO_EXPECTS_LN_LIKELIHOOD_INSTEAD_OF_MINUS_2_LN
    value = log(value);
#else
    value = -2.*log(value);
#endif
    UQ_FATAL_TEST_MACRO((domainDirection != NULL) ||
                        (gradVector      != NULL) ||
                        (hessianMatrix   != NULL) ||
                        (hessianEffect   != NULL),
                        m_env.worldRank(),
                        "GenericScalarFunction<V,M>::gradOfLn()",
                        "INCOMPLETE CODE");
  }
  return value;
}

}  // End namespace QUESO

template class QUESO::GenericScalarFunction<QUESO::GslVector, QUESO::GslMatrix>;
