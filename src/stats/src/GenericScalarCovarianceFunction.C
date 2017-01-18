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

#include <queso/GenericScalarCovarianceFunction.h>

namespace QUESO {

// Default constructor -----------------------------
template<class V,class M>
GenericScalarCovarianceFunction<V,M>::GenericScalarCovarianceFunction(
  const char*                  prefix,
  const VectorSet<V,M>& domainSet,
  double (*covRoutinePtr)(const V& positionVector1, const V& positionVector2, const void* routineDataPtr),
  const void*                  routinesDataPtr)
  :
  BaseScalarCovarianceFunction<V,M>(prefix,domainSet),
  m_covRoutinePtr                         (covRoutinePtr),
  m_routineDataPtr                        (routinesDataPtr)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering GenericScalarCovarianceFunction<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving GenericScalarCovarianceFunction<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor ---------------------------------------
template<class V,class M>
GenericScalarCovarianceFunction<V,M>::~GenericScalarCovarianceFunction()
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering GenericScalarCovarianceFunction<V,M>::destructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving GenericScalarCovarianceFunction<V,M>::destructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Math methods -------------------------------------
template<class V,class M>
double
GenericScalarCovarianceFunction<V,M>::value(const V& positionVector1, const V& positionVector2) const
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering GenericScalarCovarianceFunction<V,M>::value()"
                            << std::endl;
  }

  queso_require_msg(m_covRoutinePtr, "m_covRoutinePtr = NULL");

  double result = m_covRoutinePtr(positionVector1, positionVector2, m_routineDataPtr);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving GenericScalarCovarianceFunction<V,M>::value()"
                            << std::endl;
  }

  return result;
}

}  // End namespace QUESO
