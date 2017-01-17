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

#include <queso/ExponentialScalarCovarianceFunction.h>

namespace QUESO {

// Default constructor -----------------------------
template<class V,class M>
ExponentialScalarCovarianceFunction<V,M>::ExponentialScalarCovarianceFunction(
  const char*                  prefix,
  const VectorSet<V,M>& basicDomainSet,
  double                       sigma,
  double                       a)
  :
  BaseScalarCovarianceFunction<V,M>(prefix,basicDomainSet),
  m_sigma(sigma),
  m_a    (a)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering ExponentialScalarCovarianceFunction<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving ExponentialScalarCovarianceFunction<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor ---------------------------------------
template<class V,class M>
ExponentialScalarCovarianceFunction<V,M>::~ExponentialScalarCovarianceFunction()
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering ExponentialScalarCovarianceFunction<V,M>::destructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving ExponentialScalarCovarianceFunction<V,M>::destructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Math methods -------------------------------------
template<class V,class M>
double
ExponentialScalarCovarianceFunction<V,M>::value(const V& domainVector1, const V& domainVector2) const
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering ExponentialScalarCovarianceFunction<V,M>::value()"
                          << std::endl;
  }

  double result = 0.;

  double exponent = -(domainVector1 - domainVector2).norm2Sq()/(m_sigma*m_sigma);

  result = m_a*std::exp(exponent);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving ExponentialScalarCovarianceFunction<V,M>::value()"
                          << std::endl;
  }

  return result;
}

}  // End namespace QUESO
