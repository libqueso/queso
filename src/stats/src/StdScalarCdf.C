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

#include <queso/StdScalarCdf.h>

namespace QUESO {

// Default constructor -----------------------------
template<class T>
StdScalarCdf<T>::StdScalarCdf(
  const BaseEnvironment& env,
  const char*                   prefix,
  const std::vector<T>&         cdfGrid,
  const std::vector<double>&    cdfValues)
  :
  BaseScalarCdf<T>(env,((std::string)(prefix)+"").c_str()),
  m_cdfGrid              (env,prefix,cdfGrid),
  m_cdfValues            (cdfValues),
  m_sampledCdfGrid       (NULL)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering StdScalarCdf<T>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  m_sampledCdfGrid = new SampledScalarCdf<T>(env,prefix,m_cdfGrid,m_cdfValues);

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving StdScalarCdf<T>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor ---------------------------------------
template<class T>
StdScalarCdf<T>::~StdScalarCdf()
{
  delete m_sampledCdfGrid;
}
// Math methods--------------------------------------
template<class T>
double
StdScalarCdf<T>::value(T paramValue) const
{
  return m_sampledCdfGrid->value(paramValue);
}
//---------------------------------------------------
template<class T>
T
StdScalarCdf<T>::inverse(double cdfValue) const
{
  return m_sampledCdfGrid->inverse(cdfValue);
}
//---------------------------------------------------
template<class T>
void
StdScalarCdf<T>::getSupport(T& minHorizontal, T& maxHorizontal) const
{
  return m_sampledCdfGrid->getSupport(minHorizontal,maxHorizontal);
}
// I/O methods---------------------------------------
template <class T>
void
StdScalarCdf<T>::print(std::ostream& os) const
{
  m_sampledCdfGrid->print(os);
  return;
}
//---------------------------------------------------
template<class T>
void
StdScalarCdf<T>::subWriteContents(
  const std::string&            varNamePrefix,
  const std::string&            fileName,
  const std::string&            fileType,
  const std::set<unsigned int>& allowedSubEnvIds) const
{
  m_sampledCdfGrid->subWriteContents(varNamePrefix,fileName,fileType,allowedSubEnvIds);
  return;
}

}  // End namespace QUESO
