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

#include <queso/ScalarCdf.h>

namespace QUESO {

// Default constructor -----------------------------
template<class T>
BaseScalarCdf<T>::BaseScalarCdf(
  const BaseEnvironment& env,
  const char*                   prefix)
  :
  m_env          (env),
  m_prefix       ((std::string)(prefix)+""),
  m_minHorizontal(-INFINITY),
  m_maxHorizontal( INFINITY)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering BaseScalarCdf<T>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving BaseScalarCdf<T>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }
}
// Destructor ---------------------------------------
template<class T>
BaseScalarCdf<T>::~BaseScalarCdf()
{
}
// Environment methods ------------------------------
template <class T>
const BaseEnvironment&
BaseScalarCdf<T>::env() const
{
  return m_env;
}
// --------------------------------------------------
template <class T>
const std::string&
BaseScalarCdf<T>::prefix() const
{
  return m_prefix;
}
// I/O methods---------------------------------------
template<class T>
void
BaseScalarCdf<T>::subWriteContents(
  const std::string&            varNamePrefix,
  const std::string&            fileName,
  const std::string&            fileType,
  const std::set<unsigned int>& allowedSubEnvIds) const
{
  std::cerr << "WARNING: BaseScalarCdf<T>::subWriteContents() being used..."
            << std::endl;
  return;
}

//---------------------------------------------------
// Method outside either class definition------------
//---------------------------------------------------

//*****************************************************
// Horizontal distance
//*****************************************************

//! It calculated the maximum horizontal distance between two CDFs.
template <class T>
double
horizontalDistance(const BaseScalarCdf<T>& cdf1,
                   const BaseScalarCdf<T>& cdf2,
                   double epsilon)
{
  double maxDistance     = 0.;
  double xForMaxDistance = 0.;

  double x1 = cdf1.inverse(epsilon*.5);
  double x2 = cdf1.inverse(1.-epsilon*.5);
  if (cdf1.env().subDisplayFile()) {
    *cdf1.env().subDisplayFile() << "In horizontalDistance()"
                                 << ", cdf1.prefix() = " << cdf1.prefix()
                                 << ", cdf2.prefix() = " << cdf2.prefix()
                                 << ", epsilon = "       << epsilon
                                 << ": x1 = "            << x1
                                 << ", x2 = "            << x2
                                 << std::endl;
  }

  //if (cdf1.env().subDisplayFile()) {
  //  *cdf1.env().subDisplayFile() << "In horizontalDistance: x1 = " << x1
  //                              << ", x2 = " << x2
  //                              << std::endl;
  //}

  double numEvaluationPoints = 1001.;
  for (double i = 0.; i < numEvaluationPoints; ++i) {
    double ratio = i/(numEvaluationPoints-1.); // IMPORTANT: Yes, '-1.'
    double x = (1.-ratio)*x1 + ratio*x2;
    double y = cdf2.inverse(cdf1.value(x));
    double d = fabs(x-y);
    if ((cdf1.env().subDisplayFile()) && (cdf1.env().displayVerbosity() >= 3)) {
      *cdf1.env().subDisplayFile() << "In horizontalDistance"
                                   << ": i = "                  << i
                                   << ", x = "                  << x
                                   << ", cdf1.value(x) = "      << cdf1.value(x)
                                   << ", y = "                  << y
                                   << ", d = "                  << d
                                   << ", currentMaxDistance = " << maxDistance
                                   << std::endl;
    }
    if (maxDistance < d) {
      maxDistance     = d;
      xForMaxDistance = x;
      if ((cdf1.env().subDisplayFile()) && (cdf1.env().displayVerbosity() >= 3)) {
        *cdf1.env().subDisplayFile() << "In horizontalDistance"
                                     << ": i = "               << i
                                     << ", NOW maxDistance = " << maxDistance
                                     << ", xForMaxDistance = " << xForMaxDistance
                                     << std::endl;
      }
    }
  }

  if (cdf1.env().subDisplayFile()) {
    *cdf1.env().subDisplayFile() << "In horizontalDistance()"
                                 << ", cdf1.prefix() = "   << cdf1.prefix()
                                 << ", cdf2.prefix() = "   << cdf2.prefix()
                                 << ", epsilon = "         << epsilon
                                 << ": maxDistance = "     << maxDistance
                                 << ", xForMaxDistance = " << xForMaxDistance
                                 << std::endl;
  }

  return maxDistance;
}

}  // End namespace QUESO
