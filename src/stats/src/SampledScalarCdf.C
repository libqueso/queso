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

#include <queso/SampledScalarCdf.h>

namespace QUESO {

// Default constructor -----------------------------
template<class T>
SampledScalarCdf<T>::SampledScalarCdf(
  const BaseEnvironment& env,
  const char*                   prefix,
  const BaseOneDGrid<T>& cdfGrid,
  const std::vector<double>&    cdfValues)
  :
  BaseScalarCdf<T>(env,((std::string)(prefix)+"").c_str()),
  m_cdfGrid              (cdfGrid  ),
  m_cdfValues            (cdfValues)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering SampledScalarCdf<T>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  //for (unsigned int i = 0; i < m_cdfValues.size(); ++i) {
  //  m_sortedCdfValues[i] = m_cdfValues[i];
  //}
  //std::sort(m_sortedCdfValues.begin(), m_sortedCdfValues.end());

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving SampledScalarCdf<T>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor ---------------------------------------
template<class T>
SampledScalarCdf<T>::~SampledScalarCdf()
{
}
// Math methods--------------------------------------
template<class T>
double
SampledScalarCdf<T>::value(T paramValue) const
{
  double result = 0.;
  if (paramValue <= m_cdfGrid[0]) {
    result = 0.;
  }
  else if (m_cdfGrid[m_cdfGrid.size()-1] <= paramValue) {
    result = 1.;
  }
  else {
    unsigned int intervalId = m_cdfGrid.findIntervalId(paramValue);
    queso_require_less_msg(intervalId, (m_cdfGrid.size()-1), "invalid intervalId");

    double intervalLen = m_cdfGrid[intervalId+1] - m_cdfGrid[intervalId];
    double ratio = (paramValue - m_cdfGrid[intervalId])/intervalLen;
#if 0
    *m_env.subDisplayFile() << "In SampledScalarCdf::value()"
                            << ": paramValue = "              << paramValue
                            << ", intervalId = "              << intervalId
                            << ", cdfGrid.size() = "          << m_cdfGrid.size()
                            << ", m_cdfGrid[intervalId] = "   << m_cdfGrid[intervalId]
                            << ", m_cdfGrid[intervalId+1] = " << m_cdfGrid[intervalId+1]
                            << ", intervalLen = "             << intervalLen
                            << ", ratio = "                   << ratio
                            << std::endl;
#endif
    queso_require_greater_equal_msg(ratio, 0., "invalid ratio");

    result = (1.-ratio)*m_cdfValues[intervalId] + ratio*m_cdfValues[intervalId+1];
  }

  return result;
}
//---------------------------------------------------
template<class T>
T
SampledScalarCdf<T>::inverse(double cdfValue) const
{
  //*m_env.subDisplayFile() << "In SampledScalarCdf::inverse(): cdfValue = " << cdfValue
  //                       << std::endl;
  queso_require_msg(!((cdfValue < 0.) || (1. < cdfValue)), "invalid cdfValue");
  double result = 0.;
  unsigned int i = 0;
  unsigned int j = m_cdfValues.size()-1;
  bool searchPosition = true;
  do {
    if (cdfValue == m_cdfValues[i]) {
      while ((0 < i) && (cdfValue == m_cdfValues[i-1])) --i;
      result = m_cdfGrid[i];
      searchPosition = false;
    }

    if (cdfValue == m_cdfValues[j]) {
      while ((0 < j) && (cdfValue == m_cdfValues[j-1])) --j;
      result = m_cdfGrid[j];
      searchPosition = false;
    }

    if ((j-i) <= 0) {
      queso_error_msg("invalid pair of values 'i' and 'j'");
    }
    else if ((j-i) == 1) {
      double ratio = (cdfValue-m_cdfValues[i])/(m_cdfValues[j]-m_cdfValues[i]);
      result = (1.-ratio)*m_cdfGrid[i] + ratio*m_cdfGrid[j];
      searchPosition = false;
    }
    else {
      unsigned int k= (unsigned int) ((i+j)*.5);
      if (cdfValue < m_cdfValues[k]) {
        j = k;
      }
      else if (cdfValue == m_cdfValues[k]) {
        while ((0 < k) && (cdfValue == m_cdfValues[k-1])) --k;
        result = m_cdfGrid[k];
        searchPosition = false;
      }
      else {
        i = k;
      }
    }
  } while (searchPosition);

  return result;
}
//---------------------------------------------------
template<class T>
void
SampledScalarCdf<T>::getSupport(T& minHorizontal, T& maxHorizontal) const
{
  if ((m_minHorizontal == -INFINITY) ||
      (m_maxHorizontal ==  INFINITY)) {
    queso_require_msg(!((m_minHorizontal != -INFINITY) || (m_maxHorizontal != INFINITY)), "unexpected values of m_minHorizontal and/or m_maxHorizontal");

    unsigned int iMax = m_cdfGrid.size();

    for (unsigned int i = 0; i < iMax; ++i) {
      if (m_cdfValues[i] > 0.) {
        if (i > 0) --i;
        m_minHorizontal = m_cdfGrid[i];
        break;
      }
    }

    queso_require_not_equal_to_msg(m_minHorizontal, -INFINITY, "unexpected value for m_minHorizontal");

    if (iMax == 1) {
      queso_require_equal_to_msg(m_cdfValues[iMax - 1], 1., "unexpected value for case 'iMax = 1'");
      m_maxHorizontal = m_cdfGrid[iMax-1];
    }
    else for (unsigned int i = 0; i < iMax; ++i) {
      if (m_cdfValues[iMax-1-i] < 1.) {
        if (i > 0) --i;
        m_maxHorizontal = m_cdfGrid[iMax-1-i];
        break;
      }
    }

    queso_require_not_equal_to_msg(m_maxHorizontal, INFINITY, "unexpected value for m_maxHorizontal");
  }

  minHorizontal = m_minHorizontal;
  maxHorizontal = m_maxHorizontal;

  return;
}
// I/O methods---------------------------------------
template <class T>
void
SampledScalarCdf<T>::print(std::ostream& os) const
{
  // Print values *of* grid points
  os << m_cdfGrid;

  // Print *cdf* values *at* grid points
  os << m_prefix << "values_sub" << m_env.subIdString() << " = zeros(" << m_cdfValues.size()
     << ","                                                            << 1
     << ");"
     << std::endl;
  os << m_prefix << "values_sub" << m_env.subIdString() << " = [";
  for (unsigned int j = 0; j < m_cdfValues.size(); ++j) {
    os << m_cdfValues[j] << " ";
  }
  os << "];"
     << std::endl;

  return;
}
//---------------------------------------------------
template<class T>
void
SampledScalarCdf<T>::subWriteContents(
  const std::string&            varNamePrefix,
  const std::string&            fileName,
  const std::string&            fileType,
  const std::set<unsigned int>& allowedSubEnvIds) const
{
  queso_require_greater_equal_msg(m_env.subRank(), 0, "unexpected subRank");

  FilePtrSetStruct filePtrSet;
  if (m_env.openOutputFile(fileName,
                           fileType, // "m or hdf"
                           allowedSubEnvIds,
                           false,
                           filePtrSet)) {

    // Grid
    *filePtrSet.ofsVar << varNamePrefix << "grid_sub" << m_env.subIdString() << " = zeros(" << m_cdfGrid.size()
                       << ","                                                               << 1
                       << ");"
                       << std::endl;
    *filePtrSet.ofsVar << varNamePrefix << "grid_sub" << m_env.subIdString() << " = [";

    unsigned int savedPrecision = filePtrSet.ofsVar->precision();
    filePtrSet.ofsVar->precision(16);
    for (unsigned int j = 0; j < m_cdfGrid.size(); ++j) {
      *filePtrSet.ofsVar << m_cdfGrid[j] << " ";
    }
    filePtrSet.ofsVar->precision(savedPrecision);

    *filePtrSet.ofsVar << "];\n";

    // Values
    *filePtrSet.ofsVar << varNamePrefix << "values_sub" << m_env.subIdString() << " = zeros(" << m_cdfValues.size()
                       << ","                                                                 << 1
                       << ");"
                       << std::endl;
    *filePtrSet.ofsVar << varNamePrefix << "values_sub" << m_env.subIdString() << " = [";

    savedPrecision = filePtrSet.ofsVar->precision();
    filePtrSet.ofsVar->precision(16);
    for (unsigned int j = 0; j < m_cdfValues.size(); ++j) {
      *filePtrSet.ofsVar << m_cdfValues[j] << " ";
    }
    filePtrSet.ofsVar->precision(savedPrecision);

    *filePtrSet.ofsVar << "];\n";

    // Close file
    m_env.closeFile(filePtrSet,fileType);
  }

  return;
}

}  // End namespace QUESO
