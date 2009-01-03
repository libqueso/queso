/* uq/examples/queso/pyramid/uqTgaStorage.h
 *
 * Copyright (C) 2008 The QUESO Team, http://queso.ices.utexas.edu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __UQ_TGA_STORAGE_H__
#define __UQ_TGA_STORAGE_H__

#include <uq1D1DFunction.h>
#include <uqTgaDefines.h>
#include <uqEnvironment.h>
#include <uqDefines.h>

template<class P_V, class P_M>
class
uqTgaStorageClass
{
public:
  uqTgaStorageClass(const std::vector<double>& times,
#ifdef QUESO_TGA_USES_OLD_COMPATIBLE_CODE
                    const std::vector<double>& temps,
#endif
                    const std::vector<double>& values,
                    const std::vector<double>* variances,
                          bool                 dataIsContinuousWithTime);
  uqTgaStorageClass();
 ~uqTgaStorageClass();

#ifdef QUESO_TGA_USES_OLD_COMPATIBLE_CODE
#else
        void                 reset         (unsigned int newSize, bool dataIsContinuousWithTime);
        void                 setInstantData(unsigned int id,
                                            double       time,
                                            double       value,
                                            double       variance);
        void                 set           (const std::vector<double>& times,
                                            const std::vector<double>& values,
                                            const std::vector<double>* variances,
                                            bool                       dataIsContinuousWithTime);
#endif

        void                 getValue (double        time,
                                       unsigned int& startingTimeId,
                                       double&       returnValue,
                                       bool*         timeWasMatchedExactly) const;
  const std::vector<double>& times    () const;
#ifdef QUESO_TGA_USES_OLD_COMPATIBLE_CODE
  const std::vector<double>& temps    () const;
#endif
  const std::vector<double>& values   () const;
  const std::vector<double>& variances() const;
        bool                 dataIsContinuousWithTime() const;

  void  printForMatlab(std::ofstream& ofs, const std::string& prefixName) const;

protected:
        std::vector<double> m_times;
#ifdef QUESO_TGA_USES_OLD_COMPATIBLE_CODE
        std::vector<double> m_temps;
#endif
        std::vector<double> m_values;
        std::vector<double> m_variances;
        bool                m_dataIsContinuousWithTime;
};

template<class P_V, class P_M>
uqTgaStorageClass<P_V,P_M>::uqTgaStorageClass(
  const std::vector<double>& times,
#ifdef QUESO_TGA_USES_OLD_COMPATIBLE_CODE
  const std::vector<double>& temps,
#endif
  const std::vector<double>& values,
  const std::vector<double>* variances,
  bool dataIsContinuousWithTime)
  :
  m_times    (times.size(), 0.),
#ifdef QUESO_TGA_USES_OLD_COMPATIBLE_CODE
  m_temps    (temps.size(),0.),
#endif
  m_values   (values.size(),0.),
  m_variances(values.size(),1.),
  m_dataIsContinuousWithTime(dataIsContinuousWithTime)
{
  unsigned int tmpSize = m_times.size();
  for (unsigned int i = 0; i < tmpSize; ++i) {
    m_times [i] = times [i];
#ifdef QUESO_TGA_USES_OLD_COMPATIBLE_CODE
    m_temps [i] = temps [i];
#endif
    m_values[i] = values[i];
    if (variances) m_variances[i] = (*variances)[i];
  }
}

template<class P_V, class P_M>
uqTgaStorageClass<P_V,P_M>::uqTgaStorageClass()
  :
  m_times    (0),
#ifdef QUESO_TGA_USES_OLD_COMPATIBLE_CODE
  m_temps    (0),
#endif
  m_values   (0),
  m_variances(0),
  m_dataIsContinuousWithTime(false)
{
}

template<class P_V, class P_M>
uqTgaStorageClass<P_V,P_M>::~uqTgaStorageClass()
{
}

#ifdef QUESO_TGA_USES_OLD_COMPATIBLE_CODE
#else
template<class P_V, class P_M>
void
uqTgaStorageClass<P_V,P_M>::reset(
  unsigned int newSize,
  bool         dataIsContinuousWithTime)
{
  m_times.resize    (newSize,0.);
#ifdef QUESO_TGA_USES_OLD_COMPATIBLE_CODE
  m_temps.resize    (newSize,0.);
#endif
  m_values.resize   (newSize,0.);
  m_variances.resize(newSize,1.);
  m_dataIsContinuousWithTime = dataIsContinuousWithTime;

  return;
}

template<class P_V, class P_M>
void
uqTgaStorageClass<P_V,P_M>::setInstantData(
  unsigned int id,
  double       time,
  double       value,
  double       variance)
{
  UQ_FATAL_TEST_MACRO(id >= m_times.size(),
                      UQ_UNAVAILABLE_RANK,
                      "uqTgaStorageClass<P_V,P_M>::setInstantData()",
                      "id is too big");

  m_times    [id] = time;
  m_values   [id] = value;
  m_variances[id] = variance;

  return;
}

template<class P_V, class P_M>
void
uqTgaStorageClass<P_V,P_M>::set(
  const std::vector<double>& times,
  const std::vector<double>& values,
  const std::vector<double>* variances,
  bool                       dataIsContinuousWithTime)
{
  unsigned int tmpSize = times.size();
  this->reset(tmpSize,dataIsContinuousWithTime);
  for (unsigned int i = 0; i < tmpSize; ++i) {
    m_times [i] = times [i];
    m_values[i] = values[i];
    if (variances) m_variances[i] = (*variances)[i];
  }

  return;
}
#endif

template<class P_V, class P_M>
void   
uqTgaStorageClass<P_V,P_M>::getValue(
  double        time,
  unsigned int& startingTimeId,
  double&       returnValue,
  bool*         timeWasMatchedExactly) const
{
  returnValue = 0.;

  unsigned int tmpSize = m_times.size();
  //std::cout << "In uqTgaStorageClass<P_V,P_M>::getValue()"
  //          << ": time = "         << time
  //          << ", tmpSize = "      << tmpSize
  //          << ", m_times[0] = "   << m_times[0]
  //          << ", m_times[max] = " << m_times[tmpSize-1]
  //          << std::endl;

  UQ_FATAL_TEST_MACRO(tmpSize == 0,
                      UQ_UNAVAILABLE_RANK,
                      "uqTgaStorageClass<P_V,P_M>::getValue()",
                      "m_times.size() = 0");

  UQ_FATAL_TEST_MACRO(time < m_times[0],
                      UQ_UNAVAILABLE_RANK,
                      "uqTgaStorageClass<P_V,P_M>::getValue()",
                      "time < m_times[0]");

  UQ_FATAL_TEST_MACRO(m_times[tmpSize-1] < time,
                      UQ_UNAVAILABLE_RANK,
                      "uqTgaStorageClass<P_V,P_M>::getValue()",
                      "m_times[max] < time");

  unsigned int i = 0;
  for (i = 0; i < tmpSize; ++i) {
    if (time <= m_times[i]) break;
  }

  if (time == m_times[i]) {
    if (timeWasMatchedExactly) *timeWasMatchedExactly = true;
    returnValue = m_values[i];
  }
  else {
    if (timeWasMatchedExactly) *timeWasMatchedExactly = false;
    //if ((9130.0 < time) && (time < 9131.0)) {
    //  std::cout << "time = " << time
    //            << "i = " << i
    //            << "time[i-1] = " << m_times[i-1]
    //            << "value[i-1] = " << m_values[i-1]
    //            << "time[i] = " << m_times[i]
    //            << "value[i] = " << m_values[i]
    //            << std::endl;
    //}
    if (m_dataIsContinuousWithTime) {
      double ratio = (time - m_times[i-1])/(m_times[i]-m_times[i-1]);
      returnValue = m_values[i-1] + ratio * (m_values[i]-m_values[i-1]);
    }
    else {
      // Leave returnValue = 0.;
    }
  }

  return;
}

template<class P_V, class P_M>
const std::vector<double>&
uqTgaStorageClass<P_V,P_M>::times() const
{
  return m_times;
}

#ifdef QUESO_TGA_USES_OLD_COMPATIBLE_CODE
template<class P_V, class P_M>
const std::vector<double>&
uqTgaStorageClass<P_V,P_M>::temps() const
{
  return m_temps;
}

#endif

template<class P_V, class P_M>
const std::vector<double>&
uqTgaStorageClass<P_V,P_M>::values() const
{
  return m_values;
}

template<class P_V, class P_M>
const std::vector<double>&
uqTgaStorageClass<P_V,P_M>::variances() const
{
  return m_variances;
}

template<class P_V, class P_M>
bool
uqTgaStorageClass<P_V,P_M>::dataIsContinuousWithTime() const
{
  return m_dataIsContinuousWithTime;
}

template<class P_V, class P_M>
void
uqTgaStorageClass<P_V,P_M>::printForMatlab(
  std::ofstream&     ofs,
  const std::string& prefixName) const
{
  unsigned int tmpSize = m_times.size();
  if (tmpSize == 0) {
    tmpSize = 1;
    ofs << "\n" << prefixName << "Time = zeros("  << tmpSize << ",1);"
        << "\n" << prefixName << "Value = zeros(" << tmpSize << ",1);";
  }
  else {
    ofs << "\n" << prefixName << "Time = zeros("  << tmpSize << ",1);"
        << "\n" << prefixName << "Value = zeros(" << tmpSize << ",1);";
    for (unsigned int i = 0; i < tmpSize; ++i) {
      ofs << "\n" << prefixName << "Time("  << i+1 << ",1) = " << m_times[i]  << ";"
          << "\n" << prefixName << "Value(" << i+1 << ",1) = " << m_values[i] << ";";
    }
  }

  return;
}
#endif // __UQ_TGA_STORAGE _H__
