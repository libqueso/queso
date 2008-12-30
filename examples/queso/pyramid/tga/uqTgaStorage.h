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
                    const std::vector<double>& temps,
                    const std::vector<double>& values,
                    const std::vector<double>* variances,
                          bool                 dataIsContinuousWithTime);
  uqTgaStorageClass(bool dataIsContinuousWithTime);
 ~uqTgaStorageClass();

        void                 resizeData    (unsigned int newSize);
        void                 setInstantData(unsigned int id,
                                            double       time,
                                            double       temp,
                                            double       value,
                                            double       variance);
        void                 set           (const std::vector<double>& times,
                                            const std::vector<double>& temps,
                                            const std::vector<double>& values,
                                            const std::vector<double>* variances);

        double               value    (double time) const;
  const std::vector<double>& times    () const;
  const std::vector<double>& temps    () const;
  const std::vector<double>& values   () const;
  const std::vector<double>& variances() const;
        bool                 dataIsContinuousWithTime() const;

  void  printForMatlab(std::ofstream& ofs, const std::string& prefixName) const;

protected:
        std::vector<double> m_times;
        std::vector<double> m_temps;
        std::vector<double> m_values;
        std::vector<double> m_variances;
        bool                m_dataIsContinuousWithTime;
};

template<class P_V, class P_M>
uqTgaStorageClass<P_V,P_M>::uqTgaStorageClass(
  const std::vector<double>& times,
  const std::vector<double>& temps,
  const std::vector<double>& values,
  const std::vector<double>* variances,
  bool dataIsContinuousWithTime)
  :
  m_times    (times.size(), 0.),
  m_temps    (temps.size(), 0.),
  m_values   (values.size(),0.),
  m_variances(values.size(),1.),
  m_dataIsContinuousWithTime(dataIsContinuousWithTime)
{
  unsigned int tmpSize = m_times.size();
  for (unsigned int i = 0; i < tmpSize; ++i) {
    m_times [i] = times [i];
    m_temps [i] = temps [i];
    m_values[i] = values[i];
    if (variances) m_variances[i] = (*variances)[i];
  }
}

template<class P_V, class P_M>
uqTgaStorageClass<P_V,P_M>::uqTgaStorageClass(
  bool dataIsContinuousWithTime)
  :
  m_times    (0),
  m_temps    (0),
  m_values   (0),
  m_variances(0),
  m_dataIsContinuousWithTime(dataIsContinuousWithTime)
{
}

template<class P_V, class P_M>
uqTgaStorageClass<P_V,P_M>::~uqTgaStorageClass()
{
}

template<class P_V, class P_M>
void
uqTgaStorageClass<P_V,P_M>::resizeData(unsigned int newSize)
{
  m_times.resize    (newSize,0.);
  m_temps.resize    (newSize,0.);
  m_values.resize   (newSize,0.);
  m_variances.resize(newSize,1.);

  return;
}

template<class P_V, class P_M>
void
uqTgaStorageClass<P_V,P_M>::setInstantData(
  unsigned int id,
  double       time,
  double       temp,
  double       value,
  double       variance)
{
  m_times    [id] = time;
  m_temps    [id] = temp;
  m_values   [id] = value;
  m_variances[id] = variance;

  return;
}

template<class P_V, class P_M>
void
uqTgaStorageClass<P_V,P_M>::set(
  const std::vector<double>& times,
  const std::vector<double>& temps,
  const std::vector<double>& values,
  const std::vector<double>* variances)
{
  unsigned int tmpSize = times.size();
  this->resizeData(tmpSize);
  for (unsigned int i = 0; i < tmpSize; ++i) {
    m_times [i] = times [i];
    m_temps [i] = temps [i];
    m_values[i] = values[i];
    if (variances) m_variances[i] = (*variances)[i];
  }

  return;
}

template<class P_V, class P_M>
double
uqTgaStorageClass<P_V,P_M>::value(double time) const
{
  double result = 0.;

  unsigned int tmpSize = m_times.size();
  //std::cout << "In uqTgaStorageClass<P_V,P_M>::value()"
  //          << ": time = "         << time
  //          << ", tmpSize = "      << tmpSize
  //          << ", m_times[0] = "   << m_times[0]
  //          << ", m_times[max] = " << m_times[tmpSize-1]
  //          << std::endl;

  UQ_FATAL_TEST_MACRO(tmpSize == 0,
                      UQ_UNAVAILABLE_RANK,
                      "uqTgaStorageClass<P_V,P_M>::value()",
                      "m_times.size() = 0");

  UQ_FATAL_TEST_MACRO(time < m_times[0],
                      UQ_UNAVAILABLE_RANK,
                      "uqTgaStorageClass<P_V,P_M>::value()",
                      "time < m_times[0]");

  UQ_FATAL_TEST_MACRO(m_times[tmpSize-1] < time,
                      UQ_UNAVAILABLE_RANK,
                      "uqTgaStorageClass<P_V,P_M>::value()",
                      "m_times[max] < time");

  unsigned int i = 0;
  for (i = 0; i < tmpSize; ++i) {
    if (time <= m_times[i]) break;
  }

  if (time == m_times[i]) {
    result = m_values[i];
  }
  else {
    //if ((9130.0 < time) && (time < 9131.0)) {
    //  std::cout << "time = " << time
    //            << "i = " << i
    //            << "time[i-1] = " << m_times[i-1]
    //            << "value[i-1] = " << m_values[i-1]
    //            << "time[i] = " << m_times[i]
    //            << "value[i] = " << m_values[i]
    //            << std::endl;
    //}
    double ratio = (time - m_times[i-1])/(m_times[i]-m_times[i-1]);
    result = m_values[i-1] + ratio * (m_values[i]-m_values[i-1]);
  }

  return result;
}

template<class P_V, class P_M>
const std::vector<double>&
uqTgaStorageClass<P_V,P_M>::times() const
{
  return m_times;
}

template<class P_V, class P_M>
const std::vector<double>&
uqTgaStorageClass<P_V,P_M>::temps() const
{
  return m_temps;
}

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
        << "\n" << prefixName << "Temp = zeros("  << tmpSize << ",1);"
        << "\n" << prefixName << "Value = zeros(" << tmpSize << ",1);";
  }
  else {
    ofs << "\n" << prefixName << "Time = zeros("  << tmpSize << ",1);"
        << "\n" << prefixName << "Temp = zeros("  << tmpSize << ",1);"
        << "\n" << prefixName << "Value = zeros(" << tmpSize << ",1);";
    for (unsigned int i = 0; i < tmpSize; ++i) {
      ofs << "\n" << prefixName << "Time("  << i+1 << ",1) = " << m_times[i]  << ";"
          << "\n" << prefixName << "Temp("  << i+1 << ",1) = " << m_temps[i]  << ";"
          << "\n" << prefixName << "Value(" << i+1 << ",1) = " << m_values[i] << ";";
    }
  }

  return;
}
#endif // __UQ_TGA_STORAGE _H__
