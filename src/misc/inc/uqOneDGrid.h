//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012 The PECOS Development Team
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

#ifndef __UQ_ONE_D_GRID_FUNCTION_H__
#define __UQ_ONE_D_GRID_FUNCTION_H__

#include <uqEnvironment.h>
#include <math.h>

//*****************************************************
// Classes to accomodate a one dimensional grid
//*****************************************************

//*****************************************************
// Base class
//*****************************************************
template<class T>
class uqBaseOneDGridClass {
public:
           uqBaseOneDGridClass(const uqBaseEnvironmentClass& env,
                               const char* prefix);
  virtual ~uqBaseOneDGridClass();

  virtual unsigned int size          ()                     const = 0;
  virtual T            operator[]    (unsigned int i)       const = 0;
  virtual unsigned int findIntervalId(const T& paramValue)  const = 0; 
          void         print         (std::ostream& ofsvar) const;

protected:
  const uqBaseEnvironmentClass& m_env;
        std::string             m_prefix;
};

template<class T>
uqBaseOneDGridClass<T>::uqBaseOneDGridClass(
  const uqBaseEnvironmentClass& env,
  const char*                   prefix)
  :
  m_env   (env),
  m_prefix((std::string)(prefix)+"grid")
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqBaseOneDGridClass<T>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqBaseOneDGridClass<T>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }
}

template<class T>
uqBaseOneDGridClass<T>::~uqBaseOneDGridClass()
{
}

template <class T>
void
uqBaseOneDGridClass<T>::print(std::ostream& os) const
{
  // Print values *of* grid points
  os << m_prefix << "_sub" << m_env.subIdString() << " = zeros(" << this->size()
     << ","                                                      << 1
     << ");"
     << std::endl;
  os << m_prefix << "_sub" << m_env.subIdString() << " = [";
  for (unsigned int j = 0; j < this->size(); ++j) {
    os << (*this)[j] << " ";
  }
  os << "];"
     << std::endl;

  return;
}

template <class T>
std::ostream& operator<< (std::ostream& os, const uqBaseOneDGridClass<T>& obj)
{
  obj.print(os);
  return os;
}

//*****************************************************
// Uniform grid class
//*****************************************************
template<class T>
class uqUniformOneDGridClass : public uqBaseOneDGridClass<T> {
public:
  uqUniformOneDGridClass(const uqBaseEnvironmentClass& env,
                         const char*               prefix,
                               unsigned int        size,
                               T                   minPosition,
                               T                   maxPosition);
 ~uqUniformOneDGridClass();

  unsigned int size          ()                    const;
  T            operator[]    (unsigned int i)      const;
  unsigned int findIntervalId(const T& paramValue) const; 

protected:
  using uqBaseOneDGridClass<T>::m_env;
  using uqBaseOneDGridClass<T>::m_prefix;

  unsigned int m_size;
  T            m_minPosition;
  T            m_maxPosition;
};

template<class T>
uqUniformOneDGridClass<T>::uqUniformOneDGridClass(
  const uqBaseEnvironmentClass& env,
  const char*               prefix,
        unsigned int        size,
        T                   minPosition,
        T                   maxPosition)
  :
  uqBaseOneDGridClass<T>(env,prefix),
  m_size       (size),
  m_minPosition(minPosition),
  m_maxPosition(maxPosition)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqUniformOneDGridClass<T>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqUniformOneDGridClass<T>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }
}

template<class T>
uqUniformOneDGridClass<T>::~uqUniformOneDGridClass()
{
}

template<class T>
unsigned int
uqUniformOneDGridClass<T>::size() const
{
  return m_size;
}

template<class T>
T
uqUniformOneDGridClass<T>::operator[](unsigned int i) const
{
  UQ_FATAL_TEST_MACRO(i >= m_size,
                      m_env.worldRank(),
                      "uqUniformOneDGridClass<V,M>::operator[]",
                      "too large i");

  T ratio = ((T) i)/(((T)m_size)-1.); // IMPORTANT: Yes, '-1.'
  T position = (1.-ratio)*m_minPosition + ratio*m_maxPosition;
  return position;
}

template<class T>
unsigned int
uqUniformOneDGridClass<T>::findIntervalId(const T& paramValue) const
{
  UQ_FATAL_TEST_MACRO((paramValue < m_minPosition) || (m_maxPosition < paramValue),
                      m_env.worldRank(),
                      "uqUniformOneDGridClass<V,M>::findIntervalId[]",
                      "paramValue is out of domain");

  T ratio = (paramValue - m_minPosition)/(m_maxPosition - m_minPosition);
  unsigned int i = (unsigned int) (ratio*(m_size-1.));
  if ((i > 0                  ) && 
      ((*this)[i] > paramValue)) {
    i--;
  }

  return i;
}

//*****************************************************
// Std grid class
//*****************************************************
template<class T>
class uqStdOneDGridClass : public uqBaseOneDGridClass<T> {
public:
  uqStdOneDGridClass(const uqBaseEnvironmentClass& env,
                     const char*                   prefix,
                     const std::vector<T>&         points);
 ~uqStdOneDGridClass();

  unsigned int size          ()                    const;
  T            operator[]    (unsigned int i)      const;
  unsigned int findIntervalId(const T& paramValue) const; 

protected:
  using uqBaseOneDGridClass<T>::m_env;
  using uqBaseOneDGridClass<T>::m_prefix;

  std::vector<T> m_points;
};

template<class T>
uqStdOneDGridClass<T>::uqStdOneDGridClass(
  const uqBaseEnvironmentClass& env,
  const char*                   prefix,
  const std::vector<T>&         points)
  :
  uqBaseOneDGridClass<T>(env,prefix),
  m_points              (points)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqStdOneDGridClass<T>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqStdOneDGridClass<T>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }
}

template<class T>
uqStdOneDGridClass<T>::~uqStdOneDGridClass()
{
}

template<class T>
unsigned int
uqStdOneDGridClass<T>::size() const
{
  return m_points.size();
}

template<class T>
T
uqStdOneDGridClass<T>::operator[](unsigned int i) const
{
  UQ_FATAL_TEST_MACRO(i >= m_points.size(),
                      m_env.worldRank(),
                      "uqStdOneDGridClass<V,M>::operator[]",
                      "too large i");

  return m_points[i];
}

template<class T>
unsigned int
uqStdOneDGridClass<T>::findIntervalId(const T& paramValue) const
{
  UQ_FATAL_TEST_MACRO((paramValue < m_points[0]) || (m_points[m_points.size()-1] < paramValue),
                      m_env.worldRank(),
                      "uqStdOneDGridClass<V,M>::findIntervalId[]",
                      "paramValue is out of domain");

  unsigned int iMax = m_points.size();
  unsigned int i = 1; // Yes, '1'
  for (i = 1; i < iMax; ++i) { // Yes, '1'
    if (paramValue < m_points[i]) {
      i--;
      break;
    }
  }

  return i;
}

#endif // __UQ_ONE_D_GRID_FUNCTION_H__
