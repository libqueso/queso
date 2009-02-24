/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of the QUESO Library (Quantification of Uncertainty
 * for Estimation, Simulation and Optimization).
 *
 * QUESO is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QUESO is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QUESO. If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * $Id$
 *
 * Brief description of this file: 
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

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
        std::string         m_prefix;
};

template<class T>
uqBaseOneDGridClass<T>::uqBaseOneDGridClass(
  const uqBaseEnvironmentClass& env,
  const char*               prefix)
  :
  m_env   (env),
  m_prefix((std::string)(prefix)+"grid")
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqBaseOneDGridClass<T>::constructor()"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqBaseOneDGridClass<T>::constructor()"
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
  os << m_prefix << " = zeros(" << this->size()
     << ","                     << 1
     << ");"
     << std::endl;
  os << m_prefix << " = [";
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
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqUniformOneDGridClass<T>::constructor()"
              << ": prefix = " << m_prefix
              << std::endl;
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqUniformOneDGridClass<T>::constructor()"
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
                      m_env.rank(),
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
                      m_env.rank(),
                      "uqUniformOneDGridClass<V,M>::findIntervalId[]",
                      "paramValue is out of domain");

  T ratio = (paramValue - m_minPosition)/(m_maxPosition - m_minPosition);
  return (unsigned int) (ratio*(m_size-1.));
}

#endif // __UQ_ONE_D_GRID_FUNCTION_H__
