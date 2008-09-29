/* uq/libs/mcmc/src/uqGslVectorRV.C
 *
 * Copyright (C) 2008 The PECOS Team, http://queso.ices.utexas.edu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <uqVectorRV.h>
#include <uqGslMatrix.h>

template<>
void
uqBaseVectorRVClass<uqGslVectorClass,uqGslMatrixClass>::createMinValues() const
{
  m_minValues = new uqGslVectorClass(m_env,m_imageSpace.map());

  for (unsigned int i = 0; i < m_components.size(); ++i) {
    if (m_components[i]) (*m_minValues)[i] = m_components[i]->minValue();
  }

  return; 
}

template<>
void
uqBaseVectorRVClass<uqGslVectorClass,uqGslMatrixClass>::createMaxValues() const
{
  m_maxValues = new uqGslVectorClass(m_env,m_imageSpace.map());

  for (unsigned int i = 0; i < m_components.size(); ++i) {
    if (m_components[i]) (*m_maxValues)[i] = m_components[i]->maxValue();
  }

  return;
}

template<>
void
uqBaseVectorRVClass<uqGslVectorClass,uqGslMatrixClass>::createExpectValues() const
{
  m_expectValues = new uqGslVectorClass(m_env,m_imageSpace.map());

  for (unsigned int i = 0; i < m_components.size(); ++i) {
    if (m_components[i]) (*m_expectValues)[i] = m_components[i]->expectValue();
  }

  return; 
}

template<>
void
uqBaseVectorRVClass<uqGslVectorClass,uqGslMatrixClass>::createStdDevValues() const
{
  m_stdDevValues = new uqGslVectorClass(m_env,m_imageSpace.map());

  for (unsigned int i = 0; i < m_components.size(); ++i) {
    if (m_components[i]) (*m_stdDevValues)[i] = m_components[i]->stdDevValue();
  }

  return;
}

#if 0
template<>
void
uqBaseVectorRVClass<uqGslVectorClass,uqGslMatrixClass>::printParameterNames(std::ostream& os, bool printHorizontally) const
{
  if (printHorizontally) { 
    for (unsigned int i = 0; i < this->dim(); ++i) {
      os << "'" << this->parameter(i).name() << "'"
         << " ";
    }
  }
  else {
    for (unsigned int i = 0; i < this->dim(); ++i) {
      os << "'" << this->parameter(i).name() << "'"
         << std::endl;
    }
  }

  return;
}
#endif
