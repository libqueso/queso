/* uq/libs/mcmc/src/uqGslParamSpace.C
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

#include <uqParamSpace.h>
#include <uqGslMatrix.h>

template<>
void
uqParamSpaceClass<uqGslVectorClass,uqGslMatrixClass>::createInitialValues() const
{
  m_initialValues = new uqGslVectorClass(m_env,*m_map);

  for (unsigned int i = 0; i < m_parameters.size(); ++i) {
    if (m_parameters[i]) (*m_initialValues)[i] = m_parameters[i]->initialValue();
  }

  return;
}

template<>
void
uqParamSpaceClass<uqGslVectorClass,uqGslMatrixClass>::createMinValues() const
{
  m_minValues = new uqGslVectorClass(m_env,*m_map);

  for (unsigned int i = 0; i < m_parameters.size(); ++i) {
    if (m_parameters[i]) (*m_minValues)[i] = m_parameters[i]->minValue();
  }

  return; 
}

template<>
void
uqParamSpaceClass<uqGslVectorClass,uqGslMatrixClass>::createMaxValues() const
{
  m_maxValues = new uqGslVectorClass(m_env,*m_map);

  for (unsigned int i = 0; i < m_parameters.size(); ++i) {
    if (m_parameters[i]) (*m_maxValues)[i] = m_parameters[i]->maxValue();
  }

  return;
}

template<>
void
uqParamSpaceClass<uqGslVectorClass,uqGslMatrixClass>::createPriorMuValues() const
{
  m_priorMuValues = new uqGslVectorClass(m_env,*m_map);

  for (unsigned int i = 0; i < m_parameters.size(); ++i) {
    if (m_parameters[i]) (*m_priorMuValues)[i] = m_parameters[i]->priorMu();
  }

  return; 
}

template<>
void
uqParamSpaceClass<uqGslVectorClass,uqGslMatrixClass>::createPriorSigmaValues() const
{
  m_priorSigmaValues = new uqGslVectorClass(m_env,*m_map);

  for (unsigned int i = 0; i < m_parameters.size(); ++i) {
    if (m_parameters[i]) (*m_priorSigmaValues)[i] = m_parameters[i]->priorSigma();
  }

  return;
}

template<>
void
uqParamSpaceClass<uqGslVectorClass,uqGslMatrixClass>::createComponentsNames() const
{
  m_componentsNames.clear();
  m_componentsNames.resize(m_parameters.size());
  for (unsigned int i = 0; i < m_parameters.size(); ++i) {
    m_componentsNames[i] = m_parameters[i]->name();
  }

  return;
}

#if 0
template<>
void
uqParamSpaceClass<uqGslVectorClass,uqGslMatrixClass>::printParameterNames(std::ostream& os, bool printHorizontally) const
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
