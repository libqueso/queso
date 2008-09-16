/* uq/libs/mcmc/src/uqGslQoISpace.C
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

#include <uqQoISpace.h>
#include <uqGslMatrix.h>

template<>
void
uqQoISpaceClass<uqGslVectorClass,uqGslMatrixClass>::createComponentsNames() const
{
  m_componentsNames.clear();
  m_componentsNames.resize(m_qois.size());
  for (unsigned int i = 0; i < m_qois.size(); ++i) {
    m_componentsNames[i] = m_qois[i]->name();
  }

  return;
}
