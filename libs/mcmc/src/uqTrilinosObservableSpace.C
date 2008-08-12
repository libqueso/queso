/* uq/libs/mcmc/src/uqTrilinosObservableSpace.C
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

#include <uqObservableSpace.h>
#include <uqTrilinosMatrix.h>

template<>
void
uqObservableSpaceClass<uqTrilinosVectorClass,uqTrilinosMatrixClass>::createNumbersOfObservations() const
{
  m_numbersOfObservations = this->newVector();//new uqTrilinosVectorClass(m_env,this->dim());

  for (unsigned int i = 0; i < m_observables.size(); ++i) {
    if (m_observables[i]) (*m_numbersOfObservations)[i] = (double) m_observables[i]->numberOfObservations();
  }

  return;
}

template<>
void
uqObservableSpaceClass<uqTrilinosVectorClass,uqTrilinosMatrixClass>::createPriorVariances() const
{
  m_priorVariances = this->newVector();//new uqTrilinosVectorClass(m_env,this->dim());

  for (unsigned int i = 0; i < m_observables.size(); ++i) {
    if (m_observables[i]) (*m_priorVariances)[i] = m_observables[i]->priorVariance();
  }

  return; 
}

template<>
void
uqObservableSpaceClass<uqTrilinosVectorClass,uqTrilinosMatrixClass>::createVarianceAccuracies() const
{
  m_varianceAccuracies = this->newVector();//new uqTrilinosVectorClass(m_env,this->dim());

  for (unsigned int i = 0; i < m_observables.size(); ++i) {
    if (m_observables[i]) (*m_varianceAccuracies)[i] = m_observables[i]->varianceAccuracy();
  }

  return;
}
