/* uq/libs/mcmc/inc/uqRandomVariable.h
 *
 * Copyright (C) 2008 The PECOS Team, http://www.ices.utexas.edu/centers/pecos
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

#ifndef __UQ_RANDOM_VARIABLE_H__
#define __UQ_RANDOM_VARIABLE_H__

#include <uqProbDensity.h>
#include <uqRealizer.h>

//*****************************************************
// Base class
//*****************************************************
template<class V, class M>
class uqRandomVariableClass {
public:
  uqRandomVariableClass(const uqProbDensity_BaseClass<V,M>* probDensity,
                        const uqRealizer_BaseClass   <V,M>* realizer);
  virtual ~uqRandomVariableClass();

  const uqProbDensity_BaseClass<V,M>& probDensity()       const;
  const uqRealizer_BaseClass   <V,M>& realizer   ()       const;
        void                          realization(V& vec) const;

protected:
  const uqProbDensity_BaseClass<V,M>* m_probDensity;
  const uqRealizer_BaseClass   <V,M>* m_realizer;
};

template<class V, class M>
uqRandomVariableClass<V,M>::uqRandomVariableClass(
  const uqProbDensity_BaseClass<V,M>* probDensity,
  const uqRealizer_BaseClass   <V,M>* realizer)
  :
  m_probDensity(probDensity),
  m_realizer   (realizer)
{
}

template<class V, class M>
uqRandomVariableClass<V,M>::~uqRandomVariableClass()
{
}

template<class V, class M>
const uqProbDensity_BaseClass<V,M>&
uqRandomVariableClass<V,M>::probDensity() const
{
  UQ_FATAL_TEST_MACRO(m_probDensity == NULL,
                      UQ_UNAVAILABLE_RANK,
                      "uqRandomVariableClass<V,M>::probDensity()",
                      "m_probDensity is NULL");

  return *m_probDensity;
}

template<class V, class M>
const uqRealizer_BaseClass<V,M>&
uqRandomVariableClass<V,M>::realizer() const
{
  UQ_FATAL_TEST_MACRO(m_realizer == NULL,
                      UQ_UNAVAILABLE_RANK,
                      "uqRandomVariableClass<V,M>::realizer()",
                      "m_realizer is NULL");

  return *m_realizer;
}

template<class V, class M>
void
uqRandomVariableClass<V,M>::realization(V& vec) const
{
  UQ_FATAL_TEST_MACRO(m_realizer == NULL,
                      UQ_UNAVAILABLE_RANK,
                      "uqRandomVariableClass<V,M>::realization()",
                      "m_realizer is NULL");

  m_realizer->nextSample(vec);

  return;
}

#endif // __UQ_RANDOM_VARIABLE_H__
