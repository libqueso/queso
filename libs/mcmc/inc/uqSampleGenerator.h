/* uq/libs/mcmc/inc/uqSampleGenerator.h
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

#ifndef __UQ_SAMPLE_GENERATOR_H__
#define __UQ_SAMPLE_GENERATOR_H__

#include <uqEnvironment.h>
#include <math.h>

//*****************************************************
// Classes to accomodate a probability density routine.
//*****************************************************

//*****************************************************
// Base class
//*****************************************************
template<class V, class M>
class uqSampleGenerator_BaseClass {
public:
  uqSampleGenerator_BaseClass(double (*routinePtr)(const V& paramValues, const void* routineDataPtr),
                                const void* routineDataPtr);
  virtual ~uqSampleGenerator_BaseClass();
  virtual void nextSample(const V& paramValues) const;

protected:
  double (*m_routinePtr)(const V& paramValues, const void* routineDataPtr);
  const void* m_routineDataPtr;
};

template<class V, class M>
uqSampleGenerator_BaseClass<V,M>::uqSampleGenerator_BaseClass(
  double (*routinePtr)(const V& paramValues, const void* routineDataPtr),
  const void* routineDataPtr)
  :
  m_routinePtr    (routinePtr),
  m_routineDataPtr(routineDataPtr)
{
}

template<class V, class M>
uqSampleGenerator_BaseClass<V,M>::~uqSampleGenerator_BaseClass()
{
}

template<class V, class M>
void
uqSampleGenerator_BaseClass<V,M>::nextSample(const V& paramValues) const
{
  return m_routinePtr(paramValues, m_routineDataPtr);
}
#endif // __UQ_SAMPLE_GENERATOR_H__
