/* uq/libs/queso/inc/uqTKRealizer.h
 *
 * Copyright (C) 2008 The QUESO Team, http://queso.ices.utexas.edu/
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

#ifndef __UQ_TK_REALIZER_H__
#define __UQ_TK_REALIZER_H__

#include <uqEnvironment.h>
#include <math.h>

//*****************************************************
// Base class
//*****************************************************
template<class V, class M>
class uqBaseTKRealizerClass {
public:
  uqBaseTKRealizerClass(double (*routinePtr)(const V& paramValues, const void* routineDataPtr),
                                const void* routineDataPtr);
  virtual ~uqBaseTKRealizerClass();
  virtual void nextSample(const V& paramValues) const;

protected:
  double (*m_routinePtr)(const V& paramValues, const void* routineDataPtr);
  const void* m_routineDataPtr;
};

template<class V, class M>
uqBaseTKRealizerClass<V,M>::uqBaseTKRealizerClass(
  double (*routinePtr)(const V& paramValues, const void* routineDataPtr),
  const void* routineDataPtr)
  :
  m_routinePtr    (routinePtr),
  m_routineDataPtr(routineDataPtr)
{
}

template<class V, class M>
uqBaseTKRealizerClass<V,M>::~uqBaseTKRealizerClass()
{
}

template<class V, class M>
void
uqBaseTKRealizerClass<V,M>::nextSample(const V& paramValues) const
{
  return m_routinePtr(paramValues, m_routineDataPtr);
}
#endif // __UQ_TK_REALIZER_H__
