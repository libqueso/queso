/* uq/libs/mcmc/inc/uqDefaultPrior.h
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

#ifndef __UQ_DEFAULT_PRIOR_H__
#define __UQ_DEFAULT_PRIOR_H__

#include <uqEnvironment.h>

//*****************************************************
// Default prior() routine provided by PECOS toolkit.
//*****************************************************
template<class V, class M>
struct
uqDefault_M2lPriorRoutine_DataType
{
  const V* paramPriorMus;
  const V* paramPriorVariances;
};

template<class V, class M>
double
uqDefault_M2lPriorRoutine(const V& paramValues, const void* functionDataPtr)
{
  if ((paramValues.env().verbosity() >= 10) && (paramValues.env().rank() == 0)) {
    std::cout << "Entering uqDefault_M2lPriorRoutine()" << std::endl;
  }

  const V& paramPriorMus       = *((uqDefault_M2lPriorRoutine_DataType<V,M> *) functionDataPtr)->paramPriorMus;
  const V& paramPriorVariances = *((uqDefault_M2lPriorRoutine_DataType<V,M> *) functionDataPtr)->paramPriorVariances;

  V diffVec(paramValues - paramPriorMus);
  double result = ((diffVec*diffVec)/paramPriorVariances).sumOfComponents();

  if ((paramValues.env().verbosity() >= 10) && (paramValues.env().rank() == 0)) {
    std::cout << "Leaving uqDefault_M2lPriorRoutine()" << std::endl;
  }

  return result;
}

#endif // __UQ_DEFAULT_PRIOR_H__
