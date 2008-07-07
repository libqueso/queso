/* uq/libs/mcmc/inc/uqPriorAndLikelihoodInterfaces.h
 *
 * Copyright (C) 2008 The PECOS Team, http://www.ices.utexas.edu/centers/pecos

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

#ifndef __UQ_PRIOR_AND_LIKELIHOOD_INTERFACES_H__
#define __UQ_PRIOR_AND_LIKELIHOOD_INTERFACES_H__

#include <uqEnvironment.h>

//*****************************************************
// Default prior() routine provided by PECOS toolkit.
//*****************************************************
template<class V, class M>
struct
uqDefault_M2lPriorFunction_DataType
{
  const V* paramPriorMus;
  const V* paramPriorSigmas;
};

template<class V, class M>
double
uqDefault_M2lPriorFunction_Routine(const V& paramValues, const void* functionDataPtr)
{
  const V& paramPriorMus    = *((uqDefault_M2lPriorFunction_DataType<V,M> *) functionDataPtr)->paramPriorMus;
  const V& paramPriorSigmas = *((uqDefault_M2lPriorFunction_DataType<V,M> *) functionDataPtr)->paramPriorSigmas;

  return ((paramValues - paramPriorMus)/paramPriorSigmas).norm2Sq();
}

//*****************************************************
// User must provide a prior() routine *if* user
// does not want to use the default prior() above.
//*****************************************************
template<class V, class M>
struct uq_M2lPriorFunction_Class {
  double operator()(const V& paramValues, const void* functionDataPtr, bool useDefault);
};

//*****************************************************
// User must provide a likelihood() routine.
//*****************************************************
template<class V, class M>
struct uq_M2lLikelihoodFunction_Class {
  void operator()(const V& paramValues, const void* functionDataPtr, V& resultValues);
};

//*****************************************************
// User must provide a stateEvolution() routine *if* user
// wants to compute predictions using uqPrediction<>().
//*****************************************************
template<class V, class M>
struct uq_StateEvolutionFunction_Class {
  void operator()(const V& paramValues, const void* functionDataPtr, std::vector<V*>& states);
};
#endif // __UQ_PRIOR_AND_LIKELIHOOD_INTERFACES_H__
