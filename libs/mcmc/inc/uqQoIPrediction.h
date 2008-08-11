/* uq/libs/mcmc/inc/uqQoIPrediction.h
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

#ifndef __UQ_QOI_PREDICTION_H__
#define __UQ_QOI_PREDICTION_H__

#include <uqEnvironment.h>

//*****************************************************
// Class to accomodate the QoI predicition function.
// User must provide a prediction() routine.
//*****************************************************
template<class V, class M>
class uq_QoIPredictionFunction_Class {
public:
  uq_QoIPredictionFunction_Class(void (*functionPtr)(const V& paramValues, const void* functionDataPtr, std::vector<V*>& predictionValues),
                                 const void* functionDataPtr);
 ~uq_QoIPredictionFunction_Class();
 void operator()(const V& paramValues, std::vector<V*>& predictionValues);

private:
  void (*m_functionPtr)(const V& paramValues, const void* functionDataPtr, std::vector<V*>& predictionValues);
  const void* m_functionDataPtr;
};

template<class V, class M>
uq_QoIPredictionFunction_Class<V,M>::uq_QoIPredictionFunction_Class(
  void (*functionPtr)(const V& paramValues, const void* functionDataPtr, std::vector<V*>& predictionValues),
  const void* functionDataPtr)
  :
  m_functionPtr(functionPtr),
  m_functionDataPtr(functionDataPtr)
{
}

template<class V, class M>
uq_QoIPredictionFunction_Class<V,M>::~uq_QoIPredictionFunction_Class()
{
}

template<class V, class M>
void
uq_QoIPredictionFunction_Class<V,M>::operator()(const V& paramValues, std::vector<V*>& predictionValues)
{
  return m_functionPtr(paramValues, m_functionDataPtr, predictionValues);
}
#endif // __UQ_QOI_PREDICTION_H__
