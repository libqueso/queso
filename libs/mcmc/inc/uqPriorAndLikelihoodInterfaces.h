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
  if ((paramValues.env().verbosity() >= 10) && (paramValues.env().rank() == 0)) {
    std::cout << "Entering uqDefault_M2lPriorFunction_Routine()" << std::endl;
  }

  const V& paramPriorMus    = *((uqDefault_M2lPriorFunction_DataType<V,M> *) functionDataPtr)->paramPriorMus;
  const V& paramPriorSigmas = *((uqDefault_M2lPriorFunction_DataType<V,M> *) functionDataPtr)->paramPriorSigmas;

  double result = ((paramValues - paramPriorMus)/paramPriorSigmas).norm2Sq();

  if ((paramValues.env().verbosity() >= 10) && (paramValues.env().rank() == 0)) {
    std::cout << "Leaving uqDefault_M2lPriorFunction_Routine()" << std::endl;
  }

  return result;
}

//*****************************************************
// Class to accomodate the prior function.
// User must provide a prior() routine *if* user
// does not want to use the default prior() above.
//*****************************************************
template<class V, class M>
class tmp_M2lPriorFunction_Class {
public:
  tmp_M2lPriorFunction_Class(double (*funcPtr)(const V& paramValues, const void* functionDataPtr));
 ~tmp_M2lPriorFunction_Class();
  double operator()(const V& paramValues, const void* functionDataPtr);

private:
  double (*m_funcPtr)(const V& paramValues, const void* functionDataPtr);

};

template<class V, class M>
tmp_M2lPriorFunction_Class<V,M>::tmp_M2lPriorFunction_Class(double (*funcPtr)(const V& paramValues, const void* functionDataPtr))
  :
  m_funcPtr(funcPtr)
{
  //if (m_funcPtr == NULL) m_funcPtr(uqDefault_M2lPriorFunction_Routine<V,M>);
}

template<class V, class M>
tmp_M2lPriorFunction_Class<V,M>::~tmp_M2lPriorFunction_Class()
{
}

template<class V, class M>
double
tmp_M2lPriorFunction_Class<V,M>::operator()(const V& paramValues, const void* functionDataPtr)
{
  return m_funcPtr(paramValues, functionDataPtr);
}

//*****************************************************
// Class to accomodate the likehood function.
// User must provide a likelihood() routine.
//*****************************************************
template<class V, class M>
class tmp_M2lLikelihoodFunction_Class {
public:
  tmp_M2lLikelihoodFunction_Class(void (*funcPtr)(const V& paramValues, const void* functionDataPtr, V& resultValues));
 ~tmp_M2lLikelihoodFunction_Class();
  void operator()(const V& paramValues, const void* functionDataPtr, V& resultValues);

private:
  void (*m_funcPtr)(const V& paramValues, const void* functionDataPtr, V& resultValues);

};

template<class V, class M>
tmp_M2lLikelihoodFunction_Class<V,M>::tmp_M2lLikelihoodFunction_Class(void (*funcPtr)(const V& paramValues, const void* functionDataPtr, V& resultValues))
  :
  m_funcPtr(funcPtr)
{
}

template<class V, class M>
tmp_M2lLikelihoodFunction_Class<V,M>::~tmp_M2lLikelihoodFunction_Class()
{
}

template<class V, class M>
void
tmp_M2lLikelihoodFunction_Class<V,M>::operator()(const V& paramValues, const void* functionDataPtr, V& resultValues)
{
  return m_funcPtr(paramValues, functionDataPtr, resultValues);
}

//*****************************************************
// Class to accomodate the QoI predicition function.
// User must provide a prediction() routine.
//*****************************************************
template<class V, class M>
class tmp_QoIPredictionFunction_Class {
public:
  tmp_QoIPredictionFunction_Class(void (*funcPtr)(const V& paramValues, const void* functionDataPtr, std::vector<V*>& predictions));
 ~tmp_QoIPredictionFunction_Class();
  void operator()(const V& paramValues, const void* functionDataPtr, std::vector<V*>& predictions);

private:
  void (*m_funcPtr)(const V& paramValues, const void* functionDataPtr, std::vector<V*>& predictions);

};

template<class V, class M>
tmp_QoIPredictionFunction_Class<V,M>::tmp_QoIPredictionFunction_Class(void (*funcPtr)(const V& paramValues, const void* functionDataPtr, std::vector<V*>& predictions))
  :
  m_funcPtr(funcPtr)
{
}

template<class V, class M>
tmp_QoIPredictionFunction_Class<V,M>::~tmp_QoIPredictionFunction_Class()
{
}

template<class V, class M>
void
tmp_QoIPredictionFunction_Class<V,M>::operator()(const V& paramValues, const void* functionDataPtr, std::vector<V*>& predictions)
{
  return m_funcPtr(paramValues, functionDataPtr, predictions);
}
#endif // __UQ_PRIOR_AND_LIKELIHOOD_INTERFACES_H__
