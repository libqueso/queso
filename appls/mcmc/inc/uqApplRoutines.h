#ifndef __UQ_APPL_ROUTINES_H__
#define __UQ_APPL_ROUTINES_H__

#include <uqEnvironment.h>

//*****************************************************
// User must provide an application routine,
// to be called by main().
//*****************************************************
template<class V, class M>
void
uqAppl(const uqEnvironmentClass& env);

//*****************************************************
// User must provide a prior() routine *if* user
// does not want to use the default prior() above.
//*****************************************************
template<class V, class M>
double
uqAppl_M2lPriorFunction_Routine(const V& paramValues, const void* functionDataPtr);

template<class V, class M>
double
uq_M2lPriorFunction_Class<V,M>::operator()(const V& paramValues, const void* functionDataPtr, bool useDefault)
{
  double value = 0.;
  if (useDefault) {
    value = uqDefault_M2lPriorFunction_Routine<V,M>(paramValues, functionDataPtr);
  }
  else {
    value = uqAppl_M2lPriorFunction_Routine<V,M>(paramValues, functionDataPtr);
  }
  return value;
}

//*****************************************************
// User must provide a likelihood() routine.
//*****************************************************
template<class V, class M>
void
uqAppl_M2lLikelihoodFunction_Routine(const V& paramValues, const void* functionDataPtr, V& resultValues);

template<class V, class M>
void
uq_M2lLikelihoodFunction_Class<V,M>::operator()(const V& paramValues, const void* functionDataPtr, V& resultValues)
{
  return uqAppl_M2lLikelihoodFunction_Routine<V,M>(paramValues, functionDataPtr, resultValues);
}

//*****************************************************
// User must provide a stateEvolution() routine *if* user
// wants to compute predictions using uqPrediction<>().
//*****************************************************
template<class V, class M>
void
uqAppl_StateEvolutionFunction_Routine(const V& paramValues, const void* functionDataPtr, std::vector<V*>& states);

template<class V, class M>
void
uq_StateEvolutionFunction_Class<V,M>::operator()(const V& paramValues, const void* functionDataPtr, std::vector<V*>& states)
{
  return uqAppl_StateEvolutionFunction_Routine<V,M>(paramValues, functionDataPtr, states);
}

#endif // __UQ_APPL_ROUTINES_H__
