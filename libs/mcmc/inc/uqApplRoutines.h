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
// Default prior() routine provided by ICES UQ library.
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
double
uqAppl_M2lPriorFunction_Routine(const V& paramValues, const void* functionDataPtr);

template<class V, class M>
struct uq_M2lPriorFunction_Class {
  double operator()(const V& paramValues, const void* functionDataPtr, bool useDefault)
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
};

//*****************************************************
// User must provide a likelihood() routine.
//*****************************************************
template<class V, class M>
void
uqAppl_M2lLikelihoodFunction_Routine(const V& paramValues, const void* functionDataPtr, V& resultValues);

template<class V, class M>
struct uq_M2lLikelihoodFunction_Class {
  void operator()(const V& paramValues, const void* functionDataPtr, V& resultValues)
  {
    return uqAppl_M2lLikelihoodFunction_Routine<V,M>(paramValues, functionDataPtr, resultValues);
  }
};
#endif // __UQ_APPL_ROUTINES_H__
