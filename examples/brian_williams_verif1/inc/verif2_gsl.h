#ifndef __VERIF2_GSL_H__
#define __VERIF2_GSL_H__

#include <uqGslMatrix.h>

void   solveSip         (const uqFullEnvironmentClass& env);

double likelihoodRoutine(const uqGslVectorClass& paramValues,
                         const uqGslVectorClass* paramDirection,
                         const void*             functionDataPtr,
                         uqGslVectorClass*       gradVector,
                         uqGslMatrixClass*       hessianMatrix,
                         uqGslVectorClass*       hessianEffect);

struct likelihoodDataStruct {
  uqGslVectorClass* aVec;    // p x 1
  double            sigmaEps;
  uqGslVectorClass* ySamples; // n x 1
};

#endif // __VERIF2_GSL_H__
