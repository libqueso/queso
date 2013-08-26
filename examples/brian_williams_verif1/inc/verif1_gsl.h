#ifndef __VERIF1_GSL_H__
#define __VERIF1_GSL_H__

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
  uqGslVectorClass* bVec;    // p x 1
  double            sigmaTotal;
  uqGslVectorClass* ySamples; // n x 1
};

#endif // __VERIF1_GSL_H__
