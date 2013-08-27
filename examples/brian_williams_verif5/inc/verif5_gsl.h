#ifndef __VERIF5_GSL_H__
#define __VERIF5_GSL_H__

#include <uqGslMatrix.h>

void   solveSip         (const uqFullEnvironmentClass& env, bool useML);

double likelihoodRoutine(const uqGslVectorClass& paramValues,
                         const uqGslVectorClass* paramDirection,
                         const void*             functionDataPtr,
                         uqGslVectorClass*       gradVector,
                         uqGslMatrixClass*       hessianMatrix,
                         uqGslVectorClass*       hessianEffect);

struct likelihoodDataStruct {
  uqGslVectorClass* bVec;     // p x 1
  double            sigmaTotal;
  uqGslVectorClass* ySamples; // n x 1
};

#endif // __VERIF5_GSL_H__
