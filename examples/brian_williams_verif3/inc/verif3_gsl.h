#ifndef __VERIF3_GSL_H__
#define __VERIF3_GSL_H__

#include <uqGslMatrix.h>

void   solveSips             (const uqFullEnvironmentClass& env);

double likelihoodRoutineForX0(const uqGslVectorClass& paramValues,
                              const uqGslVectorClass* paramDirection,
                              const void*             functionDataPtr,
                              uqGslVectorClass*       gradVector,
                              uqGslMatrixClass*       hessianMatrix,
                              uqGslVectorClass*       hessianEffect);

double likelihoodRoutineForX (const uqGslVectorClass& paramValues,
                              const uqGslVectorClass* paramDirection,
                              const void*             functionDataPtr,
                              uqGslVectorClass*       gradVector,
                              uqGslMatrixClass*       hessianMatrix,
                              uqGslVectorClass*       hessianEffect);

double routinePriorPdfForX   (const uqGslVectorClass& domainVector,
                              const uqGslVectorClass* domainDirection,
                              const void*             routineDataPtr,
                              uqGslVectorClass*       gradVector,
                              uqGslMatrixClass*       hessianMatrix,
                              uqGslVectorClass*       hessianEffect);

struct likelihoodDataStructForX0 {
  uqGslMatrixClass* zMat;     // p x n
  double            sigmaEps;
  uqGslVectorClass* ySamples; // n x 1
  uqGslMatrixClass* sigmaMatInverse; // p x p
};

struct likelihoodDataStructForX {
  uqGslMatrixClass* zMat;     // p x n
  double            sigmaEps;
  uqGslVectorClass* ySamples; // n x 1
  uqGslMatrixClass* sigmaMatInverse; // p x p
};

#endif // __VERIF3_GSL_H__
