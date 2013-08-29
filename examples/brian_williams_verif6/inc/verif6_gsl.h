#ifndef __VERIF6_GSL_H__
#define __VERIF6_GSL_H__

#include <uqGslMatrix.h>

void   solveSip         (const uqFullEnvironmentClass& env, bool useML);

double likelihoodRoutine(const uqGslVectorClass& paramValues,
                         const uqGslVectorClass* paramDirection,
                         const void*             functionDataPtr,
                         uqGslVectorClass*       gradVector,
                         uqGslMatrixClass*       hessianMatrix,
                         uqGslVectorClass*       hessianEffect);

struct likelihoodDataStruct {
  std::vector<double>* as;
  uqGslVectorClass*    bVec;   // p x 1
  std::vector<double>* sigmas;
  uqGslVectorClass*    tVec;   // n x 1
  uqGslVectorClass*    yVec;   // n x 1
};

#endif // __VERIF6_GSL_H__
