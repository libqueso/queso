//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
// 
// $Id$
//
//--------------------------------------------------------------------------
#ifndef __EX_LIKELIHOOD_H__
#define __EX_LIKELIHOOD_H__

#include <uqGslMatrix.h>

struct
likelihoodRoutine_DataType
{
  const uqGslVectorClass* meanVector;
  const uqGslMatrixClass* covMatrix;
};

double likelihoodRoutine(
  const uqGslVectorClass& paramValues,
  const uqGslVectorClass* paramDirection,
  const void*             functionDataPtr,
  uqGslVectorClass*       gradVector,
  uqGslMatrixClass*       hessianMatrix,
  uqGslVectorClass*       hessianEffect);

#endif
