//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
// 
// $Id$
//
//--------------------------------------------------------------------------
#ifndef __EX_QOI_H__
#define __EX_QOI_H__

#include <uqGslMatrix.h>
#include <EpetraExt_DistArray.h>

struct
qoiRoutine_DataType
{
  double coef1;
  double coef2;
};

void
qoiRoutine(
  const uqGslVectorClass&                        paramValues,
  const uqGslVectorClass*                        paramDirection,
  const void*                                    functionDataPtr,
        uqGslVectorClass&                        qoiValues,
        EpetraExt::DistArray<uqGslVectorClass*>* gradVectors,
        EpetraExt::DistArray<uqGslMatrixClass*>* hessianMatrices,
        EpetraExt::DistArray<uqGslVectorClass*>* hessianEffects);

#endif
