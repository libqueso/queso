//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
// 
// $Id$
//
//--------------------------------------------------------------------------

#include <example_qoi.h>

void
qoiRoutine(
  const uqGslVectorClass&                        paramValues,
  const uqGslVectorClass*                        paramDirection,
  const void*                                    functionDataPtr,
        uqGslVectorClass&                        qoiValues,
        EpetraExt::DistArray<uqGslVectorClass*>* gradVectors,
        EpetraExt::DistArray<uqGslMatrixClass*>* hessianMatrices,
        EpetraExt::DistArray<uqGslVectorClass*>* hessianEffects)
{
  // Logic just to avoid warnings from INTEL compiler
  const uqGslVectorClass* aux1 = paramDirection;
  if (aux1) {};
  EpetraExt::DistArray<uqGslVectorClass*>* aux2 = gradVectors;
  if (aux2) {};
  aux2 = hessianEffects;
  EpetraExt::DistArray<uqGslMatrixClass*>* aux3 = hessianMatrices;
  if (aux3) {};

  // Just checking: the user, at the application level, expects
  // vector 'paramValues' to have size 2 and
  // vector 'qoiValues' to have size 1.
  UQ_FATAL_TEST_MACRO(paramValues.sizeGlobal() != 2,
                      UQ_UNAVAILABLE_RANK,
                      "qoiRoutine()",
                      "paramValues vector does not have size 2");

  UQ_FATAL_TEST_MACRO(qoiValues.sizeGlobal() != 1,
                      UQ_UNAVAILABLE_RANK,
                      "qoiRoutine()",
                      "qoiValues vector does not have size 1");

  // Actual code
  //
  // This code exemplifies multiple Monte Carlo solvers, each calling this
  // qoi routine.
  //
  // In this simple example, only node 0 in each subenvironment does the job
  // even though there might be more than one node per subenvironment.
  //
  // In a more realistic situation, if the user is asking for multiple nodes per
  // subenvironment, then the model code in the qoi and likelihood routines
  // might really demand more than one node.
  //
  // Here we use 'env.subRank()' only. A realistic application might want to use
  // 'env.subComm()' or 'env.subComm().Comm()'
  const uqBaseEnvironmentClass& env = paramValues.env();
  if (env.subRank() == 0) {
    double coef1 = ((qoiRoutine_DataType *) functionDataPtr)->coef1;
    double coef2 = ((qoiRoutine_DataType *) functionDataPtr)->coef2;
    qoiValues[0] = (coef1*paramValues[0] + coef2*paramValues[1]);
  }
  else {
    qoiValues[0] = 0.;
  }

  return;
}
