//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, 
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-
// 
// $Id: example_qoi.C 37704 2013-03-08 21:10:36Z karl $
//
//--------------------------------------------------------------------------

#include <simple_sfp_example_qoi.h>

void
qoiRoutine(
  const uqGslVectorClass&                    paramValues,
  const uqGslVectorClass*                    paramDirection,
  const void*                                functionDataPtr,
        uqGslVectorClass&                    qoiValues,
        uqDistArrayClass<uqGslVectorClass*>* gradVectors,
        uqDistArrayClass<uqGslMatrixClass*>* hessianMatrices,
        uqDistArrayClass<uqGslVectorClass*>* hessianEffects)
{
  // Logic just to avoid warnings from INTEL compiler
  const uqGslVectorClass* aux1 = paramDirection;
  if (aux1) {};
  uqDistArrayClass<uqGslVectorClass*>* aux2 = gradVectors;
  if (aux2) {};
  aux2 = hessianEffects;
  uqDistArrayClass<uqGslMatrixClass*>* aux3 = hessianMatrices;
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
  // This code exemplifies multiple Monte Carlo solvers, each calling this qoi routine. 
  // In this simple example, only node 0 in each sub-environment does the job even though 
  // there might be more than one node per sub-environment.
  // In a more realistic situation, if the user is asking for multiple nodes per sub-
  // environment, then the model code in the qoi routine might really demand more than one
  // node. Here we use 'env.subRank()' only. A realistic application might want to use 
  // either 'env.subComm()' or 'env.subComm().Comm()'.
  
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
