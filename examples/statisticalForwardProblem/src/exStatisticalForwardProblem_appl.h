/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of the QUESO Library (Quantification of Uncertainty
 * for Estimation, Simulation and Optimization).
 *
 * QUESO is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QUESO is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QUESO. If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * $Id$
 *
 * Brief description of this file: 
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __EX_STATISTICAL_FORWARD_PROBLEM_APPL_H__
#define __EX_STATISTICAL_FORWARD_PROBLEM_APPL_H__

#include <exStatisticalForwardProblem_qoi.h>
#include <uqStatisticalForwardProblem.h>
#include <uqCovCond.h>

//********************************************************
// The driving routine: called by main()
//********************************************************
template<class P_V,class P_M, class Q_V, class Q_M>
void 
uqAppl(const uqBaseEnvironmentClass& env)
{
  if (env.fullRank() == 0) {
    std::cout << "Beginning run of 'exStatisticalForwardProblem_example'\n"
              << std::endl;
  }

  //******************************************************
  // Step 1 of 6: Instantiate the parameter space (witt P_V and P_M)
  // It has dimension equal to 2
  //******************************************************
  if (env.fullRank() == 0) {
    std::cout << "Executing step 1 of 6: instantiation of parameter space ...\n"
              << std::endl;
  }

  uqVectorSpaceClass<P_V,P_M> paramSpace(env,"param_",2,NULL);

  //******************************************************
  // Step 2 of 6: Instantiate the parameter domain
  //******************************************************
  if (env.fullRank() == 0) {
    std::cout << "Executing step 2 of 6: instantiation of parameter domain ...\n"
              << std::endl;
  }

  P_V paramMins(paramSpace.zeroVector());
  paramMins[0] = 3.1;
  paramMins[1] = 0.7;
  P_V paramMaxs(paramSpace.zeroVector());
  paramMaxs[0] = 101.;
  paramMaxs[1] = 99.;
  uqBoxSubsetClass<P_V,P_M> paramDomain("param_",paramSpace,paramMins,paramMaxs);

  //******************************************************
  // Step 3 of 6: Instantiate the qoi space (with Q_V and Q_M)
  // It has dimension equal to 1
  //******************************************************
  if (env.fullRank() == 0) {
    std::cout << "Executing step 3 of 6: instantiation of qoi space ...\n"
              << std::endl;
  }

  uqVectorSpaceClass<Q_V,Q_M> qoiSpace(env,"qoi_",1,NULL);

  //******************************************************
  // Step 4 of 6: Instantiate the qoi function object (data + routine), to be used by QUESO.
  //******************************************************
  if (env.fullRank() == 0) {
    std::cout << "Executing step 4 of 6: instantiation of qoi function object ...\n"
              << std::endl;
  }

  qoiRoutine_DataType<P_V,P_M,Q_V,Q_M> qoiRoutine_Data;
  qoiRoutine_Data.p1MultiplicativeFactor = 1000.;
  qoiRoutine_Data.p1ExponentFactor       = 2.;
  qoiRoutine_Data.p2MultiplicativeFactor = 3.;
  qoiRoutine_Data.p2ExponentFactor       = 1.1;

  uqGenericVectorFunctionClass<P_V,P_M,Q_V,Q_M> qoiFunctionObj("like_",
                                                               paramDomain,
                                                               qoiSpace,
                                                               qoiRoutine<P_V,P_M,Q_V,Q_M>,
                                                               (void *) &qoiRoutine_Data);

  //******************************************************
  // Step 5 of 6: Instantiate the inverse problem
  //******************************************************
  if (env.fullRank() == 0) {
    std::cout << "Executing step 5 of 6: instantiation of forward problem ...\n"
              << std::endl;
  }

  uqUniformVectorRVClass<P_V,P_M> paramRv("param_", // Extra prefix before the default "rv_" prefix
                                          paramDomain);

  uqGenericVectorRVClass<Q_V,Q_M> qoiRv("qoi_", // Extra prefix before the default "rv_" prefix
                                         qoiSpace);

  uqStatisticalForwardProblemClass<P_V,P_M,Q_V,Q_M> fp("", // No extra prefix before the default "fp_" prefix
                                                       paramRv,
                                                       qoiFunctionObj,
                                                       qoiRv);

  //******************************************************
  // Step 6 of 6: Solve the inverse problem
  //******************************************************
  if (env.fullRank() == 0) {
    std::cout << "Executing step 6 of 6: solution of forward problem ...\n"
              << std::endl;
  }

  fp.solveWithMonteCarlo();

  //******************************************************
  // Release memory before leaving routine.
  //******************************************************

  if (env.fullRank() == 0) {
    std::cout << "Finishing run of 'exStatisticalForwardProblem_example'"
              << std::endl;
  }

  return;
}
#endif // __EX_STATISTICAL_FORWARD_PROBLEM_APPL_H__
