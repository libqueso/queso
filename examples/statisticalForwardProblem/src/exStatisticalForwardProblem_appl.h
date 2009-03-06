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
#include <uqAsciiTable.h>
#include <uqCovCond.h>

//********************************************************
// The driving routine: called by main()
//********************************************************
template<class P_V,class P_M, class Q_V, class Q_M>
void 
uqAppl(const uqBaseEnvironmentClass& env)
{
  if (env.rank() == 0) {
    std::cout << "Beginning run of 'exStatisticalForwardProblem_example'\n"
              << std::endl;
  }

  UQ_FATAL_TEST_MACRO(env.isThereInputFile() == false,
                      env.rank(),
                      "uqAppl()",
                      "input file must be specified in command line, after the '-i' option");

  //******************************************************
  // Read Ascii file with information on parameters
  //******************************************************
  uqAsciiTableClass<P_V,P_M> paramsTable(env,
                                         2,    // # of rows
                                         2,    // # of cols after 'parameter name': min + max
                                         NULL, // All extra columns are of 'double' type
                                         "inputData/params.tab");

  const EpetraExt::DistArray<std::string>& paramNames = paramsTable.stringColumn(0);
  P_V                                      paramMins    (paramsTable.doubleColumn(1));
  P_V                                      paramMaxs    (paramsTable.doubleColumn(2));

  uqVectorSpaceClass<P_V,P_M> paramSpace(env,
                                         "param_", // Extra prefix before the default "space_" prefix
                                         paramsTable.numRows(),
                                         &paramNames);

  uqBoxSubsetClass<P_V,P_M> paramDomain("param_",
                                        paramSpace,
                                        paramMins,
                                        paramMaxs);

  //******************************************************
  // Read Ascii file with information on qois
  //******************************************************
  uqAsciiTableClass<Q_V,Q_M> qoisTable(env,
                                         1,    // # of rows
                                         0,    // # of cols after 'qoi name': none
                                         NULL, // All extra columns are of 'double' type
                                         "inputData/qois.tab");

  const EpetraExt::DistArray<std::string>& qoiNames = qoisTable.stringColumn(0);

  uqVectorSpaceClass<Q_V,Q_M> qoiSpace(env,
                                       "qoi_", // Extra prefix before the default "space_" prefix
                                       qoisTable.numRows(),
                                       &qoiNames);

  //******************************************************
  // Instantiate a qoi function object (data + routine), to be used by QUESO.
  //******************************************************
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
  // Deal with forward problem
  //******************************************************
  uqUniformVectorRVClass<P_V,P_M> paramRv("param_", // Extra prefix before the default "rv_" prefix
                                          paramDomain);

  uqGenericVectorRVClass<Q_V,Q_M> qoiRv("qoi_", // Extra prefix before the default "rv_" prefix
                                         qoiSpace);

  uqStatisticalForwardProblemClass<P_V,P_M,Q_V,Q_M> fp("", // No extra prefix before the default "fp_" prefix
                                                       paramRv,
                                                       qoiFunctionObj,
                                                       qoiRv);

  //******************************************************
  // Solve forward problem = set 'pdf' and 'realizer' of 'qoiRv'
  //******************************************************
  fp.solveWithMonteCarlo();

  //******************************************************
  // Release memory before leaving routine.
  //******************************************************

  if (env.rank() == 0) {
    std::cout << "Finishing run of 'exStatisticalForwardProblem_example'"
              << std::endl;
  }

  return;
}
#endif // __EX_STATISTICAL_FORWARD_PROBLEM_APPL_H__
