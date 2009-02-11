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

#ifndef __UQ_MONTE_CARLO_EX_H__
#define __UQ_MONTE_CARLO_EX_H__

#include <uqStatisticalForwardProblem.h>
#include <uqAsciiTable.h>
#include <uqCovCond.h>

//********************************************************
// The qoi routine: provided by user and called by QUESO
//********************************************************
template<class P_V,class P_M, class Q_V, class Q_M>
struct
qoiRoutine_DataType
{
  double p1MultiplicativeFactor;
  double p1ExponentFactor;
  double p2MultiplicativeFactor;
  double p2ExponentFactor;
};

template<class P_V,class P_M, class Q_V, class Q_M>
void
qoiRoutine(
  const P_V&                        paramValues,
  const P_V*                        paramDirection,
  const void*                       functionDataPtr,
        Q_V&                        qoiValues,
        EpetraExt::DistArray<P_V*>* gradVectors,
        EpetraExt::DistArray<P_M*>* hessianMatrices,
        EpetraExt::DistArray<P_V*>* hessianEffects)
{
  double a1 = ((qoiRoutine_DataType<P_V,P_M,Q_V,Q_M> *) functionDataPtr)->p1MultiplicativeFactor;
  double e1 = ((qoiRoutine_DataType<P_V,P_M,Q_V,Q_M> *) functionDataPtr)->p1ExponentFactor;
  double a2 = ((qoiRoutine_DataType<P_V,P_M,Q_V,Q_M> *) functionDataPtr)->p2MultiplicativeFactor;
  double e2 = ((qoiRoutine_DataType<P_V,P_M,Q_V,Q_M> *) functionDataPtr)->p2ExponentFactor;

  double p1 = paramValues[0];
  double p2 = paramValues[1];

  qoiValues[0] = a1*pow(p1,e1) + a2*pow(p2,e2);

  return;
}

//********************************************************
// The driving routine: called by main()
//********************************************************
template<class P_V,class P_M, class Q_V, class Q_M>
void 
uqAppl(const uqBaseEnvironmentClass& env)
{
  if (env.rank() == 0) {
    std::cout << "Beginning run of 'uqMonteCarloExample'\n"
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

  uqStatisticalForwardProblemClass<P_V,P_M,Q_V,Q_M> fp("", // No extra prefix before the default "ip_" prefix
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
    std::cout << "Finishing run of 'uqMonteCarloExample'"
              << std::endl;
  }

  return;
}
#endif // __UQ_MONTE_CARLO_EX_H__
