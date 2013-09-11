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

#ifndef EX_PHYSICS_1_LIKELIHOOD_H
#define EX_PHYSICS_1_LIKELIHOOD_H

#include <uqVectorSpace.h>
#include <uqDefines.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>

template<class P_V,class P_M>
struct
exPhysics1LikelihoodInfoStruct
{
  exPhysics1LikelihoodInfoStruct(const uqVectorSpace<P_V,P_M>& paramSpace,
                                 const std::string&                 inpName);
 ~exPhysics1LikelihoodInfoStruct();

  const uqVectorSpace <P_V,P_M>& m_paramSpace;
};

template<class P_V,class P_M>
exPhysics1LikelihoodInfoStruct<P_V,P_M>::exPhysics1LikelihoodInfoStruct(
  const uqVectorSpace<P_V,P_M>& paramSpace,
  const std::string&                 inpName)
  :
  m_paramSpace(paramSpace)
{
  if ((paramSpace.env().displayVerbosity() >= 30) && (paramSpace.env().fullRank() == 0)) {
    std::cout << "Entering exPhysics1LikelihoodInfoStruct<P_V,P_M>::constructor()"
              << inpName << "'\n"
              << std::endl;
  }

  // Read experimental data
  // Open input file on experimental data
  FILE *inp;
  inp = fopen(inpName.c_str(),"r");
  if (inp == NULL) {
    if (paramSpace.env().fullRank() == 0) {
      std::cout << "Input file " << inpName.c_str() << " was not found" << std::endl;
    }
    UQ_FATAL_TEST_MACRO((inp == NULL),
                        paramSpace.env().fullRank(),
                        "exPhysics1LikelihoodInfoStruct<P_V,P_M>::constructor()",
                        "input file not found");
  }

  // Read kinetic parameters and convert heating rate to K/s
  double tmpTemp;
  double tmpW;
  double tmpV;
  while (fscanf(inp,"%lf %lf %lf",&tmpTemp,&tmpW,&tmpV) != EOF) {
  }

  // Close input file on experimental data
  fclose(inp);

  if ((paramSpace.env().displayVerbosity() >= 30) && (paramSpace.env().fullRank() == 0)) {
    std::cout << "Leaving exPhysics1LikelihoodInfoStruct<P_V,P_M>::constructor()"
              << inpName << "'\n"
              << std::endl;
  }
}

template<class P_V,class P_M>
exPhysics1LikelihoodInfoStruct<P_V,P_M>::~exPhysics1LikelihoodInfoStruct()
{
}

// The actual (user defined) likelihood routine
template<class P_V,class P_M>
double
exPhysics1LikelihoodRoutine(
  const P_V&  paramValues,
  const P_V*  paramDirection,
  const void* functionDataPtr,
  P_V*        gradVector,
  P_M*        hessianMatrix,
  P_V*        hessianEffect)
{
  UQ_FATAL_TEST_MACRO((hessianEffect != NULL) && (paramDirection == NULL),
                      paramValues.env().fullRank(),
                      "exPhysics1LikelihoodRoutine<P_V,P_M>()",
                      "hessianEffect is being requested but no paramDirection is supplied");

  // Set all possible return values to zero
  double resultValue = 0.;
  if (gradVector   ) *gradVector    *= 0.;
  if (hessianMatrix) *hessianMatrix *= 0.;
  if (hessianEffect) *hessianEffect *= 0.;

  if ((paramValues.env().displayVerbosity() >= 30) && (paramValues.env().fullRank() == 0)) {
    std::cout << "Entering exPhysics1LikelihoodRoutine()..." << std::endl;
  }

  if ((paramValues.env().displayVerbosity() >= 10) && (paramValues.env().fullRank() == 0)) {
    std::cout << "In exPhysics1LikelihoodRoutine()"
              << ", A = " << paramValues[0]
              << ", E = " << paramValues[1]
              << std::endl;
  }

  // Decide what to compute, based on what is being requested
  bool computeLambda = false;
  bool computeWAndLambdaGradsAlso = false;

  if (gradVector) computeLambda = true;
  if (hessianMatrix) {
    computeLambda              = true;
    computeWAndLambdaGradsAlso = true;
  }
  if (hessianEffect) computeLambda = true;
  if (computeLambda) {};              // just to remove compiler warning
  if (computeWAndLambdaGradsAlso) {}; // just to remove compiler warning

  // Loop on scenarios
  const std::vector<exPhysics1LikelihoodInfoStruct<P_V,P_M> *>& vecInfo = *((const std::vector<exPhysics1LikelihoodInfoStruct<P_V,P_M> *> *) functionDataPtr);
  UQ_FATAL_TEST_MACRO(vecInfo.size() == 0,
                      paramValues.env().fullRank(),
                      "exPhysics1LikelihoodRoutine<P_V,P_M>()",
                      "vecInfo has size 0");
  for (unsigned int i = 0; i < vecInfo.size(); ++i) {
    resultValue += 0.;
  } // end loop of scenarios

  if ((paramValues.env().displayVerbosity() >= 30) && (paramValues.env().fullRank() == 0)) {
    std::cout << "Leaving exPhysics1LikelihoodRoutine()..." << std::endl;
  }

  return resultValue;
}

#endif // EX_PHYSICS_1_LIKELIHOOD_H
