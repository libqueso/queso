/* uq/examples/queso/pyramid/uqTgaLikelihood.h
 *
 * Copyright (C) 2008 The QUESO Team, http://queso.ices.utexas.edu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __UQ_TGA_LIKELIHOOD_H__
#define __UQ_TGA_LIKELIHOOD_H__

#include <uqTgaDefines.h>
#include <uqTgaMeasuredW.h>
#include <uqTgaComputableW.h>
#include <uqTgaLambda.h>
#include <uqDefines.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>

//********************************************************
// Likelihood function object for both forward problems of the validation cycle.
// A likelihood function object is provided by user and is called by the UQ library.
// This likelihood function object consists of data and routine.
//********************************************************
template<class P_V,class P_M>
struct
uqTgaLikelihoodInfoStruct
{
  uqTgaLikelihoodInfoStruct(const uqVectorSpaceClass<P_V,P_M>& paramSpace,
                            const std::string&                 inpName);
 ~uqTgaLikelihoodInfoStruct();

  const uqVectorSpaceClass<P_V,P_M>& m_paramSpace;
  uqBase1D1DFunctionClass*           m_temperatureFunctionObj;
  uqTgaMeasuredWClass<P_V,P_M>*      m_referenceW;
};

template<class P_V,class P_M>
uqTgaLikelihoodInfoStruct<P_V,P_M>::uqTgaLikelihoodInfoStruct(
  const uqVectorSpaceClass<P_V,P_M>& paramSpace,
  const std::string&                 inpName)
  :
  m_paramSpace(paramSpace)
{
  if (paramSpace.env().rank() == 0) {
    std::cout << "Entering uqTgaLikelihoodInfoStruct<P_V,P_M>::constructor()"
              << inpName << "'\n"
              << std::endl;
  }

  // Read experimental data
  // Open input file on experimental data
  FILE *inp;
  inp = fopen(inpName.c_str(),"r");

  // Read kinetic parameters and convert heating rate to K/s
  double beta;
  double initialTemp;
  unsigned int numMeasurements;
  fscanf(inp,"%lf %lf %d",&beta,&initialTemp,&numMeasurements);
  beta /= 60.;
  m_temperatureFunctionObj = new uqLinear1D1DFunctionClass(-INFINITY,INFINITY,0.,initialTemp,beta);

  std::vector<double> measuredTimes(numMeasurements,0.);
  std::vector<double> measuredTemps(numMeasurements,0.);
  std::vector<double> measuredWs   (numMeasurements,0.);
  std::vector<double> measurementVs(numMeasurements,0.);

  unsigned int whileSize = 0;
  double tmpTemp;
  double tmpW;
  double tmpV;
  while (fscanf(inp,"%lf %lf %lf",&tmpTemp,&tmpW,&tmpV) != EOF) {
    UQ_FATAL_TEST_MACRO((whileSize >= numMeasurements),
                        paramSpace.env().rank(),
                        "uqTgaLikelihoodInfoStruct<P_V,P_M>::constructor(), in uqTgaValidation.h",
                        "input file 1 has too many measurements");
    measuredTimes[whileSize] = m_temperatureFunctionObj->inverseValue(tmpTemp);
    measuredTemps[whileSize] = tmpTemp;
    measuredWs   [whileSize] = tmpW;
    measurementVs[whileSize] = tmpV;
    whileSize++;
  }
  UQ_FATAL_TEST_MACRO((whileSize != numMeasurements),
                      paramSpace.env().rank(),
                      "uqTgaLikelihoodInfoStruct<P_V,P_M>::constructor(), in uqTgaValidation.h",
                      "input file 1 has a smaller number of measurements than expected");

  // Close input file on experimental data
  fclose(inp);

  m_referenceW = new uqTgaMeasuredWClass<P_V,P_M>(paramSpace.zeroVector(),
                                                  measuredTimes,
                                                  measuredTemps,
                                                  measuredWs,
                                                  &measurementVs,
                                                  false);

  if (paramSpace.env().rank() == 0) {
    std::cout << "Leaving uqTgaLikelihoodInfoStruct<P_V,P_M>::constructor()"
              << inpName << "'\n"
              << std::endl;
  }
}

template<class P_V,class P_M>
uqTgaLikelihoodInfoStruct<P_V,P_M>::~uqTgaLikelihoodInfoStruct()
{
  delete m_referenceW;
  delete m_temperatureFunctionObj;
}

// The actual (user defined) likelihood routine
template<class P_V,class P_M>
double
tgaLikelihoodRoutine(const P_V& paramValues, const void* functionDataPtr, P_V* gradVector, P_M* hessianMatrix, P_V* hessianEffect)
{
  if ((paramValues.env().verbosity() >= 30) && (paramValues.env().rank() == 0)) {
    std::cout << "Entering tgaLikelihoodRoutine()..." << std::endl;
  }

  if ((paramValues.env().verbosity() >= 10) && (paramValues.env().rank() == 0)) {
    std::cout << "In tgaLikelihoodRoutine()"
              << ", A = " << paramValues[0]
              << ", E = " << paramValues[1]
              << std::endl;
  }

  bool computeLambda = false;
  bool computeWAndLambdaGradsAlso = false;
  bool useOdesWithDerivativeWrtTime = false; // COMPATIBILITY WITH OLD VERSION

  if (gradVector) computeLambda = true;
  if (hessianMatrix || hessianEffect) {
    computeLambda              = true;
    computeWAndLambdaGradsAlso = true;
  }
  if (computeLambda || computeWAndLambdaGradsAlso) useOdesWithDerivativeWrtTime = true;

  double resultValue = 0.;
  const std::vector<uqTgaLikelihoodInfoStruct<P_V,P_M> *>& info = *((const std::vector<uqTgaLikelihoodInfoStruct<P_V,P_M> *> *) functionDataPtr);

  /////////////////////////////////////////////////////////////////////////////
  // Loop on scenarios
  /////////////////////////////////////////////////////////////////////////////
  for (unsigned int i = 0; i < info.size(); ++i) {
    //resultValue += tgaConstraintEquation<P_V,P_M>(paramValues,*(info[i]),true,NULL,NULL,NULL);

    /////////////////////////////////////////////////////////////////////////////
    // Compute W
    /////////////////////////////////////////////////////////////////////////////
    uqTgaComputableWClass<P_V,P_M> wObj(info[i]->m_paramSpace,
                                        *(info[i]->m_temperatureFunctionObj),
                                        useOdesWithDerivativeWrtTime);
    wObj.compute(paramValues,
                 computeWAndLambdaGradsAlso,
                 info[i]->m_referenceW,
                 1900.); // COMPATIBILITY WITH OLD VERSION

    /////////////////////////////////////////////////////////////////////////////
    // Compute contribution from this scenario to the likelihood ("tmpResultValue")
    /////////////////////////////////////////////////////////////////////////////
    unsigned int numMisfits = info[i]->m_referenceW->times().size();
    std::vector<double> twiceDiffVec(numMisfits,0.);
    double tmpResultValue = 0.; // COMPATIBILITY WITH OLD VERSION
    for (unsigned int j = 0; j < numMisfits; ++j) {
      double diff = wObj.diffsForMisfit()[j];
      double ratio = diff/info[i]->m_referenceW->variances()[j];
      twiceDiffVec[j] = 2.*ratio;
      tmpResultValue += diff*diff/info[i]->m_referenceW->variances()[j]; //COMPATIBILITY WITH OLD VERSION

      if ((paramValues.env().verbosity() >= 10) && (paramValues.env().rank() == 0)) {
        std::cout << "In tgaLikelihoodRoutine()"
                  << ", misfitId = "     << j
                  << ", measuredTemp = " << info[i]->m_referenceW->temps()[j]
                  << ", measuredW = "    << info[i]->m_referenceW->ws()[j]
                  << ": computedW = "    << info[i]->m_referenceW->ws()[j] + wObj.diffsForMisfit()[j]
                  << ", misfitValue = "  << diff
                  << std::endl;
      }
    }
    resultValue += tmpResultValue;

    if ((paramValues.env().verbosity() >= 10) && (paramValues.env().rank() == 0)) {
      char stringA[64];
      char stringE[64];
      sprintf(stringA,"%12.6e",paramValues[0]);
      sprintf(stringE,"%12.6e",paramValues[1]);
      std::cout << "In tgaLikelihoodRoutine()"
                << ", A = "                              << stringA
                << ", E = "                              << stringE
                << ", beta = "                           << info[i]->m_temperatureFunctionObj->deriv(0.)
                << ": finished ode loop after "          << wObj.ws().size()
                << " iterations, with weigthedMisfit = " << tmpResultValue
                << std::endl;
    }

    /////////////////////////////////////////////////////////////////////////////
    // Compute Lambda, if necessary
    /////////////////////////////////////////////////////////////////////////////
    uqTgaLambdaClass<P_V,P_M> lambdaObj(info[i]->m_paramSpace,
                                        *(info[i]->m_temperatureFunctionObj));
    if (computeLambda) {
      lambdaObj.compute(paramValues,
                        computeWAndLambdaGradsAlso,
                        twiceDiffVec,
                        info[i]->m_referenceW->continuous(),
                        wObj);
    }

    /////////////////////////////////////////////////////////////////////////////
    // Compute contribution from this scenario to gradVector
    /////////////////////////////////////////////////////////////////////////////
    if (gradVector) {
      P_V tmpGradVector(info[i]->m_paramSpace.zeroVector());
      *gradVector += tmpGradVector;
    }

    /////////////////////////////////////////////////////////////////////////////
    // Compute contribution from this scenario to hessianMatrix and/or hessianEffect
    /////////////////////////////////////////////////////////////////////////////
    if (hessianMatrix || hessianEffect) {
      P_M* tmpMatrix = info[i]->m_paramSpace.newMatrix();

      if (hessianMatrix) *hessianMatrix += *tmpMatrix;
      if (hessianEffect) *hessianEffect += ((*tmpMatrix) * paramValues);

      delete tmpMatrix;
    }
  }

  if ((paramValues.env().verbosity() >= 30) && (paramValues.env().rank() == 0)) {
    std::cout << "Leaving tgaLikelihoodRoutine()..." << std::endl;
  }

  return resultValue;
}
#endif // __UQ_TGA_LIKELIHOOD_H__
