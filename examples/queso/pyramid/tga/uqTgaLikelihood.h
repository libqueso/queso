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
#include <uqTgaStorage.h>
#include <uqTgaComputableW.h>
#include <uqTgaLambda.h>
#include <uqTgaIntegrals.h>
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

  void setReferenceW       (const uqTgaStorageClass<P_V,P_M>& referenceW);
  void setPlottingVariables(uqTgaStorageClass<P_V,P_M>* refWPtrForPlot,
                            uqTgaStorageClass<P_V,P_M>* wPtrForPlot,
                            uqTgaStorageClass<P_V,P_M>* wAPtrForPlot,
                            uqTgaStorageClass<P_V,P_M>* wEPtrForPlot,
                            uqTgaStorageClass<P_V,P_M>* misfitPtrForPlot,
                            uqTgaStorageClass<P_V,P_M>* lambdaPtrForPlot,
                            uqTgaStorageClass<P_V,P_M>* lambdaAPtrForPlot,
                            uqTgaStorageClass<P_V,P_M>* lambdaEPtrForPlot);

  const uqVectorSpaceClass <P_V,P_M>& m_paramSpace;
  uqBase1D1DFunctionClass*            m_temperatureFunctionObj;
  uqTgaStorageClass        <P_V,P_M>* m_referenceW;

  mutable uqTgaStorageClass<P_V,P_M>* m_refWPtrForPlot;
  mutable uqTgaStorageClass<P_V,P_M>* m_wPtrForPlot;
  mutable uqTgaStorageClass<P_V,P_M>* m_wAPtrForPlot;
  mutable uqTgaStorageClass<P_V,P_M>* m_wEPtrForPlot;
  mutable uqTgaStorageClass<P_V,P_M>* m_misfitPtrForPlot;
  mutable uqTgaStorageClass<P_V,P_M>* m_lambdaPtrForPlot;
  mutable uqTgaStorageClass<P_V,P_M>* m_lambdaAPtrForPlot;
  mutable uqTgaStorageClass<P_V,P_M>* m_lambdaEPtrForPlot;
};

template<class P_V,class P_M>
uqTgaLikelihoodInfoStruct<P_V,P_M>::uqTgaLikelihoodInfoStruct(
  const uqVectorSpaceClass<P_V,P_M>& paramSpace,
  const std::string&                 inpName)
  :
  m_paramSpace            (paramSpace),
  m_temperatureFunctionObj(NULL),
  m_referenceW            (NULL),
  m_refWPtrForPlot        (NULL),
  m_wPtrForPlot           (NULL),
  m_wAPtrForPlot          (NULL),
  m_wEPtrForPlot          (NULL),
  m_misfitPtrForPlot      (NULL),
  m_lambdaPtrForPlot      (NULL),
  m_lambdaAPtrForPlot     (NULL),
  m_lambdaEPtrForPlot     (NULL)
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

  m_referenceW = new uqTgaStorageClass<P_V,P_M>(measuredTimes,
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

template<class P_V,class P_M>
void
uqTgaLikelihoodInfoStruct<P_V,P_M>::setReferenceW(const uqTgaStorageClass<P_V,P_M>& referenceW)
{
  if (m_referenceW) delete m_referenceW;
  m_referenceW = new uqTgaStorageClass<P_V,P_M>(referenceW.times(),
                                                referenceW.temps(),
                                                referenceW.values(),
                                                &referenceW.variances(),
                                                referenceW.dataIsContinuousWithTime());
  return;
}

template<class P_V,class P_M>
void
uqTgaLikelihoodInfoStruct<P_V,P_M>::setPlottingVariables(
  uqTgaStorageClass<P_V,P_M>* refWPtrForPlot,
  uqTgaStorageClass<P_V,P_M>* wPtrForPlot,
  uqTgaStorageClass<P_V,P_M>* wAPtrForPlot,
  uqTgaStorageClass<P_V,P_M>* wEPtrForPlot,
  uqTgaStorageClass<P_V,P_M>* misfitPtrForPlot,
  uqTgaStorageClass<P_V,P_M>* lambdaPtrForPlot,
  uqTgaStorageClass<P_V,P_M>* lambdaAPtrForPlot,
  uqTgaStorageClass<P_V,P_M>* lambdaEPtrForPlot)
{
  m_refWPtrForPlot    = refWPtrForPlot;
  m_wPtrForPlot       = wPtrForPlot;
  m_wAPtrForPlot      = wAPtrForPlot;
  m_wEPtrForPlot      = wEPtrForPlot;
  m_misfitPtrForPlot  = misfitPtrForPlot;
  m_lambdaPtrForPlot  = lambdaPtrForPlot;
  m_lambdaAPtrForPlot = lambdaAPtrForPlot;
  m_lambdaEPtrForPlot = lambdaEPtrForPlot;

  return;
}

// The actual (user defined) likelihood routine
template<class P_V,class P_M>
double
uqTgaLikelihoodRoutine(
  const P_V&  paramValues,
  const void* functionDataPtr,
  P_V*        gradVector,
  P_M*        hessianMatrix,
  P_V*        hessianEffect)
{
  if ((paramValues.env().verbosity() >= 30) && (paramValues.env().rank() == 0)) {
    std::cout << "Entering uqTgaLikelihoodRoutine()..." << std::endl;
  }

  if ((paramValues.env().verbosity() >= 10) && (paramValues.env().rank() == 0)) {
    std::cout << "In uqTgaLikelihoodRoutine()"
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

  // Loop on scenarios
  for (unsigned int i = 0; i < info.size(); ++i) {
    /////////////////////////////////////////////////////////////////////////////
    // Compute W
    // Compute contribution from this scenario to the likelihood ("weigthedMisfitSum")
    /////////////////////////////////////////////////////////////////////////////
    uqTgaComputableWClass<P_V,P_M> wObj(info[i]->m_paramSpace,
                                        *(info[i]->m_temperatureFunctionObj));

    double weigthedMisfitSum = 0.; // COMPATIBILITY WITH OLD VERSION
    uqTgaStorageClass<P_V,P_M> weigthedMisfitData(info[i]->m_referenceW->dataIsContinuousWithTime());
    if (useOdesWithDerivativeWrtTime) {
      wObj.computeUsingTime(paramValues,
                            computeWAndLambdaGradsAlso,
                            info[i]->m_referenceW,
                            &weigthedMisfitSum,
                            &weigthedMisfitData);
    }
    else {
      wObj.computeUsingTemp(paramValues,
                            1900., // COMPATIBILITY WITH OLD VERSION
                            info[i]->m_referenceW,
                            &weigthedMisfitSum);
    }

    resultValue += weigthedMisfitSum;

    /////////////////////////////////////////////////////////////////////////////
    // Compute Lambda, if necessary
    /////////////////////////////////////////////////////////////////////////////
    uqTgaLambdaClass<P_V,P_M> lambdaObj(info[i]->m_paramSpace,
                                        *(info[i]->m_temperatureFunctionObj));
    if (computeLambda) {
      lambdaObj.compute(paramValues,
                        computeWAndLambdaGradsAlso,
                        weigthedMisfitData,
                        wObj);
    }

    /////////////////////////////////////////////////////////////////////////////
    // Compute contribution from this scenario to gradVector
    /////////////////////////////////////////////////////////////////////////////
    if (gradVector) {
      P_V tmpGradVector(info[i]->m_paramSpace.zeroVector());

      unsigned int tmpSize = weigthedMisfitData.times().size();
      double upperIntegralLimit = weigthedMisfitData.times()[tmpSize-1];
      uqTgaLagrangianGradientWrtDesignParameters<P_V,P_M>(paramValues,
                                                          *(info[i]->m_temperatureFunctionObj),
                                                          upperIntegralLimit,                  
                                                          wObj.times(),
                                                          wObj.ws(),
                                                          lambdaObj.times(),
                                                          lambdaObj.lambdas(),
                                                          tmpGradVector);

      *gradVector += tmpGradVector;
    }

    /////////////////////////////////////////////////////////////////////////////
    // Compute contribution from this scenario to hessianMatrix and/or hessianEffect
    /////////////////////////////////////////////////////////////////////////////
    if (hessianMatrix || hessianEffect) {
      P_M* tmpMatrix = info[i]->m_paramSpace.newMatrix();

      // AQUI
      UQ_FATAL_TEST_MACRO(true,
                          UQ_UNAVAILABLE_RANK,
                          "uqTgaLikelihoodRoutine<P_V,P_M>(), hessianMatrix",
                          "INCOMPLETE CODE");

      if (hessianMatrix) *hessianMatrix += *tmpMatrix;
      if (hessianEffect) *hessianEffect += ((*tmpMatrix) * paramValues);

      delete tmpMatrix;
    }

    /////////////////////////////////////////////////////////////////////////////
    // Store computed dara, if requested and if available
    /////////////////////////////////////////////////////////////////////////////
    if (info[i]->m_refWPtrForPlot) info[i]->m_refWPtrForPlot->set(info[i]->m_referenceW->times(),
                                                                  info[i]->m_referenceW->temps(),
                                                                  info[i]->m_referenceW->values(),
                                                                  &info[i]->m_referenceW->variances());

    if (info[i]->m_wPtrForPlot) info[i]->m_wPtrForPlot->set(wObj.times(),
                                                            wObj.temps(),
                                                            wObj.ws(),
                                                            NULL);

    if (info[i]->m_wAPtrForPlot) {
      unsigned int tmpSize = wObj.times().size();
      std::vector<double> tmpVec(tmpSize,0.);
      if (computeWAndLambdaGradsAlso) for (unsigned int j = 0; j < tmpSize; ++j) {
        tmpVec[j] = (*(wObj.grads()[j]))[0];
      }
      info[i]->m_wAPtrForPlot->set(wObj.times(),
                                   wObj.temps(),
                                   tmpVec,
                                   NULL);
    }

    if (info[i]->m_wEPtrForPlot) {
      unsigned int tmpSize = wObj.times().size();
      std::vector<double> tmpVec(tmpSize,0.);
      if (computeWAndLambdaGradsAlso) for (unsigned int j = 0; j < tmpSize; ++j) {
        tmpVec[j] = (*(wObj.grads()[j]))[1];
      }
      info[i]->m_wEPtrForPlot->set(wObj.times(),
                                   wObj.temps(),
                                   tmpVec,
                                   NULL);
    }

    if (info[i]->m_misfitPtrForPlot) info[i]->m_misfitPtrForPlot->set(weigthedMisfitData.times(),
                                                                      weigthedMisfitData.temps(),
                                                                      weigthedMisfitData.values(),
                                                                      NULL);

    if (info[i]->m_lambdaPtrForPlot) info[i]->m_lambdaPtrForPlot->set(lambdaObj.times(),
                                                                      lambdaObj.temps(),
                                                                      lambdaObj.lambdas(),
                                                                      NULL);

    if (info[i]->m_lambdaAPtrForPlot) {
      unsigned int tmpSize = lambdaObj.times().size();
      std::vector<double> tmpVec(tmpSize,0.);
      if (computeWAndLambdaGradsAlso) for (unsigned int j = 0; j < tmpSize; ++j) {
        tmpVec[j] = (*(lambdaObj.grads()[j]))[0];
      }
      info[i]->m_lambdaAPtrForPlot->set(lambdaObj.times(),
                                        lambdaObj.temps(),
                                        tmpVec,
                                        NULL);
    }

    if (info[i]->m_lambdaEPtrForPlot) {
      unsigned int tmpSize = lambdaObj.times().size();
      std::vector<double> tmpVec(tmpSize,0.);
      if (computeWAndLambdaGradsAlso) for (unsigned int j = 0; j < tmpSize; ++j) {
        tmpVec[j] = (*(lambdaObj.grads()[j]))[1];
      }
      info[i]->m_lambdaEPtrForPlot->set(lambdaObj.times(),
                                        lambdaObj.temps(),
                                        tmpVec,
                                        NULL);
    }
  }

  if ((paramValues.env().verbosity() >= 30) && (paramValues.env().rank() == 0)) {
    std::cout << "Leaving uqTgaLikelihoodRoutine()..." << std::endl;
  }

  return resultValue;
}
#endif // __UQ_TGA_LIKELIHOOD_H__
