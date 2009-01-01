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
#include <uqTgaW.h>
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
  void setCheckingFlag     (bool value);

  const uqVectorSpaceClass <P_V,P_M>& m_paramSpace;
  uqBase1D1DFunctionClass*            m_temperatureFunctionObj;
  uqTgaStorageClass        <P_V,P_M>* m_referenceW;
  uqTgaWClass              <P_V,P_M>* m_wObj;
  uqTgaLambdaClass         <P_V,P_M>* m_lambdaObj;
  uqTgaStorageClass        <P_V,P_M>* m_weigthedMisfitData;

  bool                        m_checkingFlag;
  uqTgaStorageClass<P_V,P_M>* m_refWPtrForPlot;
  uqTgaStorageClass<P_V,P_M>* m_wPtrForPlot;
  uqTgaStorageClass<P_V,P_M>* m_wAPtrForPlot;
  uqTgaStorageClass<P_V,P_M>* m_wEPtrForPlot;
  uqTgaStorageClass<P_V,P_M>* m_misfitPtrForPlot;
  uqTgaStorageClass<P_V,P_M>* m_lambdaPtrForPlot;
  uqTgaStorageClass<P_V,P_M>* m_lambdaAPtrForPlot;
  uqTgaStorageClass<P_V,P_M>* m_lambdaEPtrForPlot;
};

template<class P_V,class P_M>
uqTgaLikelihoodInfoStruct<P_V,P_M>::uqTgaLikelihoodInfoStruct(
  const uqVectorSpaceClass<P_V,P_M>& paramSpace,
  const std::string&                 inpName)
  :
  m_paramSpace            (paramSpace),
  m_temperatureFunctionObj(NULL),
  m_referenceW            (NULL),
  m_wObj                  (NULL),
  m_lambdaObj             (NULL),
  m_weigthedMisfitData    (NULL),
  m_checkingFlag          (false),
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

  m_wObj               = new uqTgaWClass      <P_V,P_M>(m_paramSpace, *m_temperatureFunctionObj);
  m_lambdaObj          = new uqTgaLambdaClass <P_V,P_M>(m_paramSpace, *m_temperatureFunctionObj);
  m_weigthedMisfitData = new uqTgaStorageClass<P_V,P_M>();

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
  if (m_refWPtrForPlot   ) delete m_refWPtrForPlot;
  if (m_wPtrForPlot      ) delete m_wPtrForPlot;
  if (m_wAPtrForPlot     ) delete m_wAPtrForPlot;
  if (m_wEPtrForPlot     ) delete m_wEPtrForPlot;
  if (m_misfitPtrForPlot ) delete m_misfitPtrForPlot;
  if (m_lambdaPtrForPlot ) delete m_lambdaPtrForPlot;
  if (m_lambdaAPtrForPlot) delete m_lambdaAPtrForPlot;
  if (m_lambdaEPtrForPlot) delete m_lambdaEPtrForPlot;

  delete m_referenceW;
  delete m_lambdaObj;
  delete m_wObj;
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
uqTgaLikelihoodInfoStruct<P_V,P_M>::setCheckingFlag(bool value)
{
  m_checkingFlag = value;
  return;
}

// Just declaration of checking routine here: actual code is below
template<class P_V,class P_M>
void
uqTgaLikelihoodChecking(
  const P_V&                          paramValues,
  const uqTgaWClass      <P_V,P_M>&   wObj,
  const uqTgaLambdaClass <P_V,P_M>&   lambdaObj,
  bool                                wAndLambdaGradsAreAlsoAvaiable,
  const uqTgaStorageClass<P_V,P_M>&   weigthedMisfitData,
  const P_V*                          tmpGradVector,
  const P_M*                          tmpHessianMatrix,
  const P_V*                          tmpHessianEffect,
  uqTgaLikelihoodInfoStruct<P_V,P_M>& info);

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

#ifdef QUESO_TGA_USES_OLD_COMPATIBLE_CODE
#else
  if (gradVector) computeLambda = true;
  if (hessianMatrix || hessianEffect) {
    computeLambda              = true;
    computeWAndLambdaGradsAlso = true;
  }
#endif

  double wMaxDeltaTime      = 0.;
  double lambdaMaxDeltaTime = 0.;
  if (computeLambda) {
    wMaxDeltaTime      = .1;
    lambdaMaxDeltaTime = 5.;
  }

  double resultValue = 0.;
  const std::vector<uqTgaLikelihoodInfoStruct<P_V,P_M> *>& vecInfo = *((const std::vector<uqTgaLikelihoodInfoStruct<P_V,P_M> *> *) functionDataPtr);

  // Loop on scenarios
  for (unsigned int i = 0; i < vecInfo.size(); ++i) {
    uqTgaLikelihoodInfoStruct<P_V,P_M>& info = *(vecInfo[i]);
    uqTgaWClass      <P_V,P_M>& wObj               = *(info.m_wObj);
    uqTgaLambdaClass <P_V,P_M>& lambdaObj          = *(info.m_lambdaObj);
    uqTgaStorageClass<P_V,P_M>& weigthedMisfitData = *(info.m_weigthedMisfitData);

    /////////////////////////////////////////////////////////////////////////////
    // Compute W
    // Compute contribution from this scenario to the likelihood ("weigthedMisfitSum")
    /////////////////////////////////////////////////////////////////////////////
    double weigthedMisfitSum = 0.;
#ifdef QUESO_TGA_USES_OLD_COMPATIBLE_CODE
    wObj.computeUsingTemp(paramValues,
                          1900., // COMPATIBILITY WITH OLD VERSION
                          info.m_referenceW,
                          &weigthedMisfitSum);
#else
    wObj.compute(paramValues,
                 0.,
                 wMaxDeltaTime,
                 computeWAndLambdaGradsAlso,
                 info.m_referenceW,
                 &weigthedMisfitSum,
                 &weigthedMisfitData);
#endif

    resultValue += weigthedMisfitSum;

    /////////////////////////////////////////////////////////////////////////////
    // Compute Lambda, if necessary
    /////////////////////////////////////////////////////////////////////////////
    if (computeLambda) {
      lambdaObj.compute(paramValues,
                        lambdaMaxDeltaTime,
                        computeWAndLambdaGradsAlso,
                        weigthedMisfitData,
                        wObj);
    }

    /////////////////////////////////////////////////////////////////////////////
    // Compute contribution from this scenario to gradVector
    // Compute contribution from this scenario to hessianMatrix and/or hessianEffect
    /////////////////////////////////////////////////////////////////////////////
    P_V* tmpVector = NULL;
    if (gradVector) tmpVector = new P_V(info.m_paramSpace.zeroVector());

    P_M* tmpMatrix = NULL;
    P_V* tmpEffect = NULL;
    if (hessianMatrix || hessianEffect) {
      tmpMatrix = info.m_paramSpace.newMatrix();
      tmpEffect = new P_V(info.m_paramSpace.zeroVector());
    }

    if (computeLambda) { // Same effect as "if (gradVector || hessianMatrix || hessianEffect) {"
      unsigned int tmpSize = weigthedMisfitData.times().size();
      double lowerIntegralLimit = 0.; //weigthedMisfitData.times()[0]; // ????
      double upperIntegralLimit = weigthedMisfitData.times()[tmpSize-1];
      uqTgaIntegrals<P_V,P_M>(paramValues,
                              *(info.m_temperatureFunctionObj),
                              lowerIntegralLimit,
                              upperIntegralLimit,
                              wObj,
                              lambdaObj,
                              tmpVector,
                              tmpMatrix);
      if (hessianEffect) *tmpEffect = ((*tmpMatrix) * paramValues);

      if (gradVector)    *gradVector    += *tmpVector;
      if (hessianMatrix) *hessianMatrix += *tmpMatrix;
      if (hessianEffect) *hessianEffect += *tmpEffect;
    }

    /////////////////////////////////////////////////////////////////////////////
    // Store computed data, if requested
    // Update some vector<double>'s in 'info', in order to be plotted afterwards
    /////////////////////////////////////////////////////////////////////////////
    if (info.m_checkingFlag) uqTgaLikelihoodChecking<P_V,P_M>(paramValues,
                                                              wObj,
                                                              lambdaObj,
                                                              computeWAndLambdaGradsAlso,
                                                              weigthedMisfitData,
                                                              tmpVector,
                                                              tmpMatrix,
                                                              tmpEffect,
                                                              info);
    if (tmpVector) delete tmpVector;
    if (tmpMatrix) delete tmpMatrix;
    if (tmpEffect) delete tmpEffect;
  }

  if ((paramValues.env().verbosity() >= 30) && (paramValues.env().rank() == 0)) {
    std::cout << "Leaving uqTgaLikelihoodRoutine()..." << std::endl;
  }

  return resultValue;
}

template<class P_V,class P_M>
void
uqTgaLikelihoodChecking(
  const P_V&                          paramValues,
  const uqTgaWClass      <P_V,P_M>&   wObj,
  const uqTgaLambdaClass <P_V,P_M>&   lambdaObj,
  bool                                wAndLambdaGradsAreAlsoAvaiable,
  const uqTgaStorageClass<P_V,P_M>&   weigthedMisfitData,
  const P_V*                          tmpGradVector,
  const P_M*                          tmpHessianMatrix,
  const P_V*                          tmpHessianEffect,
  uqTgaLikelihoodInfoStruct<P_V,P_M>& info)
{
  //////////////////////////////////////////////////////////////////
  // Step 1 of 3: Store data to be printed
  //////////////////////////////////////////////////////////////////

  if (info.m_refWPtrForPlot    == NULL) info.m_refWPtrForPlot    = new uqTgaStorageClass<P_V,P_M>();
  if (info.m_wPtrForPlot       == NULL) info.m_wPtrForPlot       = new uqTgaStorageClass<P_V,P_M>();
  if (info.m_wAPtrForPlot      == NULL) info.m_wAPtrForPlot      = new uqTgaStorageClass<P_V,P_M>();
  if (info.m_wEPtrForPlot      == NULL) info.m_wEPtrForPlot      = new uqTgaStorageClass<P_V,P_M>();
  if (info.m_misfitPtrForPlot  == NULL) info.m_misfitPtrForPlot  = new uqTgaStorageClass<P_V,P_M>();
  if (info.m_lambdaPtrForPlot  == NULL) info.m_lambdaPtrForPlot  = new uqTgaStorageClass<P_V,P_M>();
  if (info.m_lambdaAPtrForPlot == NULL) info.m_lambdaAPtrForPlot = new uqTgaStorageClass<P_V,P_M>();
  if (info.m_lambdaEPtrForPlot == NULL) info.m_lambdaEPtrForPlot = new uqTgaStorageClass<P_V,P_M>();

  info.m_refWPtrForPlot->set(info.m_referenceW->times(),
                             info.m_referenceW->temps(),
                             info.m_referenceW->values(),
                             &info.m_referenceW->variances(),
                             info.m_referenceW->dataIsContinuousWithTime());

  info.m_wPtrForPlot->set(wObj.times(),
                          wObj.temps(),
                          wObj.ws(),
                          NULL,
                          true);

  unsigned int wTmpSize = wObj.times().size();
  std::vector<double> wTmpVec(wTmpSize,0.);

  if (wAndLambdaGradsAreAlsoAvaiable) for (unsigned int j = 0; j < wTmpSize; ++j) {
    wTmpVec[j] = (*(wObj.grads()[j]))[0];
  }
  info.m_wAPtrForPlot->set(wObj.times(),
                           wObj.temps(),
                           wTmpVec,
                           NULL,
                           true);

  if (wAndLambdaGradsAreAlsoAvaiable) for (unsigned int j = 0; j < wTmpSize; ++j) {
    wTmpVec[j] = (*(wObj.grads()[j]))[1];
  }
  info.m_wEPtrForPlot->set(wObj.times(),
                           wObj.temps(),
                           wTmpVec,
                           NULL,
                           true);

  info.m_misfitPtrForPlot->set(weigthedMisfitData.times(),
                               weigthedMisfitData.temps(),
                               weigthedMisfitData.values(),
                               NULL,
                               weigthedMisfitData.dataIsContinuousWithTime());

  info.m_lambdaPtrForPlot->set(lambdaObj.times(),
                               lambdaObj.temps(),
                               lambdaObj.lambdas(),
                               NULL,
                               true);

  unsigned int lambdaTmpSize = lambdaObj.times().size();
  std::vector<double> lambdaTmpVec(lambdaTmpSize,0.);

  if (wAndLambdaGradsAreAlsoAvaiable) for (unsigned int j = 0; j < lambdaTmpSize; ++j) {
    lambdaTmpVec[j] = (*(lambdaObj.grads()[j]))[0];
  }
  info.m_lambdaAPtrForPlot->set(lambdaObj.times(),
                                lambdaObj.temps(),
                                lambdaTmpVec,
                                NULL,
                                true);

  if (wAndLambdaGradsAreAlsoAvaiable) for (unsigned int j = 0; j < lambdaTmpSize; ++j) {
    lambdaTmpVec[j] = (*(lambdaObj.grads()[j]))[1];
  }
  info.m_lambdaEPtrForPlot->set(lambdaObj.times(),
                                lambdaObj.temps(),
                                lambdaTmpVec,
                                NULL,
                                true);

  double guessA = paramValues[0];
  double guessE = paramValues[1];

  double deltaA = guessA * 1.e-6;
  double deltaE = guessE * 1.e-6;

  //////////////////////////////////////////////////////////////////
  // Step 2 of 3: Compute gradient(misfit) wrt A and E, using finite differences
  //////////////////////////////////////////////////////////////////
  if (tmpGradVector) {
    P_V tmpParamValues(paramValues);

    if (tmpParamValues.env().rank() == 0) {
      std::cout << "\nIn uqTgaLikelihoodChecking()"
                << ": gradWithLM = " << *tmpGradVector
                << std::endl;
    }

    if (tmpParamValues.env().rank() == 0) {
      std::cout << "In uqTgaLikelihoodChecking()"
                << ": computing gradient(misfit) w.r.t. parameters A and E, using finite differences..."
                << std::endl;
    }

    uqTgaWClass<P_V,P_M> tmpW(info.m_paramSpace,
                              *(info.m_temperatureFunctionObj));

    tmpParamValues[0] = guessA-deltaA;
    tmpParamValues[1] = guessE;
    double valueAm = 0.;
    tmpW.compute(tmpParamValues,
                 0.,
                 .1,
                 false, // computeWAndLambdaGradsAlso
                 info.m_referenceW,
                 &valueAm,
                 NULL);
    //std::cout << "valueAm = " << valueAm << std::endl;

    tmpParamValues[0] = guessA+deltaA;
    tmpParamValues[1] = guessE;
    double valueAp = 0.;
    tmpW.compute(tmpParamValues,
                 0.,
                 .1,
                 false, // computeWAndLambdaGradsAlso
                 info.m_referenceW,
                 &valueAp,
                 NULL);
    //std::cout << "valueAp = " << valueAp << std::endl;

    tmpParamValues[0] = guessA;
    tmpParamValues[1] = guessE-deltaE;
    double valueEm = 0.;
    tmpW.compute(tmpParamValues,
                 0.,
                 .1,
                 false, // computeWAndLambdaGradsAlso
                 info.m_referenceW,
                 &valueEm,
                 NULL);
    //std::cout << "valueEm = " << valueEm << std::endl;

    tmpParamValues[0] = guessA;
    tmpParamValues[1] = guessE+deltaE;
    double valueEp = 0.;
    tmpW.compute(tmpParamValues,
                 0.,
                 .1,
                 false, // computeWAndLambdaGradsAlso
                 info.m_referenceW,
                 &valueEp,
                 NULL);
    //std::cout << "valueEp = " << valueEp << std::endl;

    P_V gradWithFD(info.m_paramSpace.zeroVector());
    gradWithFD[0] = (valueAp-valueAm)/2./deltaA;
    gradWithFD[1] = (valueEp-valueEm)/2./deltaE;
    if (tmpParamValues.env().rank() == 0) {
      std::cout << "\nIn uqTgaLikelihoodChecking()"
                << ": gradWithFD = " << gradWithFD
                << "\ngrad relative error = " << (gradWithFD - *tmpGradVector).norm2()/gradWithFD.norm2()
                << std::endl;
    }
  }

  //////////////////////////////////////////////////////////////////
  // Step 3 of 3: Compute Hessian(misfit), wrt A and E, using finite differences
  //////////////////////////////////////////////////////////////////

  if (tmpHessianMatrix) {
    P_V tmpParamValues(paramValues);
    bool savedCheckingFlag = info.m_checkingFlag;
    info.m_checkingFlag = false; // IMPORTANT

    std::vector<uqTgaLikelihoodInfoStruct<P_V,P_M>* > tmpCalLikelihoodInfoVector(1,NULL);
    tmpCalLikelihoodInfoVector[0] = &info;

    if (tmpParamValues.env().rank() == 0) {
      std::cout << "\nIn uqTgaLikelihoodChecking()"
                << ": hessianWithLM = " << *tmpHessianMatrix
                << std::endl;
    }

    // Compute \grad(misfit), using finite differences
    if (tmpParamValues.env().rank() == 0) {
      std::cout << "In uqTgaLikelihoodChecking()"
                << ": computing Hessian(misfit) w.r.t. parameters A and E, using finite differences..."
                << std::endl;
    }

    tmpParamValues[0] = guessA-deltaA;
    tmpParamValues[1] = guessE;
    P_V grad_Am(info.m_paramSpace.zeroVector());
    double tmpMisfit = uqTgaLikelihoodRoutine<P_V,P_M>(tmpParamValues,
                                                       (const void *)&tmpCalLikelihoodInfoVector,
                                                       &grad_Am,
                                                       NULL, // Hessian
                                                       NULL);
    //std::cout << "grad_Am = " << grad_Am << std::endl;

    tmpParamValues[0] = guessA+deltaA;
    tmpParamValues[1] = guessE;
    P_V grad_Ap(info.m_paramSpace.zeroVector());
    tmpMisfit = uqTgaLikelihoodRoutine<P_V,P_M>(tmpParamValues,
                                                (const void *)&tmpCalLikelihoodInfoVector,
                                                &grad_Ap,
                                                NULL, // Hessian
                                                NULL);
    //std::cout << "grad_Ap = " << grad_Ap << std::endl;

    tmpParamValues[0] = guessA;
    tmpParamValues[1] = guessE-deltaE;
    P_V grad_Em(info.m_paramSpace.zeroVector());
    tmpMisfit = uqTgaLikelihoodRoutine<P_V,P_M>(tmpParamValues,
                                                (const void *)&tmpCalLikelihoodInfoVector,
                                                &grad_Em,
                                                NULL, // Hessian
                                                NULL);
    //std::cout << "grad_Em = " << grad_Em << std::endl;

    tmpParamValues[0] = guessA;
    tmpParamValues[1] = guessE+deltaE;
    P_V grad_Ep(info.m_paramSpace.zeroVector());
    tmpMisfit = uqTgaLikelihoodRoutine<P_V,P_M>(tmpParamValues,
                                                (const void *)&tmpCalLikelihoodInfoVector,
                                                &grad_Ep,
                                                NULL, // Hessian
                                                NULL);
    //std::cout << "grad_Ep = " << grad_Ep << std::endl;

    P_V column1((.5/deltaA)*(grad_Ap-grad_Am));
    P_V column2((.5/deltaE)*(grad_Ep-grad_Em));
    P_M* hessianWithFD = info.m_paramSpace.newMatrix();
    (*hessianWithFD)(0,0)=column1[0];
    (*hessianWithFD)(1,0)=column1[1];
    (*hessianWithFD)(0,1)=column2[0];
    (*hessianWithFD)(1,1)=column2[1];
    if (tmpParamValues.env().rank() == 0) {
      std::cout << "\nIn uqTgaLikelihoodChecking()"
                << ": hessianWithFD = " << *hessianWithFD
                << "\nhessian absolute error = " << (*hessianWithFD + (-1.*(*tmpHessianMatrix)))
                << std::endl;
    }
    delete hessianWithFD;

    info.m_checkingFlag = savedCheckingFlag; // IMPORTANT
  }

  return;
}
#endif // __UQ_TGA_LIKELIHOOD_H__
