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
                            const std::string&                 inpName,
                            double*                            wMaxTimeStep,
                            double*                            lambdaMaxTimeStep,
                            unsigned int*                      integralsNumIntervals);
 ~uqTgaLikelihoodInfoStruct();

#ifdef QUESO_TGA_USES_OLD_COMPATIBLE_CODE
#else
  void changeReferenceW     (const uqBase1D1DFunctionClass& referenceW);
#endif
  void changeWeightFunctions(const uqBase1D1DFunctionClass* weightFunction,
                             const uqBase1D1DFunctionClass* tildeWeightFunction);
  void setCheckingVariables (bool performChecking, double* relativeFDStep);

  const uqVectorSpaceClass <P_V,P_M>& m_paramSpace;
  uqBase1D1DFunctionClass*            m_temperatureFunctionObj;
  bool                                m_refWIsTheInternallyCreated;
  bool                                m_weightFunctionsAreTheInternallyCreated;
#ifdef QUESO_TGA_USES_OLD_COMPATIBLE_CODE
  uqTgaStorageClass        <P_V,P_M>* m_referenceW;
#else
  const uqBase1D1DFunctionClass*      m_referenceW;
#endif
  const uqBase1D1DFunctionClass*      m_weightFunction;
  const uqBase1D1DFunctionClass*      m_tildeWeightFunction;
  double                              m_wMaxTimeStep;
  double                              m_lambdaMaxTimeStep;
  unsigned int                        m_integralsNumIntervals;

  // Internal working objects
  uqTgaWClass     <P_V,P_M>*  m_wObj;
  uqSampled1D1DFunctionClass* m_diffFunction;
  uqTgaLambdaClass<P_V,P_M>*  m_lambdaObj;

  bool                        m_performChecking;
  double                      m_relativeFDStep;
  uqSampled1D1DFunctionClass* m_refWPtrForPlot;
  uqSampled1D1DFunctionClass* m_wPtrForPlot;
  uqSampled1D1DFunctionClass* m_wAPtrForPlot;
  uqSampled1D1DFunctionClass* m_wEPtrForPlot;
  uqSampled1D1DFunctionClass* m_diffPtrForPlot;
  uqSampled1D1DFunctionClass* m_lambdaPtrForPlot;
  uqSampled1D1DFunctionClass* m_lambdaAPtrForPlot;
  uqSampled1D1DFunctionClass* m_lambdaEPtrForPlot;
};

template<class P_V,class P_M>
uqTgaLikelihoodInfoStruct<P_V,P_M>::uqTgaLikelihoodInfoStruct(
  const uqVectorSpaceClass<P_V,P_M>& paramSpace,
  const std::string&                 inpName,
  double*                            wMaxTimeStep,
  double*                            lambdaMaxTimeStep,
  unsigned int*                      integralsNumIntervals)
  :
  m_paramSpace                            (paramSpace),
  m_temperatureFunctionObj                (NULL),
  m_refWIsTheInternallyCreated            (true),
  m_weightFunctionsAreTheInternallyCreated(true),
  m_referenceW                            (NULL),
  m_weightFunction                        (NULL),
  m_tildeWeightFunction                   (NULL),
  m_wMaxTimeStep                          (0.), // IMPORTANT
  m_lambdaMaxTimeStep                     (0.), // IMPORTANT
  m_integralsNumIntervals                 (1000),
  m_wObj                                  (NULL),
  m_diffFunction                          (NULL),
  m_lambdaObj                             (NULL),
  m_performChecking                       (false),
  m_relativeFDStep                        (1.e-6),
  m_refWPtrForPlot                        (NULL),
  m_wPtrForPlot                           (NULL),
  m_wAPtrForPlot                          (NULL),
  m_wEPtrForPlot                          (NULL),
  m_diffPtrForPlot                        (NULL),
  m_lambdaPtrForPlot                      (NULL),
  m_lambdaAPtrForPlot                     (NULL),
  m_lambdaEPtrForPlot                     (NULL)
{
  if ((paramSpace.env().verbosity() >= 30) && (paramSpace.env().rank() == 0)) {
    std::cout << "Entering uqTgaLikelihoodInfoStruct<P_V,P_M>::constructor()"
              << inpName << "'\n"
              << std::endl;
  }

  if (wMaxTimeStep         ) m_wMaxTimeStep          = *wMaxTimeStep;
  if (lambdaMaxTimeStep    ) m_lambdaMaxTimeStep     = *lambdaMaxTimeStep;
  if (integralsNumIntervals) m_integralsNumIntervals = *integralsNumIntervals;

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

  m_wObj         = new uqTgaWClass     <P_V,P_M> (m_paramSpace, *m_temperatureFunctionObj);
  m_diffFunction = new uqSampled1D1DFunctionClass();
  m_lambdaObj    = new uqTgaLambdaClass<P_V,P_M> (m_paramSpace, *m_temperatureFunctionObj);

  std::vector<double> measuredTimes(numMeasurements,0.);
  std::vector<double> measuredTemps(numMeasurements,0.);
  std::vector<double> measuredWs   (numMeasurements,0.);
  std::vector<double> measurementVs(numMeasurements,0.);
  std::vector<double> vecOfDeltas  (numMeasurements,0.);
  std::vector<double> vecOfWeights (numMeasurements,0.);

  unsigned int whileSize = 0;
  double tmpTemp;
  double tmpW;
  double tmpV;
  while (fscanf(inp,"%lf %lf %lf",&tmpTemp,&tmpW,&tmpV) != EOF) {
    UQ_FATAL_TEST_MACRO((whileSize >= numMeasurements),
                        paramSpace.env().rank(),
                        "uqTgaLikelihoodInfoStruct<P_V,P_M>::constructor(), in uqTgaValidation.h",
                        "input file 1 has too many measurements");
    measuredTimes[whileSize] = (tmpTemp-initialTemp)/beta;
    measuredTemps[whileSize] = tmpTemp;
    measuredWs   [whileSize] = tmpW;
    measurementVs[whileSize] = tmpV;
    vecOfDeltas  [whileSize] = INFINITY;
    vecOfWeights [whileSize] = 1./tmpV;
    whileSize++;
  }
  UQ_FATAL_TEST_MACRO((whileSize != numMeasurements),
                      paramSpace.env().rank(),
                      "uqTgaLikelihoodInfoStruct<P_V,P_M>::constructor(), in uqTgaValidation.h",
                      "input file 1 has a smaller number of measurements than expected");

  // Close input file on experimental data
  fclose(inp);

#ifdef QUESO_TGA_USES_OLD_COMPATIBLE_CODE
  m_referenceW = new uqTgaStorageClass<P_V,P_M>(measuredTimes,
                                                measuredTemps,
                                                measuredWs,
                                                &measurementVs,
                                                false);
#else
  m_referenceW = new uqSampled1D1DFunctionClass(measuredTimes,
                                                measuredWs);

  m_weightFunction = new uqDeltaSeq1D1DFunctionClass(measuredTimes,
                                                     vecOfDeltas,
                                                     vecOfWeights);

  if ((paramSpace.env().verbosity() >= 30) && (paramSpace.env().rank() == 0)) {
    std::cout << "In uqTgaLikelihoodInfoStruct<P_V,P_M>::constructor()"
              << ": m_weightFunction.times().size() = "     << measuredTimes.size()
              << ", m_weightFunction.times()[0] = "         << measuredTimes[0]
              << ", m_weightFunction.deltaValues()[0] = "   << vecOfDeltas  [0]
              << ", m_weightFunction.intValues()[0] = "     << vecOfWeights [0]
              << ", m_weightFunction.times()[max] = "       << measuredTimes[measuredTimes.size()-1]
              << ", m_weightFunction.deltaValues()[max] = " << vecOfDeltas  [measuredTimes.size()-1]
              << ", m_weightFunction.intValues()[max] = "   << vecOfWeights [measuredTimes.size()-1]
              << std::endl;
  }

  std::vector<double> aux1(whileSize,0.);
  std::vector<double> aux2(whileSize,0.);
  std::vector<double> aux3(whileSize,0.);
  double maxTime = m_weightFunction->maxDomainValue();
  for (unsigned int j = 0; j < whileSize; ++j) {
    aux1[j] = maxTime-measuredTimes[whileSize-1-j];
    aux2[j] = vecOfDeltas          [whileSize-1-j];
    aux3[j] = vecOfWeights         [whileSize-1-j];
  }
  m_tildeWeightFunction = new uqDeltaSeq1D1DFunctionClass(aux1,
                                                          aux2,
                                                          aux3);
#endif

  if ((paramSpace.env().verbosity() >= 30) && (paramSpace.env().rank() == 0)) {
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
  if (m_diffPtrForPlot   ) delete m_diffPtrForPlot;
  if (m_lambdaPtrForPlot ) delete m_lambdaPtrForPlot;
  if (m_lambdaAPtrForPlot) delete m_lambdaAPtrForPlot;
  if (m_lambdaEPtrForPlot) delete m_lambdaEPtrForPlot;

  delete m_lambdaObj;
  delete m_diffFunction;
  delete m_wObj;
  if (m_weightFunctionsAreTheInternallyCreated) delete m_tildeWeightFunction;
  if (m_weightFunctionsAreTheInternallyCreated) delete m_weightFunction;
  if (m_refWIsTheInternallyCreated) delete m_referenceW;
  delete m_temperatureFunctionObj;
}

#ifdef QUESO_TGA_USES_OLD_COMPATIBLE_CODE
#else
template<class P_V,class P_M>
void
uqTgaLikelihoodInfoStruct<P_V,P_M>::changeReferenceW(const uqBase1D1DFunctionClass& referenceW)
{
  if (m_refWIsTheInternallyCreated) {
    delete m_referenceW;
    m_refWIsTheInternallyCreated = false;
  }
  m_referenceW = &referenceW;

  return;
}
#endif

template<class P_V,class P_M>
void
uqTgaLikelihoodInfoStruct<P_V,P_M>::changeWeightFunctions(
  const uqBase1D1DFunctionClass* weightFunction,
  const uqBase1D1DFunctionClass* tildeWeightFunction)
{
  if (m_weightFunctionsAreTheInternallyCreated) {
    delete m_weightFunction;
    delete m_tildeWeightFunction;
    m_weightFunctionsAreTheInternallyCreated = false;
  }
  m_weightFunction      = weightFunction;
  m_tildeWeightFunction = tildeWeightFunction;

  return;
}

template<class P_V,class P_M>
void
uqTgaLikelihoodInfoStruct<P_V,P_M>::setCheckingVariables(
  bool    performChecking,
  double* relativeFDStep)
{
  m_performChecking = performChecking;
  if (relativeFDStep) m_relativeFDStep = *relativeFDStep;

  return;
}

// Just declaration of checking routine here: actual code is below
template<class P_V,class P_M>
void
uqTgaLikelihoodCheckingRoutine(
  const P_V&                          paramValues,
  const P_V*                          paramDirection,
  const uqTgaWClass      <P_V,P_M>&   wObj,
  const uqTgaLambdaClass <P_V,P_M>&   lambdaObj,
  bool                                wAndLambdaGradsAreAlsoAvaiable,
  const uqSampled1D1DFunctionClass&   diffFunction,
  const P_V*                          tmpGradVector,
  const P_M*                          tmpHessianMatrix,
  const P_V*                          tmpHessianEffect,
  uqTgaLikelihoodInfoStruct<P_V,P_M>& info);

// The actual (user defined) likelihood routine
template<class P_V,class P_M>
double
uqTgaLikelihoodRoutine(
  const P_V&  paramValues,
  const P_V*  paramDirection,
  const void* functionDataPtr,
  P_V*        gradVector,
  P_M*        hessianMatrix,
  P_V*        hessianEffect)
{
  UQ_FATAL_TEST_MACRO((hessianEffect != NULL) && (paramDirection == NULL),
                      paramValues.env().rank(),
                      "uqTgaLikelihoodRoutine<P_V,P_M>()",
                      "hessianEffect is being requested but no paramDirection is supplied");

  // Set all possible return values to zero
  double resultValue = 0.;
  if (gradVector   ) *gradVector    *= 0.;
  if (hessianMatrix) *hessianMatrix *= 0.;
  if (hessianEffect) *hessianEffect *= 0.;

  if ((paramValues.env().verbosity() >= 30) && (paramValues.env().rank() == 0)) {
    std::cout << "Entering uqTgaLikelihoodRoutine()..." << std::endl;
  }

  if ((paramValues.env().verbosity() >= 10) && (paramValues.env().rank() == 0)) {
    std::cout << "In uqTgaLikelihoodRoutine()"
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

  // Loop on scenarios
  const std::vector<uqTgaLikelihoodInfoStruct<P_V,P_M> *>& vecInfo = *((const std::vector<uqTgaLikelihoodInfoStruct<P_V,P_M> *> *) functionDataPtr);
  for (unsigned int i = 0; i < vecInfo.size(); ++i) {
    uqTgaLikelihoodInfoStruct<P_V,P_M>& info = *(vecInfo[i]);
    uqTgaWClass     <P_V,P_M>&  wObj         = *(info.m_wObj);
    uqSampled1D1DFunctionClass& diffFunction = *(info.m_diffFunction);
    uqTgaLambdaClass<P_V,P_M>&  lambdaObj    = *(info.m_lambdaObj);

    if ((paramValues.env().verbosity() >= 10) && (paramValues.env().rank() == 0)) {
      if (info.m_weightFunction) {
        const uqDeltaSeq1D1DFunctionClass* deltaWeight = dynamic_cast< const uqDeltaSeq1D1DFunctionClass* >(info.m_weightFunction);
        if (deltaWeight != NULL) {
          unsigned int tmpSize = deltaWeight->domainValues().size();
          std::cout << "In uqTgaLikelihoodRoutine()"
                    << ", i = "                               << i
                    << ": deltaWeight->times().size() = "     << tmpSize
                    << ", deltaWeight->times()[0] = "         << deltaWeight->domainValues    ()[0]
                    << ", deltaWeight->deltaValues()[0] = "   << deltaWeight->imageValues     ()[0]
                    << ", deltaWeight->intValues()[0] = "     << deltaWeight->integratedValues()[0]
                    << ", deltaWeight->times()[max] = "       << deltaWeight->domainValues    ()[tmpSize-1]
                    << ", deltaWeight->deltaValues()[max] = " << deltaWeight->imageValues     ()[tmpSize-1]
                    << ", deltaWeight->intValues()[max] = "   << deltaWeight->integratedValues()[tmpSize-1]
                    << std::endl;
        }
      }
    }

    /////////////////////////////////////////////////////////////////////////////
    // Compute W
    // Compute contribution from this scenario to the likelihood ("misfitValue")
    /////////////////////////////////////////////////////////////////////////////
    double misfitValue = 0.;
#ifdef QUESO_TGA_USES_OLD_COMPATIBLE_CODE
    wObj.computeUsingTemp(paramValues,
                          1900., // COMPATIBILITY WITH OLD VERSION
                          info.m_referenceW,
                          &misfitValue);
#else
    wObj.compute(paramValues,
                 paramDirection,
                 1.,
                 0.,
                 info.m_wMaxTimeStep,
                 computeWAndLambdaGradsAlso,
                 info.m_referenceW,
                 info.m_weightFunction,
                 &misfitValue,
                 &diffFunction);
#endif

    resultValue += misfitValue;

    /////////////////////////////////////////////////////////////////////////////
    // Compute Lambda, if necessary
    /////////////////////////////////////////////////////////////////////////////
    if (computeLambda) {
      lambdaObj.compute(paramValues,
                        paramDirection,
                        info.m_lambdaMaxTimeStep,
                        computeWAndLambdaGradsAlso,
                        diffFunction,
                        info.m_tildeWeightFunction,
                        wObj);
    }

    /////////////////////////////////////////////////////////////////////////////
    // Compute contribution from this scenario to gradVector
    // Compute contribution from this scenario to hessianMatrix and/or hessianEffect
    /////////////////////////////////////////////////////////////////////////////
    P_V* tmpVector = NULL;
    if (gradVector) tmpVector = new P_V(info.m_paramSpace.zeroVector());

    P_M* tmpMatrix = NULL;
    if (hessianMatrix) tmpMatrix = info.m_paramSpace.newMatrix();

    P_V* tmpEffect = NULL;
    if (hessianEffect) tmpEffect = new P_V(info.m_paramSpace.zeroVector());

    if (computeLambda) { // Same effect as "if (gradVector || hessianMatrix || hessianEffect) {"
      double lowerIntegralLimit = 0.; //diffFunction.minDomainValue(); // ????
      double upperIntegralLimit = diffFunction.maxDomainValue(); // FIX ME
      if (info.m_weightFunction) upperIntegralLimit = info.m_weightFunction->maxDomainValue(); // FIX ME
      uqTgaIntegrals<P_V,P_M>(paramValues,
                              paramDirection,
                              *(info.m_temperatureFunctionObj),
                              lowerIntegralLimit,
                              upperIntegralLimit,
                              info.m_integralsNumIntervals,
                              wObj,
                              lambdaObj,
                              info.m_weightFunction,
                              tmpVector,
                              tmpMatrix,
                              tmpEffect);
      if (gradVector)    *gradVector    += *tmpVector;
      if (hessianMatrix) *hessianMatrix += *tmpMatrix;
      if (hessianEffect) *hessianEffect += *tmpEffect;
    }

    /////////////////////////////////////////////////////////////////////////////
    // Store computed data, if requested
    // Update some vector<double>'s in 'info', in order to be plotted afterwards
    /////////////////////////////////////////////////////////////////////////////
    if (info.m_performChecking) uqTgaLikelihoodCheckingRoutine<P_V,P_M>(paramValues,
                                                                        paramDirection,
                                                                        wObj,
                                                                        lambdaObj,
                                                                        computeWAndLambdaGradsAlso,
                                                                        diffFunction,
                                                                        tmpVector,
                                                                        tmpMatrix,
                                                                        tmpEffect,
                                                                        info);
    if (tmpVector) delete tmpVector;
    if (tmpMatrix) delete tmpMatrix;
    if (tmpEffect) delete tmpEffect;
  } // end loop of scenarios

  if ((paramValues.env().verbosity() >= 30) && (paramValues.env().rank() == 0)) {
    std::cout << "Leaving uqTgaLikelihoodRoutine()..." << std::endl;
  }

  return resultValue;
}

template<class P_V,class P_M>
void
uqTgaLikelihoodCheckingRoutine(
  const P_V&                          paramValues,
  const P_V*                          paramDirection,
  const uqTgaWClass      <P_V,P_M>&   wObj,
  const uqTgaLambdaClass <P_V,P_M>&   lambdaObj,
  bool                                wAndLambdaGradsAreAlsoAvaiable,
  const uqSampled1D1DFunctionClass&   diffFunction,
  const P_V*                          tmpGradVector,
  const P_M*                          tmpHessianMatrix,
  const P_V*                          tmpHessianEffect,
  uqTgaLikelihoodInfoStruct<P_V,P_M>& info)
{
  UQ_FATAL_TEST_MACRO((tmpHessianEffect != NULL) && (paramDirection == NULL),
                      paramValues.env().rank(),
                      "uqTgaLikelihoodCheckingRoutine<P_V,P_M>()",
                      "tmpHessianEffect is being requested but no paramDirection is supplied");

  //////////////////////////////////////////////////////////////////
  // Step 1 of 4: Store data to be printed
  //////////////////////////////////////////////////////////////////

  if (info.m_refWPtrForPlot    == NULL) info.m_refWPtrForPlot    = new uqSampled1D1DFunctionClass();
  if (info.m_wPtrForPlot       == NULL) info.m_wPtrForPlot       = new uqSampled1D1DFunctionClass();
  if (info.m_wAPtrForPlot      == NULL) info.m_wAPtrForPlot      = new uqSampled1D1DFunctionClass();
  if (info.m_wEPtrForPlot      == NULL) info.m_wEPtrForPlot      = new uqSampled1D1DFunctionClass();
  if (info.m_diffPtrForPlot    == NULL) info.m_diffPtrForPlot    = new uqSampled1D1DFunctionClass();
  if (info.m_lambdaPtrForPlot  == NULL) info.m_lambdaPtrForPlot  = new uqSampled1D1DFunctionClass();
  if (info.m_lambdaAPtrForPlot == NULL) info.m_lambdaAPtrForPlot = new uqSampled1D1DFunctionClass();
  if (info.m_lambdaEPtrForPlot == NULL) info.m_lambdaEPtrForPlot = new uqSampled1D1DFunctionClass();

#ifdef QUESO_TGA_USES_OLD_COMPATIBLE_CODE
#else
  const uqSampled1D1DFunctionClass* tmpRef = dynamic_cast< const uqSampled1D1DFunctionClass* >(info.m_referenceW);
  if (tmpRef != NULL) {
    info.m_refWPtrForPlot->set(tmpRef->domainValues(),
                               tmpRef->imageValues());
  }
  else {
    unsigned int refTmpSize = wObj.times().size();
    std::vector<double> refTmpVec(refTmpSize,0.);
    for (unsigned int j = 0; j < refTmpSize; ++j) {
      refTmpVec[j] = info.m_referenceW->value(wObj.times()[j]); // Might be slow
    }
    info.m_refWPtrForPlot->set(wObj.times(),
                               refTmpVec);
  }
#endif

  info.m_wPtrForPlot->set(wObj.times(),
                          wObj.ws());

  unsigned int wTmpSize = wObj.times().size();
  std::vector<double> wTmpVec(wTmpSize,0.);

  if (wAndLambdaGradsAreAlsoAvaiable) for (unsigned int j = 0; j < wTmpSize; ++j) {
    wTmpVec[j] = (*(wObj.grads()[j]))[0];
  }
  info.m_wAPtrForPlot->set(wObj.times(),
                           wTmpVec);

  if (wAndLambdaGradsAreAlsoAvaiable) for (unsigned int j = 0; j < wTmpSize; ++j) {
    wTmpVec[j] = (*(wObj.grads()[j]))[1];
  }
  info.m_wEPtrForPlot->set(wObj.times(),
                           wTmpVec);

  info.m_diffPtrForPlot->set(diffFunction.domainValues(),
                             diffFunction.imageValues());

  info.m_lambdaPtrForPlot->set(lambdaObj.times(),
                               lambdaObj.lambdas());

  unsigned int lambdaTmpSize = lambdaObj.times().size();
  std::vector<double> lambdaTmpVec(lambdaTmpSize,0.);

  if (wAndLambdaGradsAreAlsoAvaiable) for (unsigned int j = 0; j < lambdaTmpSize; ++j) {
    lambdaTmpVec[j] = (*(lambdaObj.grads()[j]))[0];
  }
  info.m_lambdaAPtrForPlot->set(lambdaObj.times(),
                                lambdaTmpVec);

  if (wAndLambdaGradsAreAlsoAvaiable) for (unsigned int j = 0; j < lambdaTmpSize; ++j) {
    lambdaTmpVec[j] = (*(lambdaObj.grads()[j]))[1];
  }
  info.m_lambdaEPtrForPlot->set(lambdaObj.times(),
                                lambdaTmpVec);

  double guessA = paramValues[0];
  double guessE = paramValues[1];

  if (paramValues.env().rank() == 0) {
    std::cout << "\nIn uqTgaLikelihoodCheckingRoutine()"
              << ": info.m_relativeFDStep = " << info.m_relativeFDStep
              << std::endl;
  }

  double deltaA = guessA * info.m_relativeFDStep;
  double deltaE = guessE * info.m_relativeFDStep;

  //////////////////////////////////////////////////////////////////
  // Step 2 of 4: Compute gradient(misfit) wrt A and E, using finite differences
  //////////////////////////////////////////////////////////////////
  if (tmpGradVector) {
    P_V tmpParamValues(paramValues);

    if (paramValues.env().rank() == 0) {
      std::cout << "\nIn uqTgaLikelihoodCheckingRoutine()"
                << ": gradWithLM = " << *tmpGradVector
                << std::endl;
    }

    if (paramValues.env().rank() == 0) {
      std::cout << "In uqTgaLikelihoodCheckingRoutine()"
                << ": computing gradient(misfit) w.r.t. parameters A and E, using finite differences..."
                << std::endl;
    }

    uqTgaWClass<P_V,P_M> tmpW(info.m_paramSpace,
                              *(info.m_temperatureFunctionObj));

    tmpParamValues[0] = guessA-deltaA;
    tmpParamValues[1] = guessE;
    double valueAm = 0.;
    tmpW.compute(tmpParamValues,
                 NULL,
                 1.,
                 0.,
                 info.m_wMaxTimeStep,
                 false, // computeWAndLambdaGradsAlso
#ifdef QUESO_TGA_USES_OLD_COMPATIBLE_CODE
                 NULL,
#else
                 info.m_referenceW,
#endif
                 info.m_weightFunction,
                 &valueAm,
                 NULL);
    std::cout << "valueAm = " << valueAm << std::endl;

    tmpParamValues[0] = guessA+deltaA;
    tmpParamValues[1] = guessE;
    double valueAp = 0.;
    tmpW.compute(tmpParamValues,
                 NULL,
                 1.,
                 0.,
                 info.m_wMaxTimeStep,
                 false, // computeWAndLambdaGradsAlso
#ifdef QUESO_TGA_USES_OLD_COMPATIBLE_CODE
                 NULL,
#else
                 info.m_referenceW,
#endif
                 info.m_weightFunction,
                 &valueAp,
                 NULL);
    std::cout << "valueAp = " << valueAp << std::endl;

    tmpParamValues[0] = guessA;
    tmpParamValues[1] = guessE-deltaE;
    double valueEm = 0.;
    tmpW.compute(tmpParamValues,
                 NULL,
                 1.,
                 0.,
                 info.m_wMaxTimeStep,
                 false, // computeWAndLambdaGradsAlso
#ifdef QUESO_TGA_USES_OLD_COMPATIBLE_CODE
                 NULL,
#else
                 info.m_referenceW,
#endif
                 info.m_weightFunction,
                 &valueEm,
                 NULL);
    std::cout << "valueEm = " << valueEm << std::endl;

    tmpParamValues[0] = guessA;
    tmpParamValues[1] = guessE+deltaE;
    double valueEp = 0.;
    tmpW.compute(tmpParamValues,
                 NULL,
                 1.,
                 0.,
                 info.m_wMaxTimeStep,
                 false, // computeWAndLambdaGradsAlso
#ifdef QUESO_TGA_USES_OLD_COMPATIBLE_CODE
                 NULL,
#else
                 info.m_referenceW,
#endif
                 info.m_weightFunction,
                 &valueEp,
                 NULL);
    std::cout << "valueEp = " << valueEp << std::endl;

    P_V gradWithFD(info.m_paramSpace.zeroVector());
    gradWithFD[0] = (valueAp-valueAm)/2./deltaA;
    gradWithFD[1] = (valueEp-valueEm)/2./deltaE;
    if (paramValues.env().rank() == 0) {
      std::cout << "\nIn uqTgaLikelihoodCheckingRoutine()"
                << ": gradWithFD = " << gradWithFD
                << "\ngrad relative error = " << (gradWithFD - *tmpGradVector).norm2()/gradWithFD.norm2()
                << std::endl;
    }
  }

  //////////////////////////////////////////////////////////////////
  // Step 3 of 4: Compute Hessian(misfit), wrt A and E, using finite differences
  //////////////////////////////////////////////////////////////////

  if (tmpHessianMatrix) {
    P_V tmpParamValues(paramValues);
    bool savedCheckingFlag = info.m_performChecking;
    info.m_performChecking = false; // IMPORTANT

    std::vector<uqTgaLikelihoodInfoStruct<P_V,P_M>* > tmpCalLikelihoodInfoVector(1,(uqTgaLikelihoodInfoStruct<P_V,P_M>*) NULL);
    tmpCalLikelihoodInfoVector[0] = &info;

    if (paramValues.env().rank() == 0) {
      std::cout << "\nIn uqTgaLikelihoodCheckingRoutine()"
                << ": hessianWithLM = " << *tmpHessianMatrix
                << std::endl;
    }

    // Compute \grad(misfit), using finite differences
    if (paramValues.env().rank() == 0) {
      std::cout << "In uqTgaLikelihoodCheckingRoutine()"
                << ": computing Hessian(misfit) w.r.t. parameters A and E, using finite differences..."
                << std::endl;
    }

    tmpParamValues[0] = guessA-deltaA;
    tmpParamValues[1] = guessE;
    P_V grad_Am(info.m_paramSpace.zeroVector());
    double tmpMisfit = uqTgaLikelihoodRoutine<P_V,P_M>(tmpParamValues,
                                                       NULL,
                                                       (const void *)&tmpCalLikelihoodInfoVector,
                                                       &grad_Am,
                                                       NULL, // Hessian
                                                       NULL);
    //std::cout << "grad_Am = " << grad_Am << std::endl;

    tmpParamValues[0] = guessA+deltaA;
    tmpParamValues[1] = guessE;
    P_V grad_Ap(info.m_paramSpace.zeroVector());
    tmpMisfit = uqTgaLikelihoodRoutine<P_V,P_M>(tmpParamValues,
                                                NULL,
                                                (const void *)&tmpCalLikelihoodInfoVector,
                                                &grad_Ap,
                                                NULL, // Hessian
                                                NULL);
    //std::cout << "grad_Ap = " << grad_Ap << std::endl;

    tmpParamValues[0] = guessA;
    tmpParamValues[1] = guessE-deltaE;
    P_V grad_Em(info.m_paramSpace.zeroVector());
    tmpMisfit = uqTgaLikelihoodRoutine<P_V,P_M>(tmpParamValues,
                                                NULL,
                                                (const void *)&tmpCalLikelihoodInfoVector,
                                                &grad_Em,
                                                NULL, // Hessian
                                                NULL);
    //std::cout << "grad_Em = " << grad_Em << std::endl;

    tmpParamValues[0] = guessA;
    tmpParamValues[1] = guessE+deltaE;
    P_V grad_Ep(info.m_paramSpace.zeroVector());
    tmpMisfit = uqTgaLikelihoodRoutine<P_V,P_M>(tmpParamValues,
                                                NULL,
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
    if (paramValues.env().rank() == 0) {
      std::cout << "\nIn uqTgaLikelihoodCheckingRoutine()"
                << ": hessianWithFD = " << *hessianWithFD
                << "\nhessian absolute error = " << (*hessianWithFD + (-1.*(*tmpHessianMatrix)))
                << std::endl;
    }
    delete hessianWithFD;

    info.m_performChecking = savedCheckingFlag; // IMPORTANT
  }

  //////////////////////////////////////////////////////////////////
  // Step 4 of 4: compute hessianEffect, using Hessian matrix
  //////////////////////////////////////////////////////////////////

  if (tmpHessianEffect) {
    if (paramValues.env().rank() == 0) {
      std::cout << "\nIn uqTgaLikelihoodCheckingRoutine()"
                << ": for paramDirection = "  << *paramDirection
                << ", hessianEffectWithLM = " << *tmpHessianEffect
                << std::endl;
    }

    if (paramDirection && tmpHessianMatrix) {
      P_V effectWithMatrix((*tmpHessianMatrix)*(*paramDirection));
      if (paramValues.env().rank() == 0) {
        std::cout << "\nIn uqTgaLikelihoodCheckingRoutine()"
                  << ": for paramDirection = "      << *paramDirection
                  << ", hessianEffectWithMatrix = " << effectWithMatrix
                  << "\neffect relative error = "   << (effectWithMatrix - *tmpHessianEffect).norm2()/effectWithMatrix.norm2()
                  << std::endl;
      }
    }
  }

  return;
}
#endif // __UQ_TGA_LIKELIHOOD_H__
