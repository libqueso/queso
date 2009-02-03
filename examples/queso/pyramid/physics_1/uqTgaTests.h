/* uq/examples/queso/pyramid/uqTgaTests.h
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

#ifndef __UQ_TGA_TESTS_H__
#define __UQ_TGA_TESTS_H__

#include <uqTgaTemperatureProfiles.h>

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqTgaValidationClass<P_V,P_M,Q_V,Q_M>::runTests()
{
  m_testOptions.scanOptionsValues();

  // Read parameters for temperature function (which is linear wrt time)
  // Read discrete measurements (which might be replaced by continuous measurements)
  m_calLikelihoodInfoVector.resize(1,NULL);
  m_calLikelihoodInfoVector[0] = new uqTgaLikelihoodInfoStruct<P_V,P_M>(*m_paramSpace,
                                                                        "tgaData/scenario_10_K_min.dat",
                                                                        &m_options.m_wMaxTimeStep,
                                                                        &m_options.m_lambdaMaxTimeStep,
                                                                        &m_options.m_integralsNumIntervals);

  switch (m_options.m_refTemperatureProfileId) {
    case 0:
      // Do nothing: keep the linearly increasing profile set by experimental input file
    break;

    case 1:
      m_testVars.tempFunction = new uqGeneric1D1DFunctionClass(-INFINITY,
                                                               INFINITY,
                                                               temperatureFunction1_Value,
                                                               temperatureFunction1_Deriv,
                                                               NULL);
      m_calLikelihoodInfoVector[0]->changeTempFunction(*m_testVars.tempFunction);
    break;

    case 2:
      m_testVars.tempFunction = new uqGeneric1D1DFunctionClass(-INFINITY,
                                                               INFINITY,
                                                               temperatureFunction2_Value,
                                                               temperatureFunction2_Deriv,
                                                               NULL);
      m_calLikelihoodInfoVector[0]->changeTempFunction(*m_testVars.tempFunction);
    break;

    case 3:
      m_testVars.tempFunction = new uqGeneric1D1DFunctionClass(-INFINITY,
                                                               INFINITY,
                                                               temperatureFunction3_Value,
                                                               temperatureFunction3_Deriv,
                                                               NULL);
      m_calLikelihoodInfoVector[0]->changeTempFunction(*m_testVars.tempFunction);
    break;

    case 4:
      m_testVars.tempFunction = new uqGeneric1D1DFunctionClass(-INFINITY,
                                                               INFINITY,
                                                               temperatureFunction4_Value,
                                                               temperatureFunction4_Deriv,
                                                               NULL);
      m_calLikelihoodInfoVector[0]->changeTempFunction(*m_testVars.tempFunction);
    break;

    default:
      // Do nothing: keep the linearly increasing profile set by experimental input file
    break;
  }

  // Output file
  m_testVars.ofs = new std::ofstream( (m_testOptions.m_outerPrefixName+".m").c_str(), std::ofstream::out | std::ofstream::trunc);

  if (m_testOptions.m_runTempTimeTest) runTempTimeTest();
  if ((m_options.m_refCreate              ) ||
      (m_testOptions.m_runTimingTest      ) ||
      (m_testOptions.m_runGradTest        ) ||
      (m_testOptions.m_runOptimizationTest)) {
    fabricateReferenceData();
  }
  if (m_testOptions.m_runTimingTest      ) runTimingTest      ();
  if (m_testOptions.m_runGradTest        ) runGradTest        ();
  if (m_testOptions.m_runOptimizationTest) runOptimizationTest();

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqTgaValidationClass<P_V,P_M,Q_V,Q_M>::runTempTimeTest()
{
  P_V paramValues(m_paramSpace->zeroVector());
  std::string tmpPrefixName;

  paramValues[0] = 2.6000e+11;
  paramValues[1] = 2.0000e+05;

  uqTgaWClass<P_V,P_M> tmpW(*m_paramSpace,
                            *(m_calLikelihoodInfoVector[0]->m_temperatureFunctionObj));

#ifdef QUESO_TGA_USES_OLD_COMPATIBLE_CODE
  // Compute w using derivative wrt temperature
  if (m_env.rank() == 0) {
    std::cout << "In uqTgaValidation::runTempTimeTest()"
              << ": computing w using derivative wrt temperature..."
              << std::endl;
  }

  tmpW.computeUsingTemp(paramValues,
                        822., // maximumTemp
                        NULL, // referenceW
                        NULL);

  tmpPrefixName = m_testOptions.m_outerPrefixName + "withTempW_";
  tmpW.printForMatlab(*m_testVars.ofs,tmpPrefixName);
#endif

  // Compute w using derivative wrt time
  if (m_env.rank() == 0) {
    std::cout << "In uqTgaValidation::runTempTimeTest()"
              << ": computing w using derivative wrt time..."
              << std::endl;
  }

  tmpW.compute(paramValues,
               NULL,
               1.,                           // initialValue
               0.,                           // maximumTime
               m_options.m_refMaxTimeStep, // maximum delta time
               false,                        // computeWAndLambdaGradsAlso
               NULL,                         // referenceW
               NULL,                         // weightFunction
               NULL,                         // misfitValue
               NULL);                        // diffFunction

  tmpPrefixName = m_testOptions.m_outerPrefixName + "withTimeW_";
  tmpW.printForMatlab(*m_testVars.ofs,tmpPrefixName);

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqTgaValidationClass<P_V,P_M,Q_V,Q_M>::fabricateReferenceData()
{
  //////////////////////////////////////////////////////////////////
  // Compute data for refW
  //////////////////////////////////////////////////////////////////
  std::string tmpPrefixName;
  P_V paramValues(m_paramSpace->zeroVector());

  double factor1 = m_options.m_refW1;
  double factor2 = m_options.m_refW2;
  double factor3 = m_options.m_refW3;
  double factor4 = 1. - factor1 - factor2 - factor3;

  // tmpW1
  paramValues[0] = m_options.m_refA1;
  paramValues[1] = m_options.m_refE1;
  uqTgaWClass<P_V,P_M> tmpW1(*m_paramSpace,
                             *(m_calLikelihoodInfoVector[0]->m_temperatureFunctionObj));

  tmpW1.compute(paramValues,
                NULL,
                1., //m_options.m_refW1,
                m_options.m_refMaxTime,     // maximumTime
                m_options.m_refMaxTimeStep, // maximum delta time
                false,                        // computeWAndLambdaGradsAlso
                NULL,                         // referenceW
                NULL,                         // weightFunction
                NULL,                         // misfitValue
                NULL);                        // diffFunction

  tmpPrefixName = m_testOptions.m_outerPrefixName + "tmpW1_";
  tmpW1.printForMatlab(*m_testVars.ofs,tmpPrefixName);

  // tmpW2
  paramValues[0] = m_options.m_refA2;
  paramValues[1] = m_options.m_refE2;
  uqTgaWClass<P_V,P_M> tmpW2(*m_paramSpace,
                             *(m_calLikelihoodInfoVector[0]->m_temperatureFunctionObj));

  if (factor2 > 0.) {
    tmpW2.compute(paramValues,
                  NULL,
                  1., //m_options.m_refW2,
                  m_options.m_refMaxTime,     // maximumTime
                  m_options.m_refMaxTimeStep, // maximum delta time
                  false,                        // computeWAndLambdaGradsAlso
                  NULL,                         // referenceW
                  NULL,                         // weightFunction
                  NULL,                         // misfitValue
                  NULL);                        // diffFunction

    tmpPrefixName = m_testOptions.m_outerPrefixName + "tmpW2_";
    tmpW2.printForMatlab(*m_testVars.ofs,tmpPrefixName);
  }

  // tmpW3
  paramValues[0] = m_options.m_refA3;
  paramValues[1] = m_options.m_refE3;
  uqTgaWClass<P_V,P_M> tmpW3(*m_paramSpace,
                             *(m_calLikelihoodInfoVector[0]->m_temperatureFunctionObj));

  if (factor3 > 0.) {
    tmpW3.compute(paramValues,
                  NULL,
                  1., //m_options.m_refW3,
                  m_options.m_refMaxTime,     // maximumTime
                  m_options.m_refMaxTimeStep, // maximum delta time
                  false,                        // computeWAndLambdaGradsAlso
                  NULL,                         // referenceW
                  NULL,                         // weightFunction
                  NULL,                         // misfitValue
                  NULL);                        // diffFunction

    tmpPrefixName = m_testOptions.m_outerPrefixName + "tmpW3_";
    tmpW3.printForMatlab(*m_testVars.ofs,tmpPrefixName);
  }

  // tmpW4
  paramValues[0] = m_options.m_refA4;
  paramValues[1] = m_options.m_refE4;
  uqTgaWClass<P_V,P_M> tmpW4(*m_paramSpace,
                             *(m_calLikelihoodInfoVector[0]->m_temperatureFunctionObj));

  if (factor4 > 0.) {
    tmpW4.compute(paramValues,
                  NULL,
                  1., //1.-m_options.m_refW1-m_options.m_refW2-m_options.m_refW3,
                  m_options.m_refMaxTime,     // maximumTime
                  m_options.m_refMaxTimeStep, // maximum delta time
                  false,                        // computeWAndLambdaGradsAlso
                  NULL,                         // referenceW
                  NULL,                         // weightFunction
                  NULL,                         // misfitValue
                  NULL);                        // diffFunction

    tmpPrefixName = m_testOptions.m_outerPrefixName + "tmpW4_";
    tmpW4.printForMatlab(*m_testVars.ofs,tmpPrefixName);
  }

  // tmpW1 + tmpW2 + tmpW3 + tmpW4
  unsigned int wSize = tmpW1.times().size();

  unsigned int startingId2 = 0;
  unsigned int startingId3 = 0;
  unsigned int startingId4 = 0;

  std::vector<double> values2(wSize,0.);
  std::vector<double> values3(wSize,0.);
  std::vector<double> values4(wSize,0.);
  std::vector<double> wValues(wSize,0.);
  std::vector<double> tempValues(wSize,0.);
  std::vector<double> tempDerivs(wSize,0.);

  for (unsigned int i = 0; i < wSize; ++i) {
    double time = tmpW1.times()[i];

    values2[i] = 0.;
    if (factor2 > 0.) {
    if (time <= tmpW2.times()[tmpW2.ws().size()-1]) {
      tmpW2.interpolate(time,
                        startingId2,
                        1.,
                        &values2[i],
                        NULL,
                        NULL,
                        NULL);
    }
    }

    values3[i] = 0.;
    if (factor3 > 0.) {
    if (time <= tmpW3.times()[tmpW3.ws().size()-1]) {
      tmpW3.interpolate(time,
                        startingId3,
                        1.,
                        &values3[i],
                        NULL,
                        NULL,
                        NULL);
    }
    }

    values4[i] = 0.;
    if (factor4 > 0.) {
    if (time <= tmpW4.times()[tmpW4.ws().size()-1]) {
      tmpW4.interpolate(time,
                        startingId4,
                        1.,
                        &values4[i],
                        NULL,
                        NULL,
                        NULL);
    }
    }

    wValues[i] = factor1*tmpW1.ws()[i] + factor2*values2[i] + factor3*values3[i] + factor4*values4[i];

    tempValues[i] = m_calLikelihoodInfoVector[0]->m_temperatureFunctionObj->value(time);

    tempDerivs[i] = m_calLikelihoodInfoVector[0]->m_temperatureFunctionObj->deriv(time);
  }

  // reference W
  uqSampled1D1DFunctionClass tmpW(tmpW1.times(),
                                  wValues);
  tmpPrefixName = m_testOptions.m_outerPrefixName + "tmpW_";
  tmpW.printForMatlab(*m_testVars.ofs,tmpPrefixName);

  // temperature values
  uqTgaStorageClass<P_V,P_M> tmpTemp(tmpW1.times(),
#ifdef QUESO_TGA_USES_OLD_COMPATIBLE_CODE
                                     tempValues,
#endif
                                     tempValues,
                                     NULL,
                                     false);
  tmpPrefixName = m_testOptions.m_outerPrefixName + "tmpTemp_";
  tmpTemp.printForMatlab(*m_testVars.ofs,tmpPrefixName);

  // temperature derivatives
  uqTgaStorageClass<P_V,P_M> tmpDeriv(tmpW1.times(),
#ifdef QUESO_TGA_USES_OLD_COMPATIBLE_CODE
                                      tempValues,
#endif
                                      tempDerivs,
                                      NULL,
                                      false);
  tmpPrefixName = m_testOptions.m_outerPrefixName + "tmpDerivs_";
  tmpDeriv.printForMatlab(*m_testVars.ofs,tmpPrefixName);

  if (m_env.rank() == 0) {
    std::cout << "In uqTgaValidation::fabricateReferenceData()"
              << ": tmpW.times().size() = " << wSize
              << ", tmpW.times()[0] = "     << tmpW1.times()[0]
              << ", tmpW.values()[0] = "    << wValues      [0]
              << ", tmpW.times()[max] = "   << tmpW1.times()[wSize-1]
              << ", tmpW.values()[max] = "  << wValues      [wSize-1]
              << std::endl;
  }

  //////////////////////////////////////////////////////////////////
  // Set refW
  //////////////////////////////////////////////////////////////////

  //std::vector<double> constantTmpWs(tmpW.times().size(),.5);
  if (m_options.m_refTreatDataAsContinuous) {
    m_testVars.refW = new uqSampled1D1DFunctionClass(tmpW1.times(),
                                                     wValues); // wValues or constantTmpWs
  }
  else {
    unsigned int tmpSize = tmpW1.times().size();
    double minTime = tmpW1.times()[0];
    double maxTime = tmpW1.times()[tmpSize-1];
    double refDeltaTime = (maxTime-minTime)/((double) (m_options.m_refNumDiscreteSamples-1));

    std::vector<double> vecOfTimes (m_options.m_refNumDiscreteSamples,0.);
    std::vector<double> vecOfValues(m_options.m_refNumDiscreteSamples,0.);
    //unsigned int startingId = 0;
    for (unsigned int j = 0; j < m_options.m_refNumDiscreteSamples; ++j) {
      if (j == 0) {
        vecOfTimes [j] = minTime;
        vecOfValues[j] = wValues[0];
      }
      else if (j == (m_options.m_refNumDiscreteSamples - 1)) {
        vecOfTimes [j] = maxTime;
        vecOfValues[j] = wValues[tmpSize-1];
      }
      else {
        vecOfTimes [j] = minTime + ((double) j)*refDeltaTime;
#if 1
        vecOfValues[j] = tmpW.value(vecOfTimes[j]);
#else
        tmpW.interpolate(vecOfTimes[j],
                         startingId,
                         1.,
                         &vecOfValues[j],
                         NULL,
                         NULL,
                         NULL);
#endif
      }
    }

    m_testVars.refW = new uqSampled1D1DFunctionClass(vecOfTimes,
                                                     vecOfValues);
  }

  //////////////////////////////////////////////////////////////////
  // Prepare weight functions
  //////////////////////////////////////////////////////////////////
  unsigned int refSize = m_testVars.refW->domainValues().size();
  double       maxTime = m_testVars.refW->maxDomainValue();
  std::vector<double> vecOfDeltas (refSize,INFINITY);
  std::vector<double> vecOfWeights(refSize,1.);
  std::vector<double> aux1        (refSize,0.);
  std::vector<double> aux2        (refSize,0.);
  std::vector<double> aux3        (refSize,0.);
  for (unsigned int j = 0; j < refSize; ++j) {
    aux1[j] = maxTime-m_testVars.refW->domainValues()[refSize-1-j];
    aux2[j] = vecOfDeltas                            [refSize-1-j];
    aux3[j] = vecOfWeights                           [refSize-1-j];
  }

  if (m_options.m_refTreatDataAsContinuous) {
    m_testVars.continuousWeightFunction = new uqSampled1D1DFunctionClass(m_testVars.refW->domainValues(),
                                                                         vecOfWeights);
    m_testVars.tildeContinuousWeightFunction = new uqSampled1D1DFunctionClass(aux1,
                                                                              aux3);
  }
  else {
    m_testVars.deltaWeightFunction = new uqDeltaSeq1D1DFunctionClass(m_testVars.refW->domainValues(),
                                                                     vecOfDeltas,
                                                                     vecOfWeights);
    m_testVars.tildeDeltaWeightFunction = new uqDeltaSeq1D1DFunctionClass(aux1,
                                                                          aux2,
                                                                          aux3);
  }

  //////////////////////////////////////////////////////////////////
  // Set likelihood information appropriately
  //////////////////////////////////////////////////////////////////

  // Change the reference data to be used inside the likelihood routine,
  // from the (discrete) one read from file
  // to the (discrete or continuous) one computed and saved in 'refW'.
#ifdef QUESO_TGA_USES_OLD_COMPATIBLE_CODE
#else
   m_calLikelihoodInfoVector[0]->changeReferenceW(*m_testVars.refW);
#endif

  if (m_options.m_refTreatDataAsContinuous) {
    m_calLikelihoodInfoVector[0]->changeWeightFunctions(m_testVars.continuousWeightFunction,m_testVars.tildeContinuousWeightFunction); // or NULL,NULL
  }
  else {
    m_calLikelihoodInfoVector[0]->changeWeightFunctions(m_testVars.deltaWeightFunction,m_testVars.tildeDeltaWeightFunction); // or NULL,NULL
  }

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqTgaValidationClass<P_V,P_M,Q_V,Q_M>::runTimingTest()
{
  if (m_env.rank() == 0) {
    std::cout << "Entering uqTgaValidation::runTimingTest()..."
              << std::endl;
  }

  int iRC;
  struct timeval timeval0;
  struct timeval timeval1;

  P_V paramValues(m_paramSpace->zeroVector());
  paramValues[0] = m_testOptions.m_guessA;
  paramValues[1] = m_testOptions.m_guessE;

  // Run without checking, just in order to measure time
  double guessMisfit = 0.;

  iRC = gettimeofday(&timeval0, NULL);
  guessMisfit = uqTgaLikelihoodRoutine<P_V,P_M>(paramValues,
                                                NULL,
                                                (const void *)&m_calLikelihoodInfoVector,
                                                NULL,
                                                NULL,
                                                NULL);
  iRC = gettimeofday(&timeval1, NULL);
  double total_usecs = (timeval1.tv_sec * 1.e+6 + timeval1.tv_usec) - (timeval0.tv_sec * 1.e+6 + timeval0.tv_usec);
  if (m_env.rank() == 0) {
    std::cout << "In uqTgaValidation::runTimingTest()"
              << ": total_usecs (just misfit value) = " << total_usecs
              << std::endl;
  }

  // Run without checking, just in order to measure time
  P_V gradWithLM(m_paramSpace->zeroVector());
  P_M* hessianWithLM = NULL;
  if (m_testOptions.m_computeHessian) hessianWithLM = m_paramSpace->newMatrix();

  P_V paramDirection(m_paramSpace->zeroVector());
  paramDirection[0] =  0.5;
  paramDirection[1] = -1.5;
  P_V hessianEffect(m_paramSpace->zeroVector());

  iRC = gettimeofday(&timeval0, NULL);
  guessMisfit = uqTgaLikelihoodRoutine<P_V,P_M>(paramValues,
                                                &paramDirection,
                                                (const void *)&m_calLikelihoodInfoVector,
                                                &gradWithLM,
                                                hessianWithLM,
                                                &hessianEffect);
  iRC = gettimeofday(&timeval1, NULL);
  total_usecs = (timeval1.tv_sec * 1.e+6 + timeval1.tv_usec) - (timeval0.tv_sec * 1.e+6 + timeval0.tv_usec);
  if (m_env.rank() == 0) {
    std::cout << "In uqTgaValidation::runTimingTest()"
              << ": total_usecs (with grad and hessian) = " << total_usecs
              << std::endl;
  }
  delete hessianWithLM;

  if (m_env.rank() == 0) {
    std::cout << "Leaving uqTgaValidation::runTimingTest()"
              << std::endl;
  }

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqTgaValidationClass<P_V,P_M,Q_V,Q_M>::runGradTest()
{
  if (m_env.rank() == 0) {
    std::cout << "Entering uqTgaValidation::runGradTest()..."
              << std::endl;
  }

  P_V paramValues(m_paramSpace->zeroVector());
  paramValues[0] = m_testOptions.m_guessA;
  paramValues[1] = m_testOptions.m_guessE;

  P_V gradWithLM(m_paramSpace->zeroVector());
  P_M* hessianWithLM = NULL;
  if (m_testOptions.m_computeHessian) hessianWithLM = m_paramSpace->newMatrix();
  double guessMisfit = 0.;

  P_V paramDirection(m_paramSpace->zeroVector());
  paramDirection[0] =  0.5;
  paramDirection[1] = -1.5;
  P_V hessianEffect(m_paramSpace->zeroVector());

  // Run with checking against finite differences
  m_calLikelihoodInfoVector[0]->setCheckingVariables(true,&m_testOptions.m_relativeFDStep); // IMPORTANT
  guessMisfit = uqTgaLikelihoodRoutine<P_V,P_M>(paramValues,
                                                &paramDirection,
                                                (const void *)&m_calLikelihoodInfoVector,
                                                &gradWithLM,
                                                hessianWithLM,
                                                &hessianEffect);
  m_calLikelihoodInfoVector[0]->setCheckingVariables(false,&m_testOptions.m_relativeFDStep); // IMPORTANT
  delete hessianWithLM;
  if (m_env.rank() == 0) {
    std::cout << "In uqTgaValidation::runGradTest()"
              << ": guessMisfit = " << guessMisfit
              << std::endl;
  }

  if (m_testOptions.m_writeOutput) {
    std::string tmpPrefixName;

    *m_testVars.ofs << "\n" << m_testOptions.m_outerPrefixName << "hessianIsAvailable = " << m_testOptions.m_computeHessian << ";" << std::endl;
    *m_testVars.ofs << "\n" << m_testOptions.m_outerPrefixName << "refA1 = "              << m_options.m_refA1          << ";" << std::endl;
    *m_testVars.ofs << "\n" << m_testOptions.m_outerPrefixName << "refE1 = "              << m_options.m_refE1          << ";" << std::endl;
    *m_testVars.ofs << "\n" << m_testOptions.m_outerPrefixName << "guessA = "             << m_testOptions.m_guessA         << ";" << std::endl;
    *m_testVars.ofs << "\n" << m_testOptions.m_outerPrefixName << "guessE = "             << m_testOptions.m_guessE         << ";" << std::endl;

    tmpPrefixName = m_testOptions.m_outerPrefixName + "refW_";
    m_calLikelihoodInfoVector[0]->m_refWPtrForPlot->printForMatlab(*m_testVars.ofs,tmpPrefixName);

    tmpPrefixName = m_testOptions.m_outerPrefixName + "w_";
    m_calLikelihoodInfoVector[0]->m_wPtrForPlot->printForMatlab(*m_testVars.ofs,tmpPrefixName);

    if (m_testOptions.m_computeHessian) {
      tmpPrefixName = m_testOptions.m_outerPrefixName + "wA_";
      m_calLikelihoodInfoVector[0]->m_wAPtrForPlot->printForMatlab(*m_testVars.ofs,tmpPrefixName);

      tmpPrefixName = m_testOptions.m_outerPrefixName + "wE_";
      m_calLikelihoodInfoVector[0]->m_wEPtrForPlot->printForMatlab(*m_testVars.ofs,tmpPrefixName);
    }

    tmpPrefixName = m_testOptions.m_outerPrefixName + "diff_";
    m_calLikelihoodInfoVector[0]->m_diffPtrForPlot->printForMatlab(*m_testVars.ofs,tmpPrefixName);

    tmpPrefixName = m_testOptions.m_outerPrefixName + "lambda_";
    m_calLikelihoodInfoVector[0]->m_lambdaPtrForPlot->printForMatlab(*m_testVars.ofs,tmpPrefixName);

    if (m_testOptions.m_computeHessian) {
      tmpPrefixName = m_testOptions.m_outerPrefixName + "lambdaA_";
      m_calLikelihoodInfoVector[0]->m_lambdaAPtrForPlot->printForMatlab(*m_testVars.ofs,tmpPrefixName);

      tmpPrefixName = m_testOptions.m_outerPrefixName + "lambdaE_";
      m_calLikelihoodInfoVector[0]->m_lambdaEPtrForPlot->printForMatlab(*m_testVars.ofs,tmpPrefixName);
    }
  }

  if (m_env.rank() == 0) {
    std::cout << "Leaving uqTgaValidation::runGradTest()"
              << std::endl;
  }

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqTgaValidationClass<P_V,P_M,Q_V,Q_M>::runOptimizationTest()
{
  if (m_env.rank() == 0) {
    std::cout << "Entering uqTgaValidation::runOptimizationTest()"
              << std::endl;
  }

  P_V paramValues(m_paramSpace->zeroVector());
  paramValues[0] = m_testOptions.m_guessA;
  paramValues[1] = m_testOptions.m_guessE;

  std::string tmpPrefixName;
#if 0
  // initialW
  uqTgaWClass<P_V,P_M> initialW(*m_paramSpace,
                                *(m_calLikelihoodInfoVector[0]->m_temperatureFunctionObj));

  initialW.compute(paramValues,
                   NULL,
                   1., //m_options.m_refW1,
                   0.,                           // maximumTime
                   m_options.m_refMaxTimeStep, // maximum delta time
                   false,                        // computeWAndLambdaGradsAlso
                   NULL,                         // referenceW
                   NULL,                         // weightFunction
                   NULL,                         // misfitValue
                   NULL);                        // diffFunction

  tmpPrefixName = m_testOptions.m_outerPrefixName + "initialW_";
  initialW.printForMatlab(*m_testVars.ofs,tmpPrefixName);
#endif
  double objFunctionValue = 0.;
  P_V    objFunctionGradient(m_paramSpace->zeroVector());
  P_M*   objFunctionHessian = m_paramSpace->newMatrix();
  P_V    paramsStep(m_paramSpace->zeroVector());

  // Apply Newton method
  bool newtonSucceeded = false;
  int newtonFailureReason = 1;
  unsigned int newtonLoopId = 0;
  while (newtonLoopId < m_testOptions.m_NewtonMaxIters) {
    if (m_env.rank() == 0) {
      char stringA[64];
      char stringE[64];
      sprintf(stringA,"%12.6e",paramValues[0]);
      sprintf(stringE,"%12.6e",paramValues[1]);
      std::cout << "Beggining newtonLoopId = " << newtonLoopId
                << " with params = "           << stringA << " " << stringE
                << std::endl;
    }

    // Compute objective function: value, gradient and Hessian
    objFunctionValue = uqTgaLikelihoodRoutine<P_V,P_M>(paramValues,
                                                       NULL,
                                                       (const void *)&m_calLikelihoodInfoVector,
                                                       &objFunctionGradient,
                                                       objFunctionHessian,
                                                       NULL);
    double gradNorm = objFunctionGradient.norm2();

    if (m_env.rank() == 0) {
      char stringA[64];
      char stringE[64];
      sprintf(stringA,"%12.6e",paramValues[0]);
      sprintf(stringE,"%12.6e",paramValues[1]);
      std::cout << "In_newtonLoopId = "    << newtonLoopId
                << ", with params = "      << stringA << " " << stringE
                << ": objFunctionValue = " << objFunctionValue
                << ", objFunctionGrad = "  << objFunctionGradient
                << ", gradNorm = "         << gradNorm
                << std::endl;
    }

    // Check sufficiently small gradient norm
    if (gradNorm < m_testOptions.m_NewtonAbsTol) {
      newtonSucceeded = true;
      break;
    }

    double a = (*objFunctionHessian)(0,0);
    double b = (*objFunctionHessian)(0,1);
    double c = (*objFunctionHessian)(1,0);
    double d = (*objFunctionHessian)(1,1);
    double determinant = a*d - b*c;
    double frobNorm = sqrt(a*a + b*b + c*c + d*d);
    double e2 = .5*(a + d + sqrt( (a-d)*(a-d) + 4*b*c ));
    double e1 = .5*(a + d - sqrt( (a-d)*(a-d) + 4*b*c ));

    if (m_env.rank() == 0) {
      char stringA[64];
      char stringE[64];
      sprintf(stringA,"%12.6e",paramValues[0]);
      sprintf(stringE,"%12.6e",paramValues[1]);
      std::cout << "In newtonLoopId = "    << newtonLoopId
                << ", with params = "      << stringA << " " << stringE
                << ", Hessian = \n"        << *objFunctionHessian
                << ", determinant = "      << determinant
                << ", frobNorm = "         << frobNorm
                << ", e2 = "               << e2
                << ", e1 = "               << e1
                << ", e2*e1 = "            << e2*e1
                << std::endl;
    }

    unsigned int tauId = 0;
    double tau = 0.;
    while (determinant <= 0.) {
      tauId++;
      tau += frobNorm;
      double newA = a + tau;
      double newD = d + tau;
      determinant = newA*newD - b*c;
    }

    if (tauId > 0) {
      (*objFunctionHessian)(0,0) += tau;
      (*objFunctionHessian)(1,1) += tau;

      a = (*objFunctionHessian)(0,0);
      b = (*objFunctionHessian)(0,1);
      c = (*objFunctionHessian)(1,0);
      d = (*objFunctionHessian)(1,1);
      determinant = a*d - b*c;
      frobNorm = sqrt(a*a + b*b + c*c + d*d);
      e2 = .5*(a + d + sqrt( (a-d)*(a-d) + 4*b*c ));
      e1 = .5*(a + d - sqrt( (a-d)*(a-d) + 4*b*c ));

      if (m_env.rank() == 0) {
        char stringA[64];
        char stringE[64];
        sprintf(stringA,"%12.6e",paramValues[0]);
        sprintf(stringE,"%12.6e",paramValues[1]);
        std::cout << "In newtonLoopId = "   << newtonLoopId
                  << ", with params = "     << stringA << " " << stringE
                  << ": tauId = "           << tauId
                  << ", tau = "             << tau
                  << ", new Hessian = \n"   << *objFunctionHessian
                  << ", new determinant = " << determinant
                  << ", new frobNorm = "    << frobNorm
                  << ", new e2 = "          << e2
                  << ", new e1 = "          << e1
                  << ", e2*e1 = "           << e2*e1
                  << std::endl;
      }
    }

    // Compute Newton step
    objFunctionHessian->invertMultiply(-1.*objFunctionGradient,paramsStep);
    double directionalDerivative = scalarProduct(objFunctionGradient,paramsStep);
    double paramsStepNorm = paramsStep.norm2();

    if (m_env.rank() == 0) {
      char stringA[64];
      char stringE[64];
      sprintf(stringA,"%12.6e",paramValues[0]);
      sprintf(stringE,"%12.6e",paramValues[1]);
      std::cout << "In newtonLoopId = "         << newtonLoopId
                << " with params = "            << stringA << " " << stringE
                << ": paramsStep = "            << paramsStep
                << ", directionalDerivative = " << directionalDerivative
                << std::endl;
    }

    // Perform line search
    bool lineSearchSucceeded = false;
    unsigned int lineSearchLoopId = 0;
    double alphaFactor = .5;

    double alpha0 = 1.;
    if ((paramValues[0] + paramsStep[0]) <= 0.) alpha0 = -.5*paramValues[0]/paramsStep[0]; 
    double alpha1 = 1.;
    if ((paramValues[1] + paramsStep[1]) <= 0.) alpha1 = -.5*paramValues[1]/paramsStep[1]; 

    double alpha = std::min(alpha0,alpha1);
    alpha /= alphaFactor;
    while (lineSearchLoopId < 30) {
      double c1 = 1.e-4;
      //double c2 = .9;
      P_V attemptedParamValues(m_paramSpace->zeroVector());

      alpha *= alphaFactor;
      double lineSearchThreshold = objFunctionValue + alpha * c1 * directionalDerivative;
      attemptedParamValues = paramValues + alpha*paramsStep;

      if (m_env.rank() == 0) {
        char stringA[64];
        char stringE[64];
        sprintf(stringA,"%12.6e",paramValues[0]);
        sprintf(stringE,"%12.6e",paramValues[1]);
        std::cout << "In newtonLoopId = "        << newtonLoopId
                  << " with params = "           << stringA << " " << stringE
                  << ", alpha = "                << alpha
                  << ": attemptedParamValues = " << attemptedParamValues
                  << std::endl;
      }

      double attemptedObjFunctionValue = uqTgaLikelihoodRoutine<P_V,P_M>(attemptedParamValues,
                                                                         NULL,
                                                                         (const void *)&m_calLikelihoodInfoVector,
                                                                         NULL,
                                                                         NULL,
                                                                         NULL);
      if (m_env.rank() == 0) {
        char stringA[64];
        char stringE[64];
        sprintf(stringA,"%12.6e",paramValues[0]);
        sprintf(stringE,"%12.6e",paramValues[1]);
        std::cout << "In newtonLoopId = "             << newtonLoopId
                  << " with params = "                << stringA << " " << stringE
                  << ", alpha = "                     << alpha
                  << ": attemptedObjFunctionValue = " << attemptedObjFunctionValue
                  << ", lineSearchThreshold = "       << lineSearchThreshold
                  << std::endl;
      }

      // Check sufficient decrease
      if (attemptedObjFunctionValue <= lineSearchThreshold) {
        lineSearchSucceeded = true;
        // Update position
        paramValues = attemptedParamValues;

        if (m_env.rank() == 0) {
          char stringA[64];
          char stringE[64];
          sprintf(stringA,"%12.6e",paramValues[0]);
          sprintf(stringE,"%12.6e",paramValues[1]);
          std::cout << "In_newtonLoopId = "    << newtonLoopId
                    << " with params = "       << stringA << " " << stringE
                    << ": paramsStep = "       << paramsStep
                    << ", paramsStepNorm = "   << paramsStepNorm
                    << ", successful alpha = " << alpha
                    << std::endl;
        }
        break;
      }

      lineSearchLoopId++;
    }

    if (!lineSearchSucceeded) {
      newtonFailureReason = 2;
      break;
    }

    newtonLoopId++;
  }
  delete objFunctionHessian;

  if (m_env.rank() == 0) {
    std::cout << "In uqTgaValidation::runOptimizationTest():";
    if (newtonSucceeded) {
      char stringA[64];
      char stringE[64];
      sprintf(stringA,"%12.6e",paramValues[0]);
      sprintf(stringE,"%12.6e",paramValues[1]);
      std::cout << " solution = " << stringA << " " << stringE;
    }
    else {
      std::cout << " Newton failed: reason = " << newtonFailureReason;
    }
    std::cout << std::endl;
  }

  if (newtonSucceeded) {
    // optW
    uqTgaWClass<P_V,P_M> optW(*m_paramSpace,
                              *(m_calLikelihoodInfoVector[0]->m_temperatureFunctionObj));

    optW.compute(paramValues,
                 NULL,
                 1., //m_options.m_refW1,
                 m_options.m_refMaxTime,     // maximumTime
                 m_options.m_refMaxTimeStep, // maximum delta time
                 false,                        // computeWAndLambdaGradsAlso
                 NULL,                         // referenceW
                 NULL,                         // weightFunction
                 NULL,                         // misfitValue
                 NULL);                        // diffFunction

    tmpPrefixName = m_testOptions.m_outerPrefixName + "optW_";
    optW.printForMatlab(*m_testVars.ofs,tmpPrefixName);

    // optLambda
    uqTgaLambdaClass<P_V,P_M> optLambda(*m_paramSpace,
                                        *(m_calLikelihoodInfoVector[0]->m_temperatureFunctionObj));

    optLambda.compute(paramValues,
                      NULL, // paramDirection
                      m_calLikelihoodInfoVector[0]->m_lambdaMaxTimeStep,
                      false, // computeWAndLambdaGradsAlso
                      *(m_calLikelihoodInfoVector[0]->m_diffFunction),
                      m_calLikelihoodInfoVector[0]->m_tildeWeightFunction,
                      optW);

    tmpPrefixName = m_testOptions.m_outerPrefixName + "optLambda_";
    optLambda.printForMatlab(*m_testVars.ofs,tmpPrefixName);
  }

  if (m_env.rank() == 0) {
    std::cout << "Leaving uqTgaValidation::runOptimizationTest()"
              << std::endl;
  }

  return;
}

#endif // __UQ_TGA_TESTS_H__

