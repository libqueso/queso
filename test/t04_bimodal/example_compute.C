//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011 The PECOS Development Team
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
// $Id:$
//
//--------------------------------------------------------------------------

#include <example_compute.h>
#include <example_likelihood.h>
#include <uqGslMatrix.h>
#include <uqStatisticalInverseProblem.h>
#include <uq1D1DFunction.h>

void compute(const uqFullEnvironmentClass& env) {
#if 0
  uqConstant1D1DFunctionClass w(-1.,1.,1.);

  uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass>
    paramSpace(env, "param_", 5, NULL);
  uqGslVectorClass v1(paramSpace.zeroVector());
  uqGslVectorClass v2(paramSpace.zeroVector());
  w.quadPtsWeigths<uqGslVectorClass,uqGslMatrixClass>(200,false,v1,v2);
  std::cout << "v1 = " << v1 << std::endl;
  std::cout << "v2 = " << v2 << std::endl;
#else
  ////////////////////////////////////////////////////////
  // Step 1 of 5: Instantiate the parameter space
  ////////////////////////////////////////////////////////
  uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass>
    paramSpace(env, "param_", 1, NULL);

  ////////////////////////////////////////////////////////
  // Step 2 of 5: Instantiate the parameter domain
  ////////////////////////////////////////////////////////
  uqGslVectorClass paramMins(paramSpace.zeroVector());
  paramMins.cwSet(-250.);
  uqGslVectorClass paramMaxs(paramSpace.zeroVector());
  paramMaxs.cwSet( 250.);
  uqBoxSubsetClass<uqGslVectorClass,uqGslMatrixClass>
    paramDomain("param_",paramSpace,paramMins,paramMaxs);

  ////////////////////////////////////////////////////////
  // Step 3 of 5: Instantiate the likelihood function object
  ////////////////////////////////////////////////////////
  uqGslVectorClass meanVector(paramSpace.zeroVector());
  meanVector[0] = 10.;
  uqGslMatrixClass* covMatrix = paramSpace.newMatrix();
  (*covMatrix)(0,0) = 1.;
  likelihoodRoutine_DataType likelihoodRoutine_Data;
  likelihoodRoutine_Data.meanVector = &meanVector;
  likelihoodRoutine_Data.covMatrix  = covMatrix;
  uqGenericScalarFunctionClass<uqGslVectorClass,uqGslMatrixClass>
    likelihoodFunctionObj("like_",
                          paramDomain,
                          likelihoodRoutine,
                          (void *) &likelihoodRoutine_Data,
                          true); // routine computes [-2.*ln(function)]

  ////////////////////////////////////////////////////////
  // Step 4 of 5: Instantiate the inverse problem
  ////////////////////////////////////////////////////////
  uqUniformVectorRVClass<uqGslVectorClass,uqGslMatrixClass>
    priorRv("prior_", paramDomain);
  uqGenericVectorRVClass<uqGslVectorClass,uqGslMatrixClass>
    postRv("post_", paramSpace);
  uqStatisticalInverseProblemClass<uqGslVectorClass,uqGslMatrixClass>
    ip("", NULL, priorRv, likelihoodFunctionObj, postRv);

  ////////////////////////////////////////////////////////
  // Step 5 of 5: Solve the inverse problem
  ////////////////////////////////////////////////////////
#if 0
  uqGslVectorClass paramInitials(paramSpace.zeroVector());
  paramInitials[0] = 45.;
  uqGslMatrixClass* proposalCovMatrix = paramSpace.newMatrix();
  (*proposalCovMatrix)(0,0) = 1600.;
  ip.solveWithBayesMetropolisHastings(NULL, paramInitials, proposalCovMatrix);
  delete proposalCovMatrix;
#else
  ip.solveWithBayesMLSampling();
#endif

  ////////////////////////////////////////////////////////
  // Print some statistics
  ////////////////////////////////////////////////////////
  unsigned int numPosTotal = postRv.realizer().subPeriod();
  if (env.subDisplayFile()) {
    *env.subDisplayFile() << "numPosTotal = " << numPosTotal
                          << std::endl;
  }

  uqGslVectorClass auxVec(paramSpace.zeroVector());
  unsigned int numPosSmallerThan40 = 0;
  for (unsigned int i = 0; i < numPosTotal; ++i) {
    postRv.realizer().realization(auxVec);
    if (auxVec[0] < 40.) numPosSmallerThan40++;
  }
  if (env.subDisplayFile()) {
    *env.subDisplayFile() << "numPosSmallerThan40 = " << numPosSmallerThan40
                          << ", ratio = " << ((double) numPosSmallerThan40)/((double) numPosTotal)
                          << std::endl;
  }

  uqScalarSequenceClass<double> seq1  (env,numPosSmallerThan40,"");
  uqScalarSequenceClass<double> seq2  (env,numPosTotal - numPosSmallerThan40,"");
  uqScalarSequenceClass<double> seqAll(env,numPosTotal,"");
  unsigned int i1 = 0;
  unsigned int i2 = 0;
  for (unsigned int i = 0; i < numPosTotal; ++i) {
    postRv.realizer().realization(auxVec);
    if (auxVec[0] < 40.) seq1[i1++] = auxVec[0];
    else                 seq2[i2++] = auxVec[0];
    seqAll[i] = auxVec[0];
  }

  double mean1 = seq1.subMean(0,seq1.subSequenceSize());
  if (env.subDisplayFile()) {
    *env.subDisplayFile() << "seq1.size() = "    << seq1.subSequenceSize()
                          << "\n seq1.mean() = " << mean1
                          << "\n seq1.std() = "  << sqrt(seq1.subSampleVariance(0,seq1.subSequenceSize(),mean1))
                          << std::endl;
  }

  double mean2 = seq2.subMean(0,seq2.subSequenceSize());
  if (env.subDisplayFile()) {
    *env.subDisplayFile() << "seq2.size() = "    << seq2.subSequenceSize()
                          << "\n seq2.mean() = " << mean2
                          << "\n seq2.std() = "  << sqrt(seq2.subSampleVariance(0,seq2.subSequenceSize(),mean2))
                          << std::endl;
  }

  double meanAll = seqAll.subMean(0,seqAll.subSequenceSize());
  if (env.subDisplayFile()) {
    *env.subDisplayFile() << "seqAll.size() = "    << seqAll.subSequenceSize()
                          << "\n seqAll.mean() = " << meanAll
                          << "\n seqAll.std() = "  << sqrt(seqAll.subSampleVariance(0,seqAll.subSequenceSize(),meanAll))
                          << std::endl;
  }

  ////////////////////////////////////////////////////////
  // Test if likelihood is normalized
  ////////////////////////////////////////////////////////
  unsigned int numGridPoints = 1000001;
  double xMin = paramDomain.minValues()[0];
  double xMax = paramDomain.maxValues()[0];
  double intervalSize = (xMax-xMin)/((double) numGridPoints - 1);
  double integral = 0.;
  for (unsigned int i = 0; i < numGridPoints; ++i) {
    auxVec[0] = xMin + i*intervalSize;
    integral += likelihoodFunctionObj.actualValue(auxVec,NULL,NULL,NULL,NULL);
  }
  integral *= intervalSize;
  if (env.subDisplayFile()) {
    *env.subDisplayFile() << "integral = " << integral
                          << std::endl;
  }

  // Return
  delete covMatrix;
#endif
  return;
}

#if 0
  // Test 1 of 2: finite distribution
  std::vector<double> weights(5,0.);
  weights[0] = 0.07;
  weights[1] = 0.14;
  weights[2] = 0.00;
  weights[3] = 0.67;
  weights[4] = 0.12;
  uqFiniteDistributionClass tmpFd(env,
                                  "",
                                  weights);

  std::vector<unsigned int> counts(5,0);
  unsigned int numSamples = 10000000;
  for (unsigned int i = 0; i < numSamples; ++i) {
    unsigned int index = tmpFd.sample();
    UQ_FATAL_TEST_MACRO(index >= 5,
                        env.fullRank(),
                        "compute() in example_compute.C",
                        "index sampled from finite distribtuion is too large");
    counts[index] += 1;
  }
  std::cout << "counts[0] = "    << ((double) counts[0])/((double) numSamples)
            << "\n counts[1] = " << ((double) counts[1])/((double) numSamples)
            << "\n counts[2] = " << ((double) counts[2])/((double) numSamples)
            << "\n counts[3] = " << ((double) counts[3])/((double) numSamples)
            << "\n counts[4] = " << ((double) counts[4])/((double) numSamples)
            << std::endl;
#endif
