//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2015 The PECOS Development Team
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

#include <bimodal_compute.h>
#include <bimodal_likelihood.h>
#include <queso/GslMatrix.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/1D1DFunction.h>
#include <queso/GenericScalarFunction.h>
#include <queso/GenericVectorRV.h>
#include <queso/UniformVectorRV.h>

void compute(const QUESO::FullEnvironment& env) {
  //------------------------------------------------------
  // Step 1 of 5: Instantiate the parameter space
  //------------------------------------------------------
  QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix>
    paramSpace(env, "param_", 1, NULL);

  //------------------------------------------------------
  // Step 2 of 5: Instantiate the parameter domain
  //------------------------------------------------------
  QUESO::GslVector paramMins(paramSpace.zeroVector());
  paramMins.cwSet(-250.);
  QUESO::GslVector paramMaxs(paramSpace.zeroVector());
  paramMaxs.cwSet( 250.);
  QUESO::BoxSubset<QUESO::GslVector,QUESO::GslMatrix>
    paramDomain("param_",paramSpace,paramMins,paramMaxs);

  //------------------------------------------------------
  // Step 3 of 5: Instantiate the likelihood function object
  //------------------------------------------------------
  QUESO::GslVector meanVector(paramSpace.zeroVector());
  meanVector[0] = 10.;
  QUESO::GslMatrix* covMatrix = paramSpace.newMatrix();
  (*covMatrix)(0,0) = 1.;
  likelihoodRoutine_DataType likelihoodRoutine_Data;
  likelihoodRoutine_Data.meanVector = &meanVector;
  likelihoodRoutine_Data.covMatrix  = covMatrix;
  QUESO::GenericScalarFunction<QUESO::GslVector,QUESO::GslMatrix>
    likelihoodFunctionObj("like_",
                          paramDomain,
                          likelihoodRoutine,
                          (void *) &likelihoodRoutine_Data,
                          true); // routine computes [-2.*ln(function)]

  //------------------------------------------------------
  // Step 4 of 5: Instantiate the inverse problem
  //------------------------------------------------------
  QUESO::UniformVectorRV<QUESO::GslVector,QUESO::GslMatrix>
    priorRv("prior_", paramDomain);
  QUESO::GenericVectorRV<QUESO::GslVector,QUESO::GslMatrix>
    postRv("post_", paramSpace);
  QUESO::StatisticalInverseProblem<QUESO::GslVector,QUESO::GslMatrix>
    ip("", NULL, priorRv, likelihoodFunctionObj, postRv);

  //------------------------------------------------------
  // Step 5 of 5: Solve the inverse problem
  //------------------------------------------------------
#if 0
  uqGslVector paramInitials(paramSpace.zeroVector());
  paramInitials[0] = 45.;
  uqGslMatrix* proposalCovMatrix = paramSpace.newMatrix();
  (*proposalCovMatrix)(0,0) = 1600.;
  ip.solveWithBayesMetropolisHastings(NULL, paramInitials, proposalCovMatrix);
  delete proposalCovMatrix;
#else
  ip.solveWithBayesMLSampling();
#endif

  //------------------------------------------------------
  // Print some statistics in the file 'display_sub0.txt'
  // They will be in the last part of the file.
  //------------------------------------------------------
  unsigned int numPosTotal = postRv.realizer().subPeriod();
  if (env.subDisplayFile()) {
    *env.subDisplayFile() << "numPosTotal = " << numPosTotal
                          << std::endl;
  }

  QUESO::GslVector auxVec(paramSpace.zeroVector());
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

  QUESO::ScalarSequence<double> seq1  (env,numPosSmallerThan40,"");
  QUESO::ScalarSequence<double> seq2  (env,numPosTotal - numPosSmallerThan40,"");
  QUESO::ScalarSequence<double> seqAll(env,numPosTotal,"");
  unsigned int i1 = 0;
  unsigned int i2 = 0;
  for (unsigned int i = 0; i < numPosTotal; ++i) {
    postRv.realizer().realization(auxVec);
    if (auxVec[0] < 40.) seq1[i1++] = auxVec[0];
    else                 seq2[i2++] = auxVec[0];
    seqAll[i] = auxVec[0];
  }

  double mean1 = seq1.subMeanExtra(0,seq1.subSequenceSize());
  if (env.subDisplayFile()) {
    *env.subDisplayFile() << "seq1.size() = "    << seq1.subSequenceSize()
                          << "\n seq1.mean() = " << mean1
                          << "\n seq1.std() = "  << sqrt(seq1.subSampleVarianceExtra(0,seq1.subSequenceSize(),mean1))
                          << std::endl;
  }

  double mean2 = seq2.subMeanExtra(0,seq2.subSequenceSize());
  if (env.subDisplayFile()) {
    *env.subDisplayFile() << "seq2.size() = "    << seq2.subSequenceSize()
                          << "\n seq2.mean() = " << mean2
                          << "\n seq2.std() = "  << sqrt(seq2.subSampleVarianceExtra(0,seq2.subSequenceSize(),mean2))
                          << std::endl;
  }

  double meanAll = seqAll.subMeanExtra(0,seqAll.subSequenceSize());
  if (env.subDisplayFile()) {
    *env.subDisplayFile() << "seqAll.size() = "    << seqAll.subSequenceSize()
                          << "\n seqAll.mean() = " << meanAll
                          << "\n seqAll.std() = "  << sqrt(seqAll.subSampleVarianceExtra(0,seqAll.subSequenceSize(),meanAll))
                          << std::endl;
  }

  //------------------------------------------------------
  // Test if likelihood is normalized
  //------------------------------------------------------
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

  return;
}
