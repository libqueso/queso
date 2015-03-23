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

#include <example_compute.h>
#include <example_likelihood.h>
#include <queso/GslMatrix.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/ConcatenationSubset.h>
#include <queso/ConcatenatedVectorRV.h>
#include <queso/GenericScalarFunction.h>
#include <queso/UniformVectorRV.h>
#include <sys/time.h>
#include <example_hyst.h>

void compute(const QUESO::FullEnvironment& env) {

  struct timeval timevalNow;
  gettimeofday(&timevalNow, NULL);
  std::cout << std::endl << "Beginning run of 'Hysteretic' example at "
            << ctime(&timevalNow.tv_sec);

  //------------------------------------------------------
  // Step 1 of 5: Instantiate the parameter space
  //------------------------------------------------------
  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix>
    paramSpaceA(env, "paramA_", 1, NULL);
  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix>
    paramSpaceB(env, "paramB_", 14, NULL);
  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix>
    paramSpace (env, "param_", 15, NULL);

  //------------------------------------------------------
  // Step 2 of 5: Instantiate the parameter domain
  //------------------------------------------------------
  QUESO::GslVector paramMinsA(paramSpaceA.zeroVector());
  paramMinsA.cwSet(0);
  QUESO::GslVector paramMaxsA(paramSpaceA.zeroVector());
  paramMaxsA.cwSet(5);
  QUESO::BoxSubset<QUESO::GslVector,QUESO::GslMatrix>
    paramDomainA("paramA_",paramSpaceA,paramMinsA,paramMaxsA);

  QUESO::GslVector paramMinsB(paramSpaceB.zeroVector());
  paramMinsB.cwSet(-INFINITY);
  QUESO::GslVector paramMaxsB(paramSpaceB.zeroVector());
  paramMaxsB.cwSet( INFINITY);
  QUESO::BoxSubset<QUESO::GslVector,QUESO::GslMatrix>
    paramDomainB("paramB_",paramSpaceB,paramMinsB,paramMaxsB);

  QUESO::ConcatenationSubset<QUESO::GslVector,QUESO::GslMatrix>
    paramDomain("",paramSpace,paramDomainA,paramDomainB);

  //------------------------------------------------------
  // Step 3 of 5: Instantiate the likelihood function object
  //------------------------------------------------------
  std::cout << "\tInstantiating the Likelihood; calling internally the hysteretic model"
	    << std::endl;

  likelihoodRoutine_DataType likelihoodRoutine_Data;
  likelihoodRoutine_Data.floor.resize(4,NULL);
  unsigned int numTimeSteps = 401;
  for (unsigned int i = 0; i < 4; ++i) {
    likelihoodRoutine_Data.floor[i] = new std::vector<double>(numTimeSteps,0.);
  }
  likelihoodRoutine_Data.accel.resize(numTimeSteps,0.);
  FILE *inp;
  inp = fopen("an.txt","r");
  unsigned int numObservations = 0;
  double tmpA;
  while (fscanf(inp,"%lf",&tmpA) != EOF) {
    likelihoodRoutine_Data.accel[numObservations] = tmpA;
    numObservations++;
  }

  numObservations=0;
  FILE *inp1_1;
  inp1_1=fopen("measured_data1_1.txt","r");
  while (fscanf(inp1_1,"%lf",&tmpA) != EOF) {
    (*likelihoodRoutine_Data.floor[0])[numObservations]=tmpA;
     numObservations++;
  }

  numObservations=0;
  FILE *inp1_2;
  inp1_2=fopen("measured_data1_2.txt","r");
  while (fscanf(inp1_2,"%lf",&tmpA) != EOF) {
    (*likelihoodRoutine_Data.floor[1])[numObservations]=tmpA;
    numObservations++;
  }

  numObservations=0;
  FILE *inp1_3;
  inp1_3=fopen("measured_data1_3.txt","r");
  while (fscanf(inp1_3,"%lf",&tmpA) != EOF) {
    (*likelihoodRoutine_Data.floor[2])[numObservations]=tmpA;
    numObservations++;
  }

  numObservations=0;
  FILE *inp1_4;
  inp1_4=fopen("measured_data1_4.txt","r");
  while (fscanf(inp1_4,"%lf",&tmpA) != EOF) {
    (*likelihoodRoutine_Data.floor[3])[numObservations]=tmpA;
    numObservations++;
  }

  QUESO::GenericScalarFunction<QUESO::GslVector,QUESO::GslMatrix>
    likelihoodFunctionObj("like_",
                          paramDomain,
                          likelihoodRoutine,
                          (void *) &likelihoodRoutine_Data,
                          true); // routine computes [ln(function)]

  //------------------------------------------------------
  // Step 4 of 5: Instantiate the inverse problem
  //------------------------------------------------------
  std::cout << "\tInstantiating the SIP" << std::endl;

  QUESO::UniformVectorRV<QUESO::GslVector,QUESO::GslMatrix>
    priorRvA("priorA_", paramDomainA);

  QUESO::GslVector meanVec(paramSpaceB.zeroVector());
  QUESO::GslVector diagVec(paramSpaceB.zeroVector());

  diagVec.cwSet(0.6*0.6);

  QUESO::GslMatrix covMatrix(diagVec);

  QUESO::GaussianVectorRV<QUESO::GslVector,QUESO::GslMatrix>
    priorRvB("priorB_", paramDomainB,meanVec,covMatrix);

  QUESO::ConcatenatedVectorRV<QUESO::GslVector,QUESO::GslMatrix>
    priorRv("prior_", priorRvA, priorRvB, paramDomain);

  QUESO::GenericVectorRV<QUESO::GslVector,QUESO::GslMatrix>
    postRv("post_", paramSpace);

  QUESO::StatisticalInverseProblem<QUESO::GslVector,QUESO::GslMatrix>
    ip("", NULL, priorRv, likelihoodFunctionObj, postRv);

  //------------------------------------------------------
  // Step 5 of 5: Solve the inverse problem
  //------------------------------------------------------
  std::cout << "\tSolving the SIP with Multilevel method" << std::endl;

  ip.solveWithBayesMLSampling();

  gettimeofday(&timevalNow, NULL);
  std::cout << "Ending run of 'Hysteretic' example at "
              << ctime(&timevalNow.tv_sec) << std::endl;
  return;
}

//------------------------------------------------------
//------------------------------------------------------
//------------------------------------------------------
void debug_hyst(const QUESO::FullEnvironment& env) {
  unsigned int numFloors = 4;
  unsigned int numTimeSteps = 401;

  std::vector<double> accel(numTimeSteps,0.);
  FILE *inp;
  inp = fopen("an.txt","r");
  unsigned int numObservations = 0;
  double tmpA;
  while (fscanf(inp,"%lf",&tmpA) != EOF) {
    UQ_FATAL_TEST_MACRO((numObservations >= accel.size()),
                        env.fullRank(),
                        "debug_hyst()",
                        "input file has too many lines");
    accel[numObservations] = tmpA;
    numObservations++;
  }
  UQ_FATAL_TEST_MACRO((numObservations != accel.size()),
                      env.fullRank(),
                      "debug_hyst()",
                      "input file has a smaller number of observations than expected");

  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> floorSpace(env, "floor_", numFloors, NULL);

  QUESO::GslVector kVec(floorSpace.zeroVector());
  kVec[0] = 2.20e+7;
  kVec[1] = 2.00e+7;
  kVec[2] = 1.70e+7;
  kVec[3] = 1.45e+7;

  QUESO::GslVector rVec(floorSpace.zeroVector());
  rVec[0] = 0.1;
  rVec[1] = 0.1;
  rVec[2] = 0.1;
  rVec[3] = 0.1;

  QUESO::GslVector uVec(floorSpace.zeroVector());
  uVec[0] = 0.008;
  uVec[1] = 0.008;
  uVec[2] = 0.007;
  uVec[3] = 0.007;

  double rho   = 7.959e-1 ;//0.1976;
  double gamma = 2.500e-3 ; //0.0038;

  std::vector<double> t(numTimeSteps,0.);
  QUESO::SequenceOfVectors<QUESO::GslVector,QUESO::GslMatrix> u     (floorSpace,numTimeSteps,""); // absolute displacement
  QUESO::SequenceOfVectors<QUESO::GslVector,QUESO::GslMatrix> ud    (floorSpace,numTimeSteps,""); // velocity
  QUESO::SequenceOfVectors<QUESO::GslVector,QUESO::GslMatrix> udd   (floorSpace,numTimeSteps,""); // acceleration
  QUESO::SequenceOfVectors<QUESO::GslVector,QUESO::GslMatrix> resfor(floorSpace,numTimeSteps,""); // restoring force
  QUESO::SequenceOfVectors<QUESO::GslVector,QUESO::GslMatrix> ru    (floorSpace,numTimeSteps,""); // relative displacement

  u.setPositionValues     (0,floorSpace.zeroVector());
  ud.setPositionValues    (0,floorSpace.zeroVector());
  udd.setPositionValues   (0,floorSpace.zeroVector());
  resfor.setPositionValues(0,floorSpace.zeroVector());
  ru.setPositionValues    (0,floorSpace.zeroVector());

  QUESO::GslVector massVec(floorSpace.zeroVector());
  massVec.cwSet(2.0e+4);

  hystereticModel(env,
                  massVec,
                  kVec,
                  rVec,
                  uVec,
                  rho,
                  gamma,
                  accel,
                  t, // output
                  u,
                  ud,
                  udd,
                  resfor,
                  ru);

  std::set<unsigned int> auxSet;
  auxSet.insert(0);

  // Writing some data to the file 'outputData/cpp_output.m'
  std::ofstream myFile;
  myFile.open ("outputData/cpp_output.m");

  // Write 't_cpp'
  myFile << "t_cpp = zeros(" << 1 << "," << numTimeSteps << ");\n"
          << "t_cpp = [";
  for (unsigned int j = 0; j < numTimeSteps; ++j) {
    myFile << t[j] << " ";
  }
  myFile << "];" << std::endl;

  // Write 'a_cpp'
  myFile << "a_cpp = zeros(" << 1 << "," << numTimeSteps << ");\n"
          << "a_cpp = [";
  for (unsigned int j = 0; j < numTimeSteps; ++j) {
    myFile << accel[j] << " ";
  }
  myFile << "];" << std::endl;

  QUESO::GslVector auxVec(floorSpace.zeroVector());

  // Write 'u_cpp'
  myFile << "u_cpp = zeros(" << numFloors << "," << numTimeSteps << ");\n"
          << "u_cpp = [";
  for (unsigned int i = 0; i < numFloors; ++i) {
    for (unsigned int j = 0; j < numTimeSteps; ++j) {
      u.getPositionValues(j,auxVec);
      myFile << auxVec[i] << " ";
    }
    myFile << std::endl;
  }
  myFile << "];" << std::endl;

  // Write 'ud_cpp'
  myFile << "ud_cpp = zeros(" << numFloors << "," << numTimeSteps << ");\n"
          << "ud_cpp = [";
  for (unsigned int i = 0; i < numFloors; ++i) {
    for (unsigned int j = 0; j < numTimeSteps; ++j) {
      ud.getPositionValues(j,auxVec);
      myFile << auxVec[i] << " ";
    }
    myFile << std::endl;
  }
  myFile << "];" << std::endl;

  // Write 'udd_cpp'
  myFile << "udd_cpp = zeros(" << numFloors << "," << numTimeSteps << ");\n"
          << "udd_cpp = [";
  for (unsigned int i = 0; i < numFloors; ++i) {
    for (unsigned int j = 0; j < numTimeSteps; ++j) {
      udd.getPositionValues(j,auxVec);
      myFile << auxVec[i] << " ";
    }
    myFile << std::endl;
  }
  myFile << "];" << std::endl;

  // Write 'resfor_cpp'
  myFile << "resfor_cpp = zeros(" << numFloors << "," << numTimeSteps << ");\n"
          << "resfor_cpp = [";
  for (unsigned int i = 0; i < numFloors; ++i) {
    for (unsigned int j = 0; j < numTimeSteps; ++j) {
      resfor.getPositionValues(j,auxVec);
      myFile << auxVec[i] << " ";
    }
    myFile << std::endl;
  }
  myFile << "];" << std::endl;

  // Write 'ru_cpp'
  myFile << "ru_cpp = zeros(" << numFloors << "," << numTimeSteps << ");\n"
          << "ru_cpp = [";
  for (unsigned int i = 0; i < numFloors; ++i) {
    for (unsigned int j = 0; j < numTimeSteps; ++j) {
      ru.getPositionValues(j,auxVec);
      myFile << auxVec[i] << " ";
    }
    myFile << std::endl;
  }
  myFile << "];" << std::endl;

  myFile.close();

  return;
}
