//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010 The PECOS Development Team
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
// $Id$
//
//--------------------------------------------------------------------------

#ifndef __UQ_MISCELLANEOUS_H__
#define __UQ_MISCELLANEOUS_H__

#include <gsl/gsl_rng.h>
#include <string>
#include <vector>
#include <uqEnvironment.h>

void   uqMiscReadDoublesFromString      (const std::string&   inputString,
                                         std::vector<double>& outputDoubles);
void   uqMiscReadWordsFromString        (const std::string&        inputString,
                                         std::vector<std::string>& outputWords);
//void   uqMiscExtractDoubleFromString    (std::string& inputString,
//                                         double&      outputDouble);
//void   uqMiscExtractWordFromString      (std::string& inputString,
//                                         std::string& outputWord);
int    uqMiscReadStringAndDoubleFromFile(std::ifstream& ifs,
                                         std::string&   termString,
                                         double*        termValue);
int    uqMiscReadCharsAndDoubleFromFile (std::ifstream& ifs,
                                         std::string&   termString,
                                         double*        termValue,
                                         bool&          endOfLineAchieved);
double uqMiscGammar                     (double         a,
                                         double         b,
                                         const gsl_rng* rng);
double uqMiscGetEllapsedSeconds         (struct timeval *timeval0);

double uqMiscHammingWindow              (unsigned int N, unsigned int j);

double uqMiscGaussianDensity            (double x, double mu, double sigma);

template <class T>
bool
uqMiscCheckForSameValueInAllNodes(T&                    inputValue, // Yes, 'not' const
                                  double                acceptableTreshold,
                                  const Epetra_MpiComm& comm,
                                  const char*           whereString)
{
  // Filter out those nodes that should not participate
  if (comm.MyPID() < 0) return true;

  double localValue = (double) inputValue;
  double sumValue = 0.;
  int mpiRC = MPI_Allreduce((void *) &localValue, (void *) &sumValue, (int) 1, MPI_DOUBLE, MPI_SUM, comm.Comm());
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      UQ_UNAVAILABLE_RANK,
                      whereString,
                      "failed MPI on 'sumValue' inside uqMiscCheckForSameValueInAllNodes()");

  double totalNumNodes = (double) comm.NumProc();
  double testValue = fabs(1. - localValue/(sumValue/totalNumNodes));
  unsigned int boolSum = 0;
#if 1
  unsigned int boolResult = 0;
  if (testValue > acceptableTreshold) boolResult = 1;
  mpiRC = MPI_Allreduce((void *) &boolResult, (void *) &boolSum, (int) 1, MPI_UNSIGNED, MPI_SUM, comm.Comm());
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      UQ_UNAVAILABLE_RANK,
                      whereString,
                      "failed MPI on 'boolSum' inside uqMiscCheckForSameValueInAllNodes()");

  if (boolSum > 0) { 
    comm.Barrier();
    for (int i = 0; i < comm.NumProc(); ++i) {
      if (i == comm.MyPID()) {
        std::cerr << "WARNING, "
                  << whereString
                  << ", inside uqMiscCheckForSameValueInAllNodes()"
                  << ", rank (in this communicator) = " << i
                  << ": boolSum = "       << boolSum
                  << ", localValue = "    << localValue
                  << ", sumValue = "      << sumValue
                  << ", totalNumNodes = " << totalNumNodes
                  << ", avgValue = "      << (sumValue/totalNumNodes)
                  << ", relativeTest = "  << testValue
                  << std::endl;
      }
      comm.Barrier();
    }
    comm.Barrier();

    mpiRC = MPI_Bcast((void *) &localValue, (int) 1, MPI_DOUBLE, 0, comm.Comm());
    UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                        UQ_UNAVAILABLE_RANK,
                        whereString,
                        "failed MPI on 'boolSum' inside uqMiscCheckForSameValueInAllNodes()");
    inputValue = localValue; // IMPORTANT
  }
#else
  UQ_FATAL_TEST_MACRO(testValue > acceptableTreshold,
                      UQ_UNAVAILABLE_RANK,
                      whereString,
                      "not all nodes have the same value inside uqMiscCheckForSameValueInAllNodes()");
#endif

  return (boolSum == 0);
}

template <class V>
void
uqMiscComputePositionsBetweenMinMax(V                minValues,
                                    V                maxValues,
                                    std::vector<V*>& positions)
{
  double factor = 0.5;
  switch (positions.size()) {
    case 0:
      // Do nothing
    break;

    case 1:
      positions[0] = new V((1. - factor) * minValues + factor * maxValues);
    break;

    default:
      for (unsigned int i = 0; i < positions.size(); ++i) {
        factor = ((double) i)/(((double) positions.size()) - 1.);
        positions[i] = new V((1. - factor) * minValues + factor * maxValues);
      }
    break;
  }

  return;
}

template <class V1,class V2>
void
uqMiscCheckTheParallelEnvironment(const V1& vec1, const V2& vec2)
{
  const uqBaseEnvironmentClass& env = vec1.env();

  if (env.numSubEnvironments() == (unsigned int) env.fullComm().NumProc()) {
    UQ_FATAL_TEST_MACRO(env.subRank() != 0,
                        env.worldRank(),
                        "uqMiscCheckTheParallelEnvironment<V1,V2>()",
                        "there should exist only one processor per sub environment");
    UQ_FATAL_TEST_MACRO((vec1.numOfProcsForStorage() != 1) ||
                        (vec2.numOfProcsForStorage() != 1),
                        env.worldRank(),
                        "uqMiscCheckTheParallelEnvironment<V1,V2>()",
                        "only 1 processor (per sub environment) should be necessary for the storage of a parameter vector");
  }
  else if (env.numSubEnvironments() < (unsigned int) env.fullComm().NumProc()) {
    UQ_FATAL_TEST_MACRO(env.fullComm().NumProc()%env.numSubEnvironments() != 0,
                        env.worldRank(),
                        "uqMiscCheckTheParallelEnvironment<V1,V2>()",
                        "total number of processors should be a multiple of the number of sub environments");
    unsigned int numProcsPerSubEnvironment = env.fullComm().NumProc()/env.numSubEnvironments();
    UQ_FATAL_TEST_MACRO(env.subComm().NumProc() != (int) numProcsPerSubEnvironment,
                        env.worldRank(),
                        "uqMiscCheckTheParallelEnvironment<V1,V2>()",
                        "inconsistent number of processors per sub environment");
    if ((vec1.numOfProcsForStorage() == 1) &&
        (vec2.numOfProcsForStorage() == 1)) {
      // Ok
    }
    else if ((vec1.numOfProcsForStorage() == numProcsPerSubEnvironment) &&
             (vec2.numOfProcsForStorage() == numProcsPerSubEnvironment)) {
      UQ_FATAL_TEST_MACRO(true,
                          env.worldRank(),
                          "uqMiscCheckTheParallelEnvironment<V1,V2>()",
                          "parallel vectors are not supported yet");
    }
    else {
      UQ_FATAL_TEST_MACRO(true,
                          env.worldRank(),
                          "uqMiscCheckTheParallelEnvironment<V1,V2>()",
                          "number of processors required for a vector storage should be equal to either 1 or to the number of processors in the sub environment");
    }
  }
  else {
    UQ_FATAL_TEST_MACRO(true,
                        env.worldRank(),
                        "uqMiscCheckTheParallelEnvironment<V1,V2>()",
                        "number of processors per sub environment is less than 1!");
  }

  return;
}
#endif // __UQ_MISCELLANEOUS_H__
