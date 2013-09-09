//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
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

#include <uqEnvironment.h>
#include <uqRngBase.h>
#include <string>
#include <vector>
#include <math.h>

namespace QUESO {

void         MiscReadDoublesFromString      (const std::string&        inputString,
                                               std::vector<double>&      outputDoubles);
void         MiscReadWordsFromString        (const std::string&        inputString,
                                               std::vector<std::string>& outputWords);
#if 0
void         MiscExtractDoubleFromString    (std::string&              inputString,
                                               double&                   outputDouble);
void         MiscExtractWordFromString      (std::string&              inputString,
                                               std::string&              outputWord);
#endif
int          MiscReadStringAndDoubleFromFile(std::ifstream&            ifs,
                                               std::string&              termString,
                                               double*                   termValue);
int          MiscReadCharsAndDoubleFromFile (std::ifstream&            ifs,
                                               std::string&              termString,
                                               double*                   termValue,
                                               bool&                     endOfLineAchieved);
double       MiscGammar                     (double                    a,
                                               double                    b,
                                               const RngBase*     rngObject);
double       MiscGetEllapsedSeconds         (struct timeval*           timeval0);
double       MiscHammingWindow              (unsigned int              N,
                                               unsigned int              j);
double       MiscGaussianDensity            (double                    x,
                                               double                    mu,
                                               double                    sigma);
unsigned int MiscUintDebugMessage           (unsigned int              value,
                                               const char*               message);
int          MiscIntDebugMessage            (int                       value,
                                               const char*               message);
double       MiscDoubleDebugMessage         (double                    value,
                                               const char*               message);

int          CheckFilePath                  (const char*               path);
int          GRVY_CheckDir                  (const char*               dirname);

template <class T>
bool
MiscCheckForSameValueInAllNodes(T&                    inputValue, // Yes, 'not' const
                                  double                acceptableTreshold,
                                  const MpiComm& comm,
                                  const char*           whereString)
{
  // Filter out those nodes that should not participate
  if (comm.MyPID() < 0) return true;

  double localValue = (double) inputValue;
  double sumValue = 0.;
  comm.Allreduce((void *) &localValue, (void *) &sumValue, (int) 1, RawValue_MPI_DOUBLE, RawValue_MPI_SUM,
                 whereString,
                 "failed MPI on 'sumValue' inside MiscCheckForSameValueInAllNodes()");

  double totalNumNodes = (double) comm.NumProc();
  double testValue = fabs(1. - localValue/(sumValue/totalNumNodes));
  unsigned int boolSum = 0;
#if 1
  unsigned int boolResult = 0;
  if (testValue > acceptableTreshold) boolResult = 1;
  comm.Allreduce((void *) &boolResult, (void *) &boolSum, (int) 1, RawValue_MPI_UNSIGNED, RawValue_MPI_SUM,
                 whereString,
                 "failed MPI on 'boolSum' inside MiscCheckForSameValueInAllNodes()");

  if (boolSum > 0) { 
    comm.Barrier();
    for (int i = 0; i < comm.NumProc(); ++i) {
      if (i == comm.MyPID()) {
        std::cerr << "WARNING, "
                  << whereString
                  << ", inside MiscCheckForSameValueInAllNodes()"
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

    comm.Bcast((void *) &localValue, (int) 1, RawValue_MPI_DOUBLE, 0,
               whereString,
               "failed MPI on 'boolSum' inside MiscCheckForSameValueInAllNodes()");
    inputValue = localValue; // IMPORTANT
  }
#else
  UQ_FATAL_TEST_MACRO(testValue > acceptableTreshold,
                      UQ_UNAVAILABLE_RANK,
                      whereString,
                      "not all nodes have the same value inside MiscCheckForSameValueInAllNodes()");
#endif

  return (boolSum == 0);
}

template <class V>
void
MiscComputePositionsBetweenMinMax(V                minValues,
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
MiscCheckTheParallelEnvironment(const V1& vec1, const V2& vec2)
{
  const BaseEnvironment& env = vec1.env();

  if (env.numSubEnvironments() == (unsigned int) env.fullComm().NumProc()) {
    UQ_FATAL_TEST_MACRO(env.subRank() != 0,
                        env.worldRank(),
                        "MiscCheckTheParallelEnvironment<V1,V2>()",
                        "there should exist only one processor per sub environment");
    UQ_FATAL_TEST_MACRO((vec1.numOfProcsForStorage() != 1) ||
                        (vec2.numOfProcsForStorage() != 1),
                        env.worldRank(),
                        "MiscCheckTheParallelEnvironment<V1,V2>()",
                        "only 1 processor (per sub environment) should be necessary for the storage of a parameter vector");
  }
  else if (env.numSubEnvironments() < (unsigned int) env.fullComm().NumProc()) {
    UQ_FATAL_TEST_MACRO(env.fullComm().NumProc()%env.numSubEnvironments() != 0,
                        env.worldRank(),
                        "MiscCheckTheParallelEnvironment<V1,V2>()",
                        "total number of processors should be a multiple of the number of sub environments");
    unsigned int numProcsPerSubEnvironment = env.fullComm().NumProc()/env.numSubEnvironments();
    UQ_FATAL_TEST_MACRO(env.subComm().NumProc() != (int) numProcsPerSubEnvironment,
                        env.worldRank(),
                        "MiscCheckTheParallelEnvironment<V1,V2>()",
                        "inconsistent number of processors per sub environment");
    if ((vec1.numOfProcsForStorage() == 1) &&
        (vec2.numOfProcsForStorage() == 1)) {
      // Ok
    }
    else if ((vec1.numOfProcsForStorage() == numProcsPerSubEnvironment) &&
             (vec2.numOfProcsForStorage() == numProcsPerSubEnvironment)) {
      UQ_FATAL_TEST_MACRO(true,
                          env.worldRank(),
                          "MiscCheckTheParallelEnvironment<V1,V2>()",
                          "parallel vectors are not supported yet");
    }
    else {
      UQ_FATAL_TEST_MACRO(true,
                          env.worldRank(),
                          "MiscCheckTheParallelEnvironment<V1,V2>()",
                          "number of processors required for a vector storage should be equal to either 1 or to the number of processors in the sub environment");
    }
  }
  else {
    UQ_FATAL_TEST_MACRO(true,
                        env.worldRank(),
                        "MiscCheckTheParallelEnvironment<V1,V2>()",
                        "number of processors per sub environment is less than 1!");
  }

  return;
}

}  // End namespace QUESO

#endif // __UQ_MISCELLANEOUS_H__
