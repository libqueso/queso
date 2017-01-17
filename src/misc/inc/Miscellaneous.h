//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2017 The PECOS Development Team
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

#ifndef UQ_MISCELLANEOUS_H
#define UQ_MISCELLANEOUS_H

#include <queso/Environment.h>
#include <queso/RngBase.h>
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
bool MiscCheckForSameValueInAllNodes(T & inputValue, // Yes, 'not' const
    double acceptableTreshold, const MpiComm& comm, const char * whereString);

template <class V>
void MiscComputePositionsBetweenMinMax(V minValues, V maxValues,
    std::vector<V*>& positions);

template <class V1,class V2>
void MiscCheckTheParallelEnvironment(const V1& vec1, const V2& vec2);

}  // End namespace QUESO

#endif // UQ_MISCELLANEOUS_H
