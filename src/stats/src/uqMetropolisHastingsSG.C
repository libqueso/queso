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

#include <uqMetropolisHastingsSG1.h>

uqMHRawChainInfoStruct::uqMHRawChainInfoStruct()
{
  reset();
}

uqMHRawChainInfoStruct::~uqMHRawChainInfoStruct()
{
}

uqMHRawChainInfoStruct::uqMHRawChainInfoStruct(const uqMHRawChainInfoStruct& rhs)
{
  this->copy(rhs);
}

uqMHRawChainInfoStruct&
uqMHRawChainInfoStruct::operator=(const uqMHRawChainInfoStruct& rhs)
{
  this->copy(rhs);
  return *this;
}

uqMHRawChainInfoStruct&
uqMHRawChainInfoStruct::operator+=(const uqMHRawChainInfoStruct& rhs)
{
  runTime          += rhs.runTime;
  candidateRunTime += rhs.candidateRunTime;
  targetRunTime    += rhs.targetRunTime;
  mhAlphaRunTime   += rhs.mhAlphaRunTime;
  drAlphaRunTime   += rhs.drAlphaRunTime;
  drRunTime        += rhs.drRunTime;
  amRunTime        += rhs.amRunTime;

  numTargetCalls            += rhs.numTargetCalls;
  numDRs                    += rhs.numDRs;
  numOutOfTargetSupport     += rhs.numOutOfTargetSupport;
  numOutOfTargetSupportInDR += rhs.numOutOfTargetSupportInDR;
  numRejections             += rhs.numRejections;

  return *this;
}

void
uqMHRawChainInfoStruct::reset()
{
  runTime          = 0.;
  candidateRunTime = 0.;
  targetRunTime    = 0.;
  mhAlphaRunTime   = 0.;
  drAlphaRunTime   = 0.;
  drRunTime        = 0.;
  amRunTime        = 0.;

  numTargetCalls            = 0;
  numDRs                    = 0;
  numOutOfTargetSupport     = 0;
  numOutOfTargetSupportInDR = 0;
  numRejections             = 0;
}

void
uqMHRawChainInfoStruct::copy(const uqMHRawChainInfoStruct& rhs)
{
  runTime          = rhs.runTime;
  candidateRunTime = rhs.candidateRunTime;
  targetRunTime    = rhs.targetRunTime;
  mhAlphaRunTime   = rhs.mhAlphaRunTime;
  drAlphaRunTime   = rhs.drAlphaRunTime;
  drRunTime        = rhs.drRunTime;
  amRunTime        = rhs.amRunTime;

  numTargetCalls            = rhs.numTargetCalls;
  numDRs                    = rhs.numDRs;
  numOutOfTargetSupport     = rhs.numOutOfTargetSupport;
  numOutOfTargetSupportInDR = rhs.numOutOfTargetSupportInDR;
  numRejections             = rhs.numRejections;

  return;
}

void
uqMHRawChainInfoStruct::mpiSum(const MPI_Comm& comm, uqMHRawChainInfoStruct& sumInfo) const
{
  int mpiRC = MPI_Allreduce((void *) &runTime, (void *) &sumInfo.runTime, (int) 7, MPI_DOUBLE, MPI_SUM, comm);
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      UQ_UNAVAILABLE_RANK,
                      "uqMHRawChainInfoStruct::mpiSum()",
                      "failed MPI_Allreduce() for sum of doubles");

  mpiRC = MPI_Allreduce((void *) &numTargetCalls, (void *) &sumInfo.numTargetCalls, (int) 5, MPI_UNSIGNED, MPI_SUM, comm);
  UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                      UQ_UNAVAILABLE_RANK,
                      "uqMHRawChainInfoStruct::mpiSum()",
                      "failed MPI_Allreduce() for sum of unsigned ints");

  return;
}
