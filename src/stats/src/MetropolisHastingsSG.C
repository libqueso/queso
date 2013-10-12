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
// $Id:$
//
//--------------------------------------------------------------------------

#include <queso/MetropolisHastingsSG1.h>

namespace QUESO {

// Default constructor -----------------------------
MHRawChainInfoStruct::MHRawChainInfoStruct()
{
  reset();
}
// Copy constructor----------------------------------
MHRawChainInfoStruct::MHRawChainInfoStruct(const MHRawChainInfoStruct& rhs)
{
  this->copy(rhs);
}
// Destructor ---------------------------------------
MHRawChainInfoStruct::~MHRawChainInfoStruct()
{
}
// Set methods---------------------------------------
MHRawChainInfoStruct&
MHRawChainInfoStruct::operator=(const MHRawChainInfoStruct& rhs)
{
  this->copy(rhs);
  return *this;
}
//---------------------------------------------------
MHRawChainInfoStruct&
MHRawChainInfoStruct::operator+=(const MHRawChainInfoStruct& rhs)
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
// Misc methods--------------------------------------------------
void
MHRawChainInfoStruct::reset()
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
//---------------------------------------------------
void
MHRawChainInfoStruct::copy(const MHRawChainInfoStruct& rhs)
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
//---------------------------------------------------
void
MHRawChainInfoStruct::mpiSum(const MpiComm& comm, MHRawChainInfoStruct& sumInfo) const
{
  comm.Allreduce((void *) &runTime, (void *) &sumInfo.runTime, (int) 7, RawValue_MPI_DOUBLE, RawValue_MPI_SUM,
                 "MHRawChainInfoStruct::mpiSum()",
                 "failed MPI.Allreduce() for sum of doubles");

  comm.Allreduce((void *) &numTargetCalls, (void *) &sumInfo.numTargetCalls, (int) 5, RawValue_MPI_UNSIGNED, RawValue_MPI_SUM,
                 "MHRawChainInfoStruct::mpiSum()",
                 "failed MPI.Allreduce() for sum of unsigned ints");

  return;
}

}  // End namespace QUESO
