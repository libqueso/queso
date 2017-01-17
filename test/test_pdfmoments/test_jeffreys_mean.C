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

#include <queso/Environment.h>
#include <queso/JeffreysJointPdf.h>
#include <queso/GslVector.h>

#define TOL 1e-14

int main(int argc, char ** argv)
{
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);

  QUESO::FullEnvironment env(MPI_COMM_WORLD, "", "", NULL);
#else
  QUESO::FullEnvironment env("", "", NULL);
#endif

  QUESO::VectorSpace<> paramSpace(env, "param_", 3, NULL);

  QUESO::GslVector paramMins(paramSpace.zeroVector());
  paramMins[0] = 1.0;
  paramMins[1] = 2.0;
  paramMins[2] = 3.0;

  QUESO::GslVector paramMaxs(paramSpace.zeroVector());
  paramMaxs[0] = 2.0;
  paramMaxs[1] = 4.0;
  paramMaxs[2] = 6.0;

  QUESO::BoxSubset<> paramDomain("param_", paramSpace, paramMins, paramMaxs);

  QUESO::JeffreysJointPdf<> pdf("", paramDomain);

  QUESO::GslVector mean(paramSpace.zeroVector());
  pdf.distributionMean(mean);

  const char *msg = "JeffreysJointPdf mean is incorrect";
  QUESO::GslVector real_mean = paramMins;
  real_mean += paramMaxs;
  real_mean /= 2;
  queso_require_less_equal_msg(std::abs(mean[0]-real_mean[0]), TOL, msg);
  queso_require_less_equal_msg(std::abs(mean[1]-real_mean[1]), TOL, msg);
  queso_require_less_equal_msg(std::abs(mean[2]-real_mean[2]), TOL, msg);

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return 0;
}
