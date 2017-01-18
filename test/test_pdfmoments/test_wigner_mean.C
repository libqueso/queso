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
#include <queso/WignerJointPdf.h>
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

  QUESO::VectorSpace<> paramSpace(env, "param_", 2, NULL);

  QUESO::GslVector center(paramSpace.zeroVector());
  center[0] = 2.4;
  center[1] = 1.2;

  double radius = 1.5;
  QUESO::WignerJointPdf<> pdf("", paramSpace, center, radius);

  QUESO::GslVector mean(paramSpace.zeroVector());
  pdf.distributionMean(mean);

  const char *msg = "WignerJointPdf mean is incorrect";
  queso_require_less_equal_msg(std::abs(mean[0]-center[0]), TOL, msg);
  queso_require_less_equal_msg(std::abs(mean[1]-center[1]), TOL, msg);

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return 0;
}
