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

#include <queso/Environment.h>
#include <queso/GammaJointPdf.h>
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

  QUESO::GslVector k(paramSpace.zeroVector());
  k[0] = 1.4;
  k[1] = 0.1;

  QUESO::GslVector theta(paramSpace.zeroVector());
  theta[0] = 2.0;
  theta[1] = 0.8;
  QUESO::GammaJointPdf<> pdf("", paramSpace, k, theta);

  QUESO::GslVector mean(paramSpace.zeroVector());
  pdf.distributionMean(mean);

  const char *msg = "ConcatenatedJointPdf mean is incorrect";
  double real_mean0 = k[0] * theta[0];
  double real_mean1 = k[1] * theta[1];
  queso_require_less_equal_msg(std::abs(mean[0]-real_mean0), TOL, msg);
  queso_require_less_equal_msg(std::abs(mean[1]-real_mean1), TOL, msg);

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return 0;
}
