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
#include <queso/GslVector.h>
#include <queso/LogNormalJointPdf.h>

#define TOL 1e-14

int main(int argc, char ** argv)
{
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);

  QUESO::FullEnvironment env(MPI_COMM_WORLD, "", "", NULL);
#else
  QUESO::FullEnvironment env("", "", NULL);
#endif

  unsigned int dim = 3;
  QUESO::VectorSpace<> paramSpace(env, "param_", dim, NULL);

  QUESO::GslVector paramMins(paramSpace.zeroVector());
  paramMins.cwSet(-INFINITY);

  QUESO::GslVector paramMaxs(paramSpace.zeroVector());
  paramMaxs.cwSet(INFINITY);

  QUESO::BoxSubset<> paramDomain("param_", paramSpace, paramMins, paramMaxs);

  QUESO::GslVector mean(paramSpace.zeroVector());
  QUESO::GslVector var(paramSpace.zeroVector());
  mean[0] = 2.0;
  mean[1] = 3.0;
  mean[2] = 4.0;
  var[0] = 5.0;
  var[1] = 6.0;
  var[2] = 7.0;

  QUESO::LogNormalJointPdf<> pdf("", paramDomain, mean, var);

  // This point is the mode of the log normal distribution defined above
  QUESO::GslVector point(paramSpace.zeroVector());
  point[0] = std::exp(mean[0] - var[0]);
  point[1] = std::exp(mean[1] - var[1]);
  point[2] = std::exp(mean[2] - var[2]);

  QUESO::GslVector lnGradVector(paramSpace.zeroVector());
  QUESO::GslVector gradVector(paramSpace.zeroVector());

  // Here we test the gradient of the log of the pdf and the gradient of the
  // pdf are zero at the mode, as we expect.
  pdf.lnValue(point, NULL, &lnGradVector, NULL, NULL);
  queso_require_less_equal_msg(std::abs(lnGradVector[0]), TOL, "grad log log_normal pdf values are incorrect");
  queso_require_less_equal_msg(std::abs(lnGradVector[1]), TOL, "grad log log_normal pdf values are incorrect");
  queso_require_less_equal_msg(std::abs(lnGradVector[2]), TOL, "grad log log_normal pdf values are incorrect");

  pdf.actualValue(point, NULL, &gradVector, NULL, NULL);
  queso_require_less_equal_msg(std::abs(gradVector[0]), TOL, "grad log_normal pdf values are incorrect");
  queso_require_less_equal_msg(std::abs(gradVector[1]), TOL, "grad log_normal pdf values are incorrect");
  queso_require_less_equal_msg(std::abs(gradVector[2]), TOL, "grad log_normal pdf values are incorrect");

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return 0;
}
