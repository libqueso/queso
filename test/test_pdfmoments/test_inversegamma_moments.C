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
#include <queso/InverseGammaJointPdf.h>
#include <queso/GslMatrix.h>
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

  QUESO::GslVector alpha(paramSpace.zeroVector());
  alpha[0] = 2.4;
  alpha[1] = 2.2;

  QUESO::GslVector beta(paramSpace.zeroVector());
  beta[0] = 2.0;
  beta[1] = 0.8;
  QUESO::InverseGammaJointPdf<> pdf("", paramSpace, alpha, beta);

  QUESO::GslVector mean(paramSpace.zeroVector());
  pdf.distributionMean(mean);

  const char *msg = "InverseGammaJointPdf mean is incorrect";
  double real_mean0 = beta[0] / (alpha[0] - 1);
  double real_mean1 = beta[1] / (alpha[1] - 1);
  queso_require_less_equal_msg(std::abs(mean[0]-real_mean0), TOL, msg);
  queso_require_less_equal_msg(std::abs(mean[1]-real_mean1), TOL, msg);

  QUESO::GslMatrix var(paramSpace.zeroVector());
  pdf.distributionVariance(var);

  const char *msgv = "InverseGammaJointPdf variance is incorrect";
  double real_var0 = beta[0] * beta[0] /
    (alpha[0] - 1) / (alpha[0] - 1) / (alpha[0] - 2);
  double real_var1 = beta[1] * beta[1] /
    (alpha[1] - 1) / (alpha[1] - 1) / (alpha[1] - 2);
  queso_require_less_equal_msg(std::abs(var(0,0)-real_var0), TOL, msgv);
  queso_require_less_equal_msg(std::abs(var(0,1)), TOL, msgv);
  queso_require_less_equal_msg(std::abs(var(1,0)), TOL, msgv);
  queso_require_less_equal_msg(std::abs(var(1,1)-real_var1), TOL, msgv);

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return 0;
}
