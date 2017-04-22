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
#include <queso/BetaJointPdf.h>
#include <queso/ConcatenatedJointPdf.h>
#include <queso/GammaJointPdf.h>
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

  QUESO::VectorSpace<> paramSpace1(env, "param1_", 1, NULL);

  QUESO::GslVector paramMins(paramSpace1.zeroVector());
  paramMins.cwSet(0.0);

  QUESO::GslVector paramMaxs(paramSpace1.zeroVector());
  paramMaxs.cwSet(1.0);

  QUESO::BoxSubset<> paramDomain1("param1_", paramSpace1, paramMins, paramMaxs);

  // We should test the other cases of alpha and beta
  QUESO::GslVector alpha(paramSpace1.zeroVector());
  alpha[0] = 2.0;

  QUESO::GslVector beta(paramSpace1.zeroVector());
  beta[0] = 3.0;

  QUESO::BetaJointPdf<> pdf1("", paramDomain1, alpha, beta);

  QUESO::VectorSpace<> paramSpace2(env, "param2_", 2, NULL);

  QUESO::GslVector k(paramSpace2.zeroVector());
  k[0] = 3.0;
  k[1] = 0.5;

  QUESO::GslVector theta(paramSpace2.zeroVector());
  theta[0] = 4.0;
  theta[1] = 2.5;
  QUESO::GammaJointPdf<> pdf2("", paramSpace2, k, theta);

  QUESO::VectorSpace<> paramSpace(env, "param_", 3, NULL);

  QUESO::ConcatenatedJointPdf<> pdf("", pdf1, pdf2, paramSpace);

  QUESO::GslVector mean(paramSpace.zeroVector());
  pdf.distributionMean(mean);

  const char *msg = "ConcatenatedJointPdf mean is incorrect";
  double real_mean0 = alpha[0] / (alpha[0] + beta[0]);
  double real_mean1 = k[0] * theta[0];
  double real_mean2 = k[1] * theta[1];
  queso_require_less_equal_msg(std::abs(mean[0]-real_mean0), TOL, msg);
  queso_require_less_equal_msg(std::abs(mean[1]-real_mean1), TOL, msg);
  queso_require_less_equal_msg(std::abs(mean[2]-real_mean2), TOL, msg);

  const char *msgv = "ConcatenatedJointPdf variance is incorrect";
  QUESO::GslMatrix var(paramSpace.zeroVector());
  pdf.distributionVariance(var);
  double real_var0 = alpha[0] * beta[0] / (alpha[0] + beta[0]) /
          (alpha[0] + beta[0]) / (alpha[0] + beta[0] + 1);
  double real_var1 = k[0] * theta[0] * theta[0];
  double real_var2 = k[1] * theta[1] * theta[1];
  queso_require_less_equal_msg(std::abs(var(0,0)-real_var0), TOL, msgv);
  queso_require_less_equal_msg(std::abs(var(1,1)-real_var1), TOL, msgv);
  queso_require_less_equal_msg(std::abs(var(2,2)-real_var2), TOL, msgv);

  for (unsigned int i=0; i != 3; ++i)
    for (unsigned int j=0; j != 3; ++j)
      if (i != j)
        queso_require_less_equal_msg(std::abs(var(i,j)), TOL, msgv);

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return 0;
}
