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
#include <queso/BetaJointPdf.h>

#define TOL 1e-14

int main(int argc, char ** argv)
{
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);

  QUESO::FullEnvironment env(MPI_COMM_WORLD, "", "", NULL);
#else
  QUESO::FullEnvironment env("", "", NULL);
#endif

  QUESO::VectorSpace<> paramSpace(env, "param_", 1, NULL);

  QUESO::GslVector paramMins(paramSpace.zeroVector());
  paramMins.cwSet(0.0);

  QUESO::GslVector paramMaxs(paramSpace.zeroVector());
  paramMaxs.cwSet(1.0);

  QUESO::BoxSubset<> paramDomain("param_", paramSpace, paramMins, paramMaxs);

  // We should test the other cases of alpha and beta
  QUESO::GslVector alpha(paramSpace.zeroVector());
  alpha[0] = 2.0;

  QUESO::GslVector beta(paramSpace.zeroVector());
  beta[0] = 3.0;

  QUESO::BetaJointPdf<> pdf("", paramDomain, alpha, beta);

  // This point is the mode of the beta distribution defined above
  QUESO::GslVector point(paramSpace.zeroVector());
  point[0] = (alpha[0] - 1.0) / (alpha[0] + beta[0] - 2.0);

  QUESO::GslVector lnGradVector(paramSpace.zeroVector());
  QUESO::GslVector gradVector(paramSpace.zeroVector());

  // Here we test the gradient of the log pdf and the gradient of the pdf are
  // zero, as we expect.
  pdf.lnValue(point, NULL, &lnGradVector, NULL, NULL);
  pdf.actualValue(point, NULL, &gradVector, NULL, NULL);

  queso_require_less_equal_msg(std::abs(lnGradVector[0]), TOL, "grad log beta pdf values are incorrect");
  queso_require_less_equal_msg(std::abs(gradVector[0]), TOL, "grad beta pdf values are incorrect");

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return 0;
}
