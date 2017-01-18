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
#include <queso/LogNormalJointPdf.h>
#include <queso/GslMatrix.h>
#include <queso/GslVector.h>

#include <cmath>

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

  QUESO::GslVector lawexp(paramSpace.zeroVector());
  lawexp[0] = 2.4;
  lawexp[1] = 0.4;

  QUESO::GslVector lawvar(paramSpace.zeroVector());
  lawvar[0] = 1.2;
  lawvar[1] = 0.9;
  QUESO::LogNormalJointPdf<> pdf("", paramSpace, lawexp, lawvar);

  QUESO::GslVector mean(paramSpace.zeroVector());
  pdf.distributionMean(mean);

  double realmean0 = std::exp(lawexp[0] + lawvar[0]/2);
  double realmean1 = std::exp(lawexp[1] + lawvar[1]/2);

  const char *msg = "LogNormalJointPdf mean is incorrect";
  queso_require_less_equal_msg(std::abs(mean[0]-realmean0), TOL, msg);
  queso_require_less_equal_msg(std::abs(mean[1]-realmean1), TOL, msg);

  QUESO::GslMatrix var(paramSpace.zeroVector());
  pdf.distributionVariance(var);

  double realvar0 = (std::exp(lawvar[0])-1) * std::exp(2*lawexp[0] + lawvar[0]);
  double realvar1 = (std::exp(lawvar[1])-1) * std::exp(2*lawexp[1] + lawvar[1]);

  const char *msgv = "LogNormalJointPdf variance is incorrect";
  queso_require_less_equal_msg(std::abs(var(0,0)-realvar0), TOL, msgv);
  queso_require_less_equal_msg(std::abs(var(0,1)), TOL, msgv);
  queso_require_less_equal_msg(std::abs(var(1,0)), TOL, msgv);
  queso_require_less_equal_msg(std::abs(var(1,1)-realvar1), TOL, msgv);

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return 0;
}
