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
#include <queso/GaussianJointPdf.h>
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

  QUESO::GslVector lawexp(paramSpace.zeroVector());
  lawexp[0] = 2.4;
  lawexp[1] = 0.4;

  QUESO::GslVector lawvar(paramSpace.zeroVector());
  lawvar[0] = 1.2;
  lawvar[1] = 0.9;
  QUESO::GaussianJointPdf<> pdf1("", paramSpace, lawexp, lawvar);

  QUESO::GslMatrix lawcovar(paramSpace.zeroVector());
  lawcovar(0,0) = 1.5;
  lawcovar(0,1) = lawcovar(1,0) = 1.2;
  lawcovar(1,1) = 0.5;
  QUESO::GaussianJointPdf<> pdf2("", paramSpace, lawexp, lawcovar);

  QUESO::GslVector mean1(paramSpace.zeroVector());
  pdf1.distributionMean(mean1);

  QUESO::GslVector mean2(paramSpace.zeroVector());
  pdf1.distributionMean(mean2);

  const char *msg = "GaussianJointPdf mean is incorrect";
  queso_require_less_equal_msg(std::abs(mean1[0]-lawexp[0]), TOL, msg);
  queso_require_less_equal_msg(std::abs(mean1[1]-lawexp[1]), TOL, msg);
  queso_require_less_equal_msg(std::abs(mean2[0]-lawexp[0]), TOL, msg);
  queso_require_less_equal_msg(std::abs(mean2[1]-lawexp[1]), TOL, msg);

  QUESO::GslMatrix var1(paramSpace.zeroVector());
  pdf1.distributionVariance(var1);

  QUESO::GslMatrix var2(paramSpace.zeroVector());
  pdf2.distributionVariance(var2);

  const char *msgv = "GaussianJointPdf variance is incorrect";
  queso_require_less_equal_msg(std::abs(var1(0,0)-lawvar[0]), TOL, msgv);
  queso_require_less_equal_msg(std::abs(var1(0,1)), TOL, msgv);
  queso_require_less_equal_msg(std::abs(var1(1,0)), TOL, msgv);
  queso_require_less_equal_msg(std::abs(var1(1,1)-lawvar[1]), TOL, msgv);
  queso_require_less_equal_msg(std::abs(var2(0,0)-lawcovar(0,0)), TOL, msgv);
  queso_require_less_equal_msg(std::abs(var2(0,1)-lawcovar(0,1)), TOL, msgv);
  queso_require_less_equal_msg(std::abs(var2(1,0)-lawcovar(1,0)), TOL, msgv);
  queso_require_less_equal_msg(std::abs(var2(1,1)-lawcovar(1,1)), TOL, msgv);

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return 0;
}
