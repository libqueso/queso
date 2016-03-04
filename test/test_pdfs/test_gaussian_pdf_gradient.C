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
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/GaussianJointPdf.h>

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
  QUESO::GslMatrix var(paramSpace.zeroVector());
  mean[0] = 2.0;
  mean[1] = 3.0;
  mean[2] = 4.0;
  var(0,0) = 5.0;
  var(1,1) = 6.0;
  var(2,2) = 7.0;

  // Check the zero-gradient case
  QUESO::GaussianJointPdf<> pdf("", paramDomain, mean, var);

  QUESO::GslVector lnGradVector(paramSpace.zeroVector());
  QUESO::GslVector gradVector(paramSpace.zeroVector());
  QUESO::GslVector point(mean);

  double lnPdfValue1 = pdf.lnValue(point, NULL, &lnGradVector, NULL, NULL);
  queso_require_less_equal_msg(std::abs(lnGradVector[0]), TOL, "grad log gaussian pdf values are incorrect");
  queso_require_less_equal_msg(std::abs(lnGradVector[1]), TOL, "grad log gaussian pdf values are incorrect");
  queso_require_less_equal_msg(std::abs(lnGradVector[2]), TOL, "grad log gaussian pdf values are incorrect");

  double pdfValue1 = pdf.actualValue(point, NULL, &gradVector, NULL, NULL);
  queso_require_less_equal_msg(std::abs(gradVector[0]), TOL, "grad guassian pdf values are incorrect");
  queso_require_less_equal_msg(std::abs(gradVector[1]), TOL, "grad guassian pdf values are incorrect");
  queso_require_less_equal_msg(std::abs(gradVector[2]), TOL, "grad guassian pdf values are incorrect");


  // Check the 1-gradient case (in log space)
  mean[0] = 0.0;
  mean[1] = 0.0;
  mean[2] = 0.0;

  var(0,0) = 1.0;
  var(0,1) = 0.8;
  var(0,2) = 0.7;
  var(1,0) = 0.8;
  var(1,1) = 2.0;
  var(1,2) = 0.6;
  var(2,0) = 0.7;
  var(2,1) = 0.6;
  var(2,2) = 3.0;

  point[0] = -2.5;
  point[1] = -3.4;
  point[2] = -4.3;

  QUESO::GaussianJointPdf<> pdf2("", paramDomain, mean, var);

  lnPdfValue1 = pdf2.lnValue(point, NULL, &lnGradVector, NULL, NULL);

  queso_require_less_equal_msg(std::abs(lnGradVector[0] - 1.0), TOL, "grad log gaussian pdf2 values are incorrect");
  queso_require_less_equal_msg(std::abs(lnGradVector[1] - 1.0), TOL, "grad log gaussian pdf2 values are incorrect");
  queso_require_less_equal_msg(std::abs(lnGradVector[2] - 1.0), TOL, "grad log gaussian pdf2 values are incorrect");

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return 0;
}
