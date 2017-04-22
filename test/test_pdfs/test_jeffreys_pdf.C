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
#include <queso/JeffreysJointPdf.h>

#define TOL 1e-15

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
  paramMaxs.cwSet(INFINITY);

  QUESO::BoxSubset<> paramDomain("param_", paramSpace, paramMins, paramMaxs);

  QUESO::JeffreysJointPdf<> pdf("", paramDomain);

  QUESO::GslVector point(paramSpace.zeroVector());
  point[0] = 2.8;

  double lnPdfValue1 = pdf.lnValue(point, NULL, NULL, NULL, NULL);
  double lnPdfValue2 = std::log(1.0 / point[0]);

  double pdfValue1 = pdf.actualValue(point, NULL, NULL, NULL, NULL);
  double pdfValue2 = std::exp(lnPdfValue2);

  queso_require_less_equal_msg(std::abs(lnPdfValue1 - lnPdfValue2), TOL, "log jeffreys pdf values are incorrect");
  queso_require_less_equal_msg(std::abs(pdfValue1 - pdfValue2), TOL, "jeffreys pdf values are incorrect");

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return 0;
}
