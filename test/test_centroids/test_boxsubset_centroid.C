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
#include <queso/BoxSubset.h>
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

  unsigned int dim = 4;
  QUESO::VectorSpace<> paramSpace(env, "param_", dim, NULL);

  QUESO::GslVector paramMins(paramSpace.zeroVector());
  paramMins[0] = 1e2;
  paramMins[1] = -1e5;
  paramMins[2] = 4e-3;
  paramMins[3] = 1;

  QUESO::GslVector paramMaxs(paramSpace.zeroVector());
  paramMaxs[0] = 2e2;
  paramMaxs[1] = 1e5;
  paramMaxs[2] = 6e-3;
  paramMaxs[3] = 11;

  QUESO::BoxSubset<> paramDomain("param_", paramSpace, paramMins, paramMaxs);

  QUESO::GslVector centroid(paramSpace.zeroVector());
  paramDomain.centroid(centroid);

  const char *msg = "BoxSubset centroid is incorrect";
  queso_require_less_equal_msg(std::abs(centroid[0]-1.5e2), TOL, msg);
  queso_require_less_equal_msg(std::abs(centroid[1]), TOL, msg);
  queso_require_less_equal_msg(std::abs(centroid[2]-5e-3), TOL, msg);
  queso_require_less_equal_msg(std::abs(centroid[3]-6), TOL, msg);

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return 0;
}
