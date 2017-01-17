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
#include <queso/ConcatenationSubset.h>
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

  QUESO::VectorSpace<> paramSpace1(env, "param1_", 4, NULL);
  QUESO::VectorSpace<> paramSpace2(env, "param2_", 2, NULL);

  QUESO::GslVector paramMins1(paramSpace1.zeroVector());
  paramMins1[0] = 1e2;
  paramMins1[1] = -1e5;
  paramMins1[2] = 4e-3;
  paramMins1[3] = 1;

  QUESO::GslVector paramMaxs1(paramSpace1.zeroVector());
  paramMaxs1[0] = 2e2;
  paramMaxs1[1] = 1e5;
  paramMaxs1[2] = 6e-3;
  paramMaxs1[3] = 11;

  QUESO::BoxSubset<> paramDomain1("", paramSpace1, paramMins1, paramMaxs1);

  QUESO::GslVector paramMins2(paramSpace2.zeroVector());
  paramMins2[0] = -1e5;
  paramMins2[1] = 2e-3;

  QUESO::GslVector paramMaxs2(paramSpace2.zeroVector());
  paramMaxs2[0] = 1e5;
  paramMaxs2[1] = 4e-3;

  QUESO::BoxSubset<> paramDomain2("", paramSpace2, paramMins2, paramMaxs2);

  QUESO::VectorSpace<> paramSpace(env, "param_", 6, NULL);

  QUESO::ConcatenationSubset<QUESO::GslVector,QUESO::GslMatrix>
    paramDomain("",paramSpace,paramDomain1,paramDomain2);

  QUESO::GslVector centroid(paramSpace.zeroVector());
  paramDomain.centroid(centroid);

  const char *msg = "ConcatenationSubset centroid is incorrect";
  queso_require_less_equal_msg(std::abs(centroid[0]-1.5e2), TOL, msg);
  queso_require_less_equal_msg(std::abs(centroid[1]), TOL, msg);
  queso_require_less_equal_msg(std::abs(centroid[2]-5e-3), TOL, msg);
  queso_require_less_equal_msg(std::abs(centroid[3]-6), TOL, msg);
  queso_require_less_equal_msg(std::abs(centroid[4]), TOL, msg);
  queso_require_less_equal_msg(std::abs(centroid[5]-3e-3), TOL, msg);

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return 0;
}
