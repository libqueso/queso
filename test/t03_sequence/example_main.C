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
//
// Brief description of this file:
//
// This is an example of how to define and use some basic classes and
// algorithms of QUESO, mainly 'uqGaussianVectorRV',
// 'uqVectorSequence', realizations and statistical computations.
// The code itself is in the routine 'compute(*env)'. This routine is
// called right after the initialization of the MPI environment and of
// the QUESO environment and is available in file 'example_compute.C'.
//
//--------------------------------------------------------------------------


#include <example_compute.h>

int main(int argc, char* argv[])
{
  // Initialize environment
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc,&argv);
  QUESO::FullEnvironment* env = new QUESO::FullEnvironment(MPI_COMM_WORLD,argv[1],"",NULL);
#else
  QUESO::FullEnvironment* env = new QUESO::FullEnvironment(argv[1],"",NULL);
#endif

  // Compute
  compute(*env);

  // Finalize environment
  delete env;
#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif
  return 0;
}
