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

 /*------------------------------------------------------------------
 * Brief description of this file:
 *
 * This is an example of how to use QUESO classes and algorithms in order to define and solve
 * a statistical inverse problem (SIP) and/or a statistical forward problem (SFP).
 * The SIP consists on calibrating the magnitude 'g' of acceleration gravity using
 * measurements of the time that it takes for an object in free fall to reach the ground from
 * a given height and zero initial velocity. The solution of the SIP is the posterior
 * probability density function (PDF) of 'g'.
 * The SFP consists of calculating the maximum distance traveled by an object in projectile
 * motion. The posterior PDF of 'g' from the SIP might be used as input to the SFP.
 *
 * The code consists of 7 files:
 * - 'gravity_main.C' (this file)
 * - 'gravity_compute.C' (the driving application code)
 * - 'gravity_compute.h'
 * - 'gravity_likelihood.C' (necessary for the SIP)
 * - 'gravity_likelihood.h'
 * - 'gravity_qoi.C' (necessary for the SFP)
 * - 'gravity_qoi.h'
 *-----------------------------------------------------------------*/

#include <gravity_compute.h>

int main(int argc, char* argv[])
{
  // Initialize QUESO environment
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc,&argv);
  QUESO::FullEnvironment* env =
    new QUESO::FullEnvironment(MPI_COMM_WORLD,argv[1],"",NULL);
#else
  QUESO::FullEnvironment* env =
    new QUESO::FullEnvironment(argv[1],"",NULL);
#endif

  // Call application
  computeGravityAndTraveledDistance(*env);

  // Finalize QUESO environment
  delete env;
#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return 0;
}
