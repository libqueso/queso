/*-------------------------------------------------------------------
 *
 * Copyright (C) 2012 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of the QUESO Library (Quantification of Uncertainty
 * for Estimation, Simulation and Optimization).
 *
 * QUESO is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QUESO is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QUESO. If not, see <http://www.gnu.org/licenses/>.
 *
 *-------------------------------------------------------------------
 *
 * $Id$
 */
 /*------------------------------------------------------------------
 * Brief description of this file: 
 * 
 * This is an example of how to use QUESO classes and algorithms in
 * order to define and solve a statistical inverse problem (SIP) and/
 * or a statistical forward problem (SFP).
 * 
 * The SIP consists on calibrating the magnitude 'g' of acceleration
 * gravity using measurements of the time that it takes for an object
 * in free fall to reach the ground from a given height and zero 
 * initial velocity. The solution of the SIP is the posterior 
 * probability density function (PDF) of 'g'.
 * 
 * The SFP consists of calculating the maximum distance traveled by
 * an object in projectile motion. The posterior PDF of 'g' from the 
 * SIP might be used as input to the SFP.
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
  MPI_Init(&argc,&argv);
  uqFullEnvironmentClass* env =
    new uqFullEnvironmentClass(MPI_COMM_WORLD,argv[1],"",NULL);

  // Call application
  computeGravityAndTraveledDistance(*env);

  // Finalize QUESO environment
  delete env;
  MPI_Finalize();

  return 0;
}
