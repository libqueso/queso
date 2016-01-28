/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008-2013 The PECOS Development Team
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
 *--------------------------------------------------------------------------
 *
 * $Id$
 *
 * Brief description of this file:
 *
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#include <example_compute.h>

int main(int argc, char* argv[])
{
  // Initialize environment
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc,&argv);
  QUESO::FullEnvironment* env =
    new QUESO::FullEnvironment(MPI_COMM_WORLD,argv[1],"",NULL);
#else
  QUESO::FullEnvironment* env =
    new QUESO::FullEnvironment(argv[1],"",NULL);
#endif

  // Compute
#if 1
  compute(*env);
#else
  debug_hyst(*env);
#endif

  // Finalize environment
  delete env;
#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return 0;
}
