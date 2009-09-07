/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008,2009 The PECOS Development Team
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
 * Simple toy problem to illustrate the use of the Basic QUESO
 * interface to perform a statistical inverse problem.
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#include <stdlib.h>
#include <stdio.h>
#include <queso.h>
#include <mpi.h> // Added by prudenci on 2009/Sep/06 in order to avoid INTEL warning. Maybe Karl want to add that in queso.h (a hpct file?)

double my_likelihood(double *params);

int main(int argc, char *argv[])
{

  MPI_Init(&argc,&argv);

  printf("--> Initializing QUESO Environment...\n");

  /* Initialize QUESO environment and define input file */

  QUESO_init("queso.inp");	                   

 /* Register an application likelihood function with QUESO and run a
  * statistical inversion problem using MCMC sampling */

  QUESO_statistical_inversion(my_likelihood);     

  /* Finalize the analysis and output statistics */

  QUESO_finalize();
  
  MPI_Finalize();
  return 0;

}

double my_likelihood(double *params)
{
  static double needle   = 42.;   
  static double variance = 1.*1.;  /* sigma^2 */

  return( (params[0] - needle)*(params[0] - needle)/(variance) );
}
