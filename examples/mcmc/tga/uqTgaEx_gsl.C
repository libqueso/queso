/* uq/examples/mcmc//uqTgaEx_gsl.C
 *
 * Copyright (C) 2008 The PECOS Team, http://queso.ices.utexas.edu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <uqTgaEx.h>
#include <uqApplRoutines.h>
#include <uqGslMatrix.h>

int main(int argc, char* argv[])
{
  //************************************************
  // Initialize environment
  //************************************************
  uqEnvironmentClass* env = new uqEnvironmentClass(argc,argv);

  //************************************************
  // Call application
  //************************************************
  uqAppl<uqGslVectorClass,uqGslMatrixClass>(*env);

  //************************************************
  // Finalize environment
  //************************************************
  delete env;

  return 0;
}
