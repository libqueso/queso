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

#include <queso/InterpolationSurrogateHelper.h>

// C++
#include <iostream>

int main()
{
  std::vector<unsigned int> n_points(5);

  n_points[0] = 11;
  n_points[1] = 21;
  n_points[2] = 31;
  n_points[3] = 41;
  n_points[4] = 51;

  std::vector<unsigned int> indices(5);
  indices[0] = 1;
  indices[1] = 2;
  indices[2] = 3;
  indices[3] = 4;
  indices[4] = 5;

  unsigned int global_exact =
    indices[0] +
    indices[1]*n_points[0] +
    indices[2]*n_points[0]*n_points[1] +
    indices[3]*n_points[0]*n_points[1]*n_points[2] +
    indices[4]*n_points[0]*n_points[1]*n_points[2]*n_points[3];

  int return_flag = 0;

  unsigned int global = QUESO::InterpolationSurrogateHelper::coordToGlobal( indices, n_points );

  // This is integer arithmetic so it should be exactly zero error
  if( (global - global_exact) != 0 )
    {
      std::cerr << "ERROR: mismatch in InterpolationSurrogateHelper::coordToGlobal test." << std::endl
                << "       test  = " << global << std::endl
                << "       exact = " << global_exact << std::endl;
      return_flag = 1;
    }

  std::vector<unsigned int> indices_test;
  QUESO::InterpolationSurrogateHelper::globalToCoord( global, n_points, indices_test );

  unsigned int test_indices = 0;
  for( unsigned int d = 0; d < 5; d++ )
    {
      if( indices[d] != indices_test[d] )
        test_indices = 1;
    }

  if( test_indices == 1 )
    {
      std::cerr << "ERROR: mismatch in InterpolationSurrogateHelper::globalToCoord test." << std::endl
                << "       test  = ";
      for( unsigned int d = 0; d < 5; d++ )
        {
           std::cerr << indices_test[d] << " ";
        }
      std::cerr << std::endl;

      std::cerr << "       exact = ";
      for( unsigned int d = 0; d < 5; d++ )
        {
          std::cerr << indices[d] << " ";
        }
      std::cerr << std::endl;

      return_flag = 1;
    }

  return return_flag;
}
