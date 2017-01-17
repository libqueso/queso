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

#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/BoxSubset.h>
#include <queso/LinearLagrangeInterpolationSurrogate.h>
#include <queso/InterpolationSurrogateData.h>

#include <cstdlib>
#include <limits>

double three_d_fn( double x, double y, double z );

int main(int argc, char ** argv)
{
  std::string inputFileName = "test_InterpolationSurrogate/queso_input.txt";
  const char * test_srcdir = std::getenv("srcdir");
  if (test_srcdir)
    inputFileName = test_srcdir + ('/' + inputFileName);

#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);
  QUESO::FullEnvironment env(MPI_COMM_WORLD, inputFileName, "", NULL);
#else
  QUESO::FullEnvironment env(inputFileName, "", NULL);
#endif

  int return_flag = 0;

  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix>
    paramSpace(env,"param_", 3, NULL);

  QUESO::GslVector paramMins(paramSpace.zeroVector());
  paramMins[0] = -1;
  paramMins[1] = -0.5;
  paramMins[2] = 1.1;

  QUESO::GslVector paramMaxs(paramSpace.zeroVector());
  paramMaxs[0] = 0.9;
  paramMaxs[1] = 3.14;
  paramMaxs[2] = 2.1;

  QUESO::BoxSubset<QUESO::GslVector, QUESO::GslMatrix>
    paramDomain("param_", paramSpace, paramMins, paramMaxs);

  std::vector<unsigned int> n_points(3);
  n_points[0] = 101;
  n_points[1] = 51;
  n_points[2] = 31;

  QUESO::InterpolationSurrogateData<QUESO::GslVector, QUESO::GslMatrix>
    data(paramDomain,n_points);

  std::vector<double> values(n_points[0]*n_points[1]*n_points[2]);

  double spacing_x = (paramMaxs[0] - paramMins[0])/(n_points[0]-1);
  double spacing_y = (paramMaxs[1] - paramMins[1])/(n_points[1]-1);
  double spacing_z = (paramMaxs[2] - paramMins[2])/(n_points[2]-1);

  for( unsigned int i = 0; i < n_points[0]; i++ )
    {
      for( unsigned int j = 0; j < n_points[1]; j++ )
        {
          for( unsigned int k = 0; k < n_points[2]; k++ )
            {
              unsigned int n = i + j*n_points[0] + k*n_points[0]*n_points[1];

              double x = paramMins[0] + i*spacing_x;
              double y = paramMins[1] + j*spacing_y;
              double z = paramMins[2] + k*spacing_z;

              values[n] = three_d_fn(x,y,z);
            }
        }
    }

  data.set_values( values );

  QUESO::LinearLagrangeInterpolationSurrogate<QUESO::GslVector,QUESO::GslMatrix>
    three_d_surrogate( data );

  QUESO::GslVector domainVector(paramSpace.zeroVector());
  domainVector[0] = -0.4;
  domainVector[1] = 3.0;
  domainVector[2] = 1.5;

  double test_val = three_d_surrogate.evaluate(domainVector);

  double exact_val = three_d_fn(domainVector[0],domainVector[1],domainVector[2]);

  double tol = 2.0*std::numeric_limits<double>::epsilon();

  double rel_error = (test_val - exact_val)/exact_val;

  if( std::fabs(rel_error) > tol )
    {
      std::cerr << "ERROR: Tolerance exceeded for 3D Lagrange interpolation test."
                << std::endl
                << " test_val  = " << test_val << std::endl
                << " exact_val = " << exact_val << std::endl
                << " rel_error = " << rel_error << std::endl
                << " tol       = " << tol << std::endl;

      return_flag = 1;
    }

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif
  return return_flag;
}

double three_d_fn( double x, double y, double z )
{
  return 3.0 + 2.5*x - 3.1*y + 2.71*z + 0.1*x*y + 1.2*x*z + 0.5*y*z + 2.5*x*y*z;
}
