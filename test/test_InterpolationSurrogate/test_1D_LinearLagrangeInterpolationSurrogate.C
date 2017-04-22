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

double one_d_fn( double x );

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
    paramSpace(env,"param_", 1, NULL);

  double min_val = -5;
  double max_val = 3;

  QUESO::GslVector paramMins(paramSpace.zeroVector());
  paramMins.cwSet(min_val);
  QUESO::GslVector paramMaxs(paramSpace.zeroVector());
  paramMaxs.cwSet(max_val);

  QUESO::BoxSubset<QUESO::GslVector, QUESO::GslMatrix>
    paramDomain("param_", paramSpace, paramMins, paramMaxs);

  std::vector<unsigned int> n_points(1, 101);

  QUESO::InterpolationSurrogateData<QUESO::GslVector, QUESO::GslMatrix>
    data(paramDomain,n_points);

  std::vector<double> values(n_points[0]);

  double spacing = (max_val - min_val)/(n_points[0]-1);

  for( unsigned int i = 0; i < n_points[0]; i++ )
    {
      double x = min_val + i*spacing;
      values[i] = one_d_fn(x);
    }

  data.set_values( values );

  QUESO::LinearLagrangeInterpolationSurrogate<QUESO::GslVector,QUESO::GslMatrix>
    one_d_surrogate( data );

  QUESO::GslVector domainVector(paramSpace.zeroVector());
  domainVector.cwSet(1.221);

  double test_val = one_d_surrogate.evaluate(domainVector);

  double exact_val = one_d_fn(1.221);

  double tol = 2.0*std::numeric_limits<double>::epsilon();

  double rel_error = (test_val - exact_val)/exact_val;

  if( std::fabs(rel_error) > tol )
    {
      std::cerr << "ERROR: Tolerance exceeded for 1D Lagrange interpolation test."
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

double one_d_fn( double x )
{
  return 3.0 + 2.5*x;
}
