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

#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/InterpolationSurrogateIOASCII.h>
#include <queso/LinearLagrangeInterpolationSurrogate.h>

double four_d_fn( double x, double y, double z, double a );

int main(int argc, char ** argv)
{
  if( argc < 2 )
    {
      std::cerr << "ERROR: Must run program as '4d_interp <inputfile>'." << std::endl;
      exit(1);
    }
  std::string filename = argv[1];

  MPI_Init(&argc, &argv);
  QUESO::FullEnvironment env(MPI_COMM_WORLD, filename.c_str(), "", NULL);

  // We will read in the previously computed interpolation data
  QUESO::InterpolationSurrogateIOASCII<QUESO::GslVector, QUESO::GslMatrix>
    data_reader;

  data_reader.read( "./4d_interp_data_coarse.dat", env, "param_" );

  // Grab a reference to the data built in the reader
  const QUESO::InterpolationSurrogateData<QUESO::GslVector, QUESO::GslMatrix>&
    data = data_reader.data();

  // The reader read in the data, so now we can give the data
  // to the interpolation surrogate. This object can now be used in a likelihood
  // function for example. Here, we just illustrate calling the surrogate model
  // evaluation.
  QUESO::LinearLagrangeInterpolationSurrogate<QUESO::GslVector,QUESO::GslMatrix>
    four_d_surrogate( data );

  // A point in parameter space at which we want to use the surrogate
  QUESO::GslVector domainVector(data.get_paramDomain().vectorSpace().zeroVector());
  domainVector[0] = -0.4;
  domainVector[1] = 3.0;
  domainVector[2] = 1.5;
  domainVector[3] = 1.65;

  // Evaluate the surrogate model at the given point in parameter space
  // Because the exact function is quadrilinear, our interpolated value
  // should be exact.
  double value = four_d_surrogate.evaluate(domainVector);

  for( int r = 0; r < env.fullComm().NumProc(); r++ )
    {
      if( env.fullRank() == r )
        {
          std::cout << "======================================" << std::endl
                    << "Processor: " << env.fullRank() << std::endl
                    << "Interpolated value: " << value << std::endl
                    << "Exact value: " << four_d_fn( domainVector[0],domainVector[1],domainVector[2],domainVector[3])
                    << std::endl
                    << "======================================" << std::endl;
        }

      MPI_Barrier(MPI_COMM_WORLD);
    }

  MPI_Finalize();

  return 0;
}

double four_d_fn( double x, double y, double z, double a )
{
  return 3.0 + 2.5*x - 3.1*y + 2.71*z + 3.14*a
    + 0.1*x*y + 1.2*x*z + 0.5*y*z + 0.1*x*a + 1.1*y*a + 2.1*z*a
    + 0.3*x*y*a + 0.5*x*z*a + 1.2*y*z*a + 0.9*x*y*z
    + 2.5*x*y*z*a;
}
