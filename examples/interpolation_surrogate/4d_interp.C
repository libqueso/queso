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
#include <queso/BoxSubset.h>
#include <queso/LinearLagrangeInterpolationSurrogate.h>
#include <queso/InterpolationSurrogateBuilder.h>
#include <queso/InterpolationSurrogateIOASCII.h>

double four_d_fn_1( double x, double y, double z, double a );
double four_d_fn_2( double x, double y, double z, double a );

/* We subclass InterpolationSurrogateBuilder to implement our model evaluation.
   In this case, it just calls the simple functions, but you may do much more
   complicated things. */
template<class V, class M>
class MyInterpolationBuilder : public QUESO::InterpolationSurrogateBuilder<V,M>
{
public:
  MyInterpolationBuilder( QUESO::InterpolationSurrogateDataSet<V,M>& data )
    : QUESO::InterpolationSurrogateBuilder<V,M>(data)
  {};

  virtual ~MyInterpolationBuilder(){};

  virtual void evaluate_model( const V & domainVector, std::vector<double>& values )
  { queso_assert_equal_to( domainVector.sizeGlobal(), 4);
    queso_assert_equal_to( values.size(), 2);
    values[0] = four_d_fn_1(domainVector[0],domainVector[1],domainVector[2],domainVector[3]);
    values[1] = four_d_fn_2(domainVector[0],domainVector[1],domainVector[2],domainVector[3]);
  };
};

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

  // Define the parameter space. It's 4-dimensional in this example.
  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix>
    paramSpace(env,"param_", 4, NULL);

  // Define parameter bounds
  QUESO::GslVector paramMins(paramSpace.zeroVector());
  paramMins[0] = -1;
  paramMins[1] = -0.5;
  paramMins[2] = 1.1;
  paramMins[3] = -2.1;

  QUESO::GslVector paramMaxs(paramSpace.zeroVector());
  paramMaxs[0] = 0.9;
  paramMaxs[1] = 3.14;
  paramMaxs[2] = 2.1;
  paramMaxs[3] = 4.1;

  // Define parameter domain.
  QUESO::BoxSubset<QUESO::GslVector, QUESO::GslMatrix>
    paramDomain("param_", paramSpace, paramMins, paramMaxs);

  // How many points to use in each coordinate direction of parameter space.
  // Since we're using an interpolation surrogate, it will evenly space the points
  // in each coordinate direction. These are the points at which the model
  // will be evaluated.
  std::vector<unsigned int> n_points(4);
  n_points[0] = 101;
  n_points[1] = 51;
  n_points[2] = 31;
  n_points[3] = 41;

  // Construct data object. We have two functions we're interpolating
  // at the same time, so n_datasets = 2.
  unsigned int n_datasets = 2;
  QUESO::InterpolationSurrogateDataSet<QUESO::GslVector, QUESO::GslMatrix>
    data(paramDomain,n_points,n_datasets);

  // Construct builder. The builder will add the values from the model evaluations
  // so the data object must stay alive.
  MyInterpolationBuilder<QUESO::GslVector,QUESO::GslMatrix>
    builder( data );

  // The expensive part. The builder will now evaluate the model for all the
  // desired points in parameter space. This will build both interpolants.
  builder.build_values();

  // Now that we've built the data, we write it out so we can reuse it later
  QUESO::InterpolationSurrogateIOASCII<QUESO::GslVector, QUESO::GslMatrix>
    data_writer;

  // Write each datasets separately. We'll be able to read each one in
  // individually and construct separate LinearLagrangeInterpolationSurrogate
  // objects. If we have more than 1 processor at our disposal, have the
  // second data set written out by someone other than processor 0 since
  // processor 0 will be writing out the first data set.
  data_writer.write( "4d_interp_data_1.dat", data.get_dataset(0) );

  if( env.fullComm().NumProc() > 1 )
    data_writer.write( "4d_interp_data_2.dat", data.get_dataset(1), 1 );
  else
    data_writer.write( "4d_interp_data_2.dat", data.get_dataset(1) );

  // The builder put the model values into the data, so now we can give the data
  // to the interpolation surrogate. This object can now be used in a likelihood
  // function for example. Here, we just illustrate calling the surrogate model
  // evaluation.
  QUESO::LinearLagrangeInterpolationSurrogate<QUESO::GslVector,QUESO::GslMatrix>
    four_d_surrogate_1( data.get_dataset(0) );

  QUESO::LinearLagrangeInterpolationSurrogate<QUESO::GslVector,QUESO::GslMatrix>
    four_d_surrogate_2( data.get_dataset(1) );

  // A point in parameter space at which we want to use the surrogate
  QUESO::GslVector domainVector(paramSpace.zeroVector());
  domainVector[0] = -0.4;
  domainVector[1] = 3.0;
  domainVector[2] = 1.5;
  domainVector[3] = 1.65;

  // Evaluate the surrogate model at the given point in parameter space
  // Because the exact function is quadrilinear, our interpolated value
  // should be exact.
  double value_1 = four_d_surrogate_1.evaluate(domainVector);
  double value_2 = four_d_surrogate_2.evaluate(domainVector);


  for( int r = 0; r < env.fullComm().NumProc(); r++ )
    {
      if( env.fullRank() == r )
        {
          std::cout << "======================================" << std::endl
                    << "Processor: " << env.fullRank() << std::endl
                    << "Interpolated value for model 1: " << value_1 << std::endl
                    << "Exact value for model 1: " << four_d_fn_1( domainVector[0],domainVector[1],domainVector[2],domainVector[3])
                    << std::endl
                    << "Interpolated value for model 2: " << value_2 << std::endl
                    << "Exact value for model 2: " << four_d_fn_2( domainVector[0],domainVector[1],domainVector[2],domainVector[3])
                    << std::endl
                    << "======================================" << std::endl;
        }

      MPI_Barrier(MPI_COMM_WORLD);
    }

  MPI_Finalize();

  return 0;
}

double four_d_fn_1( double x, double y, double z, double a )
{
  return 3.0 + 2.5*x - 3.1*y + 2.71*z + 3.14*a
    + 0.1*x*y + 1.2*x*z + 0.5*y*z + 0.1*x*a + 1.1*y*a + 2.1*z*a
    + 0.3*x*y*a + 0.5*x*z*a + 1.2*y*z*a + 0.9*x*y*z
    + 2.5*x*y*z*a;
}

double four_d_fn_2( double x, double y, double z, double a )
{
  return 2.0 + 1.5*x - 2.1*y + 1.71*z + 2.14*a
    + 0.01*x*y + 0.2*x*z + 0.05*y*z + 0.01*x*a + 2.1*y*a + 3.1*z*a
    + 1.3*x*y*a + 1.5*x*z*a + 2.2*y*z*a + 1.9*x*y*z
    + 3.5*x*y*z*a;
}
