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
#include <queso/InterpolationSurrogateBuilder.h>
#include <queso/InterpolationSurrogateDataSet.h>
#include <queso/InterpolationSurrogateIOASCII.h>

#include <cstdlib>
#include <limits>

double four_d_fn_1( double x, double y, double z, double a );
double four_d_fn_2( double x, double y, double z, double a );

int test_val(double test_val, double exact_val, double tol, const std::string& test_name);

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
    queso_assert_equal_to( values.size(), 2 );
    values[0] = four_d_fn_1(domainVector[0],domainVector[1],domainVector[2],domainVector[3]);
    values[1] = four_d_fn_2(domainVector[0],domainVector[1],domainVector[2],domainVector[3]);
  };
};

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

  std::string vs_prefix = "param_";

  // Filename for writing/reading surrogate data
  std::string filename1 = "test_write_InterpolationSurrogateBuilder_1.dat";
  std::string filename2 = "test_write_InterpolationSurrogateBuilder_2.dat";

  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix>
      paramSpace(env,vs_prefix.c_str(), 4, NULL);

  // Point at which we will test the surrogate evaluation
  QUESO::GslVector domainVector(paramSpace.zeroVector());
  domainVector[0] = -0.4;
  domainVector[1] = 3.0;
  domainVector[2] = 1.5;
  domainVector[3] = 1.65;

  double exact_val_1 = four_d_fn_1(domainVector[0],domainVector[1],domainVector[2],domainVector[3]);
  double exact_val_2 = four_d_fn_2(domainVector[0],domainVector[1],domainVector[2],domainVector[3]);

  double tol = 2.0*std::numeric_limits<double>::epsilon();

  // First test surrogate build directly from the computed values
  {
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

    QUESO::BoxSubset<QUESO::GslVector, QUESO::GslMatrix>
      paramDomain("param_", paramSpace, paramMins, paramMaxs);

    std::vector<unsigned int> n_points(4);
    n_points[0] = 11;
    n_points[1] = 51;
    n_points[2] = 31;
    n_points[3] = 41;

    // One dataset for each of the two functions
    const unsigned int n_datasets = 2;

    QUESO::InterpolationSurrogateDataSet<QUESO::GslVector, QUESO::GslMatrix>
      data(paramDomain,n_points,n_datasets);

    MyInterpolationBuilder<QUESO::GslVector,QUESO::GslMatrix>
      builder( data );

    builder.build_values();

    QUESO::LinearLagrangeInterpolationSurrogate<QUESO::GslVector,QUESO::GslMatrix>
      four_d_surrogate_1( data.get_dataset(0) );

    QUESO::LinearLagrangeInterpolationSurrogate<QUESO::GslVector,QUESO::GslMatrix>
      four_d_surrogate_2( data.get_dataset(1) );

    double test_val_1 = four_d_surrogate_1.evaluate(domainVector);
    double test_val_2 = four_d_surrogate_2.evaluate(domainVector);

    return_flag  = return_flag ||
      test_val( test_val_1, exact_val_1, tol, "test_build_1" ) ||
      test_val( test_val_2, exact_val_2, tol, "test_build_2" );

    // Write the output to test reading next
    QUESO::InterpolationSurrogateIOASCII<QUESO::GslVector,QUESO::GslMatrix>
      data_writer;

    data_writer.write( filename1, data.get_dataset(0) );
    data_writer.write( filename2, data.get_dataset(1) );
  }

  // Now read the data and test
  {
    QUESO::InterpolationSurrogateIOASCII<QUESO::GslVector,QUESO::GslMatrix>
      data_reader_1, data_reader_2;

    data_reader_1.read( filename1, env, vs_prefix.c_str() );
    data_reader_2.read( filename2, env, vs_prefix.c_str() );


    // Build a new surrogate
    QUESO::LinearLagrangeInterpolationSurrogate<QUESO::GslVector,QUESO::GslMatrix>
      four_d_surrogate_1( data_reader_1.data() );

    QUESO::LinearLagrangeInterpolationSurrogate<QUESO::GslVector,QUESO::GslMatrix>
      four_d_surrogate_2( data_reader_2.data() );

    double test_val_1 = four_d_surrogate_1.evaluate(domainVector);
    double test_val_2 = four_d_surrogate_2.evaluate(domainVector);

    return_flag  = return_flag ||
      test_val( test_val_1, exact_val_1, tol, "test_read_1" ) ||
      test_val( test_val_2, exact_val_2, tol, "test_read_2" );
  }

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif
  return return_flag;
}

int test_val( double test_val, double exact_val, double tol, const std::string& test_name )
{
  int return_flag = 0;

  double rel_error = (test_val - exact_val)/exact_val;

  if( std::fabs(rel_error) > tol )
    {
      std::cerr << "ERROR: Tolerance exceeded for "+test_name
                << std::endl
                << " test_val  = " << test_val << std::endl
                << " exact_val = " << exact_val << std::endl
                << " rel_error = " << rel_error << std::endl
                << " tol       = " << tol << std::endl;

      return_flag = 1;
    }

  return return_flag;
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
