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
#include <queso/InterpolationSurrogateData.h>

double four_d_fn( double x, double y, double z, double a );

template<class V, class M>
class MyInterpolationBuilder : public QUESO::InterpolationSurrogateBuilder<V,M>
{
public:
  MyInterpolationBuilder( QUESO::InterpolationSurrogateData<V,M>& data )
    : QUESO::InterpolationSurrogateBuilder<V,M>(data)
  {};

  virtual ~MyInterpolationBuilder(){};

  virtual double evaluate_model( const V & domainVector ) const
  { queso_assert_equal_to( domainVector.sizeGlobal(), 4);
    return four_d_fn(domainVector[0],domainVector[1],domainVector[2],domainVector[3]); };
};

int main(int argc, char ** argv)
{
  MPI_Init(&argc, &argv);
  QUESO::FullEnvironment env(MPI_COMM_WORLD, "test_InterpolationSurrogate/queso_input.txt", "", NULL);

  int return_flag = 0;

  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix>
    paramSpace(env,"param_", 4, NULL);

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

  QUESO::InterpolationSurrogateData<QUESO::GslVector, QUESO::GslMatrix>
    data(paramDomain,n_points);

  MyInterpolationBuilder<QUESO::GslVector,QUESO::GslMatrix>
    builder( data );

  builder.build_values();

  QUESO::LinearLagrangeInterpolationSurrogate<QUESO::GslVector,QUESO::GslMatrix>
    four_d_surrogate( data );

  QUESO::GslVector domainVector(paramSpace.zeroVector());
  domainVector[0] = -0.4;
  domainVector[1] = 3.0;
  domainVector[2] = 1.5;
  domainVector[3] = 1.65;

  double test_val = four_d_surrogate.evaluate(domainVector);

  double exact_val = four_d_fn(domainVector[0],domainVector[1],domainVector[2],domainVector[3]);

  double tol = 2.0*std::numeric_limits<double>::epsilon();

  double rel_error = (test_val - exact_val)/exact_val;

  if( std::fabs(rel_error) > tol )
    {
      std::cerr << "ERROR: Tolerance exceeded for 4D Lagrange interpolation test."
                << std::endl
                << " test_val  = " << test_val << std::endl
                << " exact_val = " << exact_val << std::endl
                << " rel_error = " << rel_error << std::endl
                << " tol       = " << tol << std::endl;

      return_flag = 1;
    }

  return return_flag;
}

double four_d_fn( double x, double y, double z, double a )
{
  return 3.0 + 2.5*x - 3.1*y + 2.71*z + 3.14*a
    + 0.1*x*y + 1.2*x*z + 0.5*y*z + 0.1*x*a + 1.1*y*a + 2.1*z*a
    + 0.3*x*y*a + 0.5*x*z*a + 1.2*y*z*a + 0.9*x*y*z
    + 2.5*x*y*z*a;
}
