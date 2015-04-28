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

template<class V, class M>
class MyInterpolationBuilder : public QUESO::InterpolationSurrogateBuilder<V,M>
{
public:
  MyInterpolationBuilder( const QUESO::BoxSubset<V,M> & domain,
                          const std::vector<unsigned int>& n_points )
    : QUESO::InterpolationSurrogateBuilder<V,M>(domain,n_points)
  {};

  virtual ~MyInterpolationBuilder(){};

  virtual double evaluate_model( const V & domainVector ) const
  { queso_not_implemented();
    return 0.0; };
};

int main(int argc, char ** argv)
{
  MPI_Init(&argc, &argv);
  QUESO::FullEnvironment env(MPI_COMM_WORLD, "test_InterpolationSurrogate/queso_input.txt", "", NULL);

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

  MyInterpolationBuilder<QUESO::GslVector,QUESO::GslMatrix>
    builder( paramDomain, n_points );

  QUESO::LinearLagrangeInterpolationSurrogate<QUESO::GslVector,QUESO::GslMatrix>
    one_d_surrogate( builder );

  return return_flag;
}
