//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012 The PECOS Development Team
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
// 
// $Id: $
//
//--------------------------------------------------------------------------

#include <memory>
#include <vector>

#include <boost/shared_ptr.hpp>

#include <queso/uqInfiniteDimensionalMeasureBase.h>
#include <queso/uqInfiniteDimensionalGaussian.h>
#include <queso/Environment.h>
#include <queso/uqFunctionBase.h>
#include <queso/uqOperatorBase.h>

using namespace std;

namespace QUESO {

uqInfiniteDimensionalGaussian::uqInfiniteDimensionalGaussian(
    const uqFullEnvironmentClass & env,
    const uqFunctionBase & mean,
    const uqOperatorBase & precision,
    double alpha,
    double beta)
  : uqInfiniteDimensionalMeasureBase(),
    mean(mean),
    precision(precision),
    env(env),
    alpha(alpha),
    beta(beta)
{
  this->coeffs = new vector<double>(this->precision.get_num_converged(), 0.0);
}

uqInfiniteDimensionalGaussian::~uqInfiniteDimensionalGaussian()
{
  delete this->coeffs;
}

boost::shared_ptr<uqFunctionBase> uqInfiniteDimensionalGaussian::draw() const
{
  unsigned int i;

  for (i = 0; i < this->precision.get_num_converged(); i++) {
    (*this->coeffs)[i] = env.rngObject()->gaussianSample(this->beta);
  }

  boost::shared_ptr<uqFunctionBase> f(this->precision.inverse_kl_transform(*this->coeffs, this->alpha));
  return f;
}

double uqInfiniteDimensionalGaussian::get_kl_coefficient(unsigned int i) const
{
  // This is code repetition, but I'm not quite sure this belongs
  // in the operator class, because it's useful in the measure
  return (*this->coeffs)[i] / pow(this->precision.get_eigenvalue(i), this->alpha / 2.0);
}

}  // End namespace QUESO
