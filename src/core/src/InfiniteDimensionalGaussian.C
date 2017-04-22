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

#include <memory>
#include <vector>
#include <cmath>

#include <queso/InfiniteDimensionalMeasureBase.h>
#include <queso/InfiniteDimensionalGaussian.h>
#include <queso/Environment.h>
#include <queso/FunctionBase.h>
#include <queso/OperatorBase.h>
#include <queso/RngBase.h>

namespace QUESO {

InfiniteDimensionalGaussian::InfiniteDimensionalGaussian(
    const FullEnvironment& env,
    const FunctionBase & mean,
    const OperatorBase & precision,
    double alpha,
    double beta)
  : InfiniteDimensionalMeasureBase(),
    mean(mean),
    precision(precision),
    env(env),
    alpha(alpha),
    beta(beta)
{
  this->coeffs.resize(this->precision.get_num_converged(), 0.0);
}

InfiniteDimensionalGaussian::~InfiniteDimensionalGaussian()
{
}

SharedPtr<FunctionBase>::Type InfiniteDimensionalGaussian::draw()
{
  unsigned int i;

  for (i = 0; i < this->precision.get_num_converged(); i++) {
    (this->coeffs)[i] = env.rngObject()->gaussianSample(this->beta);
  }

  SharedPtr<FunctionBase>::Type f(this->precision.inverse_kl_transform(this->coeffs, this->alpha));

  // Add on the mean
  f->add(1.0, mean);

  return f;
}

double InfiniteDimensionalGaussian::get_kl_coefficient(unsigned int i) const
{
  // This is code repetition, but I'm not quite sure this belongs
  // in the operator class, because it's useful in the measure
  return (this->coeffs)[i] / std::pow(this->precision.get_eigenvalue(i), this->alpha / 2.0);
}

}  // End namespace QUESO
