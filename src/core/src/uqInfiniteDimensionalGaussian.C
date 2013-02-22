//-----------------------------------------------------------------------bl-
//
//-----------------------------------------------------------------------el-
// 
// $Id: $
//
//--------------------------------------------------------------------------

#include <memory>
#include <vector>

#include <uqInfiniteDimensionalMeasureBase.h>
#include <uqInfiniteDimensionalGaussian.h>
#include <uqEnvironment.h>
#include <uqFunctionBase.h>
#include <uqOperatorBase.h>

using namespace std;

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

std::auto_ptr<uqFunctionBase> uqInfiniteDimensionalGaussian::draw() const
{
  unsigned int i;

  for (i = 0; i < this->precision.get_num_converged(); i++) {
    (*this->coeffs)[i] = env.rngObject()->gaussianSample(this->beta);
  }

  return this->precision.inverse_kl_transform(*this->coeffs, this->alpha);
}

double uqInfiniteDimensionalGaussian::get_kl_coefficient(unsigned int i) const
{
  // This is code repetition, but I'm not quite sure this belongs
  // in the operator class, because it's useful in the measure
  return (*this->coeffs)[i] / pow(this->precision.get_eigenvalue(i), this->alpha / 2.0);
}
