//-----------------------------------------------------------------------bl-
//
//-----------------------------------------------------------------------el-
// 
// $Id: $
//
//--------------------------------------------------------------------------

#include <memory>

#include <uqInfiniteDimensionalMeasureBase.h>
#include <uqInfiniteDimensionalGaussian.h>
#include <uqEnvironment.h>
#include <uqFunctionBase.h>
#include <uqOperatorBase.h>

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
}

uqInfiniteDimensionalGaussian::~uqInfiniteDimensionalGaussian()
{
  return;
}

std::auto_ptr<uqFunctionBase> uqInfiniteDimensionalGaussian::draw() const
{
  unsigned int i;
  std::vector<double> coeffs(precision.get_num_converged(), 0.0);

  for (i = 0; i < precision.get_num_converged(); i++) {
    // Probably a better way to do this, using env.rngObject() perhaps?
    coeffs[i] = env.rngObject()->gaussianSample(this->beta);
  }

  return precision.inverse_kl_transform(coeffs, this->alpha);
}
