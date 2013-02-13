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
    const uqOperatorBase & precision)
  : uqInfiniteDimensionalMeasureBase(),
    mean(mean),
    precision(precision),
    env(env)
{
  // this->_param_space = new uqVectorSpaceClass<uqGslVectorClass, uqGslMatrixClass>(
  //     env, "inf_measure_space", precision.get_num_converged(), NULL);
  // this->_mins = new uqGslVectorClass(this->_param_space->zeroVector());
  // this->_maxs = new uqGslVectorClass(this->_param_space->zeroVector());
  // this->_mins->cwSet(-INFINITY);
  // this->_maxs->cwSet(INFINITY);
  // this->_param_domain = new uqBoxSubsetClass<uqGslVectorClass, uqGslMatrixClass>(
  //     "inf_measure_domain", *this->_param_space, *this->_mins, *this->_maxs);
  // this->_mean = new uqGslVectorClass(this->_param_space->zeroVector());
  // this->_vars = new uqGslVectorClass(this->_param_space->zeroVector());
  // this->_vars->cwSet(1.0);
  // this->_coeffs_rv = new uqGaussianVectorRVClass<uqGslVectorClass, uqGslMatrixClass>(
  //     "inf_measure", *this->_param_domain, *this->_mean, *this->_vars);
  // this->_coeffs = new uqGslVectorClass(*this->_param_space->zeroVector());
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
    // coeffs[i] = env.rngObject()->gaussianSample(1.0);
    coeffs[i] = gsl_ran_gaussian(env.rng(), 1.0);
  }

  return precision.inverse_kl_transform(coeffs);
}
