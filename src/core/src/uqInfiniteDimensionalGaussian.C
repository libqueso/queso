//-----------------------------------------------------------------------bl-
//
//-----------------------------------------------------------------------el-
// 
// $Id: $
//
//--------------------------------------------------------------------------

#include <uqFunctionBase.h>
#include <uqOperatorBase.h>

uqInfiniteDimensionalGaussian::uqInfiniteDimensionalGaussian(const uqFunctionBase & mean, const uqOperatorBase & precision)
  : uqInfiniteDimensionalMeasureBase(),
    mean(mean),
    precision(precision)
{
  return;
}

uqInfiniteDimensionalGaussian::~uqInfiniteDimensionalGaussian()
{
  return;
}

void uqInfiniteDimensionalGaussian::draw()
{
  return;
}
