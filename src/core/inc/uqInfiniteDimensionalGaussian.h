//-----------------------------------------------------------------------bl-
//
//-----------------------------------------------------------------------el-
// 
// $Id: $
//
//--------------------------------------------------------------------------

#ifndef __QUESO_INFINITEDIMENSIONALGAUSSIAN__
#define __QUESO_INFINITEDIMENSIONALGAUSSIAN__

#include <memory>
#include <uqGslVector.h>
#include <uqGslMatrix.h>
#include <uqVectorSpace.h>
#include <uqVectorSubset.h>
#include <uqVectorRV.h>
#include <uqInfiniteDimensionalMeasureBase.h>
#include <uqFunctionBase.h>

class uqFullEnvironmentClass;
class uqOperatorBase;

/*!
 * \file uqInfiniteDimensionalGaussian.h
 * \brief Class defining infinite dimensional Gaussian measures
 */

class uqInfiniteDimensionalGaussian : public uqInfiniteDimensionalMeasureBase {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Construct a Gaussian with mean \c mean, precision operator \c precision.
  /*!
   * The argument \c alpha refers to the power to which the \c precision
   * will be taken. And \c beta is the multiplicative coefficient of
   * \c preicision. The complete measure is therefore:
   *
   * N(\c mean, \c beta^2 * pow(\c precision, - \c alpha))
   *
   * It is expected that \c mean and \c precision will live longer than \c this
   */
  uqInfiniteDimensionalGaussian(const uqFullEnvironmentClass & env,
      const uqFunctionBase &mean, const uqOperatorBase &precision,
      double alpha, double beta);

  //! Destructor
  ~uqInfiniteDimensionalGaussian();
  //@}

  //! Draw from the measure, and store the result in the public member variable
  /*!
   * This updates the public memeber variable current draw
   */
  virtual std::auto_ptr<uqFunctionBase> draw() const;

private:
  // Mean
  const uqFunctionBase & mean;

  // Precision -- I suppose you saw that one coming.
  const uqOperatorBase & precision;

  // QUESO environment
  const uqFullEnvironmentClass & env;

  // Fractional power
  double alpha;

  // Multiplicative constant
  double beta;
};

#endif // __QUESO_INFINITEDIMENSIONALGAUSSIAN__
