//-----------------------------------------------------------------------bl-
//
//-----------------------------------------------------------------------el-
// 
// $Id: $
//
//--------------------------------------------------------------------------

#ifndef __QUESO_INFINITEDIMENSIONALGAUSSIAN__
#define __QUESO_INFINITEDIMENSIONALGAUSSIAN__

/*!
 * \file uqInfiniteDimensionalGaussian.h
 * \brief Class defining infinite dimensional Gaussian measures
 */

class uqOperatorBase;
class uqFunctionBase;

class uqInfiniteDimensionalGaussian {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Construct a Gaussian with mean \c mean and precision operator \c precision.
  //! It is expected that \c mean and \c precision will live longer than \c this.
  uqInfiniteDimensionalGaussian(const uqFunctionBase &mean, const uqOperatorBase &precision);

  //! Destructor
  ~uqInfiniteDimensionalGaussian();
  //@}

  //! Draw from the measure, and store the result in the public member varaible. This
  //! updates the public memeber variable current draw
  virtual void draw();

private:
  // Mean
  const uqFunctionBase & mean;

  // Precision -- I suppose you saw that one coming.
  const uqOperatorBase & precision;
};

#endif // __QUESO_INFINITEDIMENSIONALGAUSSIAN__
