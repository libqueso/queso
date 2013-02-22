//-----------------------------------------------------------------------bl-
//
//-----------------------------------------------------------------------el-
// 
// $Id: $
//
//--------------------------------------------------------------------------

#ifndef __QUESO_INFINITEDIMENSIONALMEASURE_BASE__
#define __QUESO_INFINITEDIMENSIONALMEASURE_BASE__

#include <memory>
#include <uqFunctionBase.h>

/*!
 * \file uqInfiniteDimensionalMeasureBase.h
 * \brief Abstract base class for infinite dimensional measures
 */

class uqFunctionBase;

class uqInfiniteDimensionalMeasureBase {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor
  uqInfiniteDimensionalMeasureBase();

  //! Destructor
  virtual ~uqInfiniteDimensionalMeasureBase();
  //@}

  //! Draw from the measure, and store the result in the public member varaible. This
  //! updates the public memeber variable current draw
  virtual std::auto_ptr<uqFunctionBase> draw() const = 0;

  //! Return the coefficient \c i of the KL expansion of the current draw
  /*!
   * You need to make a draw before you call this
   */
  virtual double get_kl_coefficient(unsigned int i) const = 0;
};

#endif // __QUESO_INFINITEDIMENSIONALMEASURE_BASE__
