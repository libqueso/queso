//-----------------------------------------------------------------------bl-
//
//-----------------------------------------------------------------------el-

#ifndef __QUESO_INFINITEDIMENSIONALMEASURE_BASE__
#define __QUESO_INFINITEDIMENSIONALMEASURE_BASE__

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
  ~uqInfiniteDimensionalMeasureBase();
  //@}

  //! Draw from the measure, and store the result in the public member varaible. This
  //! updates the public memeber variable current draw
  virtual void draw() = 0;

  //! Stores the most recent draw from the measure
  uqFunctionBase *currentDraw;
};

#endif // __QUESO_INFINITEDIMENSIONALMEASURE_BASE__
