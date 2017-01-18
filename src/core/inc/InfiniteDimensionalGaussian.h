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

#ifndef QUESO_INFINITEDIMENSIONALGAUSSIAN_H
#define QUESO_INFINITEDIMENSIONALGAUSSIAN_H

#include <queso/SharedPtr.h>
#include <queso/InfiniteDimensionalMeasureBase.h>
#include <queso/FunctionBase.h>

namespace QUESO {

class FullEnvironment;
class OperatorBase;

/*!
 * \file InfiniteDimensionalGaussian.h
 * \brief Class defining infinite dimensional Gaussian measures
 *
 * \class InfiniteDimensionalGaussian
 * \brief Class defining infinite dimensional Gaussian measures
 */

class InfiniteDimensionalGaussian : public InfiniteDimensionalMeasureBase {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Construct a Gaussian with mean \c mean, precision operator \c precision.  \c alpha controls regularity of draws and \c beta^2 is a multiplicative coefficient of \c precision.
  /*!
   * The argument \c alpha refers to the power to which the \c precision
   * will be taken. And \c beta is the multiplicative coefficient of
   * \c preicision. The complete measure is therefore:
   *
   * N(\c mean, \c beta^2 * \c pow(\c precision, \c alpha))
   *
   * It is expected that \c mean and \c precision will live longer than \c this
   */
  InfiniteDimensionalGaussian(const FullEnvironment& env,
      const FunctionBase &mean, const OperatorBase &precision,
      double alpha, double beta);

  //! Destructor
  ~InfiniteDimensionalGaussian();
  //@}

  //! Draw from the measure, and then return a shared pointer to the draw
  virtual SharedPtr<FunctionBase>::Type draw();

  //! Return coefficient \c i of the KL expansion of the current draw.  Must be called after draw()
  virtual double get_kl_coefficient(unsigned int i) const;

private:
  // Mean
  const FunctionBase & mean;

  // Precision -- I suppose you saw that one coming.
  const OperatorBase & precision;

  // QUESO environment
  const FullEnvironment& env;

  // Fractional power
  double alpha;

  // Multiplicative constant
  double beta;

  // The coefficients for the KL expansion
  // (not multiplied by the evals)
  std::vector<double> coeffs;
};

}  // End namespace QUESO

#endif // QUESO_INFINITEDIMENSIONALGAUSSIAN_H
