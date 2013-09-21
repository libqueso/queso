//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012 The PECOS Development Team
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
// 
// $Id: $
//
//--------------------------------------------------------------------------

#ifndef __QUESO_INFINITEDIMENSIONALGAUSSIAN__
#define __QUESO_INFINITEDIMENSIONALGAUSSIAN__

#include <boost/shared_ptr.hpp>

#include <queso/Environment.h>
#include <queso/uqOperatorBase.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/VectorSpace.h>
#include <queso/VectorSubset.h>
#include <queso/VectorRV.h>
#include <queso/uqInfiniteDimensionalMeasureBase.h>
#include <queso/uqFunctionBase.h>

namespace QUESO {

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
  virtual boost::shared_ptr<uqFunctionBase> draw() const;

  //! Return coefficient \c i of the KL expansion of the current draw
  /*!
   * You need to make a draw before you call this
   */
  virtual double get_kl_coefficient(unsigned int i) const;

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

  // The coefficients for the KL expansion
  // (not multiplied by the evals)
  std::vector<double> * coeffs;
};

}  // End namespace QUESO

#endif // __QUESO_INFINITEDIMENSIONALGAUSSIAN__
