//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
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

#ifndef QUESO_INFINITE_LIKELIHOOD_H
#define QUESO_INFINITE_LIKELIHOOD_H

namespace QUESO {

class FunctionBase;

/*!
 * \file InfiniteDimensionalLikelihoodBase.h
 * \brief Abstract class representing the likelihood
 *
 * \class InfiniteDimensionalLikelihoodBase
 * \brief Abstract class representing the likelihood.  Users must subclass this.
 *
 * This is the class you must subclass to evaluate your forward problem.  The
 * evaluation should be implemented in the evaluate() method.
 */

class InfiniteDimensionalLikelihoodBase
{
public:
  //! Construct a likelihood functional with \c obs_stddev as the observational error standard deviation
  InfiniteDimensionalLikelihoodBase(double obs_stddev);

  //! Destructor
  virtual ~InfiniteDimensionalLikelihoodBase();

  //! Set the observation standard deviation.
  virtual void set_obs_stddev(double stddev);

  //! Return the observation standard deviation.
  virtual double obs_stddev() const;

  //! Evaluate the likelihood functional at \c flow.  Subclasses must implement this method.
  virtual double evaluate(FunctionBase & flow) = 0;

private:
  double _obs_stddev;
};

}  // End namespace QUESO

#endif // QUESO_INFINITE_LIKELIHOOD_H
