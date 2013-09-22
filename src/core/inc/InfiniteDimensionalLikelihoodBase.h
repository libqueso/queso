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
// 
// $Id: $
//
//--------------------------------------------------------------------------

#ifndef __QUESO_INFINITE_LIKELIHOOD_H
#define __QUESO_INFINITE_LIKELIHOOD_H

// #include <libmesh/system.h>

namespace QUESO {

class uqFunctionBase;

// namespace libMesh {
//   class EquationSystems;
// }

// class ForwardSolver;

/*!
 * \brief Abstract class representing the likelihood
 */
class uqInfiniteDimensionalLikelihoodBase
{
public:
  /*!
   * Constructor
   */
  uqInfiniteDimensionalLikelihoodBase(double obs_stddev);

  /*!
   * Destructor
   */
  virtual ~uqInfiniteDimensionalLikelihoodBase();

  /*!
   * Set the observation standard deviation.  Default is 1.
   */
  virtual void set_obs_stddev(double stddev);

  /*!
   * Get the observation standard deviation.
   */
  virtual double obs_stddev() const;

  virtual double evaluate(uqFunctionBase & flow) = 0;

// protected:
//   ForwardSolver & fwd_solver;
//   libMesh::EquationSystems * equation_systems;

private:
  double _obs_stddev;
};

}  // End namespace QUESO

#endif // __QUESO_INFINITE_LIKELIHOOD_H
