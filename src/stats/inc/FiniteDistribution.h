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

#ifndef UQ_FINITE_DISTRIBUTION_H
#define UQ_FINITE_DISTRIBUTION_H

#include <map>
#include <vector>
#include <queso/Environment.h>

namespace QUESO {

/*! \file FiniteDistribution.h
 * \brief A templated class for a finite distribution.
 *
 * \class FiniteDistribution
 * \brief A templated class for a finite distribution.
 *
 * Unordered, discrete distribution, whose weights must be nonnegative, and are treated as unnormalized
 * probabilities.\n
 *
 * TODO: Describe me better!*/

class FiniteDistribution {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Constructor.
  FiniteDistribution(const BaseEnvironment& env,
                            const char*                   prefix,
                            const std::vector<double>&    inpWeights);
  //! Virtual destructor
  virtual ~FiniteDistribution();
  //@}

  //! @name Misc methods
  //@{
  //! Environment; access to protected attribute m_env.
  const BaseEnvironment& env    () const;
  //@}

  //! @name Statistical methods
  //@{
  //! Weights.
  const std::vector<double>&    weights() const;

  //! Samples.
  unsigned int            sample () const;
 //@}

protected:
  const BaseEnvironment& m_env;
        std::string             m_prefix;
	std::vector<double>     m_weights;

	std::map<double,unsigned int> m_map;
};

}  // End namespace QUESO

#endif // UQ_FINITE_DISTRIBUTION_H
