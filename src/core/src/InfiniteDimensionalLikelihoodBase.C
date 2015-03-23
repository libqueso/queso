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

#include <queso/InfiniteDimensionalLikelihoodBase.h>

namespace QUESO {

InfiniteDimensionalLikelihoodBase::InfiniteDimensionalLikelihoodBase(
    double obs_stddev)
  : _obs_stddev(obs_stddev)
{
}

InfiniteDimensionalLikelihoodBase::~InfiniteDimensionalLikelihoodBase()
{
}

void InfiniteDimensionalLikelihoodBase::set_obs_stddev(double stddev)
{
  this->_obs_stddev = stddev;
}

double InfiniteDimensionalLikelihoodBase::obs_stddev() const
{
  return this->_obs_stddev;
}

}  // End namespace QUESO
