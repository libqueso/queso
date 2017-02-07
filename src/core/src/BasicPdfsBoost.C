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

#include <queso/BasicPdfsBoost.h>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/beta.hpp>

namespace QUESO {

//! Constructor ---------------------------
BasicPdfsBoost::BasicPdfsBoost(int worldRank)
  :
  BasicPdfsBase(worldRank)
{
}

// Destructor ---------------------------------------
BasicPdfsBoost::~BasicPdfsBoost()
{
  //this function does nothing
}

// --------------------------------------------------
double
BasicPdfsBoost::betaPdfActualValue(double x, double alpha, double beta) const
{
  boost::math::beta_distribution<> beta_dist(alpha, beta);
  return boost::math::pdf(beta_dist, x);
}

// --------------------------------------------------
double
BasicPdfsBoost::gammaPdfActualValue(double x, double a, double b) const
{
  boost::math::gamma_distribution<> gamma_dist(a, b);
  return boost::math::pdf(gamma_dist, x);
}

}  // End namespace QUESO
