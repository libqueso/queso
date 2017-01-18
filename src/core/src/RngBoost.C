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

#include <queso/RngBoost.h>

#include <boost/random.hpp>
#include <boost/math/distributions.hpp>

namespace QUESO {

//! Constructor with seed ---------------------------
RngBoost::RngBoost(int seed, int worldRank)
  :
  RngBase(seed,worldRank)
  //m_rng(seed)
{
  resetSeed(seed);
//TODO Find a suitable test for here; Kemelli todo

//                       m_worldRank,
//                       "RngBoos::constructor()",
//                       "Kemelli todo: boost rng");
}

// Destructor ---------------------------------------
RngBoost::~RngBoost()
{
  //this function does nothing
}

// Sampling methods ---------------------------------
void
RngBoost::resetSeed(int newSeed)
{
  m_rng.seed(newSeed);

  return;
}

// --------------------------------------------------
double
RngBoost::uniformSample() const
{
  static boost::uniform_01<boost::mt19937> zeroone(m_rng);
  return zeroone();
}

// --------------------------------------------------
double
RngBoost::gaussianSample(double stdDev) const
{
  double mean = 0.; //it will be added conveniently later
  static boost::uniform_01<boost::mt19937> zeroone(m_rng);
  boost::math::normal_distribution<double>  gaussian_dist(mean, stdDev);
  return quantile(gaussian_dist, zeroone());
}

// --------------------------------------------------
double
RngBoost::betaSample(double alpha, double beta) const
{
  static boost::uniform_01<boost::mt19937> zeroone(m_rng);
  boost::math::beta_distribution<double> beta_dist(alpha, beta);
  return quantile(beta_dist, zeroone());
}

// --------------------------------------------------
double
RngBoost::gammaSample(double a, double b) const
{
  static boost::uniform_01<boost::mt19937> zeroone(m_rng);
  boost::math::gamma_distribution<double>  gamma_dist(a,b);
  return quantile(gamma_dist, zeroone());
}

}  // End namespace QUESO
