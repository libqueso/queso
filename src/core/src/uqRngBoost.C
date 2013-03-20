//-----------------------------------------------------------------------bl-
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
// $Id$
//
//--------------------------------------------------------------------------

#include <uqRngBoost.h>
#include <mpi.h>

uqRngBoostClass::uqRngBoostClass()
  :
  uqRngBaseClass()
{
  UQ_FATAL_TEST_MACRO(true,
                      m_worldRank,
                      "uqRngBoostClass::constructor(), default",
                      "should not be used by user");
}

uqRngBoostClass::uqRngBoostClass(int seed, int worldRank)
  :
  uqRngBaseClass(seed,worldRank)
  //m_rng(seed)
{
  resetSeed(seed);
//TODO Find a suitable test for here; Kemelli todo
//   UQ_FATAL_TEST_MACRO(true,
//                       m_worldRank,
//                       "uqRngBoosClass::constructor()",
//                       "Kemelli todo: boost rng");
}

uqRngBoostClass::~uqRngBoostClass()
{
  //this function does nothing
}

  
void
uqRngBoostClass::resetSeed(int newSeed)
{
  m_rng.seed(newSeed);

  return;
}

double
uqRngBoostClass::uniformSample() const
{
  static boost::uniform_01<boost::mt19937> zeroone(m_rng);
  return zeroone();
}

double
uqRngBoostClass::gaussianSample(double stdDev) const
{
  double mean = 0.; //it will be added conveniently later
  static boost::uniform_01<boost::mt19937> zeroone(m_rng);
  boost::math::normal_distribution<double>  gaussian_dist(mean, stdDev);
  return quantile(gaussian_dist, zeroone());  
}

double
uqRngBoostClass::betaSample(double alpha, double beta) const
{
  static boost::uniform_01<boost::mt19937> zeroone(m_rng); 
  boost::math::beta_distribution<double> beta_dist(alpha, beta); 
  return quantile(beta_dist, zeroone());
}

double
uqRngBoostClass::gammaSample(double a, double b) const
{
  static boost::uniform_01<boost::mt19937> zeroone(m_rng);
  boost::math::gamma_distribution<double>  gamma_dist(a,b);
  return quantile(gamma_dist, zeroone());
}

