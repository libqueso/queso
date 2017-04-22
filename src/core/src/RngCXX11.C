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

#include <queso/RngCXX11.h>

#ifdef QUESO_HAVE_CXX11

namespace QUESO {

RngCXX11::RngCXX11(int seed, int worldRank)
  :
    RngBase(seed,worldRank),
    m_rng((unsigned int)m_seed)  // m_seed set in the base class.  It is always
                                 // positive
{
}

RngCXX11::~RngCXX11()
{
}

void
RngCXX11::resetSeed(int newSeed)
{
  RngBase::resetSeed(newSeed);
  m_rng.seed(m_seed);  // m_seed was set by the previous line
}

double
RngCXX11::uniformSample() const
{
  std::uniform_real_distribution<double> d(0.0, 1.0);
  return d(m_rng);
}

double
RngCXX11::gaussianSample(double stdDev) const
{
  std::normal_distribution<double> d(0.0, stdDev);
  return d(m_rng);
}

double
RngCXX11::betaSample(double alpha, double beta) const
{
  std::gamma_distribution<double> d1(alpha, 1.0);
  std::gamma_distribution<double> d2(beta, 1.0);

  double x = d1(m_rng);  // x ~ \Gamma(alpha, 1)
  double y = d2(m_rng);  // y ~ \Gamma(beta, 1)

  // x / (x + y) ~ Beta(alpha, beta)
  // See https://en.wikipedia.org/wiki/Beta_distribution#Generating_beta-distributed_random_variates
  return x / (x + y);
}

double
RngCXX11::gammaSample(double a, double b) const
{
  std::gamma_distribution<double> d(a, b);
  return d(m_rng);
}

}  // End namespace QUESO

#endif  // QUESO_HAVE_CXX11
