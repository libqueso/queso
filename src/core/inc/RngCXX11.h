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

#ifndef QUESO_RNGCXX11_H
#define QUESO_RNGCXX11_H

#include <queso/Defines.h>

#ifdef QUESO_HAVE_CXX11

#include <queso/RngBase.h>
#include <random>

/*!
 * \file RngCXX11.h
 * \brief C++11 Random Number Generation class.
 */

/*!
 * \class RngCXX11
 * \brief Class for random number generation using std::random from C++11.
 *
 * This class is  class of random number generator, using std::random.
 */
namespace QUESO {

class RngCXX11 : public RngBase
{
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Constructor with seed.
  /*
   * Uses std::random with a generator of type: mt19937.  The widely used
   * Mersenne Twister random number generator.
   */
  RngCXX11(int seed, int worldRank);

  //! Destructor
  ~RngCXX11();
  //@}

  //! @name Sampling methods
  //@{
  //! Resets the seed with value \c newSeed.
  void resetSeed(int newSeed);

  //! Samples a value from a uniform distribution. Support: [0,1) or [a,b).
  /*!
   * This function samples from continuous uniform distribution on the range
   * [0,1). It is possible to scale this distribution so the support is defined
   * by the two parameters, a and b, which are its minimum and maximum values.
   * Support: -infinity < a < x< b< infinity.
   */
  double uniformSample() const;

  //! Samples a value from a Gaussian distribution with standard deviation
  //! given by \c stdDev.  Support:  (-infinity, infinity).
  /*!
   * The parameter mu (mean or expectation of the distribution) in this
   * Gaussian sample is set to zero, and thus, needs to be provided in an
   * alternative way (e.g., in the form of a sum. The parameter stdDev is its
   * standard deviation; its variance is therefore stdDev^2.  A random variable
   * with a Gaussian distribution is said to be normally distributed and is
   * called a normal deviate. Support: (-infinity, infinity).
   */
  double gaussianSample(double stdDev) const;

  //! Samples a value from a Beta distribution. Support: [0,1]
  /*!
   * The Beta Distribution is a continuous probability distribution; it has two
   * free parameters, which are labeled alpha and beta. The beta distribution
   * is used as a prior distribution for binomial proportions in Bayesian
   * analysis  (Evans et al.  2000, p. 34). Support (domain): [0,1].
   *
   * Since std::random doesn't have a beta distribution implementation, we
   * leverage a novel transformation involving two Gamma distributed variates
   * instead.  For details, see https://en.wikipedia.org/wiki/Beta_distribution#Generating_beta-distributed_random_variates
   */
  double betaSample(double alpha, double beta) const;

  //! Samples a value from a Gamma distribution. Support: [0,infinity).
  /*!
   * The Gamma Distribution is a continuous probability distribution; it has
   * two free parameters, which may be labeled: a shape parameter a and an
   * inverse scale parameter b, called a rate parameter.
   * Support (domain): [0, infinity).
   */
  double gammaSample(double a, double b) const;

private:
  //! Default Constructor: it should not be used.
  RngCXX11();

  //! The internal random number generator.  Mersenne Twister generator.
  mutable std::mt19937 m_rng;
};

}  // End namespace QUESO

#endif  // QUESO_HAVE_CXX11

#endif  // QUESO_RNGCXX11_H
