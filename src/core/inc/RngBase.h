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

#ifndef UQ_RNG_BASE_H
#define UQ_RNG_BASE_H

#include <queso/Defines.h>
#include <iostream>

namespace QUESO {

/*! \file RngBase.h
    \brief Random Number Generation class.
*/

/*! \class RngBase
    \brief Class for random number generation (base class for either GSL or Boost RNG).

    This class is  a “virtual” class of generic random number generator, in order
    to accommodate either GSL or Boost RNG.
*/


class RngBase
{
public:

  //! @name Constructor/Destructor methods
  //@{
  //! Constructor with seed.
  RngBase(int seed, int worldRank);

  //! Virtual destructor.
  virtual ~RngBase();
  //@}


    //! @name Sampling methods
  //@{
  //! Sets the seed.
          int    seed          () const;

  //! Resets the seed with value \c newSeed.
  virtual void   resetSeed     (int newSeed);

  //! Samples a value from a uniform distribution.
  virtual double uniformSample ()                          const = 0;

  //! Samples a value from a Gaussian distribution with standard deviation given by \c stdDev.
  virtual double gaussianSample(double stdDev)             const = 0;

  //! Samples a value from a Beta distribution.
  virtual double betaSample    (double alpha, double beta) const = 0;

  //! Samples a value from a Gamma distribution.
  virtual double gammaSample   (double a, double b)        const = 0;

  //@}
protected:
  //! Seed.
          int m_seed;

  //! Rank of processor.
          int m_worldRank;

private:
  //! Default Constructor: it should not be used.
  RngBase();

  //! Reset seed.
  void privateResetSeed();
};

}  // End namespace QUESO

#endif // UQ_RNG_BASE_H
