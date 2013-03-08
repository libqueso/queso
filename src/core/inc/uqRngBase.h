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

#ifndef __UQ_RNG_BASE_H__
#define __UQ_RNG_BASE_H__

#include <uqDefines.h>
#include <iostream>

class uqRngBaseClass
{
public:
           uqRngBaseClass();
           uqRngBaseClass(int seed, int worldRank);
  virtual ~uqRngBaseClass();

          int    seed          () const;
  virtual void   resetSeed     (int newSeed);
  virtual double uniformSample ()                          const = 0;
  virtual double gaussianSample(double stdDev)             const = 0;
  virtual double betaSample    (double alpha, double beta) const = 0;
  virtual double gammaSample   (double a, double b)        const = 0;

protected:
          int m_seed;
          int m_worldRank;

private:
          void privateResetSeed();
};

#endif // __UQ_RNG_BASE_H__
