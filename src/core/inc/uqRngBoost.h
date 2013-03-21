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

#ifndef __UQ_RNG_BOOST_H__
#define __UQ_RNG_BOOST_H__

#include <uqRngBase.h>
#include <boost/random.hpp>
#include <boost/math/distributions.hpp> 

class uqRngBoostClass : public uqRngBaseClass
{
public:
  
  uqRngBoostClass();
  uqRngBoostClass(int seed, int worldRank);
 ~uqRngBoostClass();

  void     resetSeed     (int newSeed);
  double   uniformSample ()                          const;
  double   gaussianSample(double stdDev)             const;
  double   betaSample    (double alpha, double beta) const;
  double   gammaSample   (double a, double b)        const;
	
private:  
  boost::mt19937 m_rng; // it cannot be static, as it is not initialized yet
};

#endif // __UQ_RNG_BOOST_H__

