/* uq/libs/mcmc/inc/uqCovCond.h
 *
 * Copyright (C) 2008 The PECOS Team, http://www.ices.utexas.edu/centers/pecos
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __UQ_COV_COND_H__
#define __UQ_COV_COND_H__

#include <iostream>

template <class V, class M>
void
uqCovCond(
  double   condNumber,
  const V& direction,
  M&       covMatrix,
  M&       precMatrix)
{
  //std::cout << "Entering uqCovCond()"
  //          << std::endl;

  V v1(direction);
  unsigned int size1 = v1.size();
  //std::cout << "In uqCovCond(), v1 contents are:"
  //          << std::endl
  //          << v1
  //          << std::endl;

  V v2(direction.env(),condNumber,1.0,size1); // MATLAB linspace
  v2.cwInvert();
  v2.sort();
  //std::cout << "In uqCovCond(), v2 contents are:"
  //          << std::endl
  //          << v2
  //          << std::endl;

  double v1Norm2 = v1.norm2();
  if (v1[0] >=0) v1[0] += v1Norm2;
  else           v1[0] -= v1Norm2;
  double v1Norm2Sq = v1.norm2Sq();

  M Z(direction.env(),size1,1.0);
  Z -= (2./v1Norm2Sq) * matrixProduct(v1,v1);
  //std::cout << "In uqCovCond(), Z contents are:"
  //          << std::endl
  //          << Z
  //          << std::endl;

  M Zt(Z.transpose());
  covMatrix  = Z * diagScaling(v2,   Zt);
  precMatrix = Z * diagScaling(1./v2,Zt);

  //std::cout << "Leaving uqCovCond()"
  //          << std::endl;
}
#endif // __UQ_COV_COND_H__
