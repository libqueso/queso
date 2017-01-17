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

#include <queso/CovCond.h>

namespace QUESO {

template <class V, class M>
void CovCond(double condNumber, const V & direction, M & covMatrix,
    M & precMatrix)
{
  //std::cout << "Entering CovCond()"
  //          << std::endl;

  V v1(direction);
  //std::cout << "In CovCond(), v1 contents are:"
  //          << std::endl
  //          << v1
  //          << std::endl;

  V v2(direction,condNumber,1.0); // MATLAB linspace
  v2.cwInvert();
  v2.sort();
  //std::cout << "In CovCond(), v2 contents are:"
  //          << std::endl
  //          << v2
  //          << std::endl;

  double v1Norm2 = v1.norm2();
  if (v1[0] >=0) v1[0] += v1Norm2;
  else           v1[0] -= v1Norm2;
  double v1Norm2Sq = v1.norm2Sq();

  M Z(direction,1.0);
  Z -= (2./v1Norm2Sq) * matrixProduct(v1,v1);
  //std::cout << "In CovCond(), Z contents are:"
  //          << std::endl
  //          << Z
  //          << std::endl;

  M Zt(Z.transpose());
  covMatrix  = Z * leftDiagScaling(v2,   Zt);
  precMatrix = Z * leftDiagScaling(1./v2,Zt);

  //std::cout << "Leaving CovCond()"
  //          << std::endl;
}

}  // End namespace QUESO
