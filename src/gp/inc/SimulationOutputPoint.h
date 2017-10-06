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

#ifndef UQ_SIMULATION_OUTPUT_POINT_H
#define UQ_SIMULATION_OUTPUT_POINT_H

#include <vector>

namespace QUESO {

/*!
 * \file SimulationOutputPoint.h
 * \brief This class represents a single point in space/time.
 *
 * \class SimulationOutputPoint
 * \brief This class represents a single point in space/time.
 */

class SimulationOutputPoint
{
public:
  // Initialize to (x,y,z,t) = 0
  SimulationOutputPoint()
  { _val[0] = _val[1] = _val[2] = _val[3] = 0;}

  // Getters
  double x() const { return _val[0]; }
  double y() const { return _val[1]; }
  double z() const { return _val[2]; }
  double t() const { return _val[3]; }

  // Get by index rather than coordinate name
  double val(unsigned int i) const { return _val[i]; }

  // Setters
  double & x() { return _val[0]; }
  double & y() { return _val[1]; }
  double & z() { return _val[2]; }
  double & t() { return _val[3]; }

private:
  double _val[4];  // hardcoding x,y,z,t
};


}  // End namespace QUESO

#endif // UQ_SIMULATION_OUTPUT_POINT_H
