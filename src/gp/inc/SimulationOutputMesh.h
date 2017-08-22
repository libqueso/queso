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

#ifndef UQ_SIMULATION_OUTPUT_MESH_H
#define UQ_SIMULATION_OUTPUT_MESH_H

#include <vector>

#include <queso/SimulationOutputPoint.h>

namespace QUESO {

// Forward declarations
class GslVector;

template <class V = GslVector>
class SimulationOutputMesh
{
public:
  SimulationOutputMesh(std::size_t first_solution_index = 0);

  virtual ~SimulationOutputMesh();

  std::size_t first_solution_index() { return _first_solution_index; }
  virtual std::size_t n_outputs() = 0;

  virtual double interpolateOutput(const V & solutionVector,
                                   const SimulationOutputPoint & outputPoint) const = 0;

  virtual void interpolateOutputs(const V & solutionVector,
                                  const std::vector<SimulationOutputPoint> & outputPoints,
                                  V & interpolatedValues);

protected:

  std::size_t _first_solution_index;
};


}  // End namespace QUESO

#endif // UQ_SIMULATION_OUTPUT_MESH_H
