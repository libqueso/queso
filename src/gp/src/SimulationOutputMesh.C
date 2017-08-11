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

#include <queso/SimulationOutputMesh.h>
#include <queso/GslVector.h>

namespace QUESO {

template <class V>
SimulationOutputMesh<V>::SimulationOutputMesh
  (std::size_t first_solution_index) :
  _first_solution_index(first_solution_index)
{
}



template <class V>
SimulationOutputMesh<V>::~SimulationOutputMesh()
{
}



template <class V>
void
SimulationOutputMesh<V>::interpolateOutputs
  (const V & solutionVector,
   const std::vector<SimulationOutputPoint> & outputPoints,
   V & interpolatedValues)
{
  for (unsigned int i=0; i != solutionVector.sizeLocal(); ++i)
    interpolatedValues[i] =
      this->interpolateOutput(solutionVector, outputPoints[i]);
}



}  // End namespace QUESO

template class QUESO::SimulationOutputMesh<QUESO::GslVector>;
