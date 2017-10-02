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

#include <queso/SimulationOutputPoint.h>

#include <queso/SharedPtr.h>

#include <vector>

namespace QUESO {

// Forward declarations
class GslVector;

class GPMSAOptions;

/*!
 * \file SimulationOutputMesh.h
 * \brief This abstract base class defines the interface which
 * gaussian process emulation requires of a simulation mesh
 *
 * \class SimulationOutputMesh
 * \brief This abstract base class defines the interface which
 * gaussian process emulation requires of a simulation mesh
 */


template <class V = GslVector>
class SimulationOutputMesh
{
public:
  //! Construct, and associate data with solution vector index
  SimulationOutputMesh(std::size_t first_solution_index = 0);

  virtual ~SimulationOutputMesh();

  //! Get the solution vector index at which this mesh's data begins
  std::size_t first_solution_index() const { return _first_solution_index; }

  //! Get the number of indices associated with data on this mesh
  virtual std::size_t n_outputs() const = 0;

  //! Return the value of simulation data at a particular point
  /*!
   * Using coefficients from \c solutionVector for data discretized
   * on this mesh, return the value which the solution would take at
   * the point \c outputPoint (which may not be a mesh node)
   */
  virtual double interpolateOutput(const V & solutionVector,
                                   const SimulationOutputPoint & outputPoint) const = 0;

  //! Return the value of simulation data at particular points
  /*!
   * Using coefficients from \c solutionVector for data discretized
   * on this mesh, return the value which the solution would take at
   * the points \c outputPoints (which may not be mesh nodes)
   */
  virtual void interpolateOutputs(const V & solutionVector,
                                  const std::vector<SimulationOutputPoint> & outputPoints,
                                  V & interpolatedValues);

  //! Fill solution vectors with a discrepancy basis
  /*!
   * Based on gaussian discrepancy basis generation options given by
   * \c opt (for the \c mesh_number corresponding to this mesh), 
   * fill \c bases with the coefficients of the generated discrepancy
   * basis functions. 
   */
  virtual void generateDiscrepancyBases(const GPMSAOptions & opt,
                                        unsigned int mesh_number,
                                        std::vector<typename SharedPtr<V>::Type> bases) const = 0;

protected:

  // The solution vector index at which this mesh's data begins
  std::size_t _first_solution_index;
};


}  // End namespace QUESO

#endif // UQ_SIMULATION_OUTPUT_MESH_H
