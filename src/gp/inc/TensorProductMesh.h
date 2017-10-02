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

#ifndef UQ_TENSOR_PRODUCT_MESH_H
#define UQ_TENSOR_PRODUCT_MESH_H

#include <queso/SimulationOutputMesh.h>

namespace QUESO {

/*!
 * \file TensorProductMesh.h
 * \brief This class implements the SimulationOutputMesh interface for
 * meshes defined as tensor products of 1-D meshes
 *
 * \class TensorProductMesh
 * \brief This class implements the SimulationOutputMesh interface for
 * meshes defined as tensor products of 1-D meshes
 */

template <class V = GslVector>
class TensorProductMesh : public SimulationOutputMesh<V>
{
public:
  // Construct an empty mesh
  TensorProductMesh();

  static const unsigned int max_coordinates = 4;

  // Set coordinates for 1-D meshes defining this product mesh.
  void set_x_coordinates(std::vector<double> _coord_vals);
  void set_y_coordinates(std::vector<double> _coord_vals);
  void set_z_coordinates(std::vector<double> _coord_vals);
  void set_t_coordinates(std::vector<double> _coord_vals);

  /**
   * Set the order in which indices are expected to be traversed in
   * the corresponding solution data.  The default is x,y,z,t.  This
   * corresponds to an "order" of "0,1,2,3", and it means that the
   * solution first iterates over all values in x, then increments in
   * y before repeating, then so on for all values in y, then
   * increments in z before repeating, then so on for all values in z
   * before incrementing t.
   *
   * To choose a different ordering, supply a different permutation.
   * For instance, a user whose data increments in t first, then x,
   * could supply order = "3,0,1,2".
   */
  void set_coordinate_order(std::vector<unsigned int> order);

  virtual ~TensorProductMesh();

  //! Return the product of the number of points in each dimension
  virtual std::size_t n_outputs() const;

  virtual double interpolateOutput(const V & solutionVector,
                                   const SimulationOutputPoint & outputPoint) const;

  virtual void generateDiscrepancyBases(const GPMSAOptions & opt,
                                        unsigned int mesh_number,
                                        std::vector<typename SharedPtr<V>::Type> bases) const;

private:
  // Helper function to implement set_x,y,z,t_coordinates
  void set_coordinates(unsigned int i, std::vector<double> & _coord_vals);

  // Coordinates for each of the definitional 1-D meshes
  std::vector<double> _coordinate_vals[max_coordinates];

  // Order permuting index incidence in solution data
  std::vector<unsigned int> _order;
};


}  // End namespace QUESO

#endif // UQ_TENSOR_PRODUCT_MESH_H
