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

template <class V = GslVector>
class TensorProductMesh : public SimulationOutputMesh<V>
{
public:
  TensorProductMesh();

  void set_x_coordinates(std::vector<double> _coord_vals);
  void set_y_coordinates(std::vector<double> _coord_vals);
  void set_z_coordinates(std::vector<double> _coord_vals);
  void set_t_coordinates(std::vector<double> _coord_vals);

  virtual ~TensorProductMesh();

  virtual double interpolateOutput(const V & solutionVector,
                                   const SimulationOutputPoint & outputPoint) const;

private:
  void set_coordinates(unsigned int i, std::vector<double> & _coord_vals);

  std::vector<double> _coordinate_vals[4];
};


}  // End namespace QUESO

#endif // UQ_TENSOR_PRODUCT_MESH_H
