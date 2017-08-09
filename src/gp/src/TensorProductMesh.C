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

#include <queso/TensorProductMesh.h>
#include <queso/GslVector.h>

namespace QUESO {

template <class V>
TensorProductMesh<V>::TensorProductMesh()
{
}



template <class V>
TensorProductMesh<V>::~TensorProductMesh()
{
}



template <class V>
void
TensorProductMesh<V>::set_x_coordinates
  (std::vector<double> _coord_vals)
{
  this->set_coordinates(0, _coord_vals);
}



template <class V>
void
TensorProductMesh<V>::set_y_coordinates
  (std::vector<double> _coord_vals)
{
  this->set_coordinates(1, _coord_vals);
}



template <class V>
void
TensorProductMesh<V>::set_z_coordinates
  (std::vector<double> _coord_vals)
{
  this->set_coordinates(2, _coord_vals);
}



template <class V>
void
TensorProductMesh<V>::set_t_coordinates
  (std::vector<double> _coord_vals)
{
  this->set_coordinates(3, _coord_vals);
}



template <class V>
void
TensorProductMesh<V>::set_coordinates
  (unsigned int i,
   std::vector<double> & _coord_vals)
{
  // We passed that argument by value, so it should be safe to steal.
  _coordinate_vals[i].swap(_coord_vals);

  // Assert that coordinates are sorted
}


template <class V>
double
TensorProductMesh<V>::interpolateOutput
  (const V & solutionVector,
   const SimulationOutputPoint & outputPoint) const
{
  std::vector<unsigned int> indices(4, 0);

  for (unsigned int dim = 0; dim != 4; ++dim)
    {
      if (!_coordinate_vals[dim].empty())
        {
          std::vector<double>::const_iterator lb =
            std::lower_bound(_coordinate_vals[dim].begin(),
                             _coordinate_vals[dim].end(), outputPoint.val(dim));
          queso_assert_not_equal_to(lb, _coordinate_vals[dim].end());
          indices[dim] = lb - _coordinate_vals[dim].begin();
        }
    }
}



}  // End namespace QUESO

template class QUESO::TensorProductMesh<QUESO::GslVector>;
