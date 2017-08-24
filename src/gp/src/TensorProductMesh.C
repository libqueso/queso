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
#include <queso/asserts.h>

namespace QUESO {

template <class V>
TensorProductMesh<V>::TensorProductMesh()
{
  _order.resize(max_coordinates);
  for (unsigned int i=0; i != max_coordinates; ++i)
    _order[i] = i;
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

  // Check for user error
  // We won't check for empty coordinate_vals, because the user might
  // be disabling a previously-enabled coordinate for some reason.
  for (unsigned int v=1; v < _coordinate_vals[i].size(); ++v)
    queso_require(_coordinate_vals[i][v-1] < _coordinate_vals[i][v]);
}



template <class V>
void
TensorProductMesh<V>::set_coordinate_order
  (std::vector<unsigned int> order)
{
  // We passed that argument by value, so it should be safe to steal.
  _order.swap(order);

#ifndef NDEBUG
  for (std::size_t i=0; i != _order.size(); ++i)
    {
      queso_assert_less(_order[i], max_coordinates);
      for (std::size_t j=i+1; j != _order.size(); ++j)
        queso_assert_not_equal_to(_order[i], _order[j]);
    }
#endif
}



template <class V>
std::size_t
TensorProductMesh<V>::n_outputs() const
{
  std::size_t n_outputs = 1;
  for (unsigned int dim = 0; dim != max_coordinates; ++dim)
    if (!_coordinate_vals[dim].empty())
      n_outputs *= _coordinate_vals[dim].size();

  return n_outputs;
}


template <class V>
void
TensorProductMesh<V>::generateDiscrepancyBases
  (const GPMSAOptions & opt,
   std::vector<typename SharedPtr<V>::Type> bases) const
{
  queso_not_implemented(); // FIXME
}


template <class V>
double
TensorProductMesh<V>::interpolateOutput
  (const V & solutionVector,
   const SimulationOutputPoint & outputPoint) const
{
  std::vector<std::size_t> indices(max_coordinates, 0);
  std::vector<double> lb_fraction(max_coordinates, 0);

  unsigned int num_coordinates_used = 0;

  // Find where we are in the space of each coordinate
  for (unsigned int dim = 0; dim != max_coordinates; ++dim)
    {
      if (_coordinate_vals[dim].empty())
        continue;

      std::vector<double>::const_iterator ub =
        std::lower_bound(_coordinate_vals[dim].begin(),
                         _coordinate_vals[dim].end(), outputPoint.val(dim));
      queso_assert(ub != _coordinate_vals[dim].end());
      queso_assert(ub != _coordinate_vals[dim].begin());
      indices[dim] = ub - _coordinate_vals[dim].begin() - 1;

      std::vector<double>::const_iterator lb = ub;
      lb--;
      if (ub == _coordinate_vals[dim].end())
        {
          // Hopefully we're just at the endpoint, not
          // extrapolating past it
          lb_fraction[dim] = 0;
        }
      else
        {
          lb_fraction[dim] =
            (outputPoint.val(dim) - *lb)/(*ub - *lb);
        }

      num_coordinates_used++;
    }

  const unsigned int num_points = 1 << num_coordinates_used;

  std::vector<unsigned int> point_indices(num_coordinates_used, 0);

  double interpolated_val = 0;
  for (unsigned int p = 0; p != num_points; ++p)
    {
      std::size_t solution_index = this->_first_solution_index;
      std::size_t stride = 1;
      double coefficient = 1;
      for (unsigned int predim = 0; predim != num_coordinates_used; ++predim)
        {
          const unsigned int raw_dim = _order[predim];
          const bool on_larger_point = (1 << predim) & p;
          point_indices[predim] = indices[raw_dim] + on_larger_point;
          solution_index += stride * point_indices[predim];
          stride *= _coordinate_vals[raw_dim].size();
          coefficient *= on_larger_point ?
                         (1.0 - lb_fraction[raw_dim]) :
                         lb_fraction[raw_dim];
        }

      interpolated_val += solutionVector[solution_index] * coefficient;
    }

  return interpolated_val;
}



}  // End namespace QUESO

template class QUESO::TensorProductMesh<QUESO::GslVector>;
