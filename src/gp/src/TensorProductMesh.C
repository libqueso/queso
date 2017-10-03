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
#include <queso/GPMSAOptions.h>

#include <algorithm> // std::sort
#include <cmath> // std::exp

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
   unsigned int mesh_number,
   std::vector<typename SharedPtr<V>::Type> bases) const
{
  bases.clear(); // We're filling this, not appending to it

  double kernel_spacing[max_coordinates];
  kernel_spacing[0] = opt.m_gaussianDiscrepancyDistanceX[mesh_number];
  kernel_spacing[1] = opt.m_gaussianDiscrepancyDistanceY[mesh_number];
  kernel_spacing[2] = opt.m_gaussianDiscrepancyDistanceZ[mesh_number];
  kernel_spacing[3] = opt.m_gaussianDiscrepancyDistanceT[mesh_number];

  double kernel_periodic[max_coordinates];
  kernel_periodic[0] = opt.m_gaussianDiscrepancyPeriodicX[mesh_number];
  kernel_periodic[1] = opt.m_gaussianDiscrepancyPeriodicY[mesh_number];
  kernel_periodic[2] = opt.m_gaussianDiscrepancyPeriodicZ[mesh_number];
  kernel_periodic[3] = opt.m_gaussianDiscrepancyPeriodicT[mesh_number];

  // "This could get very complicated - perhaps we do not want to
  // attack this issue now"
  for (unsigned int i=0; i != 4; ++i)
    if (kernel_periodic[i])
      queso_not_implemented();

  std::vector<double> peaks_1D[max_coordinates];

  const double support_threshold = opt.m_gaussianDiscrepancySupportThreshold;

  unsigned int num_coordinates_used = 0;
  unsigned int num_grid_points = 1;
  unsigned int num_discrepancy_bases = 1;

  // Figure out where peaks of Gaussian bases should be located, in
  // each dimension.
  for (unsigned int dim = 0; dim != max_coordinates; ++dim)
    {
      if (_coordinate_vals[dim].empty())
        continue;

      // We'll want to make this optional eventually
      const double stddev = kernel_spacing[dim];

      num_coordinates_used++;
      num_grid_points *= _coordinate_vals[dim].size();

      const double central_peak =
        (_coordinate_vals[dim].front() + _coordinate_vals[dim].back()) / 2;

      double peak = central_peak;

      do
        {
          peaks_1D[dim].push_back(peak);

          peak -= kernel_spacing[dim];
          const double boundary_distance =
            peak - _coordinate_vals[dim].front();
          if (boundary_distance < 0)
            {
              const double boundary_value = std::exp(-(boundary_distance * boundary_distance)/2/stddev/stddev);

              if (boundary_value < support_threshold)
                break;
            }
        } while (true);

      peak = central_peak;
      do
        {
          peak += kernel_spacing[dim];
          const double boundary_distance =
            _coordinate_vals[dim].back() - peak;
          if (boundary_distance < 0)
            {
              const double boundary_value = std::exp(-(boundary_distance * boundary_distance)/2/stddev/stddev);

              if (boundary_value < support_threshold)
                break;
            }
          peaks_1D[dim].push_back(peak);
        } while (true);

      num_discrepancy_bases *= peaks_1D[dim].size();
      std::sort(peaks_1D[dim].begin(), peaks_1D[dim].end());
    }

  // Start constructing bases
  double gaussian_peak_index[max_coordinates] = {}; // 0-initialize

  for (unsigned int basis = 0; basis != num_discrepancy_bases; ++basis)
    {
      double mesh_point_index[max_coordinates] = {}; // 0-initialize

      for (unsigned int ptindex = 0; ptindex != num_grid_points; ++ptindex)
        {
          double gaussian_value = 1;
          std::size_t solution_index = this->_first_solution_index;
          std::size_t stride = 1;

          unsigned int nonemptypredim = 0;
          for (unsigned int predim = 0; predim != num_coordinates_used; ++predim)
            {
              unsigned int raw_dim = _order[nonemptypredim];
              while (_coordinate_vals[raw_dim].empty())
                {
                  ++nonemptypredim;
                  raw_dim = _order[nonemptypredim];
                }
              const double stddev = kernel_spacing[raw_dim];
              const double peak = gaussian_peak_index[raw_dim];
              const double point = mesh_point_index[raw_dim];
              const double distance = peak - point;
              gaussian_value *= std::exp(-(distance * distance)/2/stddev/stddev);

              solution_index += stride * mesh_point_index[raw_dim];
              stride *= _coordinate_vals[raw_dim].size();
            }

          nonemptypredim = 0;
          for (unsigned int predim = 0; predim != num_coordinates_used; ++predim)
            {
              unsigned int raw_dim = _order[nonemptypredim];
              while (_coordinate_vals[raw_dim].empty())
                {
                  ++nonemptypredim;
                  raw_dim = _order[nonemptypredim];
                }
              if (++mesh_point_index[raw_dim] == _coordinate_vals[raw_dim].size())
                mesh_point_index[raw_dim] = 0;
              else
                break;
            }
        }

      unsigned int nonemptypredim = 0;
      for (unsigned int predim = 0; predim != num_coordinates_used; ++predim)
        {
          unsigned int raw_dim = _order[nonemptypredim];
          while (_coordinate_vals[raw_dim].empty())
            {
              ++nonemptypredim;
              raw_dim = _order[nonemptypredim];
            }
          if (++gaussian_peak_index[raw_dim] == peaks_1D[raw_dim].size())
            gaussian_peak_index[raw_dim] = 0;
          else
            break;
        }
    }

  // Finally we need to normalize the bases we've just generated.
  // This is based on:
  // D_ij := d_j(tau_i,phi_i) for tau_i, phi_i from simulation data
  // Dnorm := D / max(max(sqrt(D*D')))
  double maxnorm = 0;
  const unsigned int end_solution_index =
    this->_first_solution_index + num_grid_points;
  for (unsigned int i=this->_first_solution_index;
       i != end_solution_index; ++i)
    {
      double normval = 0;
      for (unsigned int b=0; b != bases.size(); ++b)
        {
          const V & basis = *bases[b];
          normval += basis[i]*basis[i];
        }
      maxnorm = std::max(maxnorm, normval);
    }

  const double invmaxnorm = 1./maxnorm;

  for (unsigned int b=0; b != bases.size(); ++b)
    {
      V & basis = *bases[b];
      basis *= invmaxnorm;
    }
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

      if (ub == _coordinate_vals[dim].begin())
        {
          queso_assert_equal_to(outputPoint.val(dim), *ub);
          ++ub;
        }

      indices[dim] = ub - _coordinate_vals[dim].begin() - 1;

      std::vector<double>::const_iterator lb = ub;
      lb--;
      if (ub == _coordinate_vals[dim].end())
        {
          // Hopefully we're just at the endpoint, not
          // extrapolating past it...
          queso_assert_less
            ((outputPoint.val(dim) - *lb)/
             (*lb - *_coordinate_vals[dim].begin()),
             1e-10);
          lb_fraction[dim] = 0;
        }
      else
        {
          lb_fraction[dim] =
            (outputPoint.val(dim) - *lb)/(*ub - *lb);
        }

      num_coordinates_used++;
    }

  // Interpolate from 2 points in 1D, 4 in 2D, 8 in 3D...
  const unsigned int num_points = 1 << num_coordinates_used;

  double interpolated_val = 0;
  for (unsigned int p = 0; p != num_points; ++p)
    {
      std::size_t solution_index = this->_first_solution_index;
      std::size_t stride = 1;
      double coefficient = 1;

      unsigned int nonemptypredim = 0;
      for (unsigned int predim = 0; predim != num_coordinates_used; ++predim)
        {
          unsigned int raw_dim = _order[nonemptypredim];
          while (_coordinate_vals[raw_dim].empty())
            {
              ++nonemptypredim;
              raw_dim = _order[nonemptypredim];
            }

          // Iterate in binary order over neighboring points, e.g.:
          // gridpoint_000(i) <= output_point(i) in all 3 indices i,
          // gridpoint_001(0) >= output_point(0),
          // ...
          // gridpoint_101(i) >= output_point(i) for i in 0,2
          const bool on_larger_point = (1 << predim) & p;

          // Get the index of the point we're looking at in the
          // dimension we're looking at
          const unsigned int point_index = indices[raw_dim] + on_larger_point;

          // Skip ahead in the solution: at the end of this for loop
          // solution_index will correspond to the index for
          // neighboring point p
          solution_index += stride * point_index;
          stride *= _coordinate_vals[raw_dim].size();

          // Evaluate the (bi,tri,)linear "shape function" value at p
          // by starting at 1 but multiplying by the 1-D linear shape
          // function value for each dimension.  If we're looking at
          // the larger point then the shape function goes from 0 on
          // the "left" to 1 on the "right"; otherwise the opposite.
          coefficient *= on_larger_point ?
                         lb_fraction[raw_dim] :
                         (1.0 - lb_fraction[raw_dim]);
        }

      // solution index and coefficient are finally correct here, so
      // we can add the basis function times coefficient for point p.
      interpolated_val += solutionVector[solution_index] * coefficient;

      // When *this* loop is done we'll have
      // solution(outputPoint) = sum_p coefficient_p * shapefunction_p
    }

  return interpolated_val;
}



}  // End namespace QUESO

template class QUESO::TensorProductMesh<QUESO::GslVector>;
