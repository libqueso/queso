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

#ifndef UQ_LINEAR_LAGRANGE_INTERPOLATION_SURROGATE_H
#define UQ_LINEAR_LAGRANGE_INTERPOLATION_SURROGATE_H

// QUESO
#include <queso/InterpolationSurrogateBase.h>

// C++
#include <vector>
#include <cmath>

namespace QUESO
{
  class GslVector;
  class GslMatrix;

  //! Linear Lagrange interpolation surrogate
  /*! Aribritary dimension linear Lagrange interpolant. Uses a structured
      grid to interpolate values given by data passed to the constructor. */
  template<class V = GslVector, class M = GslMatrix>
  class LinearLagrangeInterpolationSurrogate : public InterpolationSurrogateBase<V,M>
  {
  public:

    //! Constructor
    /*! The data object should be already fully populated when constructing
        this object. The user can directly set values or use the
        InterpolationSurrogateBuilder object to populate the data. */
    LinearLagrangeInterpolationSurrogate(const InterpolationSurrogateData<V,M>& data);

    virtual ~LinearLagrangeInterpolationSurrogate(){};

    //! Evaluates value of the interpolant for the given domainVector
    virtual double evaluate(const V & domainVector) const;

    //! The number of coeffs for interpolating
    unsigned int n_coeffs() const
    { return std::pow( 2, this->m_data.dim() ); };

  protected:

    //! Helper function to get lower bound indices for each dimension
    /*! domainVector contains the current value at which want want to interpolate.
        We use that to figure out the bounding values in each dimension. Then
        we populate indices with the lower index. */
    void compute_interval_indices(const V & domainVector,
                                  std::vector<unsigned int>& indices) const;

    //! Helper function to populate bounding values for the intervals in each dimension
    /*! By convention, it is assumed that the indices vector contains the lower
        bound index. See compute_interval_indices. */
    void compute_interval_values( const std::vector<unsigned int>& indices,
                                  std::vector<double>& x_min,
                                  std::vector<double>& x_max,
                                  std::vector<double>& values ) const;

    //! Evaluate multidimensional linear Lagrange interpolant
    /*! This is just a tensor product of one-dimensional interpolants */
    double eval_interpolant( const std::vector<double>& x_min,
                             const std::vector<double>& x_max,
                             const std::vector<double>& values,
                             const V & domainVector ) const;

    //! Convert local indices to single "global" index
    /*! This is for handling the local ordering of an arbitrary dimensional
        linear Lagrange elements. The function will look like
        \f$ \sum_{ijkl...} f_{ijkl...} N_{ijkl...} \f$ so we want a
        convenient way to map \f$ ijkl... \f$ to a single integer. Then, we
        can order all our data this way and easily index into it. We can
        do this since each of \f$ ijkl... \f$ is either 0 or 1. */
    unsigned int coordsToSingle( const std::vector<unsigned int>& indices ) const;

    //! Inverse of map computed in coordsToSingle
    void singleToCoords( unsigned int global, std::vector<unsigned int>& indices ) const;

    //! Compute multidimensional Lagrange polynomial as tensor product of 1D polynomials
    /*! \f$ N_{ijkl...} (\xi, \eta, \zeta,...) = N_i(\eta) N_j(\eta) N_k(\zeta) ... \f$
      The indices vector contains the i,j,k,... values */
    double tensor_product_lagrange( const std::vector<double>& x_min,
                                    const std::vector<double>& x_max,
                                    const std::vector<unsigned int>& indices,
                                    const V & domainVector ) const;

    //! Evaluate a 1-D Lagrange polynomial at point x
    /*! index tells us which of the two Lagrange functions to use.
      index must be 0 or 1. */
    double lagrange_poly( double x0, double x1, double x, unsigned int index ) const;

  private:

    LinearLagrangeInterpolationSurrogate();

  };

} // end namespace QUESO

#endif // UQ_LINEAR_LAGRANGE_INTERPOLATION_SURROGATE_H
