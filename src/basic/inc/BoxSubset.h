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

#ifndef UQ_BOX_SUBSET_H
#define UQ_BOX_SUBSET_H

#include <queso/VectorSpace.h>
#include <queso/VectorSubset.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*!
 * \class BoxSubset
 * \brief Class representing a subset of a vector space shaped like a hypercube
 *
 * This class is determined by a collection of upper and lower limits of the hypercube.
 * (line segment in \f$ R \f$, rectangle in \f$ R^2 \f$, and so on). */

template <class V = GslVector, class M = GslMatrix>
class BoxSubset : public VectorSubset<V,M> {
public:
  //! @name Constructor/Destructor methods
  //@{

  //! Shaped, default constructor.
  /*! Construct a subspace of \c vectorSpace, with min and max values given by the vectors \c minValues
   * and \c maxValues, respectively. It checks for possible inconsistencies between the values stored in
   * \c minValues and \c maxValues, and calculates the volume of the box subset, assigning it to m_volume. */
  BoxSubset(const char*                    prefix,
                   const VectorSpace<V,M>& vectorSpace,
                   const V&                       minValues,
                   const V&                       maxValues);

  //! Destructor
  ~BoxSubset();
  //@}

  //! @name Mathematical methods.
  //@{
  //! Checks whether this box subset contains vector \c vec.
  /*! It checks if both statements are true: 1) all components in \c vec are larger than
   * m_minValues, and 2) all all components in \c vec are smaller than m_maxValues. */
  bool contains (const V& vec)     const;

  //! Returns the centroid of this box subset in the vector \c vec.
  void centroid (V& vec)     const;

  //! Returns the moments of inertia of this box subset in the matrix \c mat.
  void moments (M & mat)     const;

  //! Vector of the minimum values of the box subset.
  const V&   minValues()                 const;

  //! Vector of the maximum values of the box subset.
  const V&   maxValues()                 const;
  //@}

  //! Prints the volume, the minimum and the maximum values of \c this.
  void print    (std::ostream& os) const;

protected:
  using VectorSet   <V,M>::m_env;
  using VectorSet   <V,M>::m_prefix;
  using VectorSet   <V,M>::m_volume;
  using VectorSubset<V,M>::m_vectorSpace;

  //! Vector of templated type \c V to store the minimum values of the box subset class.
  V m_minValues;

  //! Vector of templated type \c V to store the maximum values of the box subset class.
  V m_maxValues;
};

}  // End namespace QUESO

#endif // UQ_BOX_SUBSET_H
