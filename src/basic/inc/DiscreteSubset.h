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

#ifndef UQ_DISCRETE_SUBSET_H
#define UQ_DISCRETE_SUBSET_H

#include <queso/VectorSpace.h>
#include <queso/VectorSubset.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*! \class DiscreteSubset
 * \brief A templated class representing the discrete vector subsets.
 *
 * This class is used to represent the a discrete vector subset. Here the
 * notion of volume does not apply, instead there is the total number of
 * elements (vectors) that belongs to the subset.
 */

template <class V = GslVector, class M = GslMatrix>
class DiscreteSubset : public VectorSubset<V,M> {
public:
  //! @name Constructor/Destructor methods.
  //@{
  //! Default Constructor
  /*! It constructs a class object given the prefix, vector space to which it
   * belongs  and its number of elements.*/
  DiscreteSubset(const char*                    prefix,
                        const VectorSpace<V,M>& vectorSpace,
                        const std::vector<V*>&         elements);
  //! Destructor
  ~DiscreteSubset();
  //@}

  //! @name Mathematical methods.
  //@{
  //! Checks whether this discrete subset contains vector \c vec. TODO: incomplete code.
  bool contains (const V& vec)     const;
  //@}

  //! @name I/O methods.
  //@{
  //! Prints nothing.
  void print    (std::ostream& os) const;
  //@}
protected:
  using VectorSet   <V,M>::m_env;
  using VectorSet   <V,M>::m_prefix;
  using VectorSet   <V,M>::m_volume;
  using VectorSubset<V,M>::m_vectorSpace;

  //! Number of elements in the discrete vector subset.
  std::vector<V*> m_elements;
};

}  // End namespace QUESO

#endif // UQ_DISCRETE_SUBSET_H
