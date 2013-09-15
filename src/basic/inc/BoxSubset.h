//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
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
// 
// $Id$
//
//--------------------------------------------------------------------------

#ifndef UQ_BOX_SUBSET_H
#define UQ_BOX_SUBSET_H

#include <queso/VectorSpace.h>

namespace QUESO {

/*!
 * \class BoxSubset
 * \brief Class representing a subset of a vector space shaped like a hypercube
 *
 * This class is determined by a collection of upper and lower limits of the hypercube.
 * (line segment in \f$ R \f$, rectangle in \f$ R^2 \f$, and so on). */

template<class V, class M>
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

// Default, shaped constructor ----------------------
template<class V, class M>
BoxSubset<V,M>::BoxSubset(
  const char*                    prefix,
  const VectorSpace<V,M>& vectorSpace,
  const V&                       minValues,
  const V&                       maxValues)
  :
  VectorSubset<V,M>(prefix,vectorSpace,0.),
  m_minValues(minValues),
  m_maxValues(maxValues)
{
  UQ_FATAL_TEST_MACRO(minValues.sizeLocal() != maxValues.sizeLocal(),
                      m_env.worldRank(),
                      "BoxSubset<V,M>::BoxSubset()",
                      "vectors 'minValues' and 'maxValues' should have the same size");
  UQ_FATAL_TEST_MACRO(minValues.sizeLocal() != vectorSpace.dimLocal(),
                      m_env.worldRank(),
                      "BoxSubset<V,M>::BoxSubset()",
                      "sizes of vectors 'minValues' and 'maxValues' should be equal to dimension of the vector space");
  for (unsigned int i = 0; i < m_vectorSpace->dimLocal(); ++i) {
    UQ_FATAL_TEST_MACRO(minValues[i] > maxValues[i],
                        m_env.worldRank(),
                        "BoxSubset<V,M>::BoxSubset()",
                        "it should happen minValue <= maxValue for all dimensions");
  }

  m_volume = 1.;
  for (unsigned int i = 0; i < m_vectorSpace->dimLocal(); ++i) {
    m_volume *= (m_maxValues[i] - m_minValues[i]);
  }
}
// Destructor ---------------------------------------
template<class V, class M>
BoxSubset<V,M>::~BoxSubset()
{
}
// Math methods -------------------------------------
template<class V, class M>
bool
BoxSubset<V,M>::contains(const V& vec) const
{
  // prudenci, 2012-09-26: allow boundary values because of 'beta' realizer, which can generate a sample with boundary value '1'
  //return (!vec.atLeastOneComponentSmallerOrEqualThan(m_minValues) &&
  //        !vec.atLeastOneComponentBiggerOrEqualThan (m_maxValues));
  return (!vec.atLeastOneComponentSmallerThan(m_minValues) &&
          !vec.atLeastOneComponentBiggerThan (m_maxValues));
}
// --------------------------------------------------
template<class V, class M>
const V&
BoxSubset<V,M>::minValues() const
{
  return m_minValues;
}
// --------------------------------------------------
template<class V, class M>
const V&
BoxSubset<V,M>::maxValues() const
{
  return m_maxValues;
}
// I/O method ---------------------------------------
template <class V, class M>
void
BoxSubset<V,M>::print(std::ostream& os) const
{
  os << "In BoxSubset<V,M>::print()"
     << ": m_minValues = " << m_minValues
     << ", m_maxValues = " << m_maxValues
     << ", m_volume = "    << m_volume
     << std::endl;

  return;
}

}  // End namespace QUESO

#endif // UQ_BOX_SUBSET_H
