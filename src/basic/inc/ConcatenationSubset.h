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

#ifndef UQ_CONCATENATION_SUBSET_H
#define UQ_CONCATENATION_SUBSET_H

#include <queso/VectorSpace.h>
#include <queso/VectorSubset.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*! \class ConcatenationSubset
 * \brief A templated class representing the concatenation of two vector subsets.
 *
 * This class is used to represent the concatenation of two subsets. It allows concatenation
 * of both two subsets as well as a collection of subsets into one.
 */
template <class V = GslVector, class M = GslMatrix>
class ConcatenationSubset : public VectorSubset<V,M> {
public:
   //! @name Constructor/Destructor methods.
  //@{
  //! Constructor - two sets only
  /*! It concatenates set1 and set2 into a new set, m_set, of volume given by
   * set1.volume()*set2.volume(). */
  ConcatenationSubset(const char*                    prefix,
                             const VectorSpace<V,M>& vectorSpace,
                             const VectorSet<V,M>&   set1,
                             const VectorSet<V,M>&   set2);

  //! Constructor - collection of sets
  /*! It concatenates a collection of n subsets (using the std::vector to represent such collection)
   * into a new set, m_set, of volume given by sets[0].volume()*sets[1].volume()*...*sets[n-1].volume(). */
  ConcatenationSubset(const char*                                       prefix,
                             const VectorSpace<V,M>&                    vectorSpace,
                             double                                            volume,
                             const std::vector<const VectorSet<V,M>* >& sets);
  //! Destructor
  ~ConcatenationSubset();
  //@}

  //! @name Mathematical methods.
  //@{
  //! Determines whether each one of the subsets m_sets (class' private attributes) contains vector \c vec.
  bool contains(const V& vec)     const;

  //! Returns the set centroid in the vector \c vec.
  virtual       void                     centroid   (V& vec)     const;

  //! Returns the set moments of inertia in the matrix \c mat.
  virtual       void                     moments (M & mat)     const;
  //@}

  //! @name I/O methods.
  //@{
  //! Prints the subsets (via protected attribute m_sets).
  void print   (std::ostream& os) const;
  //@}
protected:
  using VectorSet   <V,M>::m_env;
  using VectorSet   <V,M>::m_prefix;
  using VectorSet   <V,M>::m_volume;
  using VectorSubset<V,M>::m_vectorSpace;

  std::vector<const VectorSet<V,M>* > m_sets;
};

}  // End namespace QUESO

#endif // UQ_CONCATENATION_SUBSET_H
