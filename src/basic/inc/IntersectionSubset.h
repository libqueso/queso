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

#ifndef UQ_INTERSECTION_SUBSET_H
#define UQ_INTERSECTION_SUBSET_H

#include <queso/VectorSpace.h>
#include <queso/VectorSubset.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*! \class IntersectionSubset
 * \brief A templated class representing the intersection of two vector sets.
 *
 * This class is used to determine if a vector  belongs to the intersection of
 * two vector sets. It is useful for handling a posterior PDF, since its domain
 * is the intersection of the domain of the prior PDF with the domain of the
 * likelihood function.*/

template <class V = GslVector, class M = GslMatrix>
class IntersectionSubset : public VectorSubset<V,M> {
public:
  //! @name Constructor/Destructor methods.
  //@{
  //! Default, shaped constructor.
  /*! Creates the class for the intersection of two vector sets, given a vector space, its volume and the
   * sets.*/
  IntersectionSubset(const char*                    prefix,
                            const VectorSpace<V,M>& vectorSpace,
                                  double                   volume,
                            const VectorSet<V,M>&   set1,
                            const VectorSet<V,M>&   set2);

  //! Destructor.
 ~IntersectionSubset();
  //@}

   //! @name Mathematical methods.
  //@{
 //! Determines whether both sets m_set1 and m_set2 (class' private attributes) contain vector \c vec.
  bool contains(const V& vec)     const;

  //! Returns the set centroid in the vector \c vec.  Not implemented.
  void centroid (V& vec)     const;
  //@}

  //! Returns the set moments of inertia in the matrix \c mat.  Not implemented.
  void moments (M& mat)     const;
  //@}

  //! @name I/O methods.
  //@{
  //! Prints both subsets (via protected attributes m_set1 and m_set2).
  void print   (std::ostream& os) const;
  //@}

protected:
  using VectorSet   <V,M>::m_env;
  using VectorSet   <V,M>::m_prefix;
  using VectorSet   <V,M>::m_volume;
  using VectorSubset<V,M>::m_vectorSpace;

  //! Vector set: m_set1.
  /*! We seek the intersection of vectors set m_set1 and m_set2.*/
  const VectorSet<V,M>& m_set1;

  //! Vector set: m_set2.
  const VectorSet<V,M>& m_set2;
};

}  // End namespace QUESO

#endif // UQ_INTERSECTION_SUBSET_H
