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

#ifndef UQ_VECTOR_SUBSET_H
#define UQ_VECTOR_SUBSET_H

#include <queso/VectorSet.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*! \file VectorSubset.h
 * \brief A templated class for handling subsets.
 *
 * \class VectorSubset
 * \brief A templated class for handling subsets.
 *
 * This class specifies a subset. The most common example of a subset is a subset of a vector space,
 * present, for instance, in the definition of a scalar function: \f$ \pi: B \subset R^n \rightarrow R \f$.
 * \f$ B \f$ is a subset of the set \f$ R^n \f$, which is also a vector space. */

template <class V = GslVector, class M = GslMatrix>
class VectorSubset : public VectorSet<V,M>
{
public:
  //! @name Constructor/Destructor methods.
  //@{
  //! Default Constructor
  /*! It should not be used by the user.*/
  VectorSubset();

  //! Shaped constructor (with volume).
  VectorSubset(const char* prefix, const VectorSpace<V,M>& vectorSpace, double volume);

  //! Destructor.
  virtual ~VectorSubset();
  //@}

    //! @name Mathematical methods.
  //@{
  //!  Vector space to which \c this set belongs to. See template specialization.
  const VectorSpace<V,M>& vectorSpace()                 const;

  //! Returns whether \c this contains vector \c vec. See template specialization.
  virtual        bool                     contains   (const V& vec)     const = 0;
  //@}

  //! @name I/O methods.
  //@{
  //! Prints nothing.
  virtual        void                     print      (std::ostream& os) const;
  //@}

protected:
  using VectorSet<V,M>::m_env;
  using VectorSet<V,M>::m_prefix;

  const VectorSpace<V,M>* m_vectorSpace;
};

}  // End namespace QUESO

#endif // UQ_VECTOR_SUBSET_H
