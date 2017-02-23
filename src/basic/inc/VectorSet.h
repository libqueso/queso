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

#ifndef UQ_VECTOR_SET_H
#define UQ_VECTOR_SET_H

#include <queso/Environment.h>
#include <queso/Defines.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*! \file VectorSet.h
 * \brief A templated class for handling sets.
 *
 * \class VectorSet
 * \brief A templated class for handling sets.
 *
 * This class allows the mathematical definition of a scalar function such as:
 * \f$ \pi: B \subset R^n \rightarrow R \f$, since it requires the specification
 * of the domain \f$ B \f$, which is a subset of the vector space \f$ R^n \f$,
 * which is itself a set.*/


template <class V, class M>
class VectorSpace;

template <class V = GslVector, class M = GslMatrix>
class VectorSet
{
public:
  //! @name Constructor/Destructor methods.
  //@{
  //! Shaped constructor.
  /*! Creates a vector set given an environment, a identifying prefix and a volume.*/
  VectorSet(const BaseEnvironment& env, const char* prefix, double volume);

  //! Virtual destructor.
  virtual ~VectorSet();
  //@}

  //! @name Environment methods
  //@{
  //! Environment.  Access to private attribute m_env.
  const BaseEnvironment&  env        ()                 const;

  //! Access to private attribute m_prefix.
  const std::string&             prefix     ()                 const;
  //@}

  //! @name Mathematical methods.
  //@{
  //! Set volume; access to private attribute m_volume.
  double                   volume     ()                 const;

  //! Vector space to which \c this set belongs to. See template specialization.
  virtual const VectorSpace<V,M>& vectorSpace()                 const = 0;

  //! Checks whether a set contains vector \c vec. See template specialization.
  virtual       bool                     contains   (const V& vec)     const = 0;

  //! Returns the set centroid in the vector \c vec. See template specialization.
  virtual       void                     centroid   (V& vec)     const = 0;

  //! Returns the set moments of inertia in the matrix \c mat. See template specialization.
  virtual       void                     moments    (M& mat)     const = 0;
  //@}

  //! @name I/O methods.
  //@{
  //! Prints nothing.
  virtual       void                     print      (std::ostream& os) const;
  friend std::ostream & operator<<(std::ostream & os,
      const VectorSet<V, M> & obj) {
    obj.print(os);
    return os;
  }
  //@}

protected:
  const BaseEnvironment& m_env;
        std::string             m_prefix;
        double                  m_volume;

private:
  //! Default Constructor
  /*! It should not be used by the user.*/
  VectorSet();
};

}  // End namespace QUESO

#endif // UQ_VECTOR_SET_H
