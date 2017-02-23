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

#ifndef UQ_VECTOR_SPACE_H
#define UQ_VECTOR_SPACE_H

#include <queso/DistArray.h>
#include <queso/Map.h>
#include <queso/VectorSet.h>
#include <cmath>

namespace QUESO {

class GslVector;
class GslMatrix;

/*!
 * \class VectorSpace
 * \brief A class representing a vector space.
 *
 * Template classes \c V and \c M are to represent a vector class and a matrix class
 * respectively. Currently (as of version 0.46.0) QUESO has matrix and vector classes
 * implemented using either GSL or Trilinos-Teuchos libraries. */

template <class V = GslVector, class M = GslMatrix>
class VectorSpace : public VectorSet<V,M>
{
public:
  //! @name Constructor/Destructor methods
  //@{

  //! Shaped constructor.
  /*! Construct a vector space with QUESO environment \c env and of dimension \c dimGlobalValue.*/
  VectorSpace(const BaseEnvironment&   env,
                     const char*                     prefix,
                     unsigned int                    dimGlobalValue,
                     const std::vector<std::string>* componentsNamesVec);

  //! Destructor
  ~VectorSpace();
  //@}


  //! @name Attribute methods
  //@{
  //! Environment.
  const BaseEnvironment&  env                     () const;

  //! Map.
  const Map&              map                     () const;

  //! Returns total number of processes.
  unsigned int                   numOfProcsForStorage    () const;


  unsigned int                   dimLocal                () const;
  unsigned int                   dimGlobal               () const;
  unsigned int                   globalIdOfFirstComponent() const;
  //@}

  //! @name Mathematical methods
  //@{
  //! Returns a vector filled with zeros
  const V&                       zeroVector              () const;

  //! Creates an empty vector of size given by Map& map. See template specialization.
  V*                             newVector               () const; // See template specialization

  //! Creates a vector of size given by Map& map and all values give by \c value. See template specialization
  V*                             newVector               (double value) const; // See template specialization

  //! Creates vector as a copy of another.
  V*                             newVector               (const V& v) const;

  //! Creates an empty matrix of size given by Map& map. See template specialization.
  M*                             newMatrix               () const; // See template specialization

  //! Creates a diagonal matrix with the elements and size of vector \c v.
  M*                             newDiagMatrix           (const V& v) const;

  //! Creates a diagonal matrix with the elements \c diagValue and size given by Map& map. See template specialization.
  M*                             newDiagMatrix           (double diagValue) const; // See template specialization

  //! Creates a diagonal matrix conditionally to values from vector \c varVec, guaranteeing that its values are neither 0, NAN nor INFINITY.
  /*! If varVec[i] is either 0, NAN or INFINITY, then this method tries to assign the value (*auxVec)[i])^2
   * to matrix(i,i). Case (*auxVec)[i])^2 is either NAN or INFINITY, then matrix(i,i)=1.*/
  M*                             newProposalMatrix       (const V* varVec, const V* auxVec) const;

  //! Accessor method to \c this. Vector space to which \c this vector set belongs to.
  /*! It is virtual in the base class 'VectorSet'*/
  const VectorSpace<V,M>&       vectorSpace             () const; //

  //! Whether \this vector contains vector \c vec.
  bool                           contains                (const V& vec) const;

  //! The (INFINITY/nonexistent) centroid of the space
  void                           centroid                (V& vec) const;

  //! The (INFINITY/nonexistent) matrix of second moments of the space
  void                           moments                 (M& vec) const;

  //! Access to private attribute m_componentsNamesArray, which is an instance of DistArray.
  const DistArray<std::string>* componentsNamesArray    () const;

  //! Access to private attribute m_componentsNamesVec.
  const std::vector<std::string>*      componentsNamesVec      () const;

  //! Returns the local component names.
  const std::string&                   localComponentName      (unsigned int localComponentId) const;
  //@}

  //! @name I/O methods
  //@{
  //! Prints the local component names.
  void                           printComponentsNames    (std::ostream& os, bool printHorizontally) const;

  //! Prints only a message.
  void                           print                   (std::ostream& os) const;
  //@}
protected:
  //! Creates a new map. See template specialization.
  Map*                    newMap                  (); // See template specialization

  using VectorSet<V,M>::m_env;
  using VectorSet<V,M>::m_prefix;
  using VectorSet<V,M>::m_volume;

  //! Global dimension.
  unsigned int                   m_dimGlobal;

  //! Map.
  const Map*              m_map;

  //! Local dimension (number of elements owned by the calling processor.).
  unsigned int                   m_dimLocal;

  //! Array of strings of the type DistArray to store the names of the components
  DistArray<std::string>* m_componentsNamesArray;

  //! Vector of strings of the type DistArray to store the names of the components
  DistArray<std::string>* m_componentsNamesVec;

  //! Empty string for the components names.
  std::string                    m_emptyComponentName;

  //! A vector of all elements equal to zero.
  V*                             m_zeroVector;

private:
  //! Default constructor
  VectorSpace();

  //! Copy constructor.
  VectorSpace(const VectorSpace<V,M>&  aux);
};

}  // End namespace QUESO

#endif // UQ_VECTOR_SPACE_H
