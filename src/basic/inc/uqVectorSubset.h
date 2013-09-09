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

#ifndef __UQ_VECTOR_SUBSET_H__
#define __UQ_VECTOR_SUBSET_H__

#include <uqVectorSpace.h>

namespace QUESO {

/*! \file uqVectorSubsetClass.h
 * \brief A templated class for handling subsets.
 * 
 * \class uqVectorSubsetClass
 * \brief A templated class for handling subsets.
 *
 * This class specifies a subset. The most common example of a subset is a subset of a vector space,
 * present, for instance, in the definition of a scalar function: \f$ \pi: B \subset R^n \rightarrow R \f$.
 * \f$ B \f$ is a subset of the set \f$ R^n \f$, which is also a vector space. */


//*****************************************************
// Base class
//*****************************************************
template <class V, class M>
class uqVectorSubsetClass : public uqVectorSetClass<V,M>
{
public:
  //! @name Constructor/Destructor methods.
  //@{ 
  //! Default Constructor
  /*! It should not be used by the user.*/
  uqVectorSubsetClass();
  
  //! Shaped constructor (with volume).
  uqVectorSubsetClass(const char* prefix, const uqVectorSpaceClass<V,M>& vectorSpace, double volume);
  
  //! Destructor.
  virtual ~uqVectorSubsetClass();
  //@}
  
    //! @name Mathematical methods.
  //@{
  //!  Vector space to which \c this set belongs to. See template specialization.
  const uqVectorSpaceClass<V,M>& vectorSpace()                 const;
  
  //! Returns whether \c this contains vector \c vec. See template specialization.
  virtual        bool                     contains   (const V& vec)     const = 0;
  //@}
  
  //! @name I/O methods.
  //@{
  //! Prints nothing.
  virtual        void                     print      (std::ostream& os) const;
  //@}
  
protected:
  using uqVectorSetClass<V,M>::m_env;
  using uqVectorSetClass<V,M>::m_prefix;

  const uqVectorSpaceClass<V,M>* m_vectorSpace;
};

// Default constructor-------------------------------
template <class V, class M>
uqVectorSubsetClass<V,M>::uqVectorSubsetClass()
  :
  uqVectorSetClass<V,M>(),
  m_vectorSpace        (NULL)
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.worldRank(),
                      "uqVectorSubsetClass<V,M>::constructor(), default",
                      "should not be used by user");
}
// Shaped constructor--------------------------------
template <class V, class M>
uqVectorSubsetClass<V,M>::uqVectorSubsetClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& vectorSpace,
  double                         volume)
  :
  uqVectorSetClass<V,M>(vectorSpace.env(),prefix,volume),
  m_vectorSpace        (&vectorSpace)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqVectorSubsetClass<V,M>::constructor()"
              << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqVectorSubsetClass<V,M>::constructor()"
              << std::endl;
  }
}
// Destructor ---------------------------------------
template <class V, class M>
uqVectorSubsetClass<V,M>::~uqVectorSubsetClass()
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering uqVectorSubsetClass<V,M>::destructor()"
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving uqVectorSubsetClass<V,M>::destructor()"
                            << std::endl;
  }
}
// Math methods -------------------------------------
template <class V, class M>
const uqVectorSpaceClass<V,M>&
uqVectorSubsetClass<V,M>::vectorSpace() const
{
  return *m_vectorSpace;
}
// I/O methods---------------------------------------
template <class V, class M>
void
uqVectorSubsetClass<V,M>::print(std::ostream& os) const
{
  os << "In uqVectorSubsetClass<V,M>::print()"
     << ": nothing to be printed" << std::endl;
  return;
}

//*****************************************************
// Box class
//*****************************************************

/*!
 * \class uqBoxSubsetClass
 * \brief Class representing a subset of a vector space shaped like a hypercube
 *
 * This class is determined by a collection of upper and lower limits of the hypercube.
 * (line segment in \f$ R \f$, rectangle in \f$ R^2 \f$, and so on). */

template<class V, class M>
class uqBoxSubsetClass : public uqVectorSubsetClass<V,M> {
public:
  //! @name Constructor/Destructor methods
  //@{

  //! Shaped, default constructor. 
  /*! Construct a subspace of \c vectorSpace, with min and max values given by the vectors \c minValues 
   * and \c maxValues, respectively. It checks for possible inconsistencies between the values stored in
   * \c minValues and \c maxValues, and calculates the volume of the box subset, assigning it to m_volume. */
  uqBoxSubsetClass(const char*                    prefix,
                   const uqVectorSpaceClass<V,M>& vectorSpace,
                   const V&                       minValues,
                   const V&                       maxValues);

  //! Destructor
  ~uqBoxSubsetClass();
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
  using uqVectorSetClass   <V,M>::m_env;
  using uqVectorSetClass   <V,M>::m_prefix;
  using uqVectorSetClass   <V,M>::m_volume;
  using uqVectorSubsetClass<V,M>::m_vectorSpace;
 
  //! Vector of templated type \c V to store the minimum values of the box subset class.
  V m_minValues;
  
  //! Vector of templated type \c V to store the maximum values of the box subset class.
  V m_maxValues;
};

// Default, shaped constructor ----------------------
template<class V, class M>
uqBoxSubsetClass<V,M>::uqBoxSubsetClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& vectorSpace,
  const V&                       minValues,
  const V&                       maxValues)
  :
  uqVectorSubsetClass<V,M>(prefix,vectorSpace,0.),
  m_minValues(minValues),
  m_maxValues(maxValues)
{
  UQ_FATAL_TEST_MACRO(minValues.sizeLocal() != maxValues.sizeLocal(),
                      m_env.worldRank(),
                      "uqBoxSubsetClass<V,M>::uqBoxSubsetClass()",
                      "vectors 'minValues' and 'maxValues' should have the same size");
  UQ_FATAL_TEST_MACRO(minValues.sizeLocal() != vectorSpace.dimLocal(),
                      m_env.worldRank(),
                      "uqBoxSubsetClass<V,M>::uqBoxSubsetClass()",
                      "sizes of vectors 'minValues' and 'maxValues' should be equal to dimension of the vector space");
  for (unsigned int i = 0; i < m_vectorSpace->dimLocal(); ++i) {
    UQ_FATAL_TEST_MACRO(minValues[i] > maxValues[i],
                        m_env.worldRank(),
                        "uqBoxSubsetClass<V,M>::uqBoxSubsetClass()",
                        "it should happen minValue <= maxValue for all dimensions");
  }

  m_volume = 1.;
  for (unsigned int i = 0; i < m_vectorSpace->dimLocal(); ++i) {
    m_volume *= (m_maxValues[i] - m_minValues[i]);
  }
}
// Destructor ---------------------------------------
template<class V, class M>
uqBoxSubsetClass<V,M>::~uqBoxSubsetClass()
{
}
// Math methods -------------------------------------
template<class V, class M>
bool
uqBoxSubsetClass<V,M>::contains(const V& vec) const
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
uqBoxSubsetClass<V,M>::minValues() const
{
  return m_minValues;
}
// --------------------------------------------------
template<class V, class M>
const V&
uqBoxSubsetClass<V,M>::maxValues() const
{
  return m_maxValues;
}
// I/O method ---------------------------------------
template <class V, class M>
void
uqBoxSubsetClass<V,M>::print(std::ostream& os) const
{
  os << "In uqBoxSubsetClass<V,M>::print()"
     << ": m_minValues = " << m_minValues
     << ", m_maxValues = " << m_maxValues
     << ", m_volume = "    << m_volume
     << std::endl;

  return;
}


//*****************************************************
// Intersection class
//*****************************************************

/*! \class uqIntersectionSubsetClass
 * \brief A templated class representing the intersection of two vector sets.
 *
 * This class is used to determine if a vector  belongs to the intersection of
 * two vector sets. It is useful for handling a posterior PDF, since its domain 
 * is the intersection of the domain of the prior PDF with the domain of the 
 * likelihood function.*/

template<class V, class M>
class uqIntersectionSubsetClass : public uqVectorSubsetClass<V,M> {
public:
  //! @name Constructor/Destructor methods.
  //@{ 
  //! Default, shaped constructor.
  /*! Creates the class for the intersection of two vector sets, given a vector space, its volume and the
   * sets.*/
  uqIntersectionSubsetClass(const char*                    prefix,
                            const uqVectorSpaceClass<V,M>& vectorSpace,
                                  double                   volume,
                            const uqVectorSetClass<V,M>&   set1,
                            const uqVectorSetClass<V,M>&   set2);
  
  //! Destructor.
 ~uqIntersectionSubsetClass();
  //@}
 
   //! @name Mathematical methods.
  //@{ 
 //! Determines whether both sets m_set1 and m_set2 (class' private attributes) contain vector \c vec.
  bool contains(const V& vec)     const;
  //@}
  
  //! @name I/O methods.
  //@{ 
  //! Prints both subsets (via protected attributes m_set1 and m_set2).
  void print   (std::ostream& os) const;
  //@}
  
protected:
  using uqVectorSetClass   <V,M>::m_env;
  using uqVectorSetClass   <V,M>::m_prefix;
  using uqVectorSetClass   <V,M>::m_volume;
  using uqVectorSubsetClass<V,M>::m_vectorSpace;
 
  //! Vector set: m_set1. 
  /*! We seek the intersection of vectors set m_set1 and m_set2.*/
  const uqVectorSetClass<V,M>& m_set1;
  
  //! Vector set: m_set2.
  const uqVectorSetClass<V,M>& m_set2;
};
// --------------------------------------------------
// Constructor/Destructor methods -------------------
// Default, shaped constructor ----------------------
template<class V, class M>
uqIntersectionSubsetClass<V,M>::uqIntersectionSubsetClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& vectorSpace,
        double                   volume,
  const uqVectorSetClass<V,M>&   set1,
  const uqVectorSetClass<V,M>&   set2)
  :
  uqVectorSubsetClass<V,M>(prefix,vectorSpace,volume),
  m_set1                  (set1),
  m_set2                  (set2)
{
}
// Destructor --------------------------------------------
template<class V, class M>
uqIntersectionSubsetClass<V,M>::~uqIntersectionSubsetClass()
{
}
// Mathematical methods-----------------------------------
template<class V, class M>
bool
uqIntersectionSubsetClass<V,M>::contains(const V& vec) const
{
  return (m_set1.contains(vec) && m_set2.contains(vec));
}
// I/O methods--------------------------------------------
template <class V, class M>
void
uqIntersectionSubsetClass<V,M>::print(std::ostream& os) const
{
  os << "In uqIntersectionSubsetClass<V,M>::print()"
     << ": m_set1 = " << m_set1
     << ", m_set2 = " << m_set2
     << std::endl;

  return;
}

//*****************************************************
// Concatenation class
//*****************************************************

/*! \class uqConcatenationSubsetClass
 * \brief A templated class representing the concatenation of two vector subsets.
 *
 * This class is used to represent the concatenation of two subsets. It allows concatenation
 * of both two subsets as well as a collection of subsets into one.
 */
template<class V, class M>
class uqConcatenationSubsetClass : public uqVectorSubsetClass<V,M> {
public:
   //! @name Constructor/Destructor methods.
  //@{ 
  //! Constructor - two sets only
  /*! It concatenates set1 and set2 into a new set, m_set, of volume given by
   * set1.volume()*set2.volume(). */
  uqConcatenationSubsetClass(const char*                    prefix,
                             const uqVectorSpaceClass<V,M>& vectorSpace,
                             const uqVectorSetClass<V,M>&   set1,
                             const uqVectorSetClass<V,M>&   set2);
  
  //! Constructor - collection of sets
  /*! It concatenates a collection of n subsets (using the std::vector to represent such collection)
   * into a new set, m_set, of volume given by sets[0].volume()*sets[1].volume()*...*sets[n-1].volume(). */
  uqConcatenationSubsetClass(const char*                                       prefix,
                             const uqVectorSpaceClass<V,M>&                    vectorSpace,
                             double                                            volume,
                             const std::vector<const uqVectorSetClass<V,M>* >& sets);
  //! Destructor
  ~uqConcatenationSubsetClass();
  //@}
  
  //! @name Mathematical methods.
  //@{ 
  //! Determines whether each one of the subsets m_sets (class' private attributes) contains vector \c vec.
  bool contains(const V& vec)     const;
  //@}
  
  //! @name I/O methods.
  //@{ 
  //! Prints the subsets (via protected attribute m_sets).
  void print   (std::ostream& os) const;
  //@}
protected:
  using uqVectorSetClass   <V,M>::m_env;
  using uqVectorSetClass   <V,M>::m_prefix;
  using uqVectorSetClass   <V,M>::m_volume;
  using uqVectorSubsetClass<V,M>::m_vectorSpace;

  std::vector<const uqVectorSetClass<V,M>* > m_sets;
};

// --------------------------------------------------
// Constructor/Destructor methods -------------------
// Default, shaped constructor ----------------------
template<class V, class M>
uqConcatenationSubsetClass<V,M>::uqConcatenationSubsetClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& vectorSpace,
  const uqVectorSetClass<V,M>&   set1,
  const uqVectorSetClass<V,M>&   set2)
  :
  uqVectorSubsetClass<V,M>(prefix,vectorSpace,set1.volume()*set2.volume()),
  m_sets                  (2,(const uqVectorSetClass<V,M>*) NULL)
{
  m_sets[0] = &set1;
  m_sets[1] = &set2;
}
// Default, shaped constructor ----------------------
template<class V, class M>
uqConcatenationSubsetClass<V,M>::uqConcatenationSubsetClass(
  const char*                                       prefix,
  const uqVectorSpaceClass<V,M>&                    vectorSpace,
  double                                            volume,
  const std::vector<const uqVectorSetClass<V,M>* >& sets)
  :
  uqVectorSubsetClass<V,M>(prefix,vectorSpace,volume),
  m_sets                  (sets.size(),(const uqVectorSetClass<V,M>*) NULL)
{
  for (unsigned int i = 0; i < m_sets.size(); ++i) {
    m_sets[i] = sets[i];
  }
}
// Destructor ---------------------------------------
template<class V, class M>
uqConcatenationSubsetClass<V,M>::~uqConcatenationSubsetClass()
{
}
// Mathematical methods ------------------------------
template<class V, class M>
bool
uqConcatenationSubsetClass<V,M>::contains(const V& vec) const
{
  bool result = true;

  std::vector<V*> vecs(m_sets.size(),(V*) NULL);
  for (unsigned int i = 0; i < vecs.size(); ++i) {
    vecs[i] = new V(m_sets[i]->vectorSpace().zeroVector());
  }

  unsigned int cummulativeSize = 0;
  for (unsigned int i = 0; i < vecs.size(); ++i) {
    vec.cwExtract(cummulativeSize,*(vecs[i]));
    cummulativeSize += vecs[i]->sizeLocal();
  }

  UQ_FATAL_TEST_MACRO(vec.sizeLocal() != cummulativeSize,
                      m_env.worldRank(),
                      "uqConcatenationSubsetClass<V,M>::contains()",
                      "incompatible vector sizes");

  for (unsigned int i = 0; i < m_sets.size(); ++i) {
    result = result && m_sets[i]->contains(*(vecs[i]));
  }
  for (unsigned int i = 0; i < vecs.size(); ++i) {
    delete vecs[i];
  }

  return (result);
}
// I/O methods --------------------------------------
template <class V, class M>
void
uqConcatenationSubsetClass<V,M>::print(std::ostream& os) const
{
  os << "In uqConcatenationSubsetClass<V,M>::print()"
     << ": m_sets.size() = " << m_sets.size()
     << std::endl;
  for (unsigned int i = 0; i < m_sets.size(); ++i) {
    os << "m_sets[" << i << "] = " << *(m_sets[i]);
    if (i < (m_sets.size()-1)) {
      os << ", ";
    }
  }
  os << std::endl;

  return;
}

//*****************************************************
// Discrete class
//*****************************************************
/*! \class uqDiscreteSubsetClass
 * \brief A templated class representing the discrete vector subsets.
 *
 * This class is used to represent the a discrete vector subset. Here the 
 * notion of volume does not apply, instead there is the total number of 
 * elements (vectors) that belongs to the subset.
 */

template<class V, class M>
class uqDiscreteSubsetClass : public uqVectorSubsetClass<V,M> {
public:
  //! @name Constructor/Destructor methods.
  //@{ 
  //! Default Constructor
  /*! It constructs a class object given the prefix, vector space to which it
   * belongs  and its number of elements.*/
  uqDiscreteSubsetClass(const char*                    prefix,
                        const uqVectorSpaceClass<V,M>& vectorSpace,
                        const std::vector<V*>&         elements);
  //! Destructor
  ~uqDiscreteSubsetClass();
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
  using uqVectorSetClass   <V,M>::m_env;
  using uqVectorSetClass   <V,M>::m_prefix;
  using uqVectorSetClass   <V,M>::m_volume;
  using uqVectorSubsetClass<V,M>::m_vectorSpace;
 
  //! Number of elements in the discrete vector subset.
  std::vector<V*> m_elements;
};

// --------------------------------------------------
// Constructor/Destructor methods -------------------
// Default constructor-------------------------------
template<class V, class M>
uqDiscreteSubsetClass<V,M>::uqDiscreteSubsetClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& vectorSpace,
  const std::vector<V*>&         elements)
  :
  uqVectorSubsetClass<V,M>(prefix,vectorSpace,0.),
  m_elements(elements.size(),NULL)
{
  m_volume = 0.;
  UQ_FATAL_TEST_MACRO(true,
                      m_env.worldRank(),
                      "uqDiscreteSubsetClass<V,M>::contains()",
                      "incomplete code");
}
// Destructor --------------------------------------------
template<class V, class M>
uqDiscreteSubsetClass<V,M>::~uqDiscreteSubsetClass()
{
}
// Mathematical methods-----------------------------------
template<class V, class M>
bool
uqDiscreteSubsetClass<V,M>::contains(const V& vec) const
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.worldRank(),
                      "uqDiscreteSubsetClass<V,M>::contains()",
                      "incomplete code");

  return false;
}
// I/O methods--------------------------------------------
template <class V, class M>
void
uqDiscreteSubsetClass<V,M>::print(std::ostream& os) const
{
  os << "In uqDiscreteSubsetClass<V,M>::print()"
     << ": nothing to print"
     << std::endl;

  return;
}

}  // End namespace QUESO

#endif // __UQ_VECTOR_SUBSET_H__
