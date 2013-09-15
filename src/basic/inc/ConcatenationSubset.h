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

#ifndef UQ_CONCATENATION_SUBSET_H
#define UQ_CONCATENATION_SUBSET_H

#include <queso/VectorSpace.h>

namespace QUESO {

/*! \class ConcatenationSubset
 * \brief A templated class representing the concatenation of two vector subsets.
 *
 * This class is used to represent the concatenation of two subsets. It allows concatenation
 * of both two subsets as well as a collection of subsets into one.
 */
template<class V, class M>
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

// --------------------------------------------------
// Constructor/Destructor methods -------------------
// Default, shaped constructor ----------------------
template<class V, class M>
ConcatenationSubset<V,M>::ConcatenationSubset(
  const char*                    prefix,
  const VectorSpace<V,M>& vectorSpace,
  const VectorSet<V,M>&   set1,
  const VectorSet<V,M>&   set2)
  :
  VectorSubset<V,M>(prefix,vectorSpace,set1.volume()*set2.volume()),
  m_sets                  (2,(const VectorSet<V,M>*) NULL)
{
  m_sets[0] = &set1;
  m_sets[1] = &set2;
}
// Default, shaped constructor ----------------------
template<class V, class M>
ConcatenationSubset<V,M>::ConcatenationSubset(
  const char*                                       prefix,
  const VectorSpace<V,M>&                    vectorSpace,
  double                                            volume,
  const std::vector<const VectorSet<V,M>* >& sets)
  :
  VectorSubset<V,M>(prefix,vectorSpace,volume),
  m_sets                  (sets.size(),(const VectorSet<V,M>*) NULL)
{
  for (unsigned int i = 0; i < m_sets.size(); ++i) {
    m_sets[i] = sets[i];
  }
}
// Destructor ---------------------------------------
template<class V, class M>
ConcatenationSubset<V,M>::~ConcatenationSubset()
{
}
// Mathematical methods ------------------------------
template<class V, class M>
bool
ConcatenationSubset<V,M>::contains(const V& vec) const
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
                      "ConcatenationSubset<V,M>::contains()",
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
ConcatenationSubset<V,M>::print(std::ostream& os) const
{
  os << "In ConcatenationSubset<V,M>::print()"
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

}  // End namespace QUESO

#endif // UQ_CONCATENATION_SUBSET_H
