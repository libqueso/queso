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

#ifndef __UQ_VECTOR_SPACE_H__
#define __UQ_VECTOR_SPACE_H__

#include <uqDistArray.h>
#include <uqMap.h>
#include <uqVectorSet.h>
#include <cmath>
//#include <math.h>

/*!
 * \class uqVectorSpaceClass
 * \brief A class representing a vector space.
 *
 * Template classes \c V and \c M are to represent a vector class and a matrix class 
 * respectively. Currently (as of version 0.46.0) QUESO has matrix and vector classes 
 * implemented using either GSL or Trilinos-Teuchos libraries. */

template <class V, class M>
class uqVectorSpaceClass : public uqVectorSetClass<V,M>
{
public:
  //! @name Constructor/Destructor methods
  //@{

  //! Default constructor
  uqVectorSpaceClass();

  //! Shaped constructor.
  /*! Construct a vector space with QUESO environment \c env and of dimension \c dimGlobalValue.*/
  uqVectorSpaceClass(const uqBaseEnvironmentClass&   env,
                     const char*                     prefix,
                     unsigned int                    dimGlobalValue,
                     const std::vector<std::string>* componentsNamesVec);
  
  //! Copy constructor.
  uqVectorSpaceClass(const uqVectorSpaceClass<V,M>&  aux);

  //! Destructor
  ~uqVectorSpaceClass();
  //@}

  
  //! @name Attribute methods
  //@{
  //! Environment.
  const uqBaseEnvironmentClass&  env                     () const;
  
  //! Map.
  const uqMapClass&              map                     () const;
  
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
  
  //! Creates an empty vector of size given by uqMapClass& map. See template specialization.
  V*                             newVector               () const; // See template specialization
        
  //! Creates a vector of size given by uqMapClass& map and all values give by \c value. See template specialization
  V*                             newVector               (double value) const; // See template specialization
        
  //! Creates vector as a copy of another.
  V*                             newVector               (const V& v) const;
	
  //! Creates an empty matrix of size given by uqMapClass& map. See template specialization.	
  M*                             newMatrix               () const; // See template specialization
  
  //! Creates a diagonal matrix with the elements and size of vector \c v.
  M*                             newDiagMatrix           (const V& v) const;
  
  //! Creates a diagonal matrix with the elements \c diagValue and size given by uqMapClass& map. See template specialization.	
  M*                             newDiagMatrix           (double diagValue) const; // See template specialization
  
  //! Creates a diagonal matrix conditionally to values from vector \c varVec, guaranteeing that its values are neither 0, NAN nor INFINITY.
  /*! If varVec[i] is either 0, NAN or INFINITY, then this method tries to assign the value (*auxVec)[i])^2
   * to matrix(i,i). Case (*auxVec)[i])^2 is either NAN or INFINITY, then matrix(i,i)=1.*/
  M*                             newProposalMatrix       (const V* varVec, const V* auxVec) const;

  //! Accessor method to \c this. Vector space to which \c this vector set belongs to.
  /*! It is virtual in the base class 'uqVectorSetClass'*/
  const uqVectorSpaceClass<V,M>&       vectorSpace             () const; // 
  
  //! Whether \this vector contains vector \c vec.
  bool                           contains                (const V& vec) const;

  //! Access to private attribute m_componentsNamesArray, which is an instance of uqDistArrayClass.
  const uqDistArrayClass<std::string>* componentsNamesArray    () const;
  
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
  uqMapClass*                    newMap                  (); // See template specialization

  using uqVectorSetClass<V,M>::m_env;
  using uqVectorSetClass<V,M>::m_prefix;
  using uqVectorSetClass<V,M>::m_volume;

  //! Global dimension.
  unsigned int                   m_dimGlobal;
  
  //! Map.
  const uqMapClass*              m_map;
  
  //! Local dimension (number of elements owned by the calling processor.).
  unsigned int                   m_dimLocal;
  
  //! Array of strings of the type uqDistArrayClass to store the names of the components
  uqDistArrayClass<std::string>* m_componentsNamesArray;
  
  //! Vector of strings of the type uqDistArrayClass to store the names of the components
  uqDistArrayClass<std::string>* m_componentsNamesVec;
  
  //! Empty string for the components names.
  std::string                    m_emptyComponentName;

  //! A vector of all elements equal to zero. 
  V*                             m_zeroVector;
};

//--------------------------------------------------------
// Constructor/Destructor methods ------------------------
// Default constructor -----------------------------------
template <class V, class M>
uqVectorSpaceClass<V,M>::uqVectorSpaceClass()
  :
  uqVectorSetClass<V,M>()
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.worldRank(),
                      "uqVectorSpaceClass<V,M>::constructor(), default",
                      "should not be used by user");
}
// Shaped constructor -----------------------------------
template <class V, class M>
uqVectorSpaceClass<V,M>::uqVectorSpaceClass(
  const uqBaseEnvironmentClass&   env,
  const char*                     prefix,
        unsigned int              dimGlobalValue,
  const std::vector<std::string>* componentsNamesVec)
  :
  uqVectorSetClass<V,M> (env,((std::string)(prefix) + "space_").c_str(),INFINITY),
  m_dimGlobal           (dimGlobalValue),
  m_map                 (newMap()),
  m_dimLocal            (m_map->NumMyElements()),
  m_componentsNamesArray(NULL),
  m_componentsNamesVec  (NULL),
  m_emptyComponentName  (""),
  m_zeroVector          (new V(m_env,*m_map))
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqVectorSpaceClass<V,M>::constructor(1)"
                            << ", with m_prefix = "                << m_prefix
                            << "\n  m_zeroVector->sizeGlobal() = " << m_zeroVector->sizeGlobal()
                            << "\n  m_dimGlobal                = " << m_dimGlobal
                            << "\n  m_zeroVector->sizeLocal()  = " << m_zeroVector->sizeLocal()
                            << "\n  m_dimLocal                 = " << m_dimLocal
                            << "\n  m_map->NumGlobalElements() = " << m_map->NumGlobalElements()
                            << "\n  componentsNamesVec         = " << componentsNamesVec
                            << std::endl;
  }
 
  if (m_zeroVector->sizeGlobal() != m_dimGlobal) {
    std::cerr << "In uqVectorSpaceClass<V,M>::constructor(1)"
              << ", with m_prefix = " << m_prefix
              << ": m_zeroVector->sizeGlobal() = " << m_zeroVector->sizeGlobal()
              << ", m_dimGlobal = "                << m_dimGlobal
              << std::endl;
  }
  UQ_FATAL_TEST_MACRO((m_zeroVector->sizeGlobal() != m_dimGlobal),
                      m_env.worldRank(),
                      "uqVectorSpaceClass<V,M>::constructor(1)",
                      "global size of 'm_zeroVector' is not equal to m_dimGlobal");

  if (m_zeroVector->sizeLocal() != m_dimLocal) {
    std::cerr << "In uqVectorSpaceClass<V,M>::constructor(1)"
              << ", with m_prefix = " << m_prefix
              << ": m_zeroVector->sizeLocal() = " << m_zeroVector->sizeLocal()
              << ", m_dimLocal = "                << m_dimLocal
              << std::endl;
  }
  UQ_FATAL_TEST_MACRO((m_zeroVector->sizeLocal() != m_dimLocal),
                      m_env.worldRank(),
                      "uqVectorSpaceClass<V,M>::constructor(1)",
                      "local size of 'm_zeroVector' is not equal to m_dimLocal");

  if (componentsNamesVec != NULL) {
    UQ_FATAL_TEST_MACRO((componentsNamesVec->size() != (size_t) m_dimGlobal),
                        m_env.worldRank(),
                        "uqVectorSpaceClass<V,M>::constructor(1)",
                        "global size of 'componentsNames' is not equal to m_dimGlobal");

    m_componentsNamesArray = new uqDistArrayClass<std::string>(*m_map,1);
    unsigned int myFirstId = this->globalIdOfFirstComponent();
    for (unsigned int i = 0; i < m_dimLocal; ++i) {
      (*m_componentsNamesArray)(i,0) = (*componentsNamesVec)[myFirstId+i];
    }

    UQ_FATAL_TEST_MACRO((m_componentsNamesArray->GlobalLength() != (int) m_dimGlobal),
                        m_env.worldRank(),
                        "uqVectorSpaceClass<V,M>::constructor(1)",
                        "global size of 'm_componentsNamesArray' is not equal to m_dimGlobal");
    UQ_FATAL_TEST_MACRO((m_componentsNamesArray->MyLength() != (int) m_dimLocal),
                        m_env.worldRank(),
                        "uqVectorSpaceClass<V,M>::constructor(1)",
                        "local size of 'm_componentsNamesArray' is not equal to m_dimLocal");
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqVectorSpaceClass<V,M>::constructor(1)"
                            << ", with m_prefix = " << m_prefix
                            << std::endl;
  }
}
// Copy constructor --------------------------------------
template <class V, class M>
uqVectorSpaceClass<V,M>::uqVectorSpaceClass(const uqVectorSpaceClass<V,M>& aux)
  :
  uqVectorSetClass<V,M> (aux.env(),((std::string)(aux.m_prefix)).c_str(),INFINITY),
  m_dimGlobal           (aux.m_dimGlobal),
  m_map                 (newMap()),
  m_dimLocal            (m_map->NumMyElements()),
  m_componentsNamesArray(NULL),
  m_componentsNamesVec  (NULL),
  m_emptyComponentName  (""),
  m_zeroVector          (new V(m_env,*m_map))
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqVectorSpaceClass<V,M>::constructor(2)"
                            << ": aux.m_componentsNamesArray = " << aux.m_componentsNamesArray
                            << ", aux.m_componentsNamesVec = "   << aux.m_componentsNamesVec
                            << std::endl;
  }

  if (aux.m_componentsNamesArray != NULL) {
    m_componentsNamesArray = new uqDistArrayClass<std::string>(*(aux.m_componentsNamesArray));
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqVectorSpaceClass<V,M>::constructor(2)"
                            << std::endl;
  }
}
// Destructor --------------------------------------------
template <class V, class M>
uqVectorSpaceClass<V,M>::~uqVectorSpaceClass()
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqVectorSpaceClass<V,M>::destructor()"
                            << std::endl;
  }

  if (m_zeroVector           != NULL) delete m_zeroVector;
  if (m_componentsNamesVec   != NULL) delete m_componentsNamesVec;
  if (m_componentsNamesArray != NULL) delete m_componentsNamesArray;
  if (m_map                  != NULL) delete m_map;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqVectorSpaceClass<V,M>::destructor()"
                            << std::endl;
  }
}

// -------------------------------------------------------
// Attribute methods--------------------------------------
template <class V, class M>
const uqBaseEnvironmentClass&
uqVectorSpaceClass<V,M>::env() const
{
  return m_env;
}
// -------------------------------------------------------
template <class V, class M>
const uqMapClass&
uqVectorSpaceClass<V,M>::map() const
{
  UQ_FATAL_TEST_MACRO(m_map == NULL,
                      m_env.worldRank(),
                      "uqVectorSpaceClass<V,M>::map()",
                      "m_map is still NULL");
  return *m_map;
}
// -------------------------------------------------------
template <class V, class M>
unsigned int
uqVectorSpaceClass<V,M>::numOfProcsForStorage() const
{
  return m_map->Comm().NumProc();
}
// -------------------------------------------------------
template <class V, class M>
unsigned int
uqVectorSpaceClass<V,M>::dimLocal() const
{
  return m_dimLocal;
}
// -------------------------------------------------------
template <class V, class M>
unsigned int
uqVectorSpaceClass<V,M>::dimGlobal() const
{
  return m_dimGlobal;
}
// -------------------------------------------------------
template <class V, class M>
unsigned int
uqVectorSpaceClass<V,M>::globalIdOfFirstComponent() const
{
  return m_map->MinMyGID();
}


template<class V, class M>
const V&
uqVectorSpaceClass<V,M>::zeroVector() const
{
  UQ_FATAL_TEST_MACRO(m_zeroVector == NULL,
                      m_env.worldRank(),
                      "uqVectorSpaceClass<V,M>::zeroVector()",
                      "m_zeroVector is still NULL");
  return *m_zeroVector;
}
// -------------------------------------------------------
template <class V, class M>
V*
uqVectorSpaceClass<V,M>::newVector(const V& v) const
{
  if (v.sizeGlobal() != m_dimGlobal) return NULL;
  if (v.sizeLocal () != m_dimLocal ) return NULL;

  return new V(v);
}

template <class V, class M>
M*
uqVectorSpaceClass<V,M>::newDiagMatrix(const V& v) const
{
  if (v.sizeGlobal() != m_dimGlobal) return NULL;
  if (v.sizeLocal () != m_dimLocal ) return NULL;

  return new M(v);
}
// -------------------------------------------------------
template <class V, class M>
M*
uqVectorSpaceClass<V,M>::newProposalMatrix(
  const V* varVec,
  const V* auxVec) const
{
  V tmpVec(*m_zeroVector);
  for (unsigned int i = 0; i < m_dimLocal; ++i) {
    double variance = INFINITY;
    if (varVec) variance = (*varVec)[i];
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqVectorSpaceClass<V,M>::newProposalMatrix()"
                              << ": i = "        << i
                              << ", variance = " << variance
                              << std::endl;
    }
    if ((variance == INFINITY) ||
        (variance == NAN     )) {
      if (auxVec) {
        tmpVec[i] = std::pow( fabs((*auxVec)[i])*0.05,2. );
        if ((tmpVec[i] == 0.      ) ||
            (tmpVec[i] == INFINITY) ||
            (tmpVec[i] == NAN     )) {
          tmpVec[i] = 1.;
        }
      }
      else {
        tmpVec[i] = 1.;
      }
    }
    else if (variance == 0.) {
      tmpVec[i] = 1.;
    }
    else {
      tmpVec[i] = variance;
    }
  }

  return newDiagMatrix(tmpVec);
}
// -------------------------------------------------------
template <class V, class M>
const uqVectorSpaceClass<V,M>&
uqVectorSpaceClass<V,M>::vectorSpace() const
{
  return *this;
}
// -------------------------------------------------------
template <class V, class M>
bool
uqVectorSpaceClass<V,M>::contains(const V& vec) const
{
  if (vec[0]) {}; // just to remove compiler warning
  return true;
}
// -------------------------------------------------------
template <class V, class M>
const uqDistArrayClass<std::string>* 
uqVectorSpaceClass<V,M>::componentsNamesArray() const
{
  return m_componentsNamesArray;
}
// -------------------------------------------------------
template <class V, class M>
const std::string&
uqVectorSpaceClass<V,M>::localComponentName(unsigned int localComponentId) const
{
  if (m_componentsNamesArray == NULL) return m_emptyComponentName;

  UQ_FATAL_TEST_MACRO(localComponentId > m_dimLocal,
                      m_env.worldRank(),
                      "uqVectorSpaceClass<V,M>::localComponentName()",
                      "localComponentId is too big");

//return (*(const_cast<uqDistArrayClass<std::string>*>(m_componentsNamesArray)))(localComponentId,0);
  return (*m_componentsNamesArray)(localComponentId,0);
}
// -------------------------------------------------------
template<class V, class M>
void
uqVectorSpaceClass<V,M>::printComponentsNames(std::ostream& os, bool printHorizontally) const
{
  if (printHorizontally) { 
    for (unsigned int i = 0; i < this->dimLocal(); ++i) {
      os << "'" << this->localComponentName(i) << "'"
         << " ";
    }
  }
  else {
    for (unsigned int i = 0; i < this->dimLocal(); ++i) {
      os << "'" << this->localComponentName(i) << "'"
         << std::endl;
    }
  }

  return;
}
// -------------------------------------------------------
template <class V, class M>
void
uqVectorSpaceClass<V,M>::print(std::ostream& os) const
{
  os << "In uqVectorSpaceClass<V,M>::print()"
     << ": nothing to be printed" << std::endl;
  return;
}
#endif // __UQ_VECTOR_SPACE_H__
