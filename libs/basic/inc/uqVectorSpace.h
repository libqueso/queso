/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of the QUESO Library (Quantification of Uncertainty
 * for Estimation, Simulation and Optimization).
 *
 * QUESO is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QUESO is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QUESO. If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * $Id$
 *
 * Brief description of this file: 
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __UQ_VECTOR_SPACE_H__
#define __UQ_VECTOR_SPACE_H__

#include <uqVectorSet.h>
#include <EpetraExt_DistArray.h>

template <class V, class M>
class uqVectorSpaceClass : public uqVectorSetClass<V,M>
{
public:
        uqVectorSpaceClass();
        uqVectorSpaceClass(const uqBaseEnvironmentClass&   env,
                           const char*                     prefix,
                           unsigned int                    dimGlobalValue,
                           const std::vector<std::string>* componentsNames);
        uqVectorSpaceClass(const uqVectorSpaceClass<V,M>&  aux);
       ~uqVectorSpaceClass();

  const uqBaseEnvironmentClass&            env                 () const;
  const Epetra_Map&                        map                 () const;
        unsigned int                       numOfProcsForStorage() const;
        unsigned int                       dimLocal            () const;
        unsigned int                       dimGlobal           () const;
        unsigned int                       globalIdOfFirstComponent() const;

  const V&                                 zeroVector          () const;
        V*                                 newVector           () const; // See template specialization
        V*                                 newVector           (double value) const; // See template specialization
        V*                                 newVector           (const V& v) const;
        M*                                 newMatrix           () const; // See template specialization
        M*                                 newDiagMatrix       (const V& v) const;
        M*                                 newDiagMatrix       (double diagValue) const; // See template specialization
        M*                                 newProposalMatrix   (const V* varVec,
                                                                const V* auxVec) const;

  const uqVectorSpaceClass<V,M>&           vectorSpace         () const; // It is virtual in the base class 'uqVectorSetClass'
        bool                               contains            (const V& vec) const;

  const EpetraExt::DistArray<std::string>* componentsNames     () const;
  const std::string&                       localComponentName  (unsigned int localComponentId) const;
        void                               printComponentsNames(std::ostream& os, bool printHorizontally) const;
        void                               print               (std::ostream& os) const;

protected:
        Epetra_Map*                        newMap              (); // See template specialization

        using uqVectorSetClass<V,M>::m_env;
        using uqVectorSetClass<V,M>::m_prefix;
        using uqVectorSetClass<V,M>::m_volume;

        unsigned int                       m_dimGlobal;
  const Epetra_Map*                        m_map;
        unsigned int                       m_dimLocal;
        EpetraExt::DistArray<std::string>* m_componentsNames;
        std::string                        m_emptyComponentName;

        V*                                 m_zeroVector;
};

template <class V, class M>
uqVectorSpaceClass<V,M>::uqVectorSpaceClass()
  :
  uqVectorSetClass<V,M>()
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.fullRank(),
                      "uqVectorSpaceClass<V,M>::constructor(), default",
                      "should not be used by user");
}

template <class V, class M>
uqVectorSpaceClass<V,M>::uqVectorSpaceClass(
  const uqBaseEnvironmentClass&   env,
  const char*                     prefix,
        unsigned int              dimGlobalValue,
  const std::vector<std::string>* componentsNames)
  :
  uqVectorSetClass<V,M>(env,((std::string)(prefix) + "space_").c_str(),INFINITY),
  m_dimGlobal          (dimGlobalValue),
  m_map                (newMap()),
  m_dimLocal           (m_map->NumMyElements()),
  m_componentsNames    (NULL),
  m_emptyComponentName (""),
  m_zeroVector         (new V(m_env,*m_map))
{

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering uqVectorSpaceClass<V,M>::constructor()"
                            << ", with m_prefix = " << m_prefix
                            << "\n  m_zeroVector->sizeGlobal() = " << m_zeroVector->sizeGlobal()
                            << "\n  m_dimGlobal                = " << m_dimGlobal
                            << "\n  m_zeroVector->sizeLocal()  = " << m_zeroVector->sizeLocal()
                            << "\n  m_dimLocal                 = " << m_dimLocal
                            << "\n  m_map->NumGlobalElements() = " << m_map->NumGlobalElements()
                            << std::endl;
  }

  if (m_zeroVector->sizeGlobal() != m_dimGlobal) {
    std::cerr << "In uqVectorSpaceClass<V,M>::constructor()"
              << ", with m_prefix = " << m_prefix
              << ": m_zeroVector->sizeGlobal() = " << m_zeroVector->sizeGlobal()
              << ", m_dimGlobal = "                << m_dimGlobal
              << std::endl;
  }
  UQ_FATAL_TEST_MACRO((m_zeroVector->sizeGlobal() != m_dimGlobal),
                      m_env.fullRank(),
                      "uqVectorSpaceClass<V,M>::constructor()",
                      "global size of 'm_zeroVector' is not equal to m_dimGlobal");

  if (m_zeroVector->sizeLocal() != m_dimLocal) {
    std::cerr << "In uqVectorSpaceClass<V,M>::constructor()"
              << ", with m_prefix = " << m_prefix
              << ": m_zeroVector->sizeLocal() = " << m_zeroVector->sizeLocal()
              << ", m_dimLocal = "                << m_dimLocal
              << std::endl;
  }
  UQ_FATAL_TEST_MACRO((m_zeroVector->sizeLocal() != m_dimLocal),
                      m_env.fullRank(),
                      "uqVectorSpaceClass<V,M>::constructor()",
                      "local size of 'm_zeroVector' is not equal to m_dimLocal");

  if (componentsNames != NULL) {
    UQ_FATAL_TEST_MACRO((componentsNames->size() != (size_t) m_dimGlobal),
                        m_env.fullRank(),
                        "uqVectorSpaceClass<V,M>::constructor()",
                        "global size of 'componentsNames' is not equal to m_dimGlobal");

    m_componentsNames = new EpetraExt::DistArray<std::string>(*m_map,1);
    unsigned int myFirstId = this->globalIdOfFirstComponent();
    //EpetraExt::DistArray<std::string>& arrayOfStrings = *m_componentsNames;
    for (unsigned int i = 0; i < m_dimLocal; ++i) {
      (*m_componentsNames)(i,0) = (*componentsNames)[myFirstId+i];
    }

    UQ_FATAL_TEST_MACRO((m_componentsNames->GlobalLength() != (int) m_dimGlobal),
                        m_env.fullRank(),
                        "uqVectorSpaceClass<V,M>::constructor()",
                        "global size of 'm_componentsNames' is not equal to m_dimGlobal");
    UQ_FATAL_TEST_MACRO((m_componentsNames->MyLength() != (int) m_dimLocal),
                        m_env.fullRank(),
                        "uqVectorSpaceClass<V,M>::constructor()",
                        "local size of 'm_componentsNames' is not equal to m_dimLocal");
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving uqVectorSpaceClass<V,M>::constructor()"
                            << ", with m_prefix = " << m_prefix
                            << std::endl;
  }
}

template <class V, class M>
uqVectorSpaceClass<V,M>::uqVectorSpaceClass(const uqVectorSpaceClass<V,M>& aux)
  :
  uqVectorSetClass<V,M>(aux.env(),((std::string)(aux.m_prefix)).c_str(),INFINITY),
  m_dimGlobal          (aux.m_dimGlobal),
  m_map                (newMap()),
  m_dimLocal           (m_map->NumMyElements()),
  m_componentsNames    (NULL),
  m_emptyComponentName (""),
  m_zeroVector         (new V(m_env,*m_map))
{
  if (aux.m_componentsNames != NULL) {
    m_componentsNames = new EpetraExt::DistArray<std::string>(*(aux.m_componentsNames));
  }
}

template <class V, class M>
uqVectorSpaceClass<V,M>::~uqVectorSpaceClass()
{
  //if (m_env.subDisplayFile()) {
  //  *m_env.subDisplayFile() << "Entering uqVectorSpaceClass<V,M>::destructor()"
  //                          << std::endl;
  //}

  if (m_zeroVector      != NULL) delete m_zeroVector;
  if (m_componentsNames != NULL) delete m_componentsNames;
  if (m_map             != NULL) delete m_map;

  //if (m_env.subDisplayFile()) {
  //  *m_env.subDisplayFile() << "Leaving uqVectorSpaceClass<V,M>::destructor()"
  //                          << std::endl;
  //}
}

template <class V, class M>
const uqVectorSpaceClass<V,M>&
uqVectorSpaceClass<V,M>::vectorSpace() const
{
  return *this;
}

template <class V, class M>
bool
uqVectorSpaceClass<V,M>::contains(const V& vec) const
{
  return true;
}

template <class V, class M>
const uqBaseEnvironmentClass&
uqVectorSpaceClass<V,M>::env() const
{
  return m_env;
}

template <class V, class M>
const Epetra_Map&
uqVectorSpaceClass<V,M>::map() const
{
  UQ_FATAL_TEST_MACRO(m_map == NULL,
                      m_env.fullRank(),
                      "uqVectorSpaceClass<V,M>::map()",
                      "m_map is still NULL");
  return *m_map;
}

template <class V, class M>
unsigned int
uqVectorSpaceClass<V,M>::numOfProcsForStorage() const
{
  return m_map->Comm().NumProc();
}

template<class V, class M>
const V&
uqVectorSpaceClass<V,M>::zeroVector() const
{
  UQ_FATAL_TEST_MACRO(m_zeroVector == NULL,
                      m_env.fullRank(),
                      "uqVectorSpaceClass<V,M>::zeroVector()",
                      "m_zeroVector is still NULL");
  return *m_zeroVector;
}

template <class V, class M>
unsigned int
uqVectorSpaceClass<V,M>::dimGlobal() const
{
  return m_dimGlobal;
}

template <class V, class M>
unsigned int
uqVectorSpaceClass<V,M>::globalIdOfFirstComponent() const
{
  return m_map->MinMyGID();
}

template <class V, class M>
unsigned int
uqVectorSpaceClass<V,M>::dimLocal() const
{
  return m_dimLocal;
}

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

template <class V, class M>
const EpetraExt::DistArray<std::string>* 
uqVectorSpaceClass<V,M>::componentsNames() const
{
  return m_componentsNames;
}

template <class V, class M>
const std::string&
uqVectorSpaceClass<V,M>::localComponentName(unsigned int localComponentId) const
{
  if (m_componentsNames == NULL) return m_emptyComponentName;

  UQ_FATAL_TEST_MACRO(localComponentId > m_dimLocal,
                      m_env.fullRank(),
                      "uqVectorSpaceClass<V,M>::localComponentName()",
                      "localComponentId is too big");

//return (*(const_cast<EpetraExt::DistArray<std::string>*>(m_componentsNames)))(localComponentId,0);
  return (*m_componentsNames)(localComponentId,0);
}

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

template <class V, class M>
void
uqVectorSpaceClass<V,M>::print(std::ostream& os) const
{
  return;
}
#endif // __UQ_VECTOR_SPACE_H__

