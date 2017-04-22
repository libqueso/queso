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

#include <queso/VectorSpace.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/TeuchosVector.h>
#include <queso/TeuchosMatrix.h>
#include <queso/DistArray.h>
#include <queso/Map.h>
#include <cmath>

namespace QUESO {

// Shaped constructor
template <class V, class M>
VectorSpace<V,M>::VectorSpace(const BaseEnvironment& env, const char* prefix,
    unsigned int dimGlobalValue,
    const std::vector<std::string>* componentsNamesVec)
  : VectorSet<V,M> (env,((std::string)(prefix) + "space_").c_str(),INFINITY),
    m_dimGlobal(dimGlobalValue),
    m_map(newMap()),
    m_dimLocal(m_map->NumMyElements()),
    m_componentsNamesArray(NULL),
    m_componentsNamesVec(NULL),
    m_emptyComponentName(""),
    m_zeroVector(new V(m_env,*m_map))
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering VectorSpace<V,M>::constructor(1)"
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
    std::cerr << "In VectorSpace<V,M>::constructor(1)"
              << ", with m_prefix = " << m_prefix
              << ": m_zeroVector->sizeGlobal() = " << m_zeroVector->sizeGlobal()
              << ", m_dimGlobal = "                << m_dimGlobal
              << std::endl;
  }
  queso_require_equal_to_msg(m_zeroVector->sizeGlobal(), m_dimGlobal, "global size of 'm_zeroVector' is not equal to m_dimGlobal");

  if (m_zeroVector->sizeLocal() != m_dimLocal) {
    std::cerr << "In VectorSpace<V,M>::constructor(1)"
              << ", with m_prefix = " << m_prefix
              << ": m_zeroVector->sizeLocal() = " << m_zeroVector->sizeLocal()
              << ", m_dimLocal = "                << m_dimLocal
              << std::endl;
  }
  queso_require_equal_to_msg(m_zeroVector->sizeLocal(), m_dimLocal, "local size of 'm_zeroVector' is not equal to m_dimLocal");

  if (componentsNamesVec != NULL) {
    queso_require_equal_to_msg(componentsNamesVec->size(), (size_t) m_dimGlobal, "global size of 'componentsNames' is not equal to m_dimGlobal");

    m_componentsNamesArray = new DistArray<std::string>(*m_map,1);
    unsigned int myFirstId = this->globalIdOfFirstComponent();
    for (unsigned int i = 0; i < m_dimLocal; ++i) {
      (*m_componentsNamesArray)(i,0) = (*componentsNamesVec)[myFirstId+i];
    }

    queso_require_equal_to_msg(m_componentsNamesArray->GlobalLength(), (int) m_dimGlobal, "global size of 'm_componentsNamesArray' is not equal to m_dimGlobal");
    queso_require_equal_to_msg(m_componentsNamesArray->MyLength(), (int) m_dimLocal, "local size of 'm_componentsNamesArray' is not equal to m_dimLocal");
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving VectorSpace<V,M>::constructor(1)"
                            << ", with m_prefix = " << m_prefix
                            << std::endl;
  }
}

// Copy constructor
template <class V, class M>
VectorSpace<V,M>::VectorSpace(const VectorSpace<V,M>& aux)
  : VectorSet<V,M>(aux.env(),((std::string)(aux.m_prefix)).c_str(),INFINITY),
    m_dimGlobal(aux.m_dimGlobal),
    m_map(newMap()),
    m_dimLocal(m_map->NumMyElements()),
    m_componentsNamesArray(NULL),
    m_componentsNamesVec(NULL),
    m_emptyComponentName(""),
    m_zeroVector(new V(m_env,*m_map))
{
}

// Destructor
template <class V, class M>
VectorSpace<V,M>::~VectorSpace()
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering VectorSpace<V,M>::destructor()"
                            << std::endl;
  }

  if (m_zeroVector           != NULL) delete m_zeroVector;
  if (m_componentsNamesVec   != NULL) delete m_componentsNamesVec;
  if (m_componentsNamesArray != NULL) delete m_componentsNamesArray;
  if (m_map                  != NULL) delete m_map;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving VectorSpace<V,M>::destructor()"
                            << std::endl;
  }
}

// Attribute methods
template <class V, class M>
const BaseEnvironment& VectorSpace<V,M>::env() const
{
  return m_env;
}

template <class V, class M>
const Map& VectorSpace<V,M>::map() const
{
  queso_require_msg(m_map, "m_map is still NULL");
  return *m_map;
}

template <class V, class M>
unsigned int VectorSpace<V,M>::numOfProcsForStorage() const
{
  return m_map->Comm().NumProc();
}

template <class V, class M>
unsigned int VectorSpace<V,M>::dimLocal() const
{
  return m_dimLocal;
}

template <class V, class M>
unsigned int VectorSpace<V,M>::dimGlobal() const
{
  return m_dimGlobal;
}

template <class V, class M>
unsigned int VectorSpace<V,M>::globalIdOfFirstComponent() const
{
  return m_map->MinMyGID();
}


template<class V, class M>
const V& VectorSpace<V,M>::zeroVector() const
{
  queso_require_msg(m_zeroVector, "m_zeroVector is still NULL");
  return *m_zeroVector;
}

template <class V, class M>
V* VectorSpace<V,M>::newVector(const V& v) const
{
  if (v.sizeGlobal() != m_dimGlobal) return NULL;
  if (v.sizeLocal () != m_dimLocal ) return NULL;

  return new V(v);
}

template <class V, class M>
M* VectorSpace<V,M>::newDiagMatrix(const V& v) const
{
  if (v.sizeGlobal() != m_dimGlobal) return NULL;
  if (v.sizeLocal () != m_dimLocal ) return NULL;

  return new M(v);
}

template <class V, class M>
M* VectorSpace<V,M>::newProposalMatrix(const V* varVec, const V* auxVec) const
{
  V tmpVec(*m_zeroVector);
  for (unsigned int i = 0; i < m_dimLocal; ++i) {
    double variance = INFINITY;
    if (varVec) variance = (*varVec)[i];
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In VectorSpace<V,M>::newProposalMatrix()"
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
const VectorSpace<V,M>& VectorSpace<V,M>::vectorSpace() const
{
  return *this;
}

template <class V, class M>
bool VectorSpace<V,M>::contains(const V& /* vec */) const
{
  return true;
}

template <class V, class M>
void VectorSpace<V,M>::centroid(V& vec) const
{
  queso_assert_equal_to (m_dimLocal, vec.sizeLocal());

  for (unsigned int i = 0; i < m_dimLocal; ++i) {
    vec[i] = INFINITY;
  }
}

template <class V, class M>
void VectorSpace<V,M>::moments(M& mat) const
{
  queso_assert_equal_to (m_dimLocal, mat.numCols());

  mat.zeroLower();
  mat.zeroUpper();
  for (unsigned int i = 0; i < m_dimLocal; ++i) {
    mat(i,i) = INFINITY;
  }
}

template <class V, class M>
const DistArray<std::string>* VectorSpace<V,M>::componentsNamesArray() const
{
  return m_componentsNamesArray;
}

template <class V, class M>
const std::string& VectorSpace<V,M>::localComponentName(
    unsigned int localComponentId) const
{
  if (m_componentsNamesArray == NULL) return m_emptyComponentName;

  queso_require_less_equal_msg(localComponentId, m_dimLocal, "localComponentId is too big");

//return (*(const_cast<DistArray<std::string>*>(m_componentsNamesArray)))(localComponentId,0);
  return (*m_componentsNamesArray)(localComponentId,0);
}

template<class V, class M>
void VectorSpace<V,M>::printComponentsNames(std::ostream& os,
    bool printHorizontally) const
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
void VectorSpace<V,M>::print(std::ostream& os) const
{
  os << "In VectorSpace<V,M>::print()"
     << ": nothing to be printed" << std::endl;
  return;
}

}  // End namespace QUESO

template class QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix>;
#ifdef QUESO_HAS_TRILINOS
template class QUESO::VectorSpace<QUESO::TeuchosVector, QUESO::TeuchosMatrix>;
#endif
