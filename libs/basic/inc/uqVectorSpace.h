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
          uqVectorSpaceClass(const uqBaseEnvironmentClass&            env,
                             const char*                              prefix,
                             unsigned int                             dimValue,
                             const EpetraExt::DistArray<std::string>* componentsNames);
         ~uqVectorSpaceClass();

  const   Epetra_Map&              map                 ()                         const;
          unsigned int             dim                 ()                         const;

  const   V&                       zeroVector          ()                         const;
          V*                       newVector           ()                         const; // See template specialization
          V*                       newVector           (double value)             const; // See template specialization
          V*                       newVector           (const V& v)               const;
          M*                       newMatrix           ()                         const; // See template specialization
          M*                       newDiagMatrix       (const V& v)               const;
          M*                       newDiagMatrix       (double diagValue)         const; // See template specialization
          M*                       newGaussianMatrix   (const V& varianceValues,
                                                        const V& initialValues)   const;

  const   uqVectorSpaceClass<V,M>& vectorSpace         ()                         const;
          bool                     contains            (const V& vec)             const;

  const   EpetraExt::DistArray<std::string>* componentsNames() const;
  const   std::string&             componentName       (unsigned int componentId) const;
          void                     printComponentsNames(std::ostream& os, bool printHorizontally) const;
          void                     print               (std::ostream& os) const;

protected:
          using uqVectorSetClass<V,M>::m_env;
          using uqVectorSetClass<V,M>::m_prefix;
          using uqVectorSetClass<V,M>::m_volume;

          unsigned int                       m_dim;
  const   EpetraExt::DistArray<std::string>* m_componentsNames;
          std::string                        m_emptyComponentName;

  const   Epetra_Map*                        m_map;
          V*                                 m_zeroVector;
};

template <class V, class M>
uqVectorSpaceClass<V,M>::uqVectorSpaceClass()
  :
  uqVectorSetClass<V,M>()
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.rank(),
                      "uqVectorSpaceClass<V,M>::constructor(), default",
                      "should not be used by user");
}

template <class V, class M>
uqVectorSpaceClass<V,M>::uqVectorSpaceClass(
  const uqBaseEnvironmentClass&            env,
  const char*                              prefix,
        unsigned int                       dimValue,
  const EpetraExt::DistArray<std::string>* componentsNames)
  :
  uqVectorSetClass<V,M>(env,((std::string)(prefix) + "space_").c_str(),INFINITY),
  m_dim                (dimValue),
  m_componentsNames    (componentsNames),
  m_emptyComponentName (""),
  m_map                (new Epetra_Map(m_dim,0,m_env.worldComm())),
  m_zeroVector         (new V(m_env,*m_map))
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqVectorSpaceClass<V,M>::constructor()"
              << std::endl;
  }

  UQ_FATAL_TEST_MACRO((m_componentsNames != NULL) && (m_componentsNames->GlobalLength() != (int) m_dim),
                      m_env.rank(),
                      "uqVectorSpaceClass<V,M>::constructor()",
                      "size of 'componentsNames' is not equal to m_dim");

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqVectorSpaceClass<V,M>::constructor()"
              << std::endl;
  }
}

template <class V, class M>
uqVectorSpaceClass<V,M>::~uqVectorSpaceClass()
{
  //std::cout << "Entering uqVectorSpaceClass<V,M>::destructor()"
  //          << std::endl;

  if (m_zeroVector != NULL) delete m_zeroVector;
  if (m_map        != NULL) delete m_map;

  //std::cout << "Leaving uqVectorSpaceClass<V,M>::destructor()"
  //          << std::endl;
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
const Epetra_Map&
uqVectorSpaceClass<V,M>::map() const
{
  UQ_FATAL_TEST_MACRO(m_map == NULL,
                      m_env.rank(),
                      "uqVectorSpaceClass<V,M>::map()",
                      "m_map is still NULL");
  return *m_map;
}

template<class V, class M>
const V&
uqVectorSpaceClass<V,M>::zeroVector() const
{
  UQ_FATAL_TEST_MACRO(m_zeroVector == NULL,
                      m_env.rank(),
                      "uqVectorSpaceClass<V,M>::zeroVector()",
                      "m_zeroVector is still NULL");
  return *m_zeroVector;
}

template <class V, class M>
unsigned int
uqVectorSpaceClass<V,M>::dim() const
{
  return m_dim;
}

template <class V, class M>
V*
uqVectorSpaceClass<V,M>::newVector(const V& v) const
{
  if (v.size() != m_dim) return NULL;

  return new V(v);
}

template <class V, class M>
M*
uqVectorSpaceClass<V,M>::newDiagMatrix(const V& v) const
{
  if (v.size() != m_dim) return NULL;

  return new M(v);
}

template <class V, class M>
M*
uqVectorSpaceClass<V,M>::newGaussianMatrix(
  const V& varianceValues,
  const V& initialValues) const
{
  V tmpVec(*m_zeroVector);
  for (unsigned int i = 0; i < m_dim; ++i) {
    double variance = varianceValues[i];
    std::cout << "In uqVectorSpaceClass<V,M>::newGaussianMatrix()"
              << ": i = "        << i
              << ", variance = " << variance
              << std::endl;
    if ((variance == INFINITY) ||
        (variance == NAN     )) {
      tmpVec[i] = pow( fabs(initialValues[i])*0.05,2. );
      if ( tmpVec[i] == 0 ) tmpVec[i] = 1.;
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
uqVectorSpaceClass<V,M>::componentName(unsigned int componentId) const
{
  if (m_componentsNames == NULL) return m_emptyComponentName;

  UQ_FATAL_TEST_MACRO(componentId > m_dim,
                      m_env.rank(),
                      "uqVectorSpaceClass<V,M>::componentName()",
                      "componentId is too big");

  return (*(const_cast<EpetraExt::DistArray<std::string>*>(m_componentsNames)))(componentId,0);
}

template<class V, class M>
void
uqVectorSpaceClass<V,M>::printComponentsNames(std::ostream& os, bool printHorizontally) const
{
  if (printHorizontally) { 
    for (unsigned int i = 0; i < this->dim(); ++i) {
      os << "'" << this->componentName(i) << "'"
         << " ";
    }
  }
  else {
    for (unsigned int i = 0; i < this->dim(); ++i) {
      os << "'" << this->componentName(i) << "'"
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

