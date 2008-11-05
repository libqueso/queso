/* uq/libs/queso/inc/uqVectorSpace.h
 *
 * Copyright (C) 2008 The QUESO Team, http://queso.ices.utexas.edu/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __UQ_VECTOR_SPACE_H__
#define __UQ_VECTOR_SPACE_H__

#include <uqEnvironment.h>
#include <uqMiscellaneous.h>
//#include <vector>
//#include <iostream>
//#include <fstream>
#include <uqDefines.h>
#include <EpetraExt_DistArray.h>

#undef UQ_VECTOR_SPACE_READS_FILE_OPTIONS

template <class V, class M>
class uqVectorSpaceClass
{
public:
          uqVectorSpaceClass();
          uqVectorSpaceClass(const uqEnvironmentClass&                env, // See template specialization
                             const char*                              prefix,
                             unsigned int                             dimValue,
                             const EpetraExt::DistArray<std::string>* componentsNames);
         ~uqVectorSpaceClass();

  const   uqEnvironmentClass&                env                 ()                         const;
  const   Epetra_Map&                        map                 ()                         const;
          unsigned int                       dim                 ()                         const;

  const   V&                                 zeroVector          ()                         const;
          V*                                 newVector           ()                         const; // See template specialization
          V*                                 newVector           (double value)             const; // See template specialization
          V*                                 newVector           (const V& v)               const;
          M*                                 newMatrix           ()                         const; // See template specialization
          M*                                 newDiagMatrix       (const V& v)               const;
          M*                                 newDiagMatrix       (double diagValue)         const; // See template specialization
          M*                                 newGaussianMatrix   (const V& varianceValues,
                                                                  const V& initialValues)   const;

  const   std::string&                       componentName       (unsigned int componentId) const;
          void                               printComponentsNames(std::ostream& os, bool printHorizontally) const;
          void                               print               (std::ostream& os) const;

protected:
#ifdef UQ_VECTOR_SPACE_READS_FILE_OPTIONS
          void                               defineMyOptions     (po::options_description& optionsDesc) const;
          void                               getMyOptionValues   (po::options_description& optionsDesc);
#endif

  const   uqEnvironmentClass&                m_env;
          std::string                        m_prefix;
          unsigned int                       m_dim;
  const   EpetraExt::DistArray<std::string>* m_componentsNames;
          std::string                        m_emptyComponentName;

  const   Epetra_Map*                        m_map;
          V*                                 m_zeroVector;

#ifdef UQ_VECTOR_SPACE_READS_FILE_OPTIONS
          po::options_description*           m_optionsDesc;
          std::string                        m_option_help;
          std::string                        m_option_dim;
#endif
};

template <class V, class M>
uqVectorSpaceClass<V,M>::uqVectorSpaceClass()
  :
  m_env(*(new uqEnvironmentClass()))
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.rank(),
                      "uqVectorSpaceClass<V,M>::constructor(), default",
                      "should not be used by user");
}

template <class V, class M>
uqVectorSpaceClass<V,M>::uqVectorSpaceClass(
  const uqEnvironmentClass& env,
  const char*               prefix,
        unsigned int        dimValue,
  const EpetraExt::DistArray<std::string>* componentsNames)
  :
  m_env               (env),
  m_prefix            ((std::string)(prefix) + "space_"),
  m_dim               (dimValue),
  m_componentsNames   (componentsNames),
  m_emptyComponentName(""),
#ifdef UQ_VECTOR_SPACE_READS_FILE_OPTIONS
  m_map               (NULL),
  m_zeroVector        (NULL)
  m_optionsDesc       (new po::options_description("Vector space options")),
  m_option_help       (m_prefix + "help"),
  m_option_dim        (m_prefix + "dim" )
#else
  m_map               (new Epetra_Map(m_dim,0,m_env.comm())),
  m_zeroVector        (new V(m_env,*m_map))
#endif
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqVectorSpaceClass<V,M>::constructor()"
              << std::endl;
  }

  UQ_FATAL_TEST_MACRO((m_componentsNames != NULL) && (m_componentsNames->GlobalLength() != (int) m_dim),
                      m_env.rank(),
                      "uqVectorSpaceClass<V,M>::constructor()",
                      "size of 'componentsNames' is not equal to m_dim");

#ifdef UQ_VECTOR_SPACE_READS_FILE_OPTIONS
  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.rank() == 0) std::cout << "After getting values of options with prefix '" << m_prefix
                                   << "', state of uqVectorSpaceClass object is:"
                                   << "\n" << *this
                                   << std::endl;

  m_map        = new Epetra_Map(m_dim,0,m_env.comm());
  m_zeroVector = new V(m_env,*m_map);
#endif
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

#ifdef UQ_VECTOR_SPACE_READS_FILE_OPTIONS
template <class V, class M>
void
uqVectorSpaceClass<V,M>::defineMyOptions(
  po::options_description& optionsDesc) const
{
  m_optionsDesc->add_options()
    (m_option_help.c_str(),                                                   "produce help message for vector space")
    (m_option_dim.c_str(),      po::value<unsigned int>()->default_value(0),  "Space dimension"                      )
  ;

  return;
}

template <class V, class M>
void
uqVectorSpaceClass<V,M>::getMyOptionValues(po::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count(m_option_help.c_str())) {
    std::cout << optionsDesc
              << std::endl;
  }

  if (m_env.allOptionsMap().count(m_option_dim.c_str())) {
    const po::variables_map& tmpMap = m_env.allOptionsMap();
    m_dim = tmpMap[m_option_dim.c_str()].as<unsigned int>();
  }

  return;
}
#endif

template <class V, class M>
const uqEnvironmentClass&
uqVectorSpaceClass<V,M>::env() const
{
  return m_env;
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
#ifdef UQ_VECTOR_SPACE_READS_FILE_OPTIONS
  os << m_option_dim << " = " << m_dim
     << std::endl;
#endif

  return;
}

template<class V, class M>
std::ostream&
operator<<(std::ostream& os, const uqVectorSpaceClass<V,M>& obj)
{
  obj.print(os);

  return os;
}

#endif // __UQ_VECTOR_SPACE_H__

