/* uq/libs/queso/inc/uqVectorSubset.h
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

#ifndef __UQ_VECTOR_SUBSET_H__
#define __UQ_VECTOR_SUBSET_H__

#include <uqVectorSet.h>
#include <EpetraExt_DistArray.h>

template <class V, class M>
class uqVectorSubsetClass : public uqVectorSetClass<V,M>
{
public:
           uqVectorSubsetClass();
           uqVectorSubsetClass(const uqBaseEnvironmentClass& env,
                              const char*                   prefix,
                              const uqVectorSpace<V,M>&     vectorSpace);
  virtual ~uqVectorSubsetClass();

           const uqVectorSpaceClass<V,M>& vectorSpace()                 const;
  virtual        bool                     contains   (const V& vec)     const = 0;
  virtual        void                     print      (std::ostream& os) const;

protected:
  using uqVectorSetClass<V,M>::m_env;
  using uqVectorSetClass<V,M>::m_prefix;

  const uqVectorSpaceClass<V,M>* m_vectorSpace;
};

template <class V, class M>
uqVectorSubsetClass<V,M>::uqVectorSubsetClass()
  :
  uqVectorSetClass<V,M>(m_env,m_prefix,0.),
  m_vectorSpace        (NULL)
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.rank(),
                      "uqVectorSubsetClass<V,M>::constructor(), default",
                      "should not be used by user");
}

template <class V, class M>
uqVectorSubsetClass<V,M>::uqVectorSubsetClass(
  const uqBaseEnvironmentClass& env,
  const char*                   prefix,
  const uqVectorSpace<V,M>&     vectorSpace)
  :
  uqVectorSetClass<V,M>(m_env,prefix,0.),
  m_vectorSpace        (&vectorSpace)
{
  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Entering uqVectorSubsetClass<V,M>::constructor()"
              << std::endl;
  }

  if ((m_env.verbosity() >= 5) && (m_env.rank() == 0)) {
    std::cout << "Leaving uqVectorSubsetClass<V,M>::constructor()"
              << std::endl;
  }
}

template <class V, class M>
uqVectorSubsetClass<V,M>::~uqVectorSubsetClass()
{
  //std::cout << "Entering uqVectorSubsetClass<V,M>::destructor()"
  //          << std::endl;

  if (m_vectorSpace != NULL) delete m_vectorSpace;

  //std::cout << "Leaving uqVectorSubsetClass<V,M>::destructor()"
  //          << std::endl;
}

template <class V, class M>
const uqVectorSpaceClass<V,M>&
uqVectorSubsetClass<V,M>::vectorSpace() const
{
  return *m_vectorSpace;
}

template <class V, class M>
void
uqVectorSubsetClass<V,M>::print(std::ostream& os) const
{
  return;
}

template<class V, class M>
std::ostream&
operator<<(std::ostream& os, const uqVectorSubsetClass<V,M>& obj)
{
  obj.print(os);

  return os;
}

#endif // __UQ_VECTOR_SUBSET_H__

